"""Isolated, identity-checked extension of the frozen Paper 3 experiment."""
import importlib.metadata
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import signal
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import revision_core as C
import revision_driver as D
import worker
from selection import AMPLITUDE, CASES, FRACTION, METHODS, SELECTION_SEED, STRESS_SEED, freeze
from support import (PACKAGE, VERSION, atomic_json, checkpoint, digest, file_sha,
                     read_json, state, task_folder, tasks_for, task_work)

HERE = Path(__file__).resolve().parent


def source_hashes():
    manifest = read_json(HERE / "source_manifest.json")
    for name, expected in manifest.items():
        path = (HERE / name).resolve()
        C.require(path.is_relative_to(HERE) and file_sha(path) == expected, "Source checksum mismatch: " + name)
    return manifest


def r_environment(args, out):
    if args.r_library is None:
        candidates = [args.project_root.parent / name / "R_library" for name in
                      ("paper3_revision_results_v1", "paper3_strengthening_results_v1", "paper3_mechanism_results_v1")]
        if args.old_results:
            candidates.insert(0, args.old_results / "R_library")
        args.r_library = next((p for p in candidates if (p / "dtangle").is_dir() and (p / "edgeR").is_dir()), out / "R_library")
    exe = shutil.which(args.rscript)
    C.require(exe is not None, "Rscript missing: specify --rscript. No research input downloads are attempted.")
    library = args.r_library.resolve()
    library.mkdir(parents=True, exist_ok=True)
    print("CHECK R", exe, "LIBRARY", library, flush=True)
    path = out / "R_versions_current.json"
    D.run_r([exe, "--vanilla", str(HERE / "bridge.R"), "--check", str(library), str(path)],
            library, out / "R_check.log", timeout=180)
    return dict(executable=exe, library=str(library), versions=read_json(path))


def load_jobs(args, out):
    if args.synthetic:
        return D.synthetic_jobs(out), []
    args.stage2_results = None
    if args.input_bundle is None:
        portable = args.project_root.parent / "paper3_mechanism_results_v1_prepared_inputs"
        if (portable / "inputs_manifest.json").exists():
            args.input_bundle = portable
    if args.external_cache is None:
        args.external_cache = args.project_root / "01_data"
    # This resolver has no network fallback. Missing or altered inputs stop the run.
    folders, gene_map, meta, counts, identities = D.resolve_all_inputs(args)
    sources, jobs = C.full_inputs(folders, gene_map)
    external = C.S.read_external(counts, meta, next(iter(sources.values()))["genes"], out / "external_inputs")
    jobs += C.S.make_external_jobs(jobs, external, out)
    jobs = [C.prepare_job(job, out, 2000) for job in jobs]
    atomic_json(out / "input_paths.json", identities)
    return jobs, identities


def execute_pending(tasks, args, cfg, out, done, work_done, total, total_work):
    started = time.monotonic()
    launch_work = 0

    def update(result):
        nonlocal done, work_done, launch_work
        done += 1
        work_done += result["work_items"]
        launch_work += result["work_items"]
        elapsed = time.monotonic() - started
        eta = (total_work - work_done) * elapsed / launch_work / 3600 if launch_work else None
        state(out, "RUNNING", workers=args.workers, signature=cfg["signature"], completed_batches=done,
              total_batches=total, completed_work=work_done, total_work=total_work, simulation_eta_hours=eta)
        print("DONE %d/%d | %.1f%% pairs | %s s | rough remaining %.2f h" %
              (done, total, 100 * work_done / total_work, round(result["elapsed_seconds"], 1), eta or 0), flush=True)

    if args.workers == 1:
        for task in tasks:
            update(worker.run_batch(task))
        return
    pool = ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("fork"))
    pending = {}
    iterator = iter(tasks)

    def enqueue():
        while len(pending) < 2 * args.workers:
            try:
                task = next(iterator)
            except StopIteration:
                break
            pending[pool.submit(worker.run_batch, task)] = task

    try:
        enqueue()
        while pending:
            completed, _ = wait(pending, timeout=30, return_when=FIRST_COMPLETED)
            for future in completed:
                del pending[future]
                update(future.result())
            enqueue()
    except BaseException:
        for future in pending:
            future.cancel()
        print("Stopping this worker pool; sealed batches remain reusable. Waiting for active workers to exit.", flush=True)
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)


def main(args):
    for name in ("workers", "repetitions", "subsets", "batch_size", "subset_batch_size", "r_timeout"):
        C.require(getattr(args, name) > 0, name + " must be positive")
    C.require(os.name != "nt" or args.workers == 1, "Use --workers 1 on Windows. Linux uses fork/shared count banks.")
    C.require(args.synthetic or (args.subsets == 50 and args.repetitions >= 200),
              "Research protocol requires 50 subsets and >=200 repetitions; --synthetic is software testing only")
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        C.require((out / "OWNER.json").exists() and read_json(out / "OWNER.json").get("package") == PACKAGE,
                  "Use a new output directory. Previous manuscript/results packages are never overwritten.")
    out.mkdir(parents=True, exist_ok=True)
    C.E.hold_output_lock(out)
    atomic_json(out / "OWNER.json", dict(package=PACKAGE, version=VERSION))
    def terminate(signum, frame):
        raise KeyboardInterrupt("Termination signal received; resume uses sealed batches only")
    if os.name != "nt":
        signal.signal(signal.SIGTERM, terminate)
    try:
        state(out, "PREFLIGHT", workers=args.workers)
        source = source_hashes()
        from run import PACKAGES
        python = {name: importlib.metadata.version(name) for name in PACKAGES}
        C.require(python == PACKAGES, "Pinned Python versions required; run with --setup")
        if args.summarize_only:
            cfg = read_json(out / "protocol.json")
            C.require(cfg["source_hashes"] == source, "Source differs from stored protocol")
        else:
            r = r_environment(args, out)
            print("Loading cached, checksum-verified count banks. No research downloads.", flush=True)
            jobs, identities = load_jobs(args, out)
            number = min(args.subsets, 3) if args.synthetic else args.subsets
            selections = freeze(jobs, number)
            units = C.plan_units(jobs, smoke=args.synthetic)
            C.require(args.synthetic or len(units) == 42, "Frozen research design requires 42 reference-target units")
            for unit in units:
                unit["conditions"] = [dict(name=name, effect=effect) for name, effect in CASES]
            settings = dict(C.SEEDS, repetitions=min(args.repetitions, 2) if args.synthetic else args.repetitions,
                batch_size=args.batch_size, subset_batch_size=args.subset_batch_size,
                subset_ids=list(jobs[0]["stress_subsets"]), random_subsets=number,
                selection_seed=SELECTION_SEED, stress_seed=STRESS_SEED, requested_fraction=FRACTION,
                amplitude=AMPLITUDE, outcomes=[dict(name=n, effect=e) for n, e in CASES],
                methods=METHODS, regression_rule="full_six", pseudocount=1.0,
                software_fixture=args.synthetic, selection_scope="within each fixed annotated reference panel",
                failed_tests="included in planned denominators; p=1 and no rejection; coverage reported separately")
            body = dict(package=PACKAGE, version=VERSION, settings=settings, units=units, source_hashes=source,
                        input_hashes={x["name"]: x["sha256"] for x in identities},
                        selection_sha256=digest(selections), R=r["versions"], python=python)
            cfg = dict(body, signature=digest(body))
            if (out / "protocol.json").exists():
                C.require(read_json(out / "protocol.json")["signature"] == cfg["signature"],
                          "Code/input/settings/R versions changed. Use a NEW output directory; workers alone may change.")
            atomic_json(out / "protocol.json", cfg)
            atomic_json(out / "selected_markers.json", selections)
            worker.JOBS = jobs
            worker.OPTIONS = dict(settings, output=str(out), r_info=r, r_timeout=args.r_timeout, signature=cfg["signature"])
            tasks = list(tasks_for(cfg))
            todo, work_done = [], 0
            print("Validating completed checkpoints before resume", flush=True)
            for task in tasks:
                receipt = checkpoint(task_folder(out, task), cfg["signature"])
                if receipt:
                    work_done += receipt["work_items"]
                else:
                    todo.append(task)
            total_work = sum(task_work(cfg, task) for task in tasks)
            done = len(tasks) - len(todo)
            print("PLAN", len(units), "units;", len(settings["subset_ids"]), "subsets including 2 controls;",
                  settings["repetitions"], "repetitions;", len(todo), "remaining batches", flush=True)
            if todo:
                # The pilot tests infrastructure/identities, not favorable FDR or power.
                state(out, "NUMERICAL_PILOT", workers=1, signature=cfg["signature"])
                first = worker.run_batch(todo.pop(0))
                done += 1
                work_done += first["work_items"]
                print("PILOT COMPLETE. Measured seconds per repetition-subset:",
                      round(first["elapsed_seconds"] / first["work_items"], 3), flush=True)
            if args.pilot_only:
                state(out, "PILOT_COMPLETE", workers=args.workers, signature=cfg["signature"],
                      completed_batches=done, total_batches=len(tasks))
                print("Re-run without --pilot-only to continue the same fixed protocol.", flush=True)
                return 0
            state(out, "RUNNING", workers=args.workers, signature=cfg["signature"])
            execute_pending(todo, args, cfg, out, done, work_done, len(tasks), total_work)
        state(out, "SUMMARIZING", workers=args.workers, signature=cfg["signature"])
        from summarize import summarize
        final = summarize(out, cfg)
        status = "SOFTWARE_FIXTURE_ONLY" if cfg["settings"]["software_fixture"] else "COMPLETE"
        receipt = dict(status=status, signature=cfg["signature"], workers=args.workers,
                       public_deposition="NOT_PERFORMED", research_inputs_downloaded=False)
        state(out, **receipt)
        shutil.copy2(out / "RUN_STATUS.json", final / "RUN_STATUS.json")
        atomic_json(final / "files_sha256.json", {p.name: file_sha(p) for p in sorted(final.iterdir())
                                                if p.is_file() and p.name != "files_sha256.json"})
        print("FINISHED", status, "\nCOLLECT ONLY:", final, flush=True)
        return 0
    except BaseException as exc:
        state(out, "INTERRUPTED_OR_FAILED", workers=args.workers, error=type(exc).__name__ + ": " + str(exc))
        raise


if __name__ == "__main__":
    from cli import parser
    raise SystemExit(main(parser().parse_args()))
