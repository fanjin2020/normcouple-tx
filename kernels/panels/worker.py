"""Small independently sealed batches; only completed checksummed batches resume."""
from pathlib import Path
import tempfile
import time

import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from extension_stats import evaluate
from selection import CASES, METHODS, stress
from support import atomic_json, checkpoint, seal, task_folder

HERE = Path(__file__).resolve().parent
JOBS, OPTIONS = [], {}


def run_batch(task):
    unit, start, block = task
    options = OPTIONS
    job = JOBS[unit["job_index"]]
    folder = task_folder(options["output"], task)
    old = checkpoint(folder, options["signature"])
    if old:
        return dict(old, resumed=True)
    folder.mkdir(parents=True, exist_ok=True)
    clock = time.monotonic()
    reps = list(range(start, min(start + options["batch_size"], options["repetitions"])))
    subsets = options["subset_ids"][block * options["subset_batch_size"]:(block + 1) * options["subset_batch_size"]]
    atomic_json(folder / "started.json", dict(started_unix=time.time(), unit=unit["unit_id"], reps=reps, subsets=subsets))
    records, diagnostics, baseline_checks = [], [], []
    with tempfile.TemporaryDirectory(prefix="working_", dir=folder) as temporary:
        work = Path(temporary)
        entries, payload = [], {}
        for rep in reps:
            base = C.draw_base(job, unit, rep, options)
            A, marker_yield, anchor = marker_properties(job, base)
            key = "r%d" % rep
            C.write_bin(work / (key + "_raw.bin"), base["draw"]["counts"][:, job["norm_genes"]].T)
            C.write_bin(work / (key + "_A.bin"), A)
            C.write_bin(work / (key + "_total.bin"), base["draw"]["mass"].sum(1))
            outcomes, r_cases, r_subsets = {}, [], []
            for name, effect in CASES:
                condition = dict(name=name, effect=effect, pattern="up", panel="frozen")
                data = C.intervene(job, base, condition, options)
                outcomes[name] = data
                filename = key + "_" + name + "_add.bin"
                C.write_bin(work / filename, data["additions"].T)
                r_cases.append(dict(name=name, additions=filename))
                if rep == 0 and name in ("null", "up05"):
                    errors = C.verify_baseline(job, unit, rep, options, data, effect)
                    baseline_checks.append(dict(rep=rep, outcome=name, errors=errors))
            for subset_id in subsets:
                additions, detail = stress(job, base, subset_id)
                filename = key + "_" + subset_id + "_stress.bin"
                C.write_bin(work / filename, additions.T)
                r_subsets.append(dict(name=subset_id, stress=filename))
                diagnostics.append(dict(rep=rep, subset_id=subset_id, **detail))
            entries.append(dict(key=key, n=unit["n"], raw=key + "_raw.bin", A=key + "_A.bin",
                                total=key + "_total.bin", marker_anchor=anchor, outcomes=r_cases, subsets=r_subsets))
            payload[rep] = outcomes, marker_yield, base["ref"]["pooled_yield"]
        C.write_json(work / "jobs.json", dict(reference_folder=job["reference_folder"], entries=entries))
        r = options["r_info"]
        D.run_r([r["executable"], "--vanilla", str(HERE / "bridge.R"), str(work), r["library"]],
                r["library"], folder / "bridge.log", options["r_timeout"])
        outputs = pd.read_csv(work / "outputs.csv", keep_default_na=False)
        expected = {(f"r{rep}_{subset}_{name}", method) for rep in reps for subset in subsets
                    for name, _ in CASES for method in METHODS if method != "oracle"}
        C.require(not outputs.duplicated(["key", "method"]).any(), "Duplicate R prediction")
        C.require(set(zip(outputs.key, outputs.method)) == expected, "Missing or unexpected R predictions")
        predictions, loaded = {}, {}
        for row in outputs.itertuples():
            path = (work / row.file).resolve()
            C.require(path.is_relative_to(work.resolve()), "Unsafe prediction path")
            if row.file not in loaded:
                loaded[row.file] = np.loadtxt(path, delimiter=",", ndmin=2)
            q = loaded[row.file]
            C.require(q.shape == (unit["n"], 7), "Prediction shape changed")
            predictions[row.key, row.method] = q, row.reason
        for rep in reps:
            outcomes, marker_yield, full_yield = payload[rep]
            oracle = {name: evaluate(data, data["p"], "full_six") for name, data in outcomes.items()}
            for subset_id in subsets:
                for name, effect in CASES:
                    data = outcomes[name]
                    for method in METHODS:
                        if method == "oracle":
                            row, arrays, fit = oracle[name]
                        else:
                            q, reason = predictions[f"r{rep}_{subset_id}_{name}", method]
                            props = C.E.convert_available_rows(q, marker_yield if method == "marker_matched" else full_yield)
                            row, arrays, fit = evaluate(data, props, "full_six", reason)
                        valid = arrays["available"].astype(bool)
                        hit = valid & (arrays["q"] <= .05)
                        null = data["truth"] == 0
                        counts = dict(planned_tests=len(null), null_tests=int(null.sum()),
                            alternative_tests=int((~null).sum()), valid_tests=int(valid.sum()),
                            covered_tests=int(arrays["covered"].sum()),
                            discoveries=int(hit.sum()), false_discoveries=int((hit & null).sum()),
                            true_discoveries=int((hit & ~null).sum()))
                        records.append(dict(rep=rep, subset_id=subset_id, outcome=name, effect=effect,
                                            method=method, **counts, **row))
        checks = pd.read_csv(work / "cache_checks.csv")
        C.require(len(checks) == len(reps) * len(subsets) * len(CASES) * 2, "Cache audit incomplete")
        C.require(checks.marker_input_identical.all(), "Direct effects changed marker input")
        checked = checks[checks.fresh_repeat_checked & checks.prediction_finite]
        C.require(checked.fresh_repeat_max_difference.lt(1e-12).all(), "Fresh prediction disagrees with cache")
        checks.to_csv(folder / "cache_checks.csv.gz", index=False, compression="gzip")
        outputs.to_csv(folder / "bridge_diagnostics.csv.gz", index=False, compression="gzip")
    frame = pd.DataFrame(records)
    C.require(len(frame) == len(reps) * len(subsets) * len(CASES) * len(METHODS), "Incomplete metric grid")
    C.require(not frame.duplicated(["rep", "subset_id", "outcome", "method"]).any(), "Duplicate metric grid")
    frame.to_csv(folder / "replicates.csv.gz", index=False, compression="gzip")
    pd.DataFrame(diagnostics).to_csv(folder / "stress_diagnostics.csv.gz", index=False, compression="gzip")
    atomic_json(folder / "baseline_checks.json", baseline_checks)
    meta = dict(unit=unit["unit_id"], start=start, block=block, repetitions=len(reps), subsets=subsets,
                work_items=len(reps) * len(subsets), elapsed_seconds=time.monotonic() - clock,
                finished_unix=time.time(), metric_rows=len(frame),
                fresh_cache_checks=int(checks.fresh_repeat_checked.sum()))
    atomic_json(folder / "batch.json", meta)
    names = ["replicates.csv.gz", "stress_diagnostics.csv.gz", "cache_checks.csv.gz",
             "bridge_diagnostics.csv.gz", "baseline_checks.json", "batch.json", "bridge.log"]
    seal(folder, names, options["signature"], **meta)
    return dict(meta, resumed=False)
