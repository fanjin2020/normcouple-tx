"""One bounded Paper 3 follow-up. Existing experiments remain read-only."""
import argparse
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import signal
import statistics
import sys
import time
from next_core import (BASE, HERE, SERVER, require, read_json, save_json, digest, sha,
                       verify_source, find_parent, parent_read_lock, verify_parent,
                       choose_python, assert_separate, write_manifest, export_review)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=SERVER/'paper1_public_evidence')
    p.add_argument('--previous-results', type=Path, help='Full bdd118b137 server run, not its review export')
    p.add_argument('--random-results', type=Path, default=SERVER/'paper3_random_marker_results_v1')
    p.add_argument('--output', type=Path, default=SERVER/'p3_next_v1')
    p.add_argument('--workers', type=int, default=36)
    p.add_argument('--repetitions', type=int, default=50, help='New repetitions per reference-target unit; frozen on first run')
    p.add_argument('--python', type=Path)
    p.add_argument('--rscript')
    p.add_argument('--r-library', type=Path, action='append', default=[])
    p.add_argument('--r-timeout', type=int, default=14400)
    p.add_argument('--cache', type=Path, help='Existing verified GSE125583 public download cache')
    for name in ('input-bundle', 'libd-dir', 'seaad-dir', 'gene-map', 'external-cache'):
        p.add_argument('--'+name, type=Path)
    for name in ('status', 'verify-status', 'export-only', 'preflight-only', 'reports-only'):
        p.add_argument('--'+name, action='store_true')
    return p


def tasks(cfg):
    settings = cfg['settings']
    for rep in range(settings['rep_start'], settings['rep_start']+settings['repetitions']):
        for unit in cfg['units']:
            yield unit, rep


def status(args):
    from support import checkpoint, lock_is_active
    if not (args.output/'protocol.json').is_file():
        print('NOT_STARTED', args.output)
        return 0
    cfg = read_json(args.output/'protocol.json')
    state = read_json(args.output/'STATE.json') if (args.output/'STATE.json').is_file() else {}
    planned = list(tasks(cfg))
    timing = []
    for unit, rep in planned:
        folder = args.output/'support_batches'/unit['unit_id']/('r%05d' % rep)
        receipt = checkpoint(folder, cfg['signature'], verify=args.verify_status)
        if receipt:
            timing.append(receipt['elapsed_seconds'])
    active = lock_is_active(args.output)
    print('STATE', state.get('status'), 'WRITER LOCK', active, 'WORKERS', state.get('workers'))
    for name in ('corrections', 'clinical_multistep', 'common_support'):
        print(name, 'COMPLETE' if checkpoint(args.output/name, cfg['signature'], verify=args.verify_status) else 'NOT_COMPLETE')
    print('NEW SUPPORT BATCHES', str(len(timing))+'/'+str(len(planned)))
    if active and len(timing) >= 3:
        eta = (len(planned)-len(timing))*statistics.median(timing)/max(1, state.get('workers', 1))/3600
        print('ROUGH SIMULATION-ONLY HOURS', round(eta, 2), '(excludes clinical checks and reporting)')
    elif not active:
        print('No active writer. The same run command resumes sealed work. Do not delete lock files.')
    else:
        print('ETA awaits at least three completed NEW batches.')
    return 0


def load_jobs(args, parent_cfg):
    import revision_core as C
    import driver as frozen
    from selection import freeze
    from matched_markers import freeze_matched
    from runtime import recorded_inputs
    recorded_inputs(args)
    args.synthetic = False
    args.old_results = args.random_results
    jobs, identities = frozen.load_jobs(args, args.output/'prepared')
    hashes = {v['name']: v['sha256'] for v in identities}
    require(hashes == parent_cfg['input_hashes'], 'Prepared banks do not match the parent experiment')
    selections = freeze(jobs, parent_cfg['settings']['random_subsets'])
    selections += freeze_matched(jobs, parent_cfg['settings']['random_subsets'], args.output/'matching')
    require(digest(selections) == parent_cfg['selection_sha256'], 'Marker sets differ from the parent experiment')
    require(C.plan_units(jobs, smoke=False) == parent_cfg['units'], 'Reference/target design changed')
    for job in jobs:
        job['budget_allocation'] = 'reference'
    save_json(args.output/'selected_markers.json', selections)
    return jobs


def make_review(output, cfg, run_status):
    from support import checkpoint
    review = output/'review'
    review.mkdir(exist_ok=True)
    stages = {}
    for name in ('corrections', 'clinical_multistep', 'common_support'):
        ready = bool(checkpoint(output/name, cfg['signature']))
        stages[name] = ready
        if ready:
            target = review/name
            target.mkdir(exist_ok=True)
            # Detailed simulation batches and clinical binary matrices stay on the server.
            for p in sorted((output/name).iterdir()):
                if p.is_file() and p.suffix in ('.csv', '.gz', '.json') and p.name != 'complete.json':
                    shutil.copy2(p, target/p.name)
            save_json(target/'server_receipt_identity.json', dict(
                receipt_sha256=sha(output/name/'complete.json'), signature=cfg['signature'],
                scope='Verified full server stage; compact review integrity is given by the outer files_sha256.json'))
    for name in ('protocol.json', 'parent_identity.json', 'runtime.json', 'R_versions.json', 'selected_markers.json'):
        if (output/name).is_file():
            shutil.copy2(output/name, review/name)
    if (output/'matching').is_dir():
        target = review/'marker_matching'
        target.mkdir(exist_ok=True)
        for p in (output/'matching').glob('*.csv'):
            shutil.copy2(p, target/p.name)
    shutil.copy2(HERE/'README.md', review/'READ_RESULTS.md')
    shutil.copy2(HERE/'source_sha256.json', review/'source_sha256.json')
    # Small, complete source snapshot: no cached data, virtualenv or analysis intermediates.
    with __import__('zipfile').ZipFile(review/'source.zip', 'w', compression=__import__('zipfile').ZIP_DEFLATED) as z:
        manifest = read_json(HERE/'source_sha256.json')
        for name in sorted([*manifest, 'source_sha256.json']):
            z.write(HERE/name, name)
    save_json(review/'STATUS.json', dict(status=run_status, modules=stages,
        original_main_batches_recomputed=0, new_support_batches=len(list(tasks(cfg))) if stages['common_support'] else 0,
        not_performed=['New biological cohort', 'Clinical truth-based FDR estimation',
            'Automatic MRN selection or FDR-control guarantee', 'MuSiC/Bisque native comparison',
            'Public-raw-to-LIBD/SEAAD-bank construction', 'Expanded perturbation directions',
            'Manuscript edits or GitHub upload'],
        diagnostic_scope='Original and reliability-qualified reports retained; no independent validation claim',
        intervention_scope='Common-support experiment changes the dose; read dose_retention.csv before interpretation'))
    write_manifest(review)
    export_review(output)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        os.environ[key] = '1'
    args = parser().parse_args(argv)
    for key in ('project_root', 'random_results', 'output'):
        setattr(args, key, getattr(args, key).expanduser().absolute())
    if args.status:
        return status(args)
    if args.export_only:
        from support import lock_is_active
        require(not lock_is_active(args.output), 'Do not export while a writer is active')
        export_review(args.output)
        return 0
    require(os.name != 'nt', 'Run scientific follow-up on the Linux server; Windows tests use test_followup.py')
    require(1 <= args.workers <= 72 and 2 <= args.repetitions <= 1000, 'Invalid workers/repetitions')
    parent = find_parent(args.previous_results, args.project_root)
    assert_separate(args.output, parent, args.random_results, HERE)
    print('PARENT READ ONLY', parent, '\nNEW OUTPUT', args.output, flush=True)
    print('NO PACKAGE INSTALLS; NO DOWNLOADS; NO ORIGINAL MAIN BATCHES RECOMPUTED', flush=True)
    selected = choose_python(args, parent, argv)
    if selected is not None:
        return selected
    code = verify_source()
    owner = dict(package='paper3_bounded_followup', source=code, parent=str(parent))
    if args.output.exists() and any(args.output.iterdir()):
        require((args.output/'OWNER.json').is_file() and read_json(args.output/'OWNER.json') == owner,
                'Output belongs to another run/source. Use a new --output; no old files changed.')
    args.output.mkdir(parents=True, exist_ok=True)
    import fcntl
    with parent_read_lock(parent), (args.output/'run.lock').open('a+b') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        save_json(args.output/'OWNER.json', owner)
        def state(name, **extra):
            save_json(args.output/'STATE.json', dict(status=name, pid=os.getpid(), workers=args.workers,
                      updated_unix=time.time(), **extra))
        def stop(signum, frame):
            raise KeyboardInterrupt('Interrupted; sealed results can resume')
        signal.signal(signal.SIGTERM, stop)
        try:
            state('VERIFYING_PARENT')
            old, identity = verify_parent(parent)
            inherited = ['paired_seed', 'external_seed', 'nuclei', 'random_subsets', 'base_subsets',
                         'subset_ids', 'support_rules', 'methods', 'regression_rule']
            settings = dict({k: old['settings'][k] for k in inherited},
                repetitions=args.repetitions, rep_start=10000, software_fixture=False,
                clinical_derivative_steps=[.01, .001, .0001, .00001],
                clinical_probe_steps=[.01, .05, .1], diagnostic_numeric_floor=1e-6,
                allocation='fixed_reference_common_positive_sample_class_support_equal_reduced_budget',
                outcomes=['null', 'up025'], scope='Common-positive-support equal-reduced-budget follow-up')
            body = dict(package='paper3_bounded_followup', version='1.0.0', source_signature=code,
                parent_identity=identity, units=old['units'], settings=settings, R=old['R'],
                python=old['python'], input_hashes=old['input_hashes'], selection_sha256=old['selection_sha256'])
            cfg = dict(body, signature=digest(body))
            if (args.output/'protocol.json').is_file():
                require(read_json(args.output/'protocol.json') == cfg, 'Frozen follow-up protocol changed; use new output')
            save_json(args.output/'protocol.json', cfg)
            save_json(args.output/'parent_identity.json', identity)
            from next_report import make_corrections
            state('CORRECTING_EXPORTS')
            make_corrections(parent, args.output, cfg['signature'])
            if args.reports_only:
                state('REPORTS_ONLY_COMPLETE')
                make_review(args.output, cfg, 'REPORTS_ONLY_COMPLETE')
                return 0
            from support import checkpoint
            if not args.preflight_only and all(checkpoint(args.output/name, cfg['signature'])
                    for name in ('clinical_multistep', 'common_support')):
                state('FOLLOWUP_COMPLETE', completed_batches=len(list(tasks(cfg))), planned_batches=len(list(tasks(cfg))))
                make_review(args.output, cfg, 'FOLLOWUP_COMPLETE')
                return 0
            from runtime import existing_r, run_r
            import revision_driver as D
            D.run_r = run_r
            recorded = read_json(parent/'runtime.json')
            args.rscript = args.rscript or recorded['rscript']
            args.r_library += [Path(p) for p in recorded.get('r_libraries', '').split(os.pathsep) if p]
            args.offline = True
            args.cache = (args.cache or args.project_root/'01_data'/'paper3_GSE125583').absolute()
            state('CHECKING_EXISTING_R')
            r = existing_r(args)
            require(r['versions'] == old['R'], 'R/package versions differ from parent. No packages installed.')
            save_json(args.output/'runtime.json', dict(python=sys.executable, rscript=r['executable'],
                r_libraries=os.environ.get('R_LIBS', ''), workers=args.workers, offline=True))
            for ref in ('LIBD', 'SEAAD'):
                require((parent/'clinical'/ref/'jobs.json').is_file(),
                        'Full clinical working inputs missing; supply full server run, not compact export')
            require((args.cache/'GSE125583_RAW.tar').is_file() and (args.cache/'GSE125583_family.soft.gz').is_file(),
                    'Public count cache missing. Use --cache EXISTING_VERIFIED_CACHE; downloads are disabled.')
            if not args.preflight_only:
                state('CLINICAL_MULTISTEP')
                from clinical_probes import run as clinical_run
                clinical_run(args, parent, cfg, r)
            state('VERIFYING_PREPARED_BANKS')
            jobs = load_jobs(args, old)
            if args.preflight_only:
                state('PREFLIGHT_COMPLETE')
                print('PREFLIGHT COMPLETE. Same command without --preflight-only starts the new work.', flush=True)
                return 0
            import common_support as S
            from support import checkpoint
            S.JOBS = jobs
            S.OPTIONS = dict(settings, output=str(args.output), signature=cfg['signature'], r_info=r, r_timeout=args.r_timeout)
            todo, timing = [], []
            planned = list(tasks(cfg))
            for task in planned:
                receipt = checkpoint(S.folder_for(args.output, task), cfg['signature'])
                if receipt:
                    timing.append(receipt['elapsed_seconds'])
                else:
                    todo.append(task)
            done = len(planned)-len(todo)
            print('FOLLOW-UP', len(cfg['units']), 'units x', args.repetitions, 'new repetitions;',
                  len(todo), 'remaining batches; workers', args.workers, flush=True)
            last_print = 0.
            def update(result):
                nonlocal done, last_print
                done += 1
                timing.append(result['elapsed_seconds'])
                eta = (len(planned)-done)*statistics.median(timing)/args.workers/3600 if len(timing) >= 3 else None
                state('RUNNING_COMMON_SUPPORT', completed_batches=done, planned_batches=len(planned), rough_eta_hours=eta)
                if time.monotonic()-last_print >= 30 or done == len(planned):
                    print('NEW BATCHES', str(done)+'/'+str(len(planned)), 'simulation ETA hours',
                          round(eta, 2) if eta is not None else 'waiting for timings', flush=True)
                    last_print = time.monotonic()
            if todo:
                state('NEW_BATCH_PILOT', completed_batches=done, planned_batches=len(planned))
                update(S.run_batch(todo.pop(0)))
            if args.workers == 1:
                for task in todo:
                    update(S.run_batch(task))
            elif todo:
                pool = mp.get_context('fork').Pool(args.workers, initializer=S.initialize_worker)
                try:
                    for result in pool.imap_unordered(S.run_batch, todo, chunksize=1):
                        update(result)
                    pool.close()
                except BaseException:
                    pool.terminate()
                    raise
                finally:
                    pool.join()
            state('SUMMARIZING_NEW_WORK')
            S.summarize(args.output, cfg)
            state('FOLLOWUP_COMPLETE', completed_batches=len(planned), planned_batches=len(planned))
            make_review(args.output, cfg, 'FOLLOWUP_COMPLETE')
            return 0
        except BaseException as exc:
            state('INTERRUPTED_OR_FAILED', error=type(exc).__name__+': '+str(exc))
            raise


if __name__ == '__main__':
    raise SystemExit(main())
