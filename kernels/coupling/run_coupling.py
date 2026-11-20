"""Bounded paired TMM controls, coefficient propagation and sensitivity validation."""
import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import signal
import sys
import time
from common import HERE, require, read_json, save_json, sha, digest, file_manifest

SERVER = Path('/home/research/git/repository')
PACKAGE = 'paper3_review_enhancement'


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=SERVER/'paper1_public_evidence')
    p.add_argument('--random-results',type=Path,default=SERVER/'paper3_random_marker_results_v1')
    p.add_argument('--phase',choices=['pilot','confirm'],default='confirm')
    p.add_argument('--output',type=Path)
    p.add_argument('--mode',choices=['all','simulation','clinical'],default='all')
    p.add_argument('--workers',type=int,default=36 if os.name!='nt' else 1)
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--r-library',type=Path,action='append',default=[])
    p.add_argument('--python',type=Path)
    p.add_argument('--cache',type=Path)
    p.add_argument('--r-timeout',type=int,default=14400)
    p.add_argument('--offline',action='store_true',help='Do not download missing public bulk objects')
    p.add_argument('--export-only',action='store_true',help='Verify and export compact results without fitting')
    for name in ('input-bundle','libd-dir','seaad-dir','gene-map','external-cache'):
        p.add_argument('--'+name,type=Path)
    for name in ('status','verify-status','synthetic','smoke','preflight-only','download-bulk','summarize-only'):
        p.add_argument('--'+name,action='store_true')
    return p


def add_engine():
    sys.path[:0] = [str(HERE/'frozen'/'engine'),str(HERE/'frozen')]


def tasks(cfg):
    for start in range(0,cfg['settings']['repetitions'],cfg['settings']['batch_size']):
        for unit in cfg['units']:
            yield unit,start


def batch_kind(start,cfg):
    return 'diagnostic' if start<cfg['settings']['diagnostic_repetitions'] else 'regular'


def remaining_hours(pending,timing,workers):
    import statistics
    seconds=0.
    for kind,number in pending.items():
        if not number:
            continue
        values=timing.get(kind,[])
        if len(values)<3:
            return None
        seconds+=number*statistics.median(values)
    return seconds/max(1,workers)/3600


def progress(args):
    path = args.output/'STATE.json'
    print(json.dumps(read_json(path),indent=2) if path.exists() else 'NOT_STARTED')
    path = args.output/'protocol.json'
    if not path.exists():
        return 0
    sys.path.insert(0,str(HERE/'frozen'))
    from support import checkpoint, lock_is_active
    cfg = read_json(path)
    done,timing,pending = 0,{'diagnostic':[],'regular':[]},{'diagnostic':0,'regular':0}
    planned = list(tasks(cfg))
    for unit,start in planned:
        folder = args.output/'batches'/unit['unit_id']/('r%05d' % start)
        receipt = checkpoint(folder,cfg['signature'],verify=args.verify_status)
        if receipt:
            done += 1
            timing[batch_kind(start,cfg)].append(receipt['elapsed_seconds'])
        else:
            pending[batch_kind(start,cfg)]+=1
    active = lock_is_active(args.output)
    print('WRITER LOCK',active,'SIMULATION BATCHES',str(done)+'/'+str(len(planned)))
    state = read_json(args.output/'STATE.json')
    hours=remaining_hours(pending,timing,state.get('workers',1))
    print('REMAINING BATCH TYPES',pending)
    if active and hours is not None:
        print('Rough simulation-only remaining hours:',round(hours,2),'(excludes clinical and summaries)')
    elif active:
        print('ETA awaits both diagnostic and regular batch timings; early diagnostic batches are heavier.')
    elif not active:
        print('No active writer. If incomplete, resume with the SAME phase/output. Do not delete run.lock.')
    return 0


def source_signature():
    manifest = read_json(HERE/'bundle_sha256.json')
    for name,expected in manifest.items():
        path = (HERE/name).resolve()
        require(path.is_relative_to(HERE) and path.is_file() and sha(path)==expected,'Source mismatch: '+name)
    return digest(manifest)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
        os.environ[key]='1'
    args = parser().parse_args(argv)
    args.output = (args.output or SERVER/('p3_review_'+args.phase)).expanduser().absolute()
    args.project_root=args.project_root.expanduser().absolute()
    args.random_results=args.random_results.expanduser().absolute()
    args.cache=args.cache or args.project_root/'01_data'/'paper3_GSE125583'
    require(not (args.offline and args.download_bulk),'Do not combine --offline and --download-bulk')
    if args.status:
        return progress(args)
    if args.export_only:
        from delivery import export_review
        return export_review(args.output)
    require(args.workers>0 and (os.name!='nt' or args.workers==1),'Positive workers required; Windows fixture uses 1')
    require(not args.output.is_relative_to(args.random_results) and not args.random_results.is_relative_to(args.output),
            'Output must be separate from existing results')
    from runtime import choose_python,existing_r,recorded_inputs,PINNED,run_r
    selected = choose_python(args,argv)
    if selected is not None:
        return selected
    code = source_signature()
    owner = dict(package=PACKAGE,source_signature=code)
    if args.output.exists() and any(args.output.iterdir()):
        require((args.output/'OWNER.json').exists() and read_json(args.output/'OWNER.json')==owner,
                'Output belongs to another source/run. Use a new directory; nothing was overwritten.')
    args.output.mkdir(parents=True,exist_ok=True)
    lock = (args.output/'run.lock').open('a+b')
    if os.name!='nt':
        import fcntl
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def stop(signum,frame):
        raise KeyboardInterrupt('Interrupted; sealed batches can resume')
    if os.name!='nt':
        signal.signal(signal.SIGTERM,stop)
    save_json(args.output/'OWNER.json',owner)
    def state(status,**extra):
        save_json(args.output/'STATE.json',dict(status=status,pid=os.getpid(),workers=args.workers,
                   updated_unix=time.time(),phase=args.phase,**extra))
    add_engine()
    import revision_core as C
    import revision_driver as D
    import driver as frozen
    import enhancement_experiment as E
    from selection import freeze
    from matched_markers import freeze_matched,variants
    from support import checkpoint
    D.run_r = run_r
    try:
        print('PACKAGE INSTALLS DISABLED; OLD FITS NOT REUSED',flush=True)
        if Path('/proc/meminfo').is_file():
            available=next((int(line.split()[1])/1024**2 for line in Path('/proc/meminfo').read_text().splitlines()
                            if line.startswith('MemAvailable:')),None)
            print('REQUESTED WORKERS',args.workers,'MemAvailable GiB',round(available,1) if available is not None else 'unknown',flush=True)
        state('PREFLIGHT')
        frozen_manifest=frozen.source_hashes()
        if args.summarize_only:
            cfg=read_json(args.output/'protocol.json')
            require(cfg['source_signature']==code,'Summary code changed')
        else:
            r=existing_r(args)
            old = None
            if not args.synthetic:
                old=read_json(args.random_results/'protocol.json')
                require(old.get('package')=='paper3_random_marker' and old['source_hashes']==frozen_manifest,
                        'Original protocol does not match the bundled generator')
                require(digest({k:v for k,v in old.items() if k!='signature'})==old['signature'],'Original protocol hash mismatch')
                require(r['versions']==old['R'],'R/package versions differ from paired experiment; choose its --rscript/--r-library')
                require(bool(old['settings'].get('software_smoke',False))==args.smoke,
                        'Parent and current software-acceptance profiles differ')
            recorded_inputs(args)
            args.old_results=args.random_results
            prepared=args.output/'prepared'
            prepared.mkdir(exist_ok=True)
            jobs,identities=frozen.load_jobs(args,prepared)
            hashes={v['name']:v['sha256'] for v in identities}
            require(args.synthetic or hashes==old['input_hashes'],'Prepared input identity differs from prior paired design')
            repetitions,random_subsets,diagnostic_reps=(10,2,2) if args.phase=='pilot' else (200,3,5)
            if args.synthetic:
                repetitions,random_subsets,diagnostic_reps=2,1,1
            elif args.smoke:
                repetitions,random_subsets,diagnostic_reps=2,3,1
            selections=freeze(jobs,random_subsets)
            selections+=freeze_matched(jobs,random_subsets,args.output/'matching')
            for job in jobs:
                job['budget_allocation']='reference'
            units=C.plan_units(jobs,smoke=args.synthetic)
            require(args.synthetic or len(units)==42,'Expected all 42 reference-target units')
            settings=dict(C.SEEDS,repetitions=repetitions,random_subsets=random_subsets,
                diagnostic_repetitions=diagnostic_reps,batch_size=1,
                subset_ids=variants(jobs[0]['stress_subsets']),base_subsets=list(jobs[0]['stress_subsets']),
                diagnostic_subsets=['zero','legacy','random_001','matched_001'],
                derivative_step=.001,probe_steps=[.01,.05,.10],
                methods=['oracle',*E.METHODS],allocation='fixed_reference_equal_sample_class_budget',
                diagnostic_outcomes=['null','up025'],regression_rule='full_six',
                software_fixture=args.synthetic,software_smoke=args.smoke,phase=args.phase,
                matching='reference-only log1p abundance, detection and specificity; within class; assignment jitter 0.1',
                support_rules=['allow_zero','positive_only'],
                infeasible_rule='report separately; no invented intervention; performance conditional on feasible designs',
                cross_engine='official dtangle plus scipy linear NNLS; no claim of all-engine equivalence',
                predictor='p90 absolute coefficient derivative / baseline HC3 SE',
                diagnostic_target='held-out high realized FDP > 0.10; NOT a true-FDR estimator',
                diagnostic_split='five folds by target-source and donor; references share donor fold',
                diagnostic_threshold='training score 90th percentile, floor 1e-6; strict greater-than')
            body=dict(package=PACKAGE,version='1.0.0',source_signature=code,settings=settings,units=units,
                python=PINNED,R=r['versions'],input_hashes=hashes,selection_sha256=digest(selections))
            cfg=dict(body,signature=digest(body))
            if (args.output/'protocol.json').exists():
                require(read_json(args.output/'protocol.json')==cfg,'Protocol changed. Use a new output; only workers/mode may change.')
            save_json(args.output/'protocol.json',cfg)
            save_json(args.output/'selected_markers.json',selections)
            save_json(args.output/'runtime.json',dict(python=sys.executable,rscript=r['executable'],
                r_libraries=os.environ.get('R_LIBS',''),input_hashes=hashes))
            if args.preflight_only:
                state('PREFLIGHT_COMPLETE')
                return 0
            if args.mode in ('all','clinical') and not args.synthetic:
                import bulk_data
                _,metadata,_,_=bulk_data.prepare(args,args.output/'clinical_input_check')
                cohort=bulk_data.cohort(metadata)[0]
                require(len(cohort)==288,'Unexpected clinical cohort before simulation')
            if args.mode in ('all','simulation'):
                E.JOBS=jobs
                E.OPTIONS=dict(settings,output=str(args.output),signature=cfg['signature'],r_info=r,r_timeout=args.r_timeout)
                todo=[]
                timing={'diagnostic':[],'regular':[]}
                pending={'diagnostic':0,'regular':0}
                planned=list(tasks(cfg))
                for task in planned:
                    receipt=checkpoint(E.folder_for(args.output,task),cfg['signature'])
                    kind=batch_kind(task[1],cfg)
                    if not receipt:
                        todo.append(task)
                        pending[kind]+=1
                    else:
                        timing[kind].append(receipt['elapsed_seconds'])
                done=len(planned)-len(todo)
                print('PLAN',len(units),'units;',repetitions,'repetitions;',len(settings['subset_ids']),
                      'marker sets;',len(todo),'remaining batches',flush=True)
                def update(result):
                    nonlocal done
                    done+=1
                    kind=batch_kind(result['start'],cfg)
                    pending[kind]-=1
                    timing[kind].append(result['elapsed_seconds'])
                    eta=remaining_hours(pending,timing,args.workers)
                    state('RUNNING_SIMULATION',completed_batches=done,total_batches=len(planned),rough_eta_hours=eta)
                    print('DONE',str(done)+'/'+str(len(planned)),'batch seconds',round(result['elapsed_seconds'],1),flush=True)
                if todo:
                    state('INFRASTRUCTURE_PILOT')
                    update(E.run_batch(todo.pop(0)))
                if args.workers==1:
                    for task in todo:
                        update(E.run_batch(task))
                elif todo:
                    state('RUNNING_SIMULATION',completed_batches=done,total_batches=len(planned))
                    pool=mp.get_context('fork').Pool(args.workers)
                    try:
                        for result in pool.imap_unordered(E.run_batch,todo,chunksize=1):
                            update(result)
                        pool.close()
                    except BaseException:
                        pool.terminate()
                        raise
                    finally:
                        pool.join()
            if args.mode in ('all','clinical') and not args.synthetic:
                state('RUNNING_CLINICAL')
                from clinical import analyze
                analyze(args,jobs,r,cfg)
        state('SUMMARIZING')
        from enhancement_report import summarize
        summary=summarize(args.output,cfg,args.mode)
        status='SOFTWARE_FIXTURE_ONLY' if cfg['settings']['software_fixture'] else (
            'PILOT_COMPLETE' if cfg['settings']['phase']=='pilot' else 'CONFIRMATION_COMPLETE')
        if not cfg['settings']['software_fixture'] and not (summary['simulation_complete'] and summary['clinical_complete']):
            status='SIMULATION_ONLY_COMPLETE' if summary['simulation_complete'] else 'CLINICAL_ONLY_COMPLETE'
        state(status,**summary)
        final=args.output/'final'
        save_json(final/'STATUS.json',dict(status=status,**summary))
        file_manifest(final)
        from delivery import export_review
        export_review(args.output)
        print('COLLECT ONLY',final,flush=True)
        return 0
    except BaseException as exc:
        state('INTERRUPTED_OR_FAILED',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        lock.close()


if __name__=='__main__':
    raise SystemExit(main())
