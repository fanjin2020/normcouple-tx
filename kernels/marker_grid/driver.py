"""Additional grid and mechanism checks, not a replacement for reviewed results."""
import argparse
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import time
from concurrent.futures import ProcessPoolExecutor,as_completed

import revision_core as C
import revision_driver as D
import extension_driver as E
import grid
from summarize import summarize

HERE=Path(__file__).resolve().parent


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('/home/research/git/repository/paper1_public_evidence'))
    p.add_argument('--output',type=Path,default=Path('/home/research/git/repository/paper3_strengthening_results_v1'))
    for name in ('libd-dir','seaad-dir','gene-map','external-cache','input-bundle','r-library','old-results',
                 'bulk-manifest','raw-counts','sample-map'):
        p.add_argument('--'+name,type=Path)
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--workers',type=int,default=36 if os.name!='nt' else 1)
    p.add_argument('--repetitions',type=int,default=200)
    p.add_argument('--batch-size',type=int,default=5)
    p.add_argument('--r-timeout',type=int,default=7200)
    p.add_argument('--only',choices=['all','grid','bulk','summarize'],default='all')
    for flag in ('offline','synthetic','smoke','pilot-only','status'):
        p.add_argument('--'+flag,action='store_true')
    for name in ('counts-url','counts-sha256','counts-provenance'):
        p.add_argument('--'+name,default=None)
    p.add_argument('--max-download-gb',type=int,default=2)
    return p


def progress(output,workers):
    path=output/'protocol.json'
    if not path.exists():
        print('No protocol yet: dependency or data preparation may still be running.'); return 0
    cfg=json.loads(path.read_text()); settings=cfg['settings']
    total=len(cfg['units'])*((settings['repetitions']+settings['batch_size']-1)//settings['batch_size'])
    times=[]; done=0; damaged=[]
    for unit in cfg['units']:
        for start in range(0,settings['repetitions'],settings['batch_size']):
            folder=output/'units'/unit['unit_id']/f'chunk_{start:05d}'
            try:
                if D.valid_checkpoint(folder,cfg['signature']):
                    done+=1; times.append(json.loads((folder/'batch.json').read_text())['elapsed_seconds'])
            except (OSError,ValueError,KeyError) as exc: damaged.append(str(exc))
    state=json.loads((output/'RUN_STATUS.json').read_text()) if (output/'RUN_STATUS.json').exists() else {}
    print('Status:',state.get('status','UNKNOWN'))
    print(f'Checksum-verified batches: {done}/{total} ({100*done/total:.1f}%)')
    print('Damaged checkpoints:',len(damaged))
    if times and done<total:
        import statistics
        median=statistics.median(times)
        active=max(1,min(workers,total-done))
        print(f'Finished batch median: {median/60:.1f} min; simulation-only estimate: {(total-done)*median/active/3600:.1f} h at {active} workers')
        print('Estimate excludes data loading, final summaries and bulk access; not a completion promise.')
    for problem in damaged[:5]: print(problem)
    return 2 if damaged else 0


def r_environment(args,out):
    E.choose_r(args)
    exe=shutil.which(str(args.rscript))
    C.require(exe is not None,'Rscript unavailable; supply --rscript')
    library=(args.r_library or out/'R_library').resolve()
    library.mkdir(parents=True,exist_ok=True)
    version=out/'R_versions.json'
    D.run_r([exe,'--vanilla',str(HERE/'grid_bridge.R'),'--check',str(library),str(version)],
            library,out/'R_check.log',timeout=180)
    return dict(executable=exe,library=str(library),versions=json.loads(version.read_text()))


def source_hashes():
    manifest=json.loads((HERE/'source_manifest.json').read_text())
    for name,digest in manifest.items(): C.require(C.sha(HERE/name)==digest,'Source was modified: '+name)
    return manifest


def main(argv=None):
    args=parser().parse_args(argv)
    if args.status: return progress(args.output,args.workers)
    C.require(args.workers>0 and args.batch_size>0,'Worker and batch counts must be positive')
    C.require(os.name!='nt' or args.workers==1,'Use 1 worker on Windows; Linux supports shared-memory fork')
    C.require(args.synthetic or args.smoke or args.repetitions>=200,'Research protocol requires at least 200 repetitions')
    out=args.output.resolve()
    if out.exists() and any(out.iterdir()):
        C.require((out/'OWNER.json').exists() and json.loads((out/'OWNER.json').read_text()).get('package')=='paper3_strengthening',
                  'Use a new output directory; old reviewed results are never overwritten')
    out.mkdir(parents=True,exist_ok=True); C.E.hold_output_lock(out)
    C.write_json(out/'OWNER.json',dict(package='paper3_strengthening',version='1.0.0'))
    C.write_json(out/'RUN_STATUS.json',dict(status='PREFLIGHT',workers=args.workers))
    source=source_hashes()
    if args.only=='summarize':
        C.require((out/'protocol.json').is_file(),'Cannot summarize without an existing protocol')
        cfg=json.loads((out/'protocol.json').read_text())
        C.require(cfg['source_hashes']==source,'Source changed since this run')
        options=dict(cfg['settings'],signature=cfg['signature'])
        final=summarize(out,cfg['units'],options)
        print('COLLECT',final,flush=True)
        C.write_json(out/'RUN_STATUS.json',dict(status='SUMMARIES_REBUILT',signature=cfg['signature'],
                     bulk_status_file=str(out/'BULK_STATUS.json')))
        shutil.copy2(out/'RUN_STATUS.json',final/'RUN_STATUS.json')
        return 0
    r_info=r_environment(args,out)
    args.stage2_results=None; args.install_r_packages=False
    if args.external_cache is None: args.external_cache=args.project_root/'01_data'
    identities=[]
    print('Preparing verified cached count banks',flush=True)
    if args.synthetic:
        jobs=D.synthetic_jobs(out)
    else:
        if args.input_bundle is None:
            portable=args.project_root.parent/'paper3_mechanism_results_v1_prepared_inputs'
            if (portable/'inputs_manifest.json').exists(): args.input_bundle=portable
        if args.only=='bulk' and args.input_bundle is None:
            folders,gene_map,identities=C.E.resolve_inputs(args)
            _,jobs=C.full_inputs(folders,gene_map)
        else:
            folders,gene_map,meta,counts,identities=E.inputs(args,out)
            sources,jobs=C.full_inputs(folders,gene_map)
            if args.only!='bulk':
                external=C.S.read_external(counts,meta,next(iter(sources.values()))['genes'],out/'external_inputs')
                jobs+=C.S.make_external_jobs(jobs,external,out)
        jobs=[C.prepare_job(j,out,2000) for j in jobs]
        C.write_json(out/'input_paths.json',identities)
    bulk=dict(status='NOT_REQUESTED')
    if args.only in ('all','bulk') and not args.synthetic:
        print('BULK: metadata-first cohort audit and authorized count validation',flush=True)
        from bulk_repair import run as bulk_run,collect
        bulk=bulk_run(args,out,jobs,r_info)
        print('BULK STATUS',bulk['status'],flush=True)
        C.write_json(out/'BULK_STATUS.json',bulk)
        if args.only=='bulk':
            final=out/'final_results'; final.mkdir(exist_ok=True)
            collect(out,final)
            C.write_json(final/'BULK_STATUS.json',bulk)
            C.write_json(out/'RUN_STATUS.json',dict(status='BULK_STAGE_FINISHED',bulk=bulk))
            return 0 if bulk['status']=='COMPLETE_RAW_COUNT_APPLICATION' else 2
    units=C.plan_units(jobs,smoke=args.synthetic)
    C.require(args.synthetic or len(units)==42,'Frozen design requires 42 donor/reference units')
    settings=dict(C.SEEDS,repetitions=2 if args.synthetic or args.smoke else args.repetitions,
                  batch_size=2 if args.synthetic or args.smoke else args.batch_size,software_fixture=args.synthetic,software_smoke=args.smoke,
                  grid=grid.stress_grid(),selection_seed=2026101401,stress_seed=2026101303,
                  regression_rule='full_six',pseudocount=1.,methods=grid.METHODS,
                  power_audit_denominator='common complete repetitions across TMM, MRN and denominator-only')
    body=dict(version='1.0.0',settings=settings,units=units,source_hashes=source,
              input_hashes={x['name']:x['sha256'] for x in identities},R=r_info['versions'],
              python={p:importlib.metadata.version(p) for p in ('numpy','pandas','scipy','h5py','matplotlib')})
    signature=C.M.digest(body)
    if (out/'protocol.json').exists():
        C.require(json.loads((out/'protocol.json').read_text())['signature']==signature,
                  'Analysis settings/code/input identity changed; use a new output directory')
    C.write_json(out/'protocol.json',dict(body,signature=signature))
    C.write_json(out/'selected_markers.json',grid.selection_manifest(jobs))
    grid.JOBS=jobs; grid.OPTIONS=dict(settings,output=str(out),r_info=r_info,r_timeout=args.r_timeout,signature=signature)
    tasks=[(u,start) for u in units for start in range(0,settings['repetitions'],settings['batch_size'])]
    C.write_json(out/'RUN_STATUS.json',dict(status='RUNNING_GRID',workers=args.workers,total_batches=len(tasks),signature=signature,bulk=bulk))
    try:
        if args.only!='summarize':
            # A numerical pilot, never selection on favorable FDR/power.
            pilot=[next(u for u in units if u['reference_source']==s) for s in sorted({u['reference_source'] for u in units})]
            for unit in pilot:
                grid.run_batch(unit,0)
                import pandas as pd
                checks=pd.read_csv(out/'units'/unit['unit_id']/'chunk_00000'/'invariance.csv')
                zero=checks[checks.grid_id.eq('zero')]
                C.require(len(zero)>0 and zero.finite.all() and zero.max_change.lt(1e-8).all(),
                          'Baseline MRN numerical pilot failed; grid expansion stopped')
            if args.pilot_only:
                C.write_json(out/'RUN_STATUS.json',dict(status='PILOT_COMPLETE',signature=signature,bulk=bulk))
                return 0
            if args.workers==1:
                for unit,start in tasks: grid.run_batch(unit,start)
            else:
                with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('fork')) as pool:
                    futures=[pool.submit(grid.run_batch,u,start) for u,start in tasks]
                    for future in as_completed(futures): future.result()
        print('Summarizing paired grid and per-gene power decomposition',flush=True)
        final=summarize(out,units,grid.OPTIONS)
        if (out/'real_bulk').exists():
            from bulk_repair import collect
            collect(out,final)
        status='SOFTWARE_FIXTURE_ONLY' if args.synthetic else ('COMPLETE_REQUESTED_STAGES' if bulk['status'] in ('NOT_REQUESTED','COMPLETE_RAW_COUNT_APPLICATION') else 'GRID_COMPLETE_BULK_INCOMPLETE')
        receipt=dict(status=status,grid='COMPLETE',bulk=bulk,signature=signature,
                     public_code_deposition='NOT_PERFORMED',new_research_results=not args.synthetic)
        C.write_json(out/'RUN_STATUS.json',receipt)
        for name in ('RUN_STATUS.json','protocol.json','selected_markers.json'):
            shutil.copy2(out/name,final/name)
        for name in ('README.md','THEORY.md'):
            shutil.copy2(HERE/name,final/name)
        print('COLLECT',final,flush=True)
        print(json.dumps(receipt,indent=2),flush=True)
        return 2 if status=='GRID_COMPLETE_BULK_INCOMPLETE' else 0
    except BaseException as exc:
        C.write_json(out/'RUN_STATUS.json',dict(status='INTERRUPTED_OR_FAILED',error=type(exc).__name__+': '+str(exc),signature=signature,bulk=bulk))
        raise


if __name__=='__main__':
    raise SystemExit(main())
