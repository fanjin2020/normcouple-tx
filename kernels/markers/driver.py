"""Independent, resumable GSE174367 disease marker anchor; never reruns a grid."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import urllib.request

import numpy as np
import pandas as pd

from anchor_core import (BASE, CLASSES, COUNTS, HASHES, META, bh, design, hc3,
                         metadata, pseudobulk, require, sha, write_json)

HERE = Path(__file__).resolve().parent
KIND = 'paper3_marker_anchor_v1'
SKIP = {'.git','.venv','R_library','node_modules','__pycache__'}
LOCK_HANDLE = None


def lock_output(out):
    global LOCK_HANDLE
    LOCK_HANDLE = (out/'run.lock').open('a+b')
    if os.name == 'nt':
        import msvcrt
        LOCK_HANDLE.seek(0)
        LOCK_HANDLE.write(b'0'); LOCK_HANDLE.flush(); LOCK_HANDLE.seek(0)
        msvcrt.locking(LOCK_HANDLE.fileno(), msvcrt.LK_NBLCK, 1)
    else:
        import fcntl
        fcntl.flock(LOCK_HANDLE, fcntl.LOCK_EX | fcntl.LOCK_NB)


def candidates(name, roots):
    seen = set()
    for root in roots:
        if root is None or not Path(root).is_dir():
            continue
        for base, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d not in SKIP and not d.startswith('.')]
            if name in files:
                p = (Path(base)/name).resolve()
                if p not in seen:
                    seen.add(p)
                    yield p


def acquire(name, args):
    explicit = args.metadata if name == META else args.counts
    choices = [explicit] if explicit else candidates(name, [args.external_cache, args.project_root, args.output/'downloads'])
    for p in choices:
        if p and Path(p).is_file() and sha(p) == HASHES[name]:
            print('REUSE', p, flush=True)
            return Path(p).resolve()
        if explicit:
            raise ValueError('Explicit public input checksum mismatch: ' + str(p))
    if args.offline:
        raise FileNotFoundError('Cached public input not found: ' + name)
    dest = args.output/'downloads'/name
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(name+'.partial')
    for attempt in range(3):
        print('DOWNLOAD', BASE+name, 'attempt', attempt+1, flush=True)
        try:
            req = urllib.request.Request(BASE+name, headers={'User-Agent': 'Paper3-public-marker-audit/1.0'})
            with urllib.request.urlopen(req, timeout=90) as source, tmp.open('wb') as target:
                total = 0
                while True:
                    block = source.read(1024*1024)
                    if not block:
                        break
                    total += len(block)
                    require(total <= 400_000_000, 'Public download exceeds fixed size limit')
                    target.write(block)
            require(sha(tmp) == HASHES[name], 'Public source checksum mismatch')
            os.replace(tmp, dest)
            return dest.resolve()
        except (OSError, ValueError):
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def r_runtime(args):
    exe = str(args.rscript) if args.rscript else shutil.which('Rscript')
    require(bool(exe), 'Rscript not found; activate the existing paper1 environment or pass --rscript')
    libs = [args.r_library] if args.r_library else []
    libs += [p/'R_library' for p in sorted(args.project_root.parent.glob('paper3_*')) if p.is_dir()]
    libs += [Path(os.environ['R_LIBS_USER'])] if os.environ.get('R_LIBS_USER') else []
    libs.append(HERE/'R_library')
    env = os.environ.copy()
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
        env[name] = '1'
    errors = []
    for lib in dict.fromkeys(p.resolve() for p in libs if p is not None):
        result = subprocess.run([exe,'--vanilla',str(HERE/'anchor_fit.R'),'--probe',str(lib)],
                                capture_output=True, text=True, encoding='utf-8', errors='replace', env=env)
        if result.returncode == 0:
            return dict(executable=exe, library=str(lib), versions=json.loads(result.stdout))
        errors.append(result.stderr[-1200:])
    raise RuntimeError('Existing R runtime needs edgeR, limma, statmod and jsonlite. '
                       'Use --r-library to select the prior successful library. No running library was changed.\n'+errors[-1])


def code_digest():
    return {p.relative_to(HERE).as_posix(): sha(p) for p in sorted(HERE.rglob('*'))
            if p.is_file() and p.suffix in ('.py','.R','.csv') and '__pycache__' not in p.parts
            and '.venv' not in p.parts and 'R_library' not in p.parts}


def validate_checkpoint(folder, config):
    receipt = folder/'complete.json'
    try:
        data = json.loads(receipt.read_text(encoding='utf-8'))
        return data['config'] == config and all((folder/k).is_file() and sha(folder/k)==v for k,v in data['files'].items())
    except (OSError, ValueError, KeyError):
        return False


def complete(folder, config, names):
    write_json(folder/'complete.json', dict(config=config, files={name:sha(folder/name) for name in names}))


def prepare_counts(args, cells, path):
    dest = args.output/'prepared'
    dest.mkdir(exist_ok=True)
    config = dict(input=HASHES[COUNTS], metadata=HASHES[META], aggregation='all_author_annotated_donor_by_class_v1')
    if validate_checkpoint(dest, config):
        print('REUSE verified pseudobulk checkpoint', flush=True)
        with np.load(dest/'pseudobulk.npz', allow_pickle=False) as f:
            return pd.Index(f['genes']), f['symbols'], f['counts'], pd.Index(f['groups'])
    print('AGGREGATE raw snRNA integer counts once', flush=True)
    genes, symbols, matrix, groups, audit = pseudobulk(path, cells)
    np.savez_compressed(dest/'pseudobulk.npz', genes=genes.to_numpy(dtype=str), symbols=symbols.astype(str),
                        counts=matrix, groups=groups.to_numpy(dtype=str))
    write_json(dest/'aggregation_audit.json', audit)
    complete(dest, config, ['pseudobulk.npz','aggregation_audit.json'])
    return genes, symbols, matrix, groups


def library_endpoints(mat, genes, samples, x, panels, cell_type, model):
    rows, point_rows, loo_rows = [], [], []
    total = mat.sum(axis=0)
    require((total > 0).all(), 'Empty pseudobulk sample')
    for reference, panel in panels.groupby('reference', sort=True):
        for scope in ('owner','all_panel'):
            subset = panel if scope=='all_panel' else panel[panel.winner_class.eq(cell_type)]
            indices = genes.get_indexer(subset.gene_id)
            ids = indices[indices>=0]
            if len(ids) == 0:
                continue
            mass = mat[ids].sum(axis=0)
            nuclei = samples.nuclei.to_numpy()
            prefix = dict(reference=reference, cell_type=cell_type, model=model, scope=scope,
                          planned_markers=len(subset), present_markers=len(ids))
            point_rows.extend(dict(**prefix, SampleID=sid, Diagnosis=diag, nuclei=int(n),
                marker_umi=int(mc), total_umi=int(tc), marker_fraction=float(mc/tc),
                marker_umi_per_nucleus=float(mc/n), total_umi_per_nucleus=float(tc/n))
                for sid,diag,n,mc,tc in zip(samples.SampleID,samples.Diagnosis,nuclei,mass,total))
            endpoints = dict(marker_fraction_log2=np.log2((mass+.5)/(total+1)),
                             marker_umi_per_nucleus_log2=np.log2((mass+.5)/nuclei),
                             total_umi_per_nucleus_log2=np.log2((total+.5)/nuclei))
            for name, y in endpoints.items():
                f = hc3(y, x)
                rows.append(dict(**prefix, endpoint=name, **f))
                for omit, sid in enumerate(samples.SampleID):
                    keep = np.arange(len(samples)) != omit
                    f = hc3(y[keep], x[keep])
                    loo_rows.append(dict(**prefix, endpoint=name, omitted_donor=sid, **f))
    return rows, point_rows, loo_rows


def run_job(job, runtime):
    folder = Path(job['folder'])
    config = job['config']
    if validate_checkpoint(folder, config):
        return job['key'], 'REUSED'
    started = time.time()
    with (folder/'R.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([runtime['executable'],'--vanilla',str(HERE/'anchor_fit.R'),str(folder),runtime['library']],
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('R job failed: ' + str(folder/'R.log'))
    complete(folder, config, ['gene_results.csv','gene_filter.csv','normalization_factors.csv','sessionInfo.txt'])
    write_json(folder/'timing.json', dict(seconds=time.time()-started))
    return job['key'], 'COMPLETE'


def compile_results(args, jobs, statuses, panels, genes, library, donor_points, loo):
    from figures import make_figures
    final = args.output/'final_results'
    final.mkdir(exist_ok=True)
    frames = []
    for status in statuses:
        ct, model = status['cell_type'], status['model']
        for normalizer in ('TMM_all','TMM_nonmarker'):
            key = ct+'__'+model+'__'+normalizer
            folder = args.output/'fits'/key
            valid = status['status']=='READY'
            result = pd.read_csv(folder/'gene_results.csv') if valid else None
            for reference, panel in panels.groupby('reference', sort=True):
                for scope in ('owner','all_panel'):
                    planned = panel if scope=='all_panel' else panel[panel.winner_class.eq(ct)]
                    t = planned.copy()
                    t['present_in_h5'] = t.gene_id.isin(genes)
                    if valid:
                        t = t.merge(result, on='gene_id', how='left', validate='one_to_one')
                        t['available'] = t.qlf_p.notna()
                        t['unavailable_reason'] = np.where(t.available,'',np.where(t.present_in_h5,'expression_filter','absent_in_h5'))
                    else:
                        t['available'] = False
                        t['unavailable_reason'] = status['reason']
                    t['cell_type'], t['model'], t['normalizer'], t['scope'] = ct, model, normalizer, scope
                    frames.append(t)
    all_results = pd.concat(frames, ignore_index=True)
    # Keep the planned panel family, including missing/untestable markers as p=1.
    for model in ('primary','extended'):
        for norm in ('TMM_all','TMM_nonmarker'):
            for scope in ('owner','all_panel'):
                take = all_results.model.eq(model)&all_results.normalizer.eq(norm)&all_results.scope.eq(scope)
                for col in ['qlf_p',*[f'treat_p_{d}' for d in (25,50,100)]]:
                    if col not in all_results:
                        all_results[col] = np.nan
                    all_results.loc[take, col+'_q_planned_family'] = bh(all_results.loc[take,col].fillna(1))
    all_results.to_csv(final/'marker_gene_associations.csv.gz', index=False, compression='gzip')
    summaries = []
    keys = ['reference','cell_type','model','normalizer','scope']
    for group, f in all_results.groupby(keys, sort=True):
        b = pd.to_numeric(f.get('log2FC',pd.Series(np.nan,index=f.index)),errors='coerce')
        for dose in (25,50,100):
            threshold = float(np.log2(1+dose/100))
            sig = f[f'treat_p_{dose}_q_planned_family'].le(.05)&f.available
            summaries.append(dict(zip(keys,group), dose_percent=dose, log2FC_threshold=threshold,
                planned=len(f), present=int(f.present_in_h5.sum()), tested=int(f.available.sum()),
                point_up=int((b>=threshold).sum()), point_down=int((b<=-threshold).sum()),
                supported_up=int((sig&b.gt(0)).sum()), supported_down=int((sig&b.lt(0)).sum()),
                point_up_fraction_tested=float((b>=threshold).sum()/f.available.sum()) if f.available.sum() else np.nan,
                point_down_fraction_tested=float((b<=-threshold).sum()/f.available.sum()) if f.available.sum() else np.nan,
                supported_up_fraction_planned=float((sig&b.gt(0)).sum()/len(f)) if f.available.any() else np.nan,
                supported_down_fraction_planned=float((sig&b.lt(0)).sum()/len(f)) if f.available.any() else np.nan))
    summary = pd.DataFrame(summaries)
    summary.to_csv(final/'marker_dose_anchor_summary.csv',index=False)
    lib = pd.DataFrame(library)
    if not lib.empty:
        lib['q_within_model_endpoint'] = np.nan
        for _, frame in lib.groupby(['model','endpoint']):
            lib.loc[frame.index,'q_within_model_endpoint'] = bh(frame.p.fillna(1))
    lib.to_csv(final/'marker_library_associations.csv',index=False)
    points = pd.DataFrame(donor_points)
    points.to_csv(final/'donor_marker_libraries.csv',index=False)
    pd.DataFrame(loo).to_csv(final/'library_leave_one_donor_out.csv',index=False)
    pd.DataFrame(statuses).to_csv(final/'design_diagnostics.csv',index=False)
    factors = []
    for job in jobs:
        t = pd.read_csv(Path(job['folder'])/'normalization_factors.csv')
        for key in ('cell_type','model','normalizer'):
            t[key] = job['config'][key]
        factors.append(t)
    pd.concat(factors,ignore_index=True).to_csv(final/'normalization_factors.csv',index=False)
    shutil.copy2(Path(jobs[0]['folder'])/'sessionInfo.txt',final/'sessionInfo.txt')
    shutil.copy2(HERE/'data'/'panel_provenance.json',final/'panel_provenance.json')
    for name in ('donors.csv','donor_cell_counts.csv','cohort_exclusions.csv','public_input_receipts.json',
                 'protocol.json','R_environment.json','METHODS.md','metadata_stage.json'):
        if (args.output/name).is_file():
            shutil.copy2(args.output/name, final/name)
    shutil.copy2(args.output/'prepared'/'aggregation_audit.json',final/'aggregation_audit.json')
    panels.to_csv(final/'fixed_marker_panels.csv',index=False)
    make_figures(summary, all_results, lib, final)
    required = len(statuses)*2
    state = dict(status='COMPLETE_WITH_DECLARED_GAPS' if any(s['status']!='READY' for s in statuses) else 'COMPLETE',
        completed_fits=len(jobs), planned_fits=required, unavailable_fits=required-len(jobs),
        independent_unit='donor within deposited broad cell type', raw_bulk_performed=False,
        simulation_rerun=False, MRN_superiority_tested=False, true_fdr_evaluated=False,
        new_external_resource=False, warning='An empirical association anchor, not confirmation of a causal marker-stress intervention.')
    write_json(final/'STATUS.json',state)
    write_json(final/'file_manifest.json',{p.name:sha(p) for p in sorted(final.iterdir()) if p.is_file() and p.name!='file_manifest.json'})
    write_json(args.output/'STATUS.json',state)
    return state


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('/home/research/git/repository/paper1_public_evidence'))
    p.add_argument('--output',type=Path)
    p.add_argument('--external-cache',type=Path)
    p.add_argument('--metadata',type=Path)
    p.add_argument('--counts',type=Path)
    p.add_argument('--rscript',type=Path)
    p.add_argument('--r-library',type=Path)
    p.add_argument('--workers',type=int,default=36)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--preflight',action='store_true',help='Metadata and R checks only; no large count download')
    p.add_argument('--status',action='store_true')
    args = p.parse_args(argv)
    args.project_root = args.project_root.resolve()
    args.output = (args.output or args.project_root.parent/'paper3_marker_anchor_results_v1').resolve()
    if args.status:
        f = args.output/'STATUS.json'
        print(f.read_text(encoding='utf-8') if f.exists() else 'NOT_STARTED')
        return 0
    require(args.workers>=1,'workers must be positive')
    args.output.mkdir(parents=True,exist_ok=True)
    require(not (args.output/'protocol.json').exists() or
            json.loads((args.output/'protocol.json').read_text(encoding='utf-8')).get('kind')==KIND,
            'Output belongs to a different analysis; use the new marker-anchor directory')
    require(not (args.output/'owner.json').exists() or
            json.loads((args.output/'owner.json').read_text())['kind']==KIND,'Different output owner')
    lock_output(args.output)
    write_json(args.output/'owner.json',dict(kind=KIND))
    try:
        write_json(args.output/'STATUS.json',dict(status='PREFLIGHT',pid=os.getpid()))
        panels = pd.read_csv(HERE/'data'/'fixed_marker_panels.csv')
        require(not panels[['reference','gene_id']].duplicated().any(),'Duplicated fixed markers')
        require(panels.groupby('reference').size().to_dict()=={'LIBD':300,'SEAAD':300},'Fixed panel identities changed')
        meta = acquire(META,args)
        cells, donors, nuclei = metadata(meta)
        donors.to_csv(args.output/'donors.csv',index=False)
        nuclei.to_csv(args.output/'donor_cell_counts.csv')
        runtime = r_runtime(args)
        write_json(args.output/'R_environment.json',runtime)
        protocol = dict(kind=KIND, source=code_digest(), public_sha256=HASHES, min_nuclei=20,
            minimum_group_donors=4, minimum_residual_df=4,
            models=dict(primary=['Age','Sex','Batch'],extended=['Age','Sex','Batch','PMI','RIN']),
            normalizers=['TMM_all','TMM_nonmarker'], dose_percent=[25,50,100],
            R_versions=runtime['versions'], multiple_testing='Both fixed reference panels and all cell types within each model/normalizer/scope; unavailable p=1',
            estimand='AD versus Control relative gene expression within author-annotated broad cell types')
        existing = args.output/'protocol.json'
        if existing.exists():
            require(json.loads(existing.read_text(encoding='utf-8'))==protocol,
                    'Protocol or dependency versions changed. Use a new output directory; do not mix results.')
        write_json(existing,protocol)
        shutil.copy2(HERE/'METHODS.md',args.output/'METHODS.md')
        stage = dict(donors=len(donors),groups=donors.Diagnosis.value_counts().to_dict(),nuclei=len(cells),
            analysis_unit='SampleID donor; no nuclei-as-replicates tests',missing_PMI=int(donors.PMI.isna().sum()))
        write_json(args.output/'metadata_stage.json',stage)
        print('METADATA',stage,flush=True)
        if args.preflight:
            write_json(args.output/'STATUS.json',dict(status='PREFLIGHT_COMPLETE_NOT_ANALYZED',**stage))
            return 0
        counts = acquire(COUNTS,args)
        write_json(args.output/'public_input_receipts.json',[dict(path=str(v),sha256=HASHES[k],url=BASE+k) for k,v in ((META,meta),(COUNTS,counts))])
        genes,symbols,matrix,groups = prepare_counts(args,cells,counts)
        jobs,statuses,exclusions,library,points,loo = [],[],[],[],[],[]
        for ct in CLASSES:
            for model in ('primary','extended'):
                m,x,audit,info = design(donors,nuclei[ct],20,model)
                info['cell_type'] = ct
                statuses.append(info)
                audit['cell_type'],audit['model'] = ct,model
                exclusions.append(audit)
                if info['status']!='READY':
                    print('UNAVAILABLE',ct,model,info['reason'],flush=True)
                    continue
                pos = groups.get_indexer(m.SampleID+'|'+ct)
                require((pos>=0).all(),'Missing pseudobulk column')
                mat = matrix[:,pos]
                lib,pts,lvs = library_endpoints(mat,genes,m,x.to_numpy(),panels,ct,model)
                library.extend(lib);points.extend(pts);loo.extend(lvs)
                for normalizer in ('TMM_all','TMM_nonmarker'):
                    key = ct+'__'+model+'__'+normalizer
                    folder = args.output/'fits'/key
                    folder.mkdir(parents=True,exist_ok=True)
                    config = dict(cell_type=ct,model=model,normalizer=normalizer,n=len(m),genes=len(genes),sample_ids=m.SampleID.tolist())
                    np.asarray(mat,dtype='<f8').ravel(order='F').tofile(folder/'counts.bin')
                    pd.DataFrame(dict(gene_id=genes,symbol=symbols,is_any_marker=genes.isin(panels.gene_id).astype(int))).to_csv(folder/'genes.csv',index=False)
                    x.to_csv(folder/'design.csv',index=False)
                    write_json(folder/'job.json',config)
                    jobs.append(dict(key=key,folder=str(folder),config=config))
        pd.concat(exclusions).to_csv(args.output/'cohort_exclusions.csv',index=False)
        done = 0
        workers = min(args.workers,len(jobs))
        require(workers>0,'No estimable cell-type designs; inspect donor metadata')
        print('RUN',len(jobs),'independent R fits; workers',workers,flush=True)
        def progress():
            write_json(args.output/'STATUS.json',dict(status='RUNNING',completed=done,total=len(jobs),workers=workers,
                unavailable_designs=sum(s['status']!='READY' for s in statuses),pid=os.getpid()))
        progress()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(run_job,job,runtime) for job in jobs]
            for future in as_completed(futures):
                key,state = future.result()
                done += 1
                print(str(done)+'/'+str(len(jobs)),key,state,flush=True)
                progress()
        result = compile_results(args,jobs,statuses,panels,genes,library,points,loo)
        print(json.dumps(result), '\nFINAL_RESULTS',args.output/'final_results',flush=True)
        return 0
    except Exception as exc:
        write_json(args.output/'STATUS.json',dict(status='FAILED',error=type(exc).__name__+': '+str(exc),
            action='Fix the reported input/environment issue and resume the same command; verified completed fits are reused.'))
        raise


if __name__=='__main__':
    raise SystemExit(main())
