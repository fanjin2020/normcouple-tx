"""One resumable post-benchmark audit. Sealed historical runs are read only."""
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import time
from types import SimpleNamespace
import zipfile
import numpy as np
import pandas as pd
from . import __version__
from .io import ROOT, atomic_json, digest, lock, read_json, receipt, separate, sha, verify_receipt
from . import imaging_audit as im
from . import probe_audit as probes
from .publication import PUB, export_published, identity_check, load_frame, cross_method_table, read


def demo(work):
    from .sensitivity import aligned,execute
    work=Path(work)
    if work.exists() and any(work.iterdir()):
        raise ValueError('Example output must be empty: '+str(work))
    work.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(50113); n=48
    ids=[f'demo{i:02}' for i in range(n)]
    x=np.tile([0.,1.],n//2); p=rng.uniform(.15,.7,n); u=rng.normal(size=n)
    y=2*p[:,None]+.3*x[:,None]+rng.normal(0,.35,(n,5))
    pd.DataFrame(y,columns=[f'g{i}' for i in range(5)]).assign(sample_id=ids).to_csv(work/'outcomes.csv',index=False)
    pd.DataFrame(dict(sample_id=ids,exposure=x,age=u)).to_csv(work/'design.csv',index=False)
    fractions=pd.concat([pd.DataFrame(dict(sample_id=ids,method=m,A=q,B=1-q)) for m,q in
        [('baseline',p),('equivalent',.1+.75*p),('exposure_aligned',p+.12*x)]])
    fractions.sample(frac=1,random_state=3).to_csv(work/'fractions.csv',index=False)
    try: aligned(pd.DataFrame({'sample_id':ids[:-1]}),ids,'demo negative control')
    except ValueError: pass
    else: raise ValueError('Missing-sample negative control failed')
    start=time.monotonic()
    execute(SimpleNamespace(outcomes=work/'outcomes.csv',design=work/'design.csv',fractions=work/'fractions.csv',
        exposure='exposure',covariates=['age'],baseline='baseline',output=work/'result'))
    rows=pd.read_csv(work/'result/design_comparison.csv').set_index('method')
    if not np.isfinite(rows.loc[['equivalent','exposure_aligned'],'beta_max']).all() or rows.loc['equivalent','beta_max']>1e-10 or rows.loc['equivalent','bh_status_changes']!=0 or rows.loc['exposure_aligned','beta_max']<1e-3:
        raise ValueError('Diagnostic example failed')
    return dict(mode='synthetic_software_check_only',family=5,samples=n,seconds=round(time.monotonic()-start,3),
        input_alignment='missing rejected; permuted rows realigned',equivalent_beta_max=float(rows.loc['equivalent','beta_max']),
        aligned_beta_max=float(rows.loc['exposure_aligned','beta_max']))


def checkpoint(root,name,key,writer):
    folder=root/name; signature=digest(key)
    if folder.exists():
        record=read_json(folder/'receipt.json')
        if record['signature']!=signature: raise ValueError('Checkpoint inputs changed: '+name+'; use a new output')
        verify_receipt(folder,record['files'])
        print('REUSE',name,flush=True); return folder
    root.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.'+name+'-',dir=root) as temp:
        scratch=Path(temp); keep=scratch/'keep'; keep.mkdir()
        writer(scratch,keep)
        files=receipt(keep,[p.relative_to(keep).as_posix() for p in keep.rglob('*') if p.is_file()])
        if not files: raise ValueError('Empty checkpoint: '+name)
        atomic_json(keep/'receipt.json',dict(signature=signature,files=files))
        os.replace(keep,folder)
    print('COMPLETE',name,flush=True); return folder


def protocol(cross,clinical,ids):
    env=read_json(clinical/'protocol.json')['environment']
    return dict(version=__version__,coefficient_identity=ids,simulation=dict(repetitions=cross['repetitions'],
        test_probe=cross['test_probe'],marker_probe=cross['marker_probe'],test_gene_count=20,
        test_law='Independent Poisson with mean n_nuclei * frozen gene scale * exposed * 0.5 at the 20 fixed indices; zero elsewhere. Does not multiply observed test counts; zero observed count alone does not force a zero addition.',
        seed='rng_for(2026111701, fixed_response_test_probe, base.key); same draws across methods and both outcome effects',
        test_seed_source='pipeline/cross_method.py',seeded_gene_definitions='simulation_probe_genes.csv (frozen 20 per reference/target); full tested families in publication/input_definitions/complete_test_gene_families.csv',
        recorded_seeds=cross['seeds'],recorded_environment=cross['environment'],recorded_svr=cross['linear_svr']),
        clinical=dict(outcome='Frozen baseline log2(CPM+1); sample IDs, adjustment set and full BH family held fixed',
        selection='Pooled baseline mean CPM, ascending with stable gene-ID tie break; fixed 10/50/90th ranked singleton and five-neighbor sets, selected without diagnosis or outcome significance',
        perturbation='Multiply only selected covariate-construction counts in AD by exp(+/-0.1); zeros remain zero. Positive-only marker control from sealed plan: first ceil(10% of mapped IDs) per class in gene ID order.',
        svr='L1-normalize reference columns and each mixture; reference globally centered/scaled by matrix-wide sample SD, each mixture independently centered/scaled across genes. Zero SD rejected. Linear nu-SVR cost 1, nu 0.25/0.5/0.75, scale=FALSE. Clip and close coefficients, select minimum mean((Xw-y)^2), then convert using marker yields.',
        historical_environment=env),
        orthogonal=dict(citation='Huuki-Myers et al., Genome Biology 26:88, doi:10.1186/s13059-025-03552-3',upstream_commit=im.COMMIT,
        design='Separate assay on the same matched tissue blocks, not an independent clinical disease cohort. Total-tissue libraries, held-out LIBD targets, six jointly defined classes, Oligo+OPC summed. The released filtered count matrix defines the assay universe. Original clinical panel: retained intersection and original anchor; no estimate below 10 mapped genes per class.',
        transfer='A separate reference-only panel, not validation of the original marker set. Four reference donors. Normalize each donor/class pseudobulk, average profiles and per-nucleus yields across donors. Among shared gene IDs, require 10 CPM in the owning class in at least three donors and positive specificity over the strongest competing class. Retain up to 25 genes per class by descending log2[(CPM_top+1)/(CPM_next+1)], gene ID tie break. Reference-derived anchor. No bulk expression values, imaging outcomes or diagnoses in selection; no fit below 10 markers per class.',
        aggregation='High/OK measured proportions on all segmented nuclei, Other excluded without renormalizing. Average sections and technical libraries within block; average block losses within donor, then weight donors equally. Class-specific common support. No clinical FDR truth.'))


def execute(args):
    repo=args.root.expanduser().resolve()
    clinical=(args.clinical or repo/'p3_complete_20261117/clinical_final').resolve()
    cross=(args.cross or repo/'p3_complete_20261117/cross_full').resolve()
    cache=(args.imaging_cache or repo/'normcouple_imaging_sources').resolve()
    work=separate(args.output or repo/'p3_evidence_20261119',clinical,cross,cache)
    separate(cache,clinical,cross)
    for p in (clinical,cross):
        if not (p/'protocol.json').is_file(): raise ValueError('Missing completed server run: '+str(p)+' (override --clinical / --cross)')
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']: os.environ[name]='1'
    libs=[Path(p).resolve() for p in args.r_library]
    if not libs:
        libs=[p for p in [repo/'paper3_revision_results_v1/R_library',repo/'paper1_public_evidence/03_inputs/313_paper3_R_library_v2'] if p.is_dir()]
    candidate=Path.home()/'miniconda3/envs/paper1/bin/Rscript'
    rscript=args.rscript or (str(candidate) if candidate.is_file() else 'Rscript')
    work.mkdir(parents=True,exist_ok=True)
    print('OUTPUT',work,'\nREUSE CLINICAL',clinical,'\nREUSE CROSS',cross,flush=True)
    with lock(work/'.writer.lock'):
        from .package import verify
        from .cross_method import _setup_runtime
        from .run import python_environment
        verify()
        python_environment()
        sign,cross_plan=probes.check_cache(clinical,cross)
        # Link the summary to sealed raw benchmark outputs, not just a copied table.
        freshly=cross_method_table(read(cross/'final/summary.csv'))
        pd.testing.assert_frame_equal(freshly,load_frame(PUB/'data','TableS16.csv'),check_dtype=False,rtol=1e-13,atol=1e-15)
        for n,path in [(17,cross/'final/paired_precision.csv'),(18,clinical/'final/method_summary.csv'),(19,clinical/'final/paired_sensitivity.csv')]:
            pd.testing.assert_frame_equal(read(path),load_frame(PUB/'data',f'TableS{n}.csv'),check_dtype=False,rtol=1e-13,atol=1e-15)
        imaging_hashes=im.sources(cache,args.allow_downloads)
        renv=subprocess.check_output([rscript,'--vanilla','-e',
            'a<-commandArgs(TRUE);.libPaths(c(a,.libPaths()));stopifnot(as.character(packageVersion("edgeR"))=="4.8.2",as.character(packageVersion("dtangle"))=="2.0.9");cat(R.version.string,as.character(packageVersion("jsonlite")))',
            *map(str,libs)],text=True).strip()
        subprocess.run([rscript,'--vanilla',str(ROOT/'tests/test_probe_engine.R'),*map(str,libs)],check=True)
        import scipy
        software=dict(manifest=sha(ROOT/'MANIFEST_SHA256.csv'),R=renv,python=platform.python_version(),
                      numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
        history=work/'checkpoints'
        c,_,engine=_setup_runtime(work)
        keys={}
        for cohort in probes.COHORTS:
            raw_hashes=[]; out_hashes=[]
            for ref in probes.REFS:
                bank=clinical/cohort/ref
                hashes=receipt(bank,['raw.bin','A.bin','receipt.json'])
                key=dict(software=software,clinical_protocol=sign,inputs=hashes,cohort=sha(bank.parent/'cohort.csv'))
                rec=read_json(bank/'receipt.json')
                if rec.get('protocol')!=sign: raise ValueError('Clinical receipt signature differs')
                verify_receipt(bank,rec['files'])
                out_hashes.append(rec['files']['outcomes.bin']); raw_hashes.append(hashes['raw.bin'])
                def run(scratch,keep):
                    frames,manifest=probes.subset_audit(bank,scratch,cohort,ref,sign,rscript,libs,engine,c.M.bh)
                    for name,frame in frames.items(): frame.to_csv(keep/(name+'.csv'),index=False)
                    atomic_json(keep/'inputs.json',manifest)
                label=cohort+'_'+ref
                keys[label]=checkpoint(history,label,key,run)
            if len(set(raw_hashes))!=1 or len(set(out_hashes))!=1: raise ValueError('Samples or response differ across references')
        def imaging(scratch,keep):
            scores,pairs,match,markers=im.compute(cache,clinical,cross,scratch,rscript,libs)
            scores.to_csv(keep/'imaging_error.csv',index=False); pairs.to_csv(keep/'imaging_blocks.csv',index=False)
            markers.to_csv(keep/'imaging_markers.csv',index=False)
            atomic_json(keep/'imaging_match.json',match)
        refs=im.reference_receipts(cross,cross_plan)
        images=checkpoint(history,'imaging',dict(software=software,sources=imaging_hashes,reference_banks=refs,
            banks=[sha(d/'inputs.json') for k,d in keys.items() if k.startswith('GSE125583')]),imaging)
        def assemble(scratch,keep):
            for name in ('dose','markers','selected_genes','probe_effects','self_effects'):
                frames=[pd.read_csv(d/(name+'.csv')) for d in keys.values()]
                nonempty=[f for f in frames if len(f)]
                if not nonempty: raise ValueError('No eligible data for '+name)
                pd.concat(nonempty,ignore_index=True).to_csv(keep/(name+'.csv'),index=False)
            for name in ('imaging_error.csv','imaging_blocks.csv','imaging_match.json','imaging_markers.csv'):
                shutil.copyfile(images/name,keep/name)
            sets=read(PUB/'input_definitions/test_gene_scales.csv')
            sets=sets[sets.panel.eq('frozen')]
            if len(sets)!=80 or not sets.groupby(['reference','target']).gene_id.nunique().eq(20).all():
                raise ValueError('Frozen 20-gene simulation probe differs')
            sets.to_csv(keep/'simulation_probe_genes.csv',index=False)
            atom=identity_check(load_frame(PUB/'data','clinical_linear_checks.csv'))
            atomic_json(keep/'protocol.json',dict(protocol(cross_plan,clinical,atom),software=software,
                imaging_source_sha256=imaging_hashes,imaging_reference_sha256=refs,clinical_inputs={k:read_json(p/'inputs.json') for k,p in keys.items()},
                demo=demo(scratch/'example')))
        final=checkpoint(work,'results',dict(software=software,stages=[sha(p/'receipt.json') for p in [*keys.values(),images]]),assemble)
        paper=export_published(work/'manuscript_replay',['pdf'],evidence=final)
        fmap=read_json(paper/'figures/figure_sources.json')
        captions=read_json(PUB/'figure_captions.json')
        if set(fmap)!=set(captions): raise ValueError('Incomplete publication figure index')
        for name,info in fmap.items():
            if info['caption']!=captions[name] or not info['sources'] or any(not (paper/'data'/s).is_file() for s in info['sources']):
                raise ValueError('Unresolved figure data: '+name)
        archive=work/'paper3_review.zip'
        with tempfile.NamedTemporaryFile(dir=work,suffix='.zip',delete=False) as f: tmp=Path(f.name)
        try:
            with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
                for path in sorted(final.iterdir()): z.write(path,'evidence/'+path.name)
                for n in (13,16,17,18,19,20,21): z.write(paper/'data'/f'TableS{n}.csv',f'manuscript/TableS{n}.csv')
                for path in [paper/'tables.json',paper/'figures/figure_sources.json']: z.write(path,'manuscript/'+path.name)
            with zipfile.ZipFile(tmp) as z:
                if z.testzip(): raise ValueError('Review archive CRC failure')
            os.replace(tmp,archive)
        finally: tmp.unlink(missing_ok=True)
        print('IMAGING PANELS:',[(p['reference'],p['panel'],p['status']) for p in read_json(final/'imaging_match.json')['panels']],flush=True)
        print('COMPLETE C1-C6; REVIEW ZIP:',archive,'\nSHA256:',sha(archive),flush=True)
