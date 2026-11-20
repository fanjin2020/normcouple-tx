"""Measured doses and bidirectional, expression-blind-to-diagnosis gene probes."""
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import t
from .io import ROOT, atomic_json, digest, read_json, sha, verify_receipt
from .sensitivity import aligned

COHORTS = ('GSE125583','GSE53697')
REFS = ('LIBD','SEAAD')
METHODS = ('dtangle_tmm','dtangle_marker')


def design(raw, genes, tests, markers):
    """Rows in the returned protocol use the 1-based indices of the R plan."""
    totals = raw.sum(axis=0)
    names = np.asarray(genes, dtype=str)
    tests = np.asarray(tests, dtype=int)
    if (raw.ndim!=2 or raw.shape[0]!=len(names) or raw.shape[1]<1 or len(tests)<10
            or not len(markers) or min(markers)<1 or max(markers)>len(names) or len(set(markers))!=len(markers)
            or len(np.unique(tests))!=len(tests) or len(set(names))!=len(names)
            or set(tests)&set(markers) or min(tests)<1 or max(tests)>len(names)
            or not np.isfinite(raw).all() or raw.min()<0 or np.any(totals<=0)):
        raise ValueError('Invalid test family, gene identities or libraries')
    mean_cpm = (raw[tests-1]/totals).mean(axis=1)*1e6
    order = tests[np.lexsort((names[tests-1],mean_cpm))]
    sets = {}
    for label,q in (('low',.10),('middle',.50),('high',.90)):
        pos = int(np.floor(q*(len(order)-1)))
        start = max(0,min(len(order)-5,pos-2))
        sets['single_'+label] = [int(order[pos])]
        sets['five_'+label] = order[start:start+5].tolist()
    jobs=[dict(id='baseline',kind='baseline',step=0.,rows=[])]
    for kind,rows in [*sets.items(),('test_family',tests.tolist())]:
        for suffix,step in [('down',-.1),('up',.1)]:
            jobs.append(dict(id=kind+'_'+suffix,kind=kind,step=step,rows=rows))
    jobs.append(dict(id='marker_up',kind='marker_control',step=.1,rows=list(markers)))
    return jobs


def stats(engine, bh, y, x, p, cov):
    fit=engine.fit_hc3(y,x,p,cov,rule='adaptive')
    good=fit['basic'] & fit['good']
    pv=np.ones(y.shape[1])
    pv[good]=2*t.sf(np.abs(fit['beta'][good]/fit['se'][good]),fit['df'])
    return fit,bh(pv),good


def check_cache(clinical,cross):
    from .publication import PUB
    expected=read_json(PUB/'published_manifest.json')['estimator_provenance']['analysis_sources']
    for key,h in expected.items():
        family,rel=key.split('/',1)
        path=(clinical if family=='clinical_estimators' else cross)/rel
        if not path.is_file() or sha(path)!=h:
            raise ValueError('Not the full submitted estimator run: '+str(path))
    protocol=read_json(clinical/'protocol.json')
    sign=digest(protocol)
    done=read_json(clinical/'complete.json')
    if done.get('protocol')!=sign: raise ValueError('Clinical protocol mismatch')
    verify_receipt(clinical,done['files'])
    cross_plan=read_json(cross/'protocol.json')
    closed=read_json(cross/'complete.json')
    if (closed.get('protocol')!=digest(cross_plan) or closed.get('datasets')!=16800
            or cross_plan['repetitions']!=200 or len(cross_plan['units'])!=42 or cross_plan['software_smoke']
            or cross_plan['seeds']['probe_seed']!=2026111701):
        raise ValueError('Cross-method run is not the completed 42-unit, 200-repetition design')
    return sign,cross_plan


def subset_audit(bank, tmp, cohort, ref, sign, rscript, libs, engine, bh):
    from .publication import PUB,load_frame
    rec=read_json(bank/'receipt.json')
    if rec.get('protocol')!=sign: raise ValueError('Unmatched clinical receipt: '+str(bank))
    verify_receipt(bank,rec['files'])
    plan=read_json(bank/'plan.json')
    ids=plan['sample_ids']; genes=np.asarray(plan['genes'])
    x=np.asarray(plan['exposure'],float); tests=np.asarray(plan['test_rows'],dtype=int)
    meta=aligned(pd.read_csv(bank.parent/'cohort.csv',dtype={'sample_id':str}),ids,'cohort')
    if not np.array_equal(meta.group.eq('AD').to_numpy(float),x): raise ValueError('Exposure differs from frozen design')
    cov=(np.column_stack([(meta.age_numeric-meta.age_numeric.mean())/meta.age_numeric.std(ddof=1),meta.sex.eq('M')])
         if cohort=='GSE125583' else meta.sequencing_batch.eq(2).to_numpy(float)[:,None])
    raw_path=bank/'raw.bin'
    if raw_path.stat().st_size!=8*len(genes)*len(ids): raise ValueError('Incorrect clinical matrix dimensions')
    raw=np.memmap(raw_path,mode='r',dtype='<f8',shape=(len(genes),len(ids)),order='F')
    raw_hash=sha(raw_path)
    if 'engine_reuse.json' in rec['files'] and raw_hash!=read_json(bank/'engine_reuse.json')['raw_sha256']:
        raise ValueError('Raw matrix differs from original estimator reuse receipt')
    y=np.memmap(bank/'outcomes.bin',mode='r',dtype='<f8',shape=(len(ids),len(tests)),order='F')
    den=pd.read_csv(bank/'outcome_denominator.csv',header=None).iloc[:,0].to_numpy(float)
    if len(den)!=len(ids) or not np.isfinite(den).all() or (den<=0).any():
        raise ValueError('Invalid sealed baseline denominator vector')
    for ii in np.array_split(np.arange(len(tests)),16):
        if not np.allclose(y[:,ii].T,np.log2(raw[tests[ii]-1]/den*1e6+1),rtol=1e-12,atol=1e-12):
            raise ValueError('Saved outcome is not the unperturbed TMM response')
    jobs=design(raw,genes,tests,plan['marker_rows'])
    if not set(plan['stress_rows']).issubset(plan['marker_rows']) or len(set(plan['stress_rows']))!=len(plan['stress_rows']):
        raise ValueError('Invalid marker positive-control subset')
    jobs[-1]['rows']=plan['stress_rows']
    raw_totals=raw.sum(axis=0)
    panel=pd.DataFrame([dict(selection=j['kind'],gene_id=genes[r-1],baseline_mean_cpm=float((raw[r-1]/raw_totals).mean()*1e6))
                        for j in jobs
                        if j['step']>0 and j['kind'].startswith(('single_','five_')) for r in j['rows']])
    atlas=pd.DataFrame(dict(row=plan['marker_rows'],cell_class=[plan['classes'][k-1] for k in plan['marker_owner']]))
    atlas['gene_id']=genes[atlas.row.to_numpy()-1]
    atlas['marker_probe']=atlas.row.isin(plan['stress_rows'])
    for _,group in atlas.groupby('cell_class'):
        selected=group.sort_values('gene_id').head(int(np.ceil(.1*len(group))))
        if set(group[group.marker_probe].row)!=set(selected.row): raise ValueError('Marker quota or sorting differs')
    atomic_json(tmp/'probes.json',jobs)
    subprocess.run([rscript,str(ROOT/'pipeline/probe_engine.R'),'clinical',str(bank),str(tmp),*map(str,libs)],check=True)
    d=pd.read_csv(tmp/'dose.csv',dtype={'sample_id':str})
    if d.duplicated(['probe','sample_id']).any() or set(d.probe)!=set(j['id'] for j in jobs):
        raise ValueError('Incomplete measured dose series')
    summaries=[]
    for j in jobs:
        dd=aligned(d[d.probe.eq(j['id'])],ids,'measured dose')
        if not np.array_equal(dd.group,x): raise ValueError('Measured dose has a different sample grouping')
        for group,subset in dd.groupby('group'):
            row=dict(probe=j['id'],group='AD' if group else 'control',samples=len(subset),genes=len(j['rows']))
            for col in ('tested_input_share','library_relative_change','effective_relative_change'):
                if not np.isfinite(subset[col]).all(): raise ValueError('Non-finite measured dose')
                for label,q in [('q25',.25),('median',.5),('q75',.75)]:
                    row[col+'_'+label]=subset[col].quantile(q)
            summaries.append(row)
    res=[]; detail=[]; classes=plan['classes']
    eligible=bool(plan['marker_methods_available'])
    old_table=load_frame(PUB/'data','TableS18.csv')
    old_table=old_table[(old_table.cohort==cohort)&(old_table.reference==ref)].set_index(['method','probe'])
    def confirm_fit(method,old_probe,fit,q,available):
        old=old_table.loc[(method,old_probe)]
        if int(available.sum())!=old.available_tests or len(q)!=old.family_size:
            raise ValueError('Reference test family/availability differs: '+method+'/'+old_probe)
        if available.any() and (int((q<=.05).sum())!=old.discoveries or
                not np.isclose(fit['exposure_residual_ss'],old.exposure_information,atol=1e-7,rtol=1e-7)):
            raise ValueError('Original coefficient design or discoveries not reproduced: '+method+'/'+old_probe)
    if eligible:
        pred=pd.read_csv(tmp/'probe_fractions.csv',dtype={'sample_id':str})
        frozen=pd.read_csv(bank/'fractions.csv',dtype={'sample_id':str})
    for method in METHODS:
        def props(frame,probe):
            mat=aligned(frame[frame.method.eq(method)&frame.probe.eq(probe)],ids,probe)[classes].to_numpy(float)
            if not np.isfinite(mat).all() or (mat<0).any() or not np.allclose(mat.sum(1),1,atol=1e-8,rtol=0):
                raise ValueError('Invalid fraction matrix')
            return mat/mat.sum(1,keepdims=True)
        if eligible:
            for new,old in [('baseline','baseline'),('test_family_up','test_counts'),('marker_up','marker_counts')]:
                if not np.allclose(props(pred,new),props(frozen,old),rtol=0,atol=1e-9):
                    raise ValueError('New fits disagree with sealed clinical estimator: '+method+'/'+old)
            p0=props(pred,'baseline'); f0,q0,g0=stats(engine,bh,y,x,p0,cov)
            confirm_fit(method,'baseline',f0,q0,g0)
        for j in jobs[1:]:
            row=dict(method=method,probe=j['id'],log_step=j['step'],tested_family=len(tests))
            res.append(row)
            if not eligible:
                row.update(available=False,reason='insufficient_marker_coverage')
                continue
            p1=props(pred,j['id']); f1,q1,g1=stats(engine,bh,y,x,p1,cov)
            if j['id'] in ('marker_up','test_family_up'):
                confirm_fit(method,'marker_counts' if j['id']=='marker_up' else 'test_counts',f1,q1,g1)
            good=g0 & g1
            if not good.any():
                row.update(available=False,reason=f0['reason'] or f1['reason'] or 'no_paired_tests')
                continue
            delta=f1['beta']-f0['beta']; direct=(f1['weights']-f0['weights'])@y
            err=float(np.max(abs(delta-direct)))
            if not np.isfinite(err) or err>1e-8*(1+np.max(abs(delta))): raise ValueError('Weight/refit mismatch')
            status=(q1<=.05)!=(q0<=.05)
            target=np.isin(tests,j['rows'])
            row.update(available=True,reason='',paired_tests=int(good.sum()),weight_check_max_error=err,
                       max_fraction_delta=float(np.abs(p1-p0).max()),bh_crossings=int(status.sum()))
            for scope,mask in [('own',target),('other',~target)]:
                ix=good & mask
                row[scope+'_paired_genes']=int(ix.sum())
                row[scope+'_beta_rms']=np.sqrt(np.mean(delta[ix]**2)) if ix.any() else np.nan
            if j['kind'].startswith(('single_','five_')):
                for k in np.flatnonzero(target):
                    detail.append(dict(method=method,probe=j['id'],gene_id=genes[tests[k]-1],
                        available=bool(good[k]),beta_change=delta[k] if good[k] else np.nan,
                        bh_crossed=bool(status[k]) if good[k] else None))
            if method=='dtangle_marker' and j['kind']!='marker_control' and (row['max_fraction_delta']>1e-10 or np.abs(delta).max()>1e-8):
                raise ValueError('Disjoint fixed-marker negative control responded')
    frames=dict(dose=pd.DataFrame(summaries),probe_effects=pd.DataFrame(res),self_effects=pd.DataFrame(detail),
                markers=atlas.drop(columns='row'),selected_genes=panel)
    for f in frames.values(): f['cohort']=cohort; f['reference']=ref
    manifest=dict(raw=raw_hash,panel_sha256=sha(bank/'A.bin'),reference_receipt=sha(bank/'receipt.json'),
        outcome_sha256=rec['files']['outcomes.bin'],cohort_sha256=sha(bank.parent/'cohort.csv'),
        marker_counts={cls:dict(mapped=int(sum(atlas.cell_class==cls)),probed=int(sum((atlas.cell_class==cls)&atlas.marker_probe))) for cls in classes})
    return frames,manifest
