"""Raw-count clinical sensitivity with explicit gene-selection policies."""
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

from .cohorts import external_cohort, expression_family, gene_mapping
from .cross_method import MUSIC_COMMIT, _reference, _setup_runtime, paired_metrics
from .inputs import stage_bulk_cache, verify_inputs
from .io import ROOT, atomic_json, digest, lock, read_json, receipt, separate, sha, verify_receipt

METHODS = ('dtangle_tmm','dtangle_marker','music_marker','music_broad','music_non_test','linear_svr_marker')
PROBES = ('baseline','test_counts','marker_counts')


def summarize_fits(s, c, y, x, cov, predictions, classes, gene_ids, sample_ids, unavailable=None):
    unavailable = unavailable or {}
    expected = {(m,p) for m in METHODS for p in PROBES}
    excluded = {(m,p) for m in unavailable for p in PROBES}
    if not set(unavailable).issubset(METHODS):
        raise ValueError('Unknown unavailable clinical method')
    fits, pvalues, props, summary, detail = {}, {}, {}, [], []
    for (method, probe), frame in predictions.groupby(['method','probe'],sort=True):
        if (method,probe) not in expected-excluded:
            raise ValueError('Unexpected clinical prediction method/probe')
        if not frame.sample_id.is_unique or set(frame.sample_id)!=set(sample_ids):
            raise ValueError('Clinical prediction sample identities differ')
        p = frame.set_index('sample_id').loc[sample_ids,classes].to_numpy(float)
        if not np.isfinite(p).all() or (p<0).any() or not np.allclose(p.sum(1),1,atol=1e-6,rtol=0):
            raise ValueError('Invalid clinical fractions')
        p /= p.sum(1,keepdims=True)
        fit = s.fit_hc3(y,x,p,cov,rule='adaptive')
        available = fit['basic'] & fit['good']
        pv = np.ones(len(gene_ids))
        pv[available] = 2*stats.t.sf(abs(fit['beta'][available]/fit['se'][available]),fit['df'])
        q = c.M.bh(pv)
        key = method, probe
        fits[key], pvalues[key], props[key] = fit, q, p
        summary.append(dict(method=method,probe=probe,samples=len(x),family_size=len(gene_ids),
            available_tests=int(available.sum()),discoveries=int((q<=.05).sum()) if available.any() else np.nan,
            composition_rank=fit['composition_rank'],residual_df=fit['df'],
            maximum_leverage=fit['leverage'],exposure_information=fit['exposure_residual_ss'],
            reason=fit['reason']))
        detail.append(pd.DataFrame(dict(method=method,probe=probe,gene_id=gene_ids,
            available=available,beta=fit['beta'],se=fit['se'],p=pv,q=q)))
    if set(fits)!=expected-excluded:
        raise ValueError('Missing clinical method/probe combination')
    for method,probe in sorted(excluded):
        summary.append(dict(method=method,probe=probe,samples=len(x),family_size=len(gene_ids),
            available_tests=0,discoveries=np.nan,composition_rank=np.nan,residual_df=np.nan,
            maximum_leverage=np.nan,exposure_information=np.nan,reason=unavailable[method]))
    comparisons = []
    for key in sorted(expected):
        method, probe = key
        baseline = ('dtangle_tmm','baseline') if probe=='baseline' else (method,'baseline')
        if key==baseline:
            continue
        if key not in fits or baseline not in fits:
            comparisons.append(dict(baseline=baseline[0],method=method,probe=probe,
                available=False,reason=unavailable.get(method,unavailable.get(baseline[0])),
                family_size=len(gene_ids),beta_rms=np.nan,bh_crossings=np.nan,
                gained=np.nan,lost=np.nan,crossing_fraction=np.nan))
            continue
        row = paired_metrics(fits[baseline],fits[key],y,props[baseline],props[key])
        infer = bool(fits[baseline]['basic'] and fits[key]['basic'] and
                     fits[baseline]['good'].any() and fits[key]['good'].any())
        a,b = pvalues[baseline]<=.05,pvalues[key]<=.05
        row.update(bh_crossings=int((a^b).sum()) if infer else np.nan,
                   gained=int((b&~a).sum()) if infer else np.nan,lost=int((a&~b).sum()) if infer else np.nan,
                   family_size=len(gene_ids),crossing_fraction=float((a^b).mean()) if infer else np.nan)
        reason='' if infer else ';'.join(sorted({f['reason'] for f in (fits[baseline],fits[key]) if f['reason']}))
        comparisons.append(dict(baseline=baseline[0],method=method,probe=probe,available=infer,reason=reason,**row))
    return pd.DataFrame(summary),pd.DataFrame(comparisons),pd.concat(detail,ignore_index=True)


def reuse_estimates(source, dest, accession, reference, protocol):
    """Reuse only sealed estimator outputs for an identical raw-input/R plan."""
    source=Path(source)
    prior=read_json(source/'protocol.json')
    old_signature=digest(prior)
    completion=read_json(source/'complete.json')
    if completion.get('protocol')!=old_signature:
        raise ValueError('Clinical cache protocol checksum differs')
    verify_receipt(source,completion['files'])
    if digest({k:v for k,v in prior.items() if k!='code'})!=digest({k:v for k,v in protocol.items() if k!='code'}):
        raise ValueError('Clinical cache scientific protocol differs')
    for name,expected in protocol['code'].items():
        if name!='pipeline/clinical_estimators.py' and prior['code'].get(name)!=expected:
            raise ValueError('Clinical cache estimator/input code differs: '+name)
    old=source/accession/reference
    record=read_json(old/'receipt.json')
    if record.get('protocol')!=old_signature:
        raise ValueError('Clinical cache reference protocol differs')
    verify_receipt(old,record['files'])
    before,after=read_json(old/'plan.json'),read_json(dest/'plan.json')
    for plan in (before,after):
        plan.pop('reference')
    if before!=after or any(sha(old/name)!=sha(dest/name) for name in ('raw.bin','A.bin')):
        raise ValueError('Clinical cache raw counts or estimator plan differs')
    for name in ('fractions.csv','engine.json','outcomes.bin','outcome_denominator.csv'):
        shutil.copy2(old/name,dest/name)
    atomic_json(dest/'engine_reuse.json',dict(source_protocol=old_signature,
        source_receipt=sha(old/'receipt.json'),raw_sha256=sha(dest/'raw.bin'),plan=digest(after)))
    print('REUSE VERIFIED CLINICAL ESTIMATES',accession,reference,flush=True)


def execute(args):
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    work = separate(args.work,args.inputs,args.bulk_cache,args.external_bulk_cache)
    if getattr(args,'estimates_from',None):
        separate(work,args.estimates_from)
    work.mkdir(parents=True,exist_ok=True)
    with lock(work/'.writer.lock'):
        verified = verify_inputs(args.inputs,args.bulk_cache)
        from .run import python_environment
        env = dict(python=python_environment(),R=json.loads(subprocess.check_output(
            [args.rscript,str(ROOT/'pipeline/cross_method.R'),'--environment',
             *map(str,args.r_library)],text=True)))
        external_raw,external_meta,external_provenance = external_cohort(
            args.external_bulk_cache,args.allow_downloads)
        c,_,s = _setup_runtime(work)
        runtime = work/'.runtime/common_support/baseline'
        sys.path.insert(0,str(runtime))
        bulk = importlib.import_module('bulk_data')
        cache = work/'.input-cache/bulk'
        stage_bulk_cache(cache,args.bulk_cache,verified)
        from types import SimpleNamespace
        raw,metadata,all_mapping,source = bulk.prepare(SimpleNamespace(cache=cache,offline=True),work/'primary_inputs')
        cohort,x,cov,_ = bulk.cohort(metadata)
        if len(cohort)!=288 or int(x.sum())!=218:
            raise ValueError('Primary clinical cohort differs from the manuscript')
        raw = raw.loc[:,cohort.sample_id]
        source['covariates'] = ['age','sex']
        sources,jobs = c.full_inputs({k:args.inputs/k for k in ('LIBD','SEAAD')},args.inputs/'gene_identity_map.csv')
        mapping = gene_mapping()
        identity = {str(p.relative_to(ROOT)):sha(p) for p in
            [ROOT/'pipeline/clinical_estimators.py',ROOT/'pipeline/clinical_estimators.R',
             ROOT/'pipeline/estimators.R',ROOT/'pipeline/cohorts.py',ROOT/'kernels/layout.json']}
        protocol = dict(schema=2,inputs=verified,external=external_provenance,environment=env,code=identity,
                        methods=METHODS,probes=PROBES,probe_log_fold=.1,
                        family='CPM>=1 in at least half of samples; union of both reference marker panels excluded',
                        outcome='Baseline all-gene TMM log2 CPM +1; unchanged across models/probes',
                        music='Official API with fixed reference moments; marker, native broad, and non-test gene policies',
                        svr='LIBSVM linear nu-SVR with nonnegative coefficient closure and marker-yield conversion',
                        marker_coverage='At least 10 mapped markers per class for marker-based methods; unavailable specifications are reported, native broad/non-test MuSiC still run')
        signature = digest(protocol)
        if (work/'protocol.json').exists() and (not args.resume or read_json(work/'protocol.json')!=json.loads(json.dumps(protocol))):
            raise ValueError('Use --resume with identical clinical inputs and code')
        atomic_json(work/'protocol.json',protocol)
        all_markers = set()
        refs = []
        for job in jobs:
            job['base'] = c.M.reference(job['source'],job['refs'])
            genes = job['source']['genes'][job['panels']['balanced']]
            all_markers.update(all_mapping.loc[all_mapping.ensembl.isin(genes),'entrez'])
            owners = c.E.official_markers(job).cell_type.to_numpy(int)
            inverse = mapping.reset_index().set_index('ensembl').entrez
            ids = [inverse.get(g,'') for g in genes]
            all_markers.update(g for g in ids if g)
            ref = _reference(job,work,verified)
            subprocess.run([args.rscript,str(ROOT/'pipeline/cross_method.R'),'--prepare',str(ref),
                            *map(str,args.r_library)],check=True)
            refs.append((job,genes,owners,ids,ref))
        cohorts = [('GSE125583',raw,cohort,x,cov,source),
                   ('GSE53697',external_raw,external_meta,external_meta.group.eq('AD').to_numpy(float),
                    external_meta.sequencing_batch.eq(2).to_numpy(float)[:,None],external_provenance)]
        summaries,comparisons,engines = [],[],[]
        for accession,counts,meta,exposure,covariates,provenance in cohorts:
            tests = expression_family(counts,all_markers)
            if accession=='GSE125583' and len(tests)!=16161:
                raise ValueError('Primary test family differs from the manuscript')
            mapped = [mapping.ensembl.get(g,'UNMAPPED:'+g) for g in counts.index]
            fixed_y = None
            for job,marker_genes,owners,ids,ref in refs:
                label = job['source']['label']
                dest = work/accession/label
                dest.mkdir(parents=True,exist_ok=True)
                marker = dest/'receipt.json'
                if marker.exists() and args.resume:
                    record = read_json(marker)
                    if record['protocol']!=signature:
                        raise ValueError('Clinical reference protocol changed')
                    verify_receipt(dest,record['files'])
                else:
                    keep = [i for i,g in enumerate(ids) if g in counts.index]
                    panel = job['panels']['balanced'][keep]
                    mark = counts.index.get_indexer([ids[i] for i in keep])
                    owner = owners[keep]
                    coverage = [int(np.sum(owner==k)) for k in range(1,8)]
                    marker_ok = min(coverage)>=10
                    unavailable = {} if marker_ok else {method:'INSUFFICIENT_MARKER_COVERAGE: '+','.join(map(str,coverage))
                        for method in ('dtangle_tmm','dtangle_marker','music_marker','linear_svr_marker')}
                    selected = []
                    for k in range(1,8):
                        candidates = np.flatnonzero(owner==k)
                        order = np.argsort(np.array(marker_genes)[keep][candidates],kind='stable')
                        selected.extend(mark[candidates[order[:max(1,int(np.ceil(.1*len(candidates))))]]])
                    a = job['base']['A'][panel]
                    pooled = job['base']['counts'][panel].sum(0)
                    full_mass = pooled/a.sum(0)
                    anchor = float(pooled.sum()/full_mass.sum())
                    full_yield = job['base']['pooled_yield']
                    c.write_bin(dest/'raw.bin',counts.to_numpy(float))
                    c.write_bin(dest/'A.bin',a)
                    atomic_json(dest/'plan.json',dict(classes=c.M.CT,genes=mapped,
                        sample_ids=meta.sample_id.tolist(),exposure=exposure.tolist(),
                        marker_rows=(mark+1).tolist(),marker_owner=owner.tolist(),
                        marker_counts_by_class=coverage,marker_methods_available=marker_ok,
                        test_rows=(tests+1).tolist(),stress_rows=(np.array(selected)+1).tolist(),
                        reference=str(ref),music_commit=MUSIC_COMMIT,anchor=anchor,
                        full_yield=full_yield.tolist(),marker_yield=(full_yield*a.sum(0)).tolist()))
                    if getattr(args,'estimates_from',None):
                        reuse_estimates(args.estimates_from,dest,accession,label,protocol)
                    else:
                        with (dest/'engine.log').open('w') as log:
                            subprocess.run([args.rscript,str(ROOT/'pipeline/clinical_estimators.R'),str(dest),
                                *map(str,args.r_library)],check=True,stdout=log,stderr=subprocess.STDOUT)
                    y = np.fromfile(dest/'outcomes.bin',dtype='<f8').reshape((len(meta),len(tests)),order='F')
                    preds = pd.read_csv(dest/'fractions.csv',dtype={'sample_id':str})
                    sm,cm,detail = summarize_fits(s,c,y,exposure,covariates,preds,c.M.CT,
                                                counts.index[tests].tolist(),meta.sample_id.tolist(),unavailable)
                    sm['minimum_markers_per_class'] = min(coverage)
                    sm['mapped_marker_total'] = sum(coverage)
                    sm.to_csv(dest/'summary.csv',index=False)
                    cm.to_csv(dest/'comparisons.csv',index=False)
                    detail.to_csv(dest/'gene_statistics.csv.gz',index=False,compression='gzip')
                    atomic_json(marker,dict(protocol=signature,files=receipt(dest,[
                        'summary.csv','comparisons.csv','gene_statistics.csv.gz','fractions.csv',
                        'engine.json','outcomes.bin','outcome_denominator.csv','plan.json',
                        *(['engine_reuse.json'] if (dest/'engine_reuse.json').exists() else [])])))
                y = np.fromfile(dest/'outcomes.bin',dtype='<f8')
                if fixed_y is None:
                    fixed_y = y
                elif not np.array_equal(y,fixed_y):
                    raise ValueError('Clinical outcomes differ across reference banks')
                summaries.append(pd.read_csv(dest/'summary.csv').assign(cohort=accession,reference=label))
                comparisons.append(pd.read_csv(dest/'comparisons.csv').assign(cohort=accession,reference=label))
                engines.append(dict(cohort=accession,reference=label,
                    reused_estimates=read_json(dest/'engine_reuse.json') if (dest/'engine_reuse.json').exists() else None,
                    **read_json(dest/'engine.json')))
                print('CLINICAL_ESTIMATORS',accession,label,len(meta),'samples',len(tests),'genes',flush=True)
            meta.to_csv(work/accession/'cohort.csv',index=False)
            atomic_json(work/accession/'source.json',provenance)
        final = work/'final'
        final.mkdir(exist_ok=True)
        pd.concat(summaries,ignore_index=True).to_csv(final/'method_summary.csv',index=False)
        pd.concat(comparisons,ignore_index=True).to_csv(final/'paired_sensitivity.csv',index=False)
        atomic_json(final/'validation.json',dict(official_parity=engines,
            same_outcomes_across_references=True,clinical_truth_available=False))
        atomic_json(work/'complete.json',dict(protocol=signature,cohorts=2,references=2,methods=6,probes=3,
            files=receipt(work,['final/method_summary.csv','final/paired_sensitivity.csv','final/validation.json'])))
        print('CLINICAL COMPLETE',final,flush=True)
