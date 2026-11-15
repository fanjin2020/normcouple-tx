"""Fresh paired computation using the frozen count generator and marker budgets."""
from pathlib import Path
import tempfile
import time
import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from extension_stats import evaluate
from selection import CASES
from reference_allocation import reference_stress
from support import atomic_json, checkpoint, seal
from coupling_math import exact_change, local_report

HERE = Path(__file__).resolve().parent
METHODS = ['tmm','tmm_exclude_factor','tmm_exclude_both',
           'marker_denominator_only','marker_matched']
INVARIANT = ['tmm_exclude_both','marker_denominator_only','marker_matched']
JOBS, OPTIONS = [], {}


def folder_for(output, task):
    unit, start = task
    return Path(output)/'batches'/unit['unit_id']/('r%05d' % start)


def exclusion_anchor(job, base, A):
    marker = job['panels']['balanced']
    class_mass = base['ref']['counts'][marker].sum(0)/A.sum(0)
    excluded = base['ref']['counts'][job['test']].sum()
    value = float(1-excluded/class_mass.sum())
    C.require(0 < value <= 1, 'Invalid non-test reference mass anchor')
    return value


def read_predictions(work, expected):
    rr = pd.read_csv(work/'outputs.csv',keep_default_na=False)
    C.require(not rr.duplicated(['key','method','step']).any(),'Duplicated R prediction')
    C.require(set(zip(rr.key,rr.method,rr.step)) == expected,'Missing/unexpected R prediction')
    pred = {}
    for row in rr.itertuples():
        path = (work/row.file).resolve()
        C.require(path.is_relative_to(work.resolve()),'Unsafe R output path')
        a = np.loadtxt(path,delimiter=',',ndmin=2)
        C.require(a.shape[1]==8,'Expected denominator plus seven fractions')
        pred[row.key,row.method,row.step] = (a[:,1:],a[:,0],row.reason,float(row.scale))
    return pred, rr


def evaluate_probe(data, p0, probes, denominators, marker, scale, h, fit, cov=None):
    if not fit['basic'] or not fit['good'].all():
        return dict(status='UNAVAILABLE_BASE_FIT',score=np.nan), []
    try:
        record, changes = local_report(data['group'],data['y'],p0,probes,h,fit['se'],cov)
        dp, dm = denominators[h/2],denominators[-h/2]
        d0 = denominators[0.]
        slope = (np.log(dp)-np.log(dm))/h
        factor = (marker*scale)/(marker*scale+d0[:,None])
        analytic = -factor*slope[:,None]/np.log(2)
        observed = (np.log2(marker/dp[:,None]*scale+1)-
                    np.log2(marker/dm[:,None]*scale+1))/h
        record.update(input_chain_max_error=float(np.max(abs(analytic-observed))),
                      log_denominator_direction_rms=float(np.sqrt(np.mean(slope**2))))
        if record.get('local_smoothness_warning',False):
            record['status']='FINITE_DIFFERENCE_SENSITIVITY_NONSMOOTH'
        for row in changes:
            row['local_smoothness_warning']=record.get('local_smoothness_warning',True)
            row['derivative_agreement']=record.get('derivative_agreement',np.nan)
        return record, changes
    except (ValueError,np.linalg.LinAlgError) as exc:
        return dict(status='UNAVAILABLE_PROBE:'+str(exc),score=np.nan), []


def run_batch(task):
    unit, start = task
    opt = OPTIONS
    job = JOBS[unit['job_index']]
    folder = folder_for(opt['output'],task)
    old = checkpoint(folder,opt['signature'])
    if old:
        return dict(old,resumed=True)
    folder.mkdir(parents=True,exist_ok=True)
    began = time.monotonic()
    reps = range(start,min(start+opt['batch_size'],opt['repetitions']))
    h = opt['derivative_step']
    allsteps = [0.,-h,h,-h/2,h/2,*opt['probe_steps']]
    records, sensitivities, approximations, invariances, stress_records, contrasts = [],[],[],[],[],[]
    with tempfile.TemporaryDirectory(prefix='working_',dir=folder) as tmp:
        work = Path(tmp)
        payload, entries, expected = {},[],set()
        axes = C.E.read_json(Path(job['reference_folder'])/'axes.json') if hasattr(C.E,'read_json') else None
        if axes is None:
            import json
            axes = json.loads((Path(job['reference_folder'])/'axes.json').read_text())
        for rep in reps:
            base = C.draw_base(job,unit,rep,opt)
            A, marker_yield, anchor = marker_properties(job,base)
            non_anchor = exclusion_anchor(job,base,A)
            afile = 'A_%d.bin' % rep
            C.write_bin(work/afile,A)
            for subset in opt['subset_ids']:
                stress, detail = reference_stress(job,base,subset)
                stress_records.append(dict(unit_id=unit['unit_id'],rep=rep,subset_id=subset,**detail))
                stressed = base['draw']['counts'][:,job['norm_genes']].T.copy()
                stressed[np.array(axes['marker_rows'])-1,:] += stress.T
                marker = stressed[np.array(axes['marker_rows'])-1,:].T
                for case,effect in CASES:
                    data = C.intervene(job,base,dict(name=case,effect=effect,pattern='up',panel='frozen'),opt)
                    if rep == 0 and case in ('null','up05'):
                        C.verify_baseline(job,unit,rep,opt,data,effect)
                    observed = stressed.copy()
                    observed[np.array(axes['test_rows'])-1,:] += data['additions'].T
                    key = 'r%d_%s_%s' % (rep,subset,case)
                    total = data['total']+stress.sum(1)
                    C.write_bin(work/(key+'.bin'),observed)
                    C.write_bin(work/(key+'_total.bin'),total)
                    steps = allsteps if rep < opt['diagnostic_repetitions'] and case in ('null','up025') else [0.]
                    entries.append(dict(key=key,n=unit['n'],raw=key+'.bin',A=afile,
                        total=key+'_total.bin',group=base['group'].tolist(),marker_anchor=anchor,
                        non_test_anchor=non_anchor,steps=steps))
                    payload[key] = (rep,subset,case,effect,data,marker_yield,base['ref']['pooled_yield'],marker,steps)
                    expected.update((key,m,t) for m in METHODS for t in steps)
        C.write_json(work/'jobs.json',dict(axes=axes,entries=entries))
        r = opt['r_info']
        D.run_r([r['executable'],'--vanilla',str(HERE/'bridge.R'),str(work),r['library']],
                r['library'],folder/'bridge.log',opt['r_timeout'])
        predictions, rr = read_predictions(work,expected)
        meta = {k:unit[k] for k in ('unit_id','reference_source','target_source','target_donor')}
        for key,(rep,subset,case,effect,data,my,fy,marker,steps) in payload.items():
            identity = dict(meta,rep=rep,subset_id=subset,outcome=case,effect=effect)
            pmap = {}
            for method in ['oracle',*METHODS]:
                if method=='oracle':
                    props, reason = data['p'],''
                else:
                    q,d,reason,scale = predictions[key,method,0.]
                    C.require(q.shape==(unit['n'],7),'Sample/column shape mismatch')
                    props = C.E.convert_available_rows(q,my if method=='marker_matched' else fy)
                    pmap[method] = props
                row, arrays, fit = evaluate(data,props,'full_six',reason)
                valid = arrays['available'].astype(bool)
                hit = valid & (arrays['q']<=.05)
                null = data['truth']==0
                row.update(planned_tests=len(null),valid_tests=int(valid.sum()),
                    discoveries=int(hit.sum()),false_discoveries=int((hit&null).sum()),
                    true_discoveries=int((hit&~null).sum()),
                    coefficient_error_rms=float(np.sqrt(np.mean((arrays['beta'][valid]-data['truth'][valid])**2))) if valid.any() else np.nan)
                records.append(dict(identity,method=method,**row))
                if method=='oracle' or len(steps)==1:
                    continue
                pp,dd = {},{0.:d}
                for t in steps:
                    if t==0.:
                        continue
                    qt,dt,_,_ = predictions[key,method,t]
                    pp[t] = C.E.convert_available_rows(qt,my if method=='marker_matched' else fy)
                    dd[t] = dt
                    if method in INVARIANT:
                        finite = np.isfinite(q).all() and np.isfinite(qt).all()
                        err = float(np.max(abs(qt-q))) if finite else np.nan
                        C.require(np.array_equal(np.isfinite(qt),np.isfinite(q)),'Invariant arm availability changed')
                        C.require(not finite or err<1e-10,'Complete exclusion failed input invariance')
                        invariances.append(dict(identity,method=method,step=t,prediction_finite=finite,max_prediction_change=err))
                score, changes = evaluate_probe(data,props,pp,dd,marker,scale,h,fit)
                sensitivities.append(dict(identity,method=method,**score,
                    fdp=row['fdp'],coefficient_error_rms=row['coefficient_error_rms'],
                    coverage=row['coverage'],availability=row['availability']))
                approximations.extend(dict(identity,method=method,**v) for v in changes)
            for method in METHODS[1:]:
                try:
                    actual,exact = exact_change(data['group'],data['y'],pmap['tmm'],pmap[method])
                    error = float(np.max(abs(actual-exact)))
                    C.require(error < 1e-7*(1+np.max(abs(actual))),'Coefficient-change identity failed')
                    contrasts.append(dict(identity,method=method,status='ESTIMABLE',
                        coefficient_change_rms=float(np.sqrt(np.mean(actual**2))),exact_identity_max_error=error))
                except (ValueError,np.linalg.LinAlgError) as exc:
                    if 'identity failed' in str(exc):
                        raise
                    contrasts.append(dict(identity,method=method,status='UNAVAILABLE:'+str(exc)))
        rr.insert(0,'unit_id',unit['unit_id'])
        rr.to_csv(folder/'r_diagnostics.csv.gz',index=False,compression='gzip')
    tables = dict(metrics=records,sensitivity=sensitivities,approximation=approximations,
                  invariance=invariances,stress=stress_records,coefficient_contrasts=contrasts)
    names = ['bridge.log','r_diagnostics.csv.gz']
    for name,rows in tables.items():
        path = folder/(name+'.csv.gz')
        pd.DataFrame(rows).to_csv(path,index=False,compression='gzip')
        names.append(path.name)
    C.require(len(records)==len(reps)*len(opt['subset_ids'])*len(CASES)*(len(METHODS)+1),'Incomplete metric grid')
    details = dict(unit=unit['unit_id'],start=start,repetitions=len(reps),
                   work_items=len(reps)*len(opt['subset_ids']),metric_rows=len(records),
                   elapsed_seconds=time.monotonic()-began,finished_unix=time.time())
    atomic_json(folder/'batch.json',details)
    seal(folder,[*names,'batch.json'],opt['signature'],**details)
    return dict(details,resumed=False)
