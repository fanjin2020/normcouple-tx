"""Matched perturbations, native linear NNLS controls and paired propagation."""
from pathlib import Path
import json
import tempfile
import time
import numpy as np
import pandas as pd
from scipy.optimize import nnls
import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from extension_stats import evaluate
from selection import CASES
from support import checkpoint,seal,atomic_json
from experiment import METHODS as DT_METHODS, read_predictions, exclusion_anchor, evaluate_probe
from coupling_math import exact_change, geometry, geometry_change
from matched_markers import stress,split_variant

HERE=Path(__file__).resolve().parent
NNLS_METHODS=['nnls_linear','nnls_marker_matched']
METHODS=[*DT_METHODS,*NNLS_METHODS]
INVARIANT=['tmm_exclude_both','marker_denominator_only','marker_matched',*NNLS_METHODS]
JOBS,OPTIONS=[],{}


def folder_for(output,task):
    unit,start=task
    return Path(output)/'batches'/unit['unit_id']/('r%05d'%start)


def nnls_prediction(A,marker,denominator):
    A,marker,denominator=np.asarray(A,float),np.asarray(marker,float),np.asarray(denominator,float)
    C.require(A.shape==(marker.shape[1],7),'Invalid NNLS reference dimensions')
    q=np.full((len(marker),7),np.nan)
    errors=[]
    for i in range(len(marker)):
        try:
            if not np.isfinite(denominator[i]) or denominator[i]<=0:
                raise ValueError('Invalid NNLS denominator')
            coef,_=nnls(A,marker[i]/denominator[i],maxiter=1000)
            if coef.sum()<=0:
                raise ValueError('Empty NNLS prediction')
            q[i]=coef/coef.sum()
        except (ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            errors.append(str(exc))
    return q,'; '.join(sorted(set(errors)))


def unavailable(identity,method,detail):
    return dict(identity,method=method,intervention_feasible=False,
        fit_status='DESIGN_INFEASIBLE:'+detail['infeasibility_reason'],
        availability=np.nan,coverage=np.nan,delivered_coverage=np.nan,fdp=np.nan,
        power=np.nan,proportion_rmse=np.nan,bias=np.nan,mean_se=np.nan,
        coefficient_error_rms=np.nan,exposure_r2=np.nan,exposure_residual_ss=np.nan,
        leverage=np.nan,planned_tests=100,valid_tests=0,discoveries=0,
        false_discoveries=0,true_discoveries=0)


def run_batch(task):
    unit,start=task
    opt=OPTIONS
    job=JOBS[unit['job_index']]
    folder=folder_for(opt['output'],task)
    old=checkpoint(folder,opt['signature'])
    if old:
        return dict(old,resumed=True)
    folder.mkdir(parents=True,exist_ok=True)
    began=time.monotonic()
    reps=range(start,min(start+opt['batch_size'],opt['repetitions']))
    h=opt['derivative_step']
    allsteps=[0.,-h,h,-h/2,h/2,*opt['probe_steps']]
    metrics,sensitivities,approximations,checks,stress_rows,contrasts=[],[],[],[],[],[]
    meta={k:unit[k] for k in ('unit_id','reference_source','target_source','target_donor')}
    with tempfile.TemporaryDirectory(prefix='working_',dir=folder) as tmp:
        work=Path(tmp)
        axes=json.loads((Path(job['reference_folder'])/'axes.json').read_text())
        payload,entries,expected={},[],set()
        for rep in reps:
            base=C.draw_base(job,unit,rep,opt)
            A,my,anchor=marker_properties(job,base)
            fy=base['ref']['pooled_yield']
            non_anchor=exclusion_anchor(job,base,A)
            afile='A_%d.bin'%rep
            C.write_bin(work/afile,A)
            budgets=set()
            for variant in opt['subset_ids']:
                additions,detail=stress(job,base,variant)
                if variant!='zero':
                    budgets.add(detail['budget_sha256'])
                stress_rows.append(dict(meta,rep=rep,subset_id=variant,**detail))
                subset,rule=split_variant(variant)
                for case,effect in CASES:
                    identity=dict(meta,rep=rep,subset_id=variant,base_subset=subset,
                                  subset_kind=subset.split('_')[0],support_rule=rule,outcome=case,effect=effect)
                    if additions is None:
                        metrics.extend(unavailable(identity,m,detail) for m in ['oracle',*METHODS])
                        if (rep<opt['diagnostic_repetitions'] and case in ('null','up025') and
                                subset in opt['diagnostic_subsets']):
                            sensitivities.extend(dict(identity,method=m,status='DESIGN_INFEASIBLE',
                                score=np.nan,fdp=np.nan,coefficient_error_rms=np.nan,
                                coverage=np.nan,availability=np.nan) for m in METHODS)
                        continue
                    data=C.intervene(job,base,dict(name=case,effect=effect,pattern='up',panel='frozen'),opt)
                    if rep==0 and case in ('null','up05'):
                        C.verify_baseline(job,unit,rep,opt,data,effect)
                    raw=base['draw']['counts'][:,job['norm_genes']].T.copy()
                    raw[np.array(axes['marker_rows'])-1,:]+=additions.T
                    marker=raw[np.array(axes['marker_rows'])-1,:].T.copy()
                    raw[np.array(axes['test_rows'])-1,:]+=data['additions'].T
                    key='r%d_%s_%s'%(rep,variant,case)
                    total=data['total']+additions.sum(1)
                    C.write_bin(work/(key+'.bin'),raw)
                    C.write_bin(work/(key+'_total.bin'),total)
                    diagnostic=(rep<opt['diagnostic_repetitions'] and case in ('null','up025') and
                                subset in opt['diagnostic_subsets'])
                    steps=allsteps if diagnostic else [0.]
                    entries.append(dict(key=key,n=unit['n'],raw=key+'.bin',A=afile,total=key+'_total.bin',
                        group=base['group'].tolist(),marker_anchor=anchor,non_test_anchor=non_anchor,steps=steps))
                    payload[key]=(identity,data,A,my,fy,marker,steps)
                    expected.update((key,m,t) for m in DT_METHODS for t in steps)
            C.require(len(budgets)==1,'Subset/support variants received different budgets')
        C.write_json(work/'jobs.json',dict(axes=axes,entries=entries))
        r=opt['r_info']
        D.run_r([r['executable'],'--vanilla',str(HERE/'bridge.R'),str(work),r['library']],
                r['library'],folder/'bridge.log',opt['r_timeout'])
        pred,rr=read_predictions(work,expected)
        for key,(identity,data,A,my,fy,marker,steps) in payload.items():
            pmap={}
            for method in ['oracle',*METHODS]:
                if method=='oracle':
                    props,reason=data['p'],''
                elif method in DT_METHODS:
                    q,d,reason,scale=pred[key,method,0.]
                    props=C.E.convert_available_rows(q,my if method=='marker_matched' else fy)
                else:
                    d=pred[key,'tmm',0.][1] if method=='nnls_linear' else marker.sum(1)
                    ref=A if method=='nnls_linear' else A/A.sum(0)
                    q,reason=nnls_prediction(ref,marker,d)
                    props=C.E.convert_available_rows(q,fy if method=='nnls_linear' else my)
                pmap[method]=props
                row,arrays,fit=evaluate(data,props,'full_six',reason)
                valid=arrays['available'].astype(bool)
                hit=valid & (arrays['q']<=.05)
                null=data['truth']==0
                row.update(intervention_feasible=True,fit_status='ESTIMABLE' if fit['basic'] else 'UNAVAILABLE',
                    planned_tests=len(null),valid_tests=int(valid.sum()),discoveries=int(hit.sum()),
                    false_discoveries=int((hit&null).sum()),true_discoveries=int((hit&~null).sum()),
                    coefficient_error_rms=float(np.sqrt(np.mean((arrays['beta'][valid]-data['truth'][valid])**2))) if valid.any() else np.nan)
                metrics.append(dict(identity,method=method,**row))
                if method=='oracle' or len(steps)==1:
                    continue
                probes,denoms={},{0.:d}
                for t in steps:
                    if t==0.:
                        continue
                    if method in DT_METHODS:
                        qt,dt,_,_=pred[key,method,t]
                        pt=C.E.convert_available_rows(qt,my if method=='marker_matched' else fy)
                    else:
                        dt=pred[key,'tmm',t][1] if method=='nnls_linear' else marker.sum(1)
                        qt,_=nnls_prediction(ref,marker,dt)
                        pt=C.E.convert_available_rows(qt,fy if method=='nnls_linear' else my)
                    probes[t],denoms[t]=pt,dt
                    if method in INVARIANT:
                        finite=np.isfinite(q).all() and np.isfinite(qt).all()
                        change=float(np.max(abs(q-qt))) if finite else np.nan
                        C.require(np.array_equal(np.isfinite(q),np.isfinite(qt)) and
                                  (not finite or change<1e-8),'Expected prediction invariance failed: '+method)
                        checks.append(dict(identity,method=method,step=t,prediction_finite=finite,max_prediction_change=change))
                if method in DT_METHODS:
                    score,changes=evaluate_probe(data,props,probes,denoms,marker,scale,h,fit)
                else:
                    from coupling_math import local_report
                    try:
                        score,changes=local_report(data['group'],data['y'],props,probes,h,fit['se'])
                    except (ValueError,np.linalg.LinAlgError) as exc:
                        score,changes=dict(status='UNAVAILABLE:'+str(exc),score=np.nan),[]
                sensitivities.append(dict(identity,method=method,**score,fdp=row['fdp'],
                    coefficient_error_rms=row['coefficient_error_rms'],coverage=row['coverage'],availability=row['availability']))
                approximations.extend(dict(identity,method=method,**v) for v in changes)
            for method in ['oracle',*METHODS[1:]]:
                try:
                    actual,exact=exact_change(data['group'],data['y'],pmap['tmm'],pmap[method])
                    err=float(np.max(abs(actual-exact)))
                    C.require(err<1e-7*(1+np.max(abs(actual))),'Coefficient identity failed')
                    contrasts.append(dict(identity,method=method,status='ESTIMABLE',
                        coefficient_change_rms=float(np.sqrt(np.mean(actual**2))),exact_identity_max_error=err,
                        **geometry_change(data['group'],pmap['tmm'],pmap[method])))
                except (ValueError,np.linalg.LinAlgError) as exc:
                    if 'identity failed' in str(exc):
                        raise
                    contrasts.append(dict(identity,method=method,status='UNAVAILABLE:'+str(exc)))
            if np.isfinite(pmap['nnls_linear']).all() and np.isfinite(pmap['nnls_marker_matched']).all():
                err=float(np.max(abs(pmap['nnls_linear']-pmap['nnls_marker_matched'])))
                checks.append(dict(identity,method='nnls_reference_yield_equivalence',step=0.,
                                   prediction_finite=True,max_prediction_change=err))
        rr.insert(0,'unit_id',unit['unit_id'])
        rr.to_csv(folder/'r_diagnostics.csv.gz',index=False,compression='gzip')
    names=['bridge.log','r_diagnostics.csv.gz']
    for name,rows in dict(metrics=metrics,sensitivity=sensitivities,approximation=approximations,
                         invariance=checks,stress=stress_rows,coefficient_contrasts=contrasts).items():
        pd.DataFrame(rows).to_csv(folder/(name+'.csv.gz'),index=False,compression='gzip')
        names.append(name+'.csv.gz')
    C.require(len(metrics)==len(reps)*len(opt['subset_ids'])*3*(len(METHODS)+1),'Incomplete planned grid')
    details=dict(unit=unit['unit_id'],start=start,repetitions=len(reps),work_items=len(reps)*len(opt['subset_ids']),
                 metric_rows=len(metrics),elapsed_seconds=time.monotonic()-began,finished_unix=time.time())
    atomic_json(folder/'batch.json',details)
    seal(folder,[*names,'batch.json'],opt['signature'],**details)
    return dict(details,resumed=False)
