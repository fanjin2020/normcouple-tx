"""Compact outputs with donor-separated diagnostics and explicit uncertainty scopes."""
from pathlib import Path
import hashlib
import shutil
import numpy as np
import pandas as pd
from scipy import stats
from common import save_json


def read_nonempty(path):
    try:
        return pd.read_csv(path,keep_default_na=False,na_values=[''])
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def auc(labels,scores):
    labels=np.asarray(labels,bool)
    if not labels.any() or labels.all():
        return np.nan
    ranks=stats.rankdata(scores)
    n=int(labels.sum()); m=len(labels)-n
    return float((ranks[labels].sum()-n*(n+1)/2)/(n*m))


def validate_diagnostic(frame):
    frame=frame.copy()
    frame['donor_group']=frame.target_source.astype(str)+'|'+frame.target_donor.astype(str)
    groups=sorted(frame.donor_group.unique(),key=lambda s:hashlib.sha256(('p3-coupling-fold|'+s).encode()).hexdigest())
    if len(groups)<5:
        return pd.DataFrame(),pd.DataFrame([dict(status='INSUFFICIENT_DONORS',donors=len(groups))])
    folds={g:i%5 for i,g in enumerate(groups)}
    frame['fold']=frame.donor_group.map(folds)
    frame['high_realized_fdp']=frame.fdp>.10
    frame['alert']=False
    frame['threshold']=np.nan
    frame['evaluated']=False
    rows=[]
    for (method,outcome),allrows in frame.groupby(['method','outcome']):
        for fold in range(5):
            train=allrows[allrows.fold.ne(fold)&np.isfinite(allrows.score)&np.isfinite(allrows.fdp)]
            test=allrows[allrows.fold.eq(fold)&np.isfinite(allrows.score)&np.isfinite(allrows.fdp)]
            if train.empty or test.empty:
                continue
            assert not set(train.donor_group)&set(test.donor_group)
            threshold=max(1e-6,float(train.score.quantile(.9)))
            alert=test.score>threshold
            frame.loc[test.index,'alert']=alert
            frame.loc[test.index,'threshold']=threshold
            frame.loc[test.index,'evaluated']=True
            labels=test.high_realized_fdp
            rows.append(dict(method=method,outcome=outcome,fold=fold,status='HELD_OUT',
                train_donors=train.donor_group.nunique(),test_donors=test.donor_group.nunique(),
                test_records=len(test),unavailable_test_records=int((allrows.fold.eq(fold)).sum()-len(test)),
                threshold=threshold,alert_rate=float(alert.mean()),
                high_realized_fdp_rate=float(labels.mean()),auc=auc(labels,test.score),
                sensitivity=float(alert[labels].mean()) if labels.any() else np.nan,
                specificity=float((~alert[~labels]).mean()) if (~labels).any() else np.nan,
                interpretation='Predicts realized FDP event in simulation; not population FDR control'))
    if not rows:
        rows=[dict(status='NO_VALID_DIAGNOSTIC_SCORES',donors=len(groups))]
    return frame,pd.DataFrame(rows)


def summarize(output,cfg,mode):
    from support import checkpoint
    from run_coupling import tasks
    from experiment import folder_for
    output=Path(output)
    final=output/'final'
    final.mkdir(exist_ok=True)
    planned=list(tasks(cfg))
    folders=[]
    for task in planned:
        folder=folder_for(output,task)
        if checkpoint(folder,cfg['signature']):
            folders.append(folder)
    simulation_complete=len(folders)==len(planned)
    if mode in ('all','simulation') and not simulation_complete:
        raise ValueError('Cannot label a partial simulation complete: %d/%d sealed batches' % (len(folders),len(planned)))
    diagnostic_records,diagnostic_valid=0,0
    if simulation_complete:
        merged={name:pd.concat([read_nonempty(f/(name+'.csv.gz')) for f in folders],ignore_index=True)
                for name in ('metrics','sensitivity','approximation','invariance','coefficient_contrasts',
                             'stress','r_diagnostics')}
        metrics=merged['metrics']
        if metrics.duplicated(['unit_id','rep','subset_id','outcome','method']).any():
            raise ValueError('Duplicated paired result keys')
        expected=len(cfg['units'])*cfg['settings']['repetitions']*len(cfg['settings']['subset_ids'])*3*6
        if len(metrics)!=expected:
            raise ValueError('Incomplete planned metric table')
        metrics['subset_kind']=np.where(metrics.subset_id.str.startswith('random_'),'random',metrics.subset_id)
        fields=['fdp','power','availability','coverage','delivered_coverage','proportion_rmse',
                'bias','mean_se','coefficient_error_rms','exposure_r2','exposure_residual_ss','leverage']
        keys=['unit_id','reference_source','target_source','target_donor','subset_kind','outcome','method']
        paired=metrics.groupby(keys+['rep'],as_index=False)[fields].mean()
        donor=paired.groupby(keys,as_index=False)[fields].mean()
        donor.to_csv(final/'donor_results.csv',index=False)
        mc=paired.groupby(keys)[fields].agg(['mean','sem','count'])
        mc.columns=['_'.join(c) for c in mc.columns]
        mc.reset_index().to_csv(final/'within_unit_monte_carlo.csv',index=False)
        summary=[]
        for key,g in donor.groupby(['reference_source','target_source','subset_kind','outcome','method']):
            row=dict(zip(['reference_source','target_source','subset_kind','outcome','method'],key))
            row['target_donors']=g.target_donor.nunique()
            for col in fields:
                row[col]=float(g[col].mean())
                row[col+'_donor_q25']=float(g[col].quantile(.25))
                row[col+'_donor_q75']=float(g[col].quantile(.75))
            summary.append(row)
        pd.DataFrame(summary).to_csv(final/'method_summary.csv',index=False)
        index=['unit_id','rep','subset_id','outcome']
        base=metrics[metrics.method.eq('tmm')].set_index(index)
        delta=[]
        for method,g in metrics[~metrics.method.eq('tmm')].groupby('method'):
            own=g.set_index(index)
            if not own.index.is_unique or set(own.index)!=set(base.index):
                raise ValueError('Unpaired method grid')
            other=base.reindex(own.index)
            d=own[['reference_source','target_source','target_donor','subset_kind']].copy()
            for field in fields:
                d[field+'_minus_tmm']=own[field]-other[field]
            d['method']=method
            delta.append(d.reset_index())
        differences=pd.concat(delta)
        deltafields=[f+'_minus_tmm' for f in fields]
        differences.groupby(keys,as_index=False)[deltafields].mean().to_csv(final/'paired_donor_differences.csv',index=False)
        sensitivity=merged['sensitivity']
        diagnostic_records=len(sensitivity)
        expected_diagnostics=len(cfg['units'])*cfg['settings']['diagnostic_repetitions']*len(cfg['settings']['subset_ids'])*2*5
        if diagnostic_records!=expected_diagnostics:
            raise ValueError('Missing planned sensitivity records')
        diagnostic_valid=int(np.isfinite(sensitivity.score).sum())
        if len(sensitivity):
            heldout,performance=validate_diagnostic(sensitivity)
            heldout.to_csv(final/'diagnostic_heldout.csv.gz',index=False,compression='gzip')
            performance.to_csv(final/'diagnostic_performance.csv',index=False)
        for name,table in merged.items():
            table.to_csv(final/(name+'.csv.gz'),index=False,compression='gzip')
    clinical_complete=bool(checkpoint(output/'clinical',cfg['signature']))
    if mode in ('all','clinical') and not cfg['settings']['software_fixture'] and not clinical_complete:
        raise ValueError('Requested clinical analysis is not sealed')
    if clinical_complete:
        dest=final/'clinical'
        dest.mkdir(exist_ok=True)
        receipt=checkpoint(output/'clinical',cfg['signature'])
        for filename in receipt['files']:
            shutil.copy2(output/'clinical'/filename,dest/filename)
    for filename in ('protocol.json','selected_markers.json','R_versions.json','runtime.json'):
        if (output/filename).exists():
            shutil.copy2(output/filename,final/filename)
    text=(Path(__file__).with_name('RESULTS_GUIDE.md')).read_text(encoding='utf-8')
    (final/'README.md').write_text(text,encoding='utf-8')
    return dict(simulation_complete=simulation_complete,clinical_complete=clinical_complete,
                simulation_batches=len(folders),planned_simulation_batches=len(planned),
                diagnostic_records=diagnostic_records,finite_diagnostic_scores=diagnostic_valid,
                original_fits_reused=False,public_deposition='NOT_PERFORMED')
