"""Equal-donor summaries with paired Monte Carlo uncertainty and explicit n."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
from theory import variance_decomposition

METRICS = ['fdp','power','availability','coverage','delivered_coverage','null_coverage',
           'null_delivered_coverage','proportion_rmse','abs_leakage','mean_se','bias',
           'exposure_r2','exposure_residual_ss','leverage','marker_library_mean','marker_mass_fraction']
CONDITIONAL_FIT = {'coverage','null_coverage','abs_leakage','mean_se','bias','leverage'}
GROUP = ['reference_source','target_source','grid_id','fraction','amplitude','selection','outcome','method']


def equal_donor(frame, keys, value='value'):
    rows, donors = [], []
    for key, f in frame.groupby(keys,dropna=False):
        attrs=dict(zip(keys,key if isinstance(key,tuple) else (key,)))
        ds=f.groupby('target_donor',dropna=False)[value].agg(['mean','var','count'])
        valid=ds['mean'].notna(); n=int(valid.sum())
        selected_mean=float(ds.loc[valid,'mean'].mean()) if n else np.nan
        mean=selected_mean if n==len(ds) else np.nan
        variances=(ds.loc[valid,'var']/ds.loc[valid,'count']).to_numpy()
        se=float(np.sqrt(variances.sum())/n) if n==len(ds) and n and np.isfinite(variances).all() else np.nan
        rows.append(dict(attrs,mean=mean,mcse=se,low=mean-1.96*se,high=mean+1.96*se,
                         planned_donors=len(ds),contributing_donors=n,contributing_repetitions=int(ds['count'].sum()),
                         contributing_donor_mean=selected_mean,
                         uncertainty='Monte Carlo conditional on fixed donor banks; not population donor uncertainty'))
        for donor,row in ds.iterrows():
            donors.append(dict(attrs,target_donor=donor,mean=row['mean'],repetitions=int(row['count'])))
    return pd.DataFrame(rows),pd.DataFrame(donors)


def power_audit(frame, unit):
    rows=[]
    methods=['tmm','marker_matched','marker_denominator_only']
    for (outcome,gene),f in frame.groupby(['outcome','gene_id'],dropna=False):
        available=f.pivot(index='rep',columns='method',values='available').reindex(columns=methods)
        common=available.fillna(0).eq(1).all(axis=1)
        reps=available.index[common]
        blocks={m:f[f.method.eq(m)&f.rep.isin(reps)].sort_values('rep') for m in methods}
        for m,d in blocks.items():
            C.require(len(d)==len(reps),'Power common-fit denominator mismatch')
            a=d.beta.to_numpy()-d.truth.to_numpy()
            row=variance_decomposition(a,d.leakage.to_numpy(),d.noise_projection.to_numpy())
            if len(d):
                row.update(truth=float(d.truth.iloc[0]),rms_se=float(np.sqrt(np.mean(d.se**2))),
                           mean_residual_variance=float(d.residual_variance.mean()),
                           mean_max_sample_hc3_contribution=float(d.max_sample_hc3_contribution.mean()),
                           mean_t=float((d.beta/d.se).mean()),
                           mean_absolute_t=float((d.beta/d.se).abs().mean()),
                           conditional_rejection=float(d.q.le(.05).mean()),
                           mean_beta=float(d.beta.mean()))
            baseline=blocks['tmm']
            if len(d):
                row.update(paired_beta_difference=float(np.mean(d.beta.to_numpy()-baseline.beta.to_numpy())),
                           paired_se_difference=float(np.mean(d.se.to_numpy()-baseline.se.to_numpy())))
            rows.append(dict(unit_id=unit['unit_id'],reference_source=unit['reference_source'],
                             target_source=unit['target_source'],target_donor=unit['target_donor'],
                             outcome=outcome,gene_id=gene,method=m,planned_repetitions=len(available),
                             denominator='common_available_repetitions_all_three_methods',**row))
    return rows


def summarize(output,units,options):
    output=Path(output); final=output/'final_results'; final.mkdir(exist_ok=True)
    all_records=[]; all_theory=[]; all_inv=[]; powers=[]
    for unit in units:
        blocks=[]
        for start in range(0,options['repetitions'],options['batch_size']):
            chunk=output/'units'/unit['unit_id']/f'chunk_{start:05d}'
            C.require(D.valid_checkpoint(chunk,options['signature']),'Incomplete or damaged grid checkpoint')
            attrs={k:unit[k] for k in ('unit_id','reference_source','target_source','target_donor')}
            for name,target in [('replicates.csv.gz',all_records),('theory_checks.csv',all_theory),('invariance.csv',all_inv)]:
                f=pd.read_csv(chunk/name,keep_default_na=False,na_values=[''])
                for k,v in attrs.items(): f[k]=v
                target.append(f)
            blocks.append(pd.read_csv(chunk/'power_gene_fits.csv.gz',keep_default_na=False,na_values=['']))
        powers.extend(power_audit(pd.concat(blocks,ignore_index=True),unit))
    records=pd.concat(all_records,ignore_index=True)
    C.require(not records.duplicated(['unit_id','rep','case_id','method']).any(),'Duplicate grid fits')
    counts=records.groupby(['unit_id','case_id','method']).rep.nunique()
    C.require(counts.eq(options['repetitions']).all(),'Unequal repetition denominators')
    records.to_csv(output/'replicate_audit.csv.gz',index=False,compression='gzip')
    # Availability and BH/FDP/power retain every planned fit; conditional metrics
    # are paired only when both fits satisfy the declared complete-fit rule.
    summaries=[]; donors=[]; pairs=[]
    base=records[records.method.eq('tmm')]
    key=['unit_id','rep','case_id']
    for metric in METRICS:
        a,b=equal_donor(records.assign(value=pd.to_numeric(records[metric])),GROUP)
        a['metric']=metric; b['metric']=metric
        summaries.append(a); donors.append(b)
        for method in ['marker_matched','marker_denominator_only','observed','oracle']:
            f=records[records.method.eq(method)].merge(base[key+[metric,'availability']].loc[:,~pd.Index(key+[metric,'availability']).duplicated()],
                 on=key,how='inner',validate='one_to_one',suffixes=('','_tmm'))
            ref=metric+'_tmm'
            C.require(ref in f,'Missing paired metric')
            common=np.isfinite(f[metric]) & np.isfinite(f[ref])
            if metric in CONDITIONAL_FIT: common &= f.availability.eq(1) & f.availability_tmm.eq(1)
            f['value']=np.where(common,f[metric]-f[ref],np.nan)
            pair,_=equal_donor(f,GROUP)
            pair['metric']=metric; pair['reference']='tmm'
            pair['denominator']='common_complete_fits' if metric in CONDITIONAL_FIT else 'finite_paired_values_all_planned_fits'
            pairs.append(pair)
    summary=pd.concat(summaries,ignore_index=True)
    paired=pd.concat(pairs,ignore_index=True)
    C.write_csv(final/'grid_summary.csv',summary)
    C.write_csv(final/'grid_donor_summary.csv',pd.concat(donors,ignore_index=True))
    C.write_csv(final/'grid_paired_vs_TMM.csv',paired)
    C.write_csv(final/'power_mechanism_gene_audit.csv',powers)
    inv=pd.concat(all_inv,ignore_index=True); theory=pd.concat(all_theory,ignore_index=True)
    C.write_csv(final/'marker_invariance_checks.csv',inv)
    C.write_csv(final/'theory_checks_summary.csv',theory.groupby(['kind','eligible'],dropna=False).agg(
        checks=('case_id','size'),max_violation=('max_violation','max')).reset_index())
    theory.to_csv(output/'theory_checks_full.csv.gz',index=False,compression='gzip')
    failures=records.groupby(GROUP+['reason'],dropna=False).size().reset_index(name='planned_fits')
    C.write_csv(final/'failure_reasons.csv',failures)
    C.write_json(final/'numerical_audit.json',dict(software_fixture=options['software_fixture'],
        units=len(units),repetitions=options['repetitions'],planned_fits=len(records),
        invariant_planned=len(inv),invariant_finite=int(inv.finite.sum()),
        invariant_unavailable=int((~inv.finite).sum()),
        invariant_max_change=float(inv.max_change.max()),
        bound_max_violation=float(theory.max_violation.max()),
        coefficient_checks_unavailable=int((theory.kind.eq('coefficient') & ~theory.eligible).sum()),
        interpretation='Finite bounds are conditional identities, not FDR or global estimator-stability guarantees.'))
    figures(final,paired,options['software_fixture'])
    return final


def figures(final,paired,fixture=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
    f=paired[paired.method.eq('marker_matched') & paired.outcome.eq('up05') & paired.selection.eq('gene_id') & paired.fraction.gt(0)]
    routes=list(f.groupby(['reference_source','target_source'],sort=True))
    fig,axes=plt.subplots(len(routes),2,figsize=(7.2,2.65*len(routes)),squeeze=False,layout='constrained')
    for i,((reference,target),data) in enumerate(routes):
        for j,metric in enumerate(['fdp','power']):
            d=data[data.metric.eq(metric)]
            table=d.pivot(index='amplitude',columns='fraction',values='mean').reindex(index=[.25,.5,1.],columns=[.05,.1,.2])*100
            vals=table.to_numpy(); lim=max(1,float(np.nanmax(np.abs(vals))))
            ax=axes[i,j]
            im=ax.imshow(vals,cmap='RdBu_r',norm=TwoSlopeNorm(0,-lim,lim),aspect='auto')
            for y in range(3):
                for x in range(3):
                    v=vals[y,x]
                    ax.text(x,y,f'{v:+.1f}' if np.isfinite(v) else 'NA',ha='center',va='center',
                            color='white' if np.isfinite(v) and abs(v)>lim*.55 else 'black')
            ax.set(xticks=range(3),xticklabels=['5%','10%','20%'],yticks=range(3),yticklabels=['0.25','0.5','1.0'],
                   xlabel='Requested marker fraction',ylabel='Added-count amplitude',
                   title=f'{reference} to {target}: '+('FDR' if metric=='fdp' else 'Power'))
            fig.colorbar(im,ax=ax,label='MRN - TMM (percentage points)',shrink=.8)
    if fixture: fig.suptitle('SOFTWARE FIXTURE ONLY')
    fig.savefig(final/'marker_stress_grid.png',dpi=300)
    fig.savefig(final/'marker_stress_grid.pdf')
    plt.close(fig)
    (final/'figure_caption.txt').write_text(
        'Paired MRN minus TMM differences, alternative outcome up05. Values are equal-donor means in percentage points. '
        'FDR is mean per-repetition false discovery proportion; power includes unavailable tests in the planned denominator. '
        'The exact selected marker counts (ceil within each class), alternative marker set and all null results are retained in the tables. '
        'Monte Carlo intervals conditional on fixed donor banks are in grid_paired_vs_TMM.csv. Blue is a negative difference, '
        'not uniformly a better result: a reduction in power is a cost. Each panel has its own symmetric color scale.\n',encoding='utf-8')
