"""Current-run summaries only; feasible-design performance and donor-separated checks."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from common import require,save_json,read_json
from report import read_nonempty,validate_diagnostic

FIELDS=['fdp','power','coverage','delivered_coverage','availability','proportion_rmse',
        'bias','mean_se','coefficient_error_rms','exposure_r2','exposure_residual_ss','leverage']
CONTEXT=['reference_source','target_source','subset_kind','support_rule','outcome','method']
DONOR=['unit_id','target_donor',*CONTEXT]
PAIR=['unit_id','rep','subset_id','outcome']


def method_differences(metrics):
    base=metrics[metrics.method.eq('tmm')].set_index(PAIR)
    require(base.index.is_unique,'Duplicated baseline keys')
    rows=[]
    for method,group in metrics[~metrics.method.eq('tmm')].groupby('method'):
        own=group.set_index(PAIR)
        require(own.index.is_unique and set(own.index)==set(base.index),'Unpaired method grid')
        other=base.reindex(own.index)
        result=own[['reference_source','target_source','target_donor','subset_kind','support_rule']].copy()
        result['method']=method
        result['common_design_feasible']=own.intervention_feasible & other.intervention_feasible
        for field in FIELDS:
            result[field+'_minus_tmm']=(own[field]-other[field]).where(result.common_design_feasible)
        rows.append(result.reset_index())
    return pd.concat(rows,ignore_index=True)


def subset_comparisons(metrics):
    rows=[]
    keys=['unit_id','rep','outcome','method']
    for (rule,subset),group in metrics[metrics.base_subset.ne('zero')].groupby(['support_rule','base_subset']):
        targets=[('predefined_same_support','legacy',rule)]
        if subset.startswith('matched_'):
            targets.append(('matched_vs_random','random_'+subset.split('_')[1],rule))
        if rule=='positive_only':
            targets.append(('positive_vs_allow_zero',subset,'allow_zero'))
        for name,baseline,other_rule in targets:
            if subset==baseline and rule==other_rule:
                continue
            a=group.set_index(keys)
            b=metrics[metrics.base_subset.eq(baseline)&metrics.support_rule.eq(other_rule)].set_index(keys)
            require(a.index.is_unique and b.index.is_unique and set(a.index)==set(b.index),'Missing subset pair')
            b=b.reindex(a.index)
            common=a.intervention_feasible & b.intervention_feasible
            out=a[['reference_source','target_source','target_donor','subset_kind']].copy()
            out['support_rule']=rule
            out['comparison']=name
            out['base_subset']=subset
            out['common_design_feasible']=common
            for field in FIELDS:
                out[field+'_difference']=(a[field]-b[field]).where(common)
            rows.append(out.reset_index())
    return pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()


def donor_tables(metrics,final):
    # Equal subset weight within donor; equal target-donor weight within setting.
    detail_keys=['unit_id','target_donor','base_subset',*CONTEXT]
    detail=metrics.groupby(detail_keys,dropna=False)[FIELDS].mean().reset_index()
    counts=metrics.groupby(detail_keys,dropna=False).agg(planned_records=('rep','size'),
               feasible_records=('intervention_feasible','sum')).reset_index()
    detail=detail.merge(counts,on=detail_keys,validate='one_to_one')
    detail.to_csv(final/'donor_subset_results.csv',index=False)
    mc=metrics.groupby(detail_keys,dropna=False)[FIELDS].agg(['mean','sem','count'])
    mc.columns=['__'.join(c) for c in mc.columns]
    mc.reset_index().to_csv(final/'within_unit_monte_carlo.csv',index=False)
    donor=detail.groupby(DONOR,dropna=False)[FIELDS].mean().reset_index()
    donor.to_csv(final/'donor_results.csv',index=False)
    summary=donor.groupby(CONTEXT,dropna=False)[FIELDS].mean().reset_index()
    counts=metrics.groupby(CONTEXT,dropna=False).agg(planned_records=('rep','size'),
             feasible_records=('intervention_feasible','sum'),target_donors=('target_donor','nunique')).reset_index()
    summary=summary.merge(counts,on=CONTEXT,validate='one_to_one')
    spread=donor.groupby(CONTEXT,dropna=False)[FIELDS].quantile([.25,.75]).unstack()
    spread.columns=[field+('__donor_q25' if q==.25 else '__donor_q75') for field,q in spread.columns]
    summary=summary.merge(spread.reset_index(),on=CONTEXT,validate='one_to_one')
    summary['design_feasibility']=summary.feasible_records/summary.planned_records
    summary['performance_scope']='conditional on feasible intervention; equal subset and donor weighting'
    summary.to_csv(final/'method_summary.csv',index=False)
    loo=[]
    for key,group in donor.groupby(CONTEXT,dropna=False):
        for excluded in group.target_donor.unique():
            keep=group[group.target_donor.ne(excluded)]
            if len(keep):
                loo.append(dict(zip(CONTEXT,key),excluded_target_donor=excluded,
                    retained_donors=keep.target_donor.nunique(),**keep[FIELDS].mean().to_dict()))
    pd.DataFrame(loo).to_csv(final/'target_donor_leave_one_out.csv',index=False)
    return summary


def direct_controls(metrics,final):
    selected=[('marker_matched','marker_denominator_only'),('marker_matched','oracle'),
              ('tmm','oracle'),('nnls_linear','oracle')]
    rows=[]
    for method,baseline in selected:
        a=metrics[metrics.method.eq(method)].set_index(PAIR)
        b=metrics[metrics.method.eq(baseline)].set_index(PAIR)
        require(a.index.is_unique and b.index.is_unique and set(a.index)==set(b.index),'Unpaired direct controls')
        b=b.reindex(a.index)
        out=a[['reference_source','target_source','target_donor','base_subset','subset_kind','support_rule']].copy()
        out['comparison']=method+'_minus_'+baseline
        out['common_feasible']=a.intervention_feasible & b.intervention_feasible
        for field in FIELDS:
            out[field+'_difference']=(a[field]-b[field]).where(out.common_feasible)
        rows.append(out.reset_index())
    frame=pd.concat(rows,ignore_index=True)
    keys=['unit_id','target_donor','reference_source','target_source','base_subset',
          'subset_kind','support_rule','outcome','comparison']
    columns=[field+'_difference' for field in FIELDS]
    stats=frame.groupby(keys,dropna=False)[columns].agg(['mean','sem','count'])
    stats.columns=['__'.join(c) for c in stats.columns]
    counts=frame.groupby(keys,dropna=False).agg(planned_pairs=('rep','size'),common_feasible_pairs=('common_feasible','sum'))
    stats.join(counts).reset_index().to_csv(final/'paired_controls.csv',index=False)


def diagnostic_details(frame,final):
    if frame.empty or 'evaluated' not in frame:
        return
    keys=['reference_source','target_source','subset_kind','support_rule','method','outcome']
    donor_rows=[]
    for key,g in frame.groupby(keys+['donor_group'],dropna=False):
        valid=g[g.evaluated & np.isfinite(g.fdp)]
        high=valid.high_realized_fdp.astype(bool)
        alert=valid.alert.astype(bool)
        donor_rows.append(dict(zip(keys+['donor_group'],key),planned_scores=len(g),
            evaluated_scores=len(valid),valid_fraction=len(valid)/len(g),
            sensitivity=float(alert[high].mean()) if high.any() else np.nan,
            false_positive_rate=float(alert[~high].mean()) if (~high).any() else np.nan,
            alert_rate=float(alert.mean()) if len(valid) else np.nan,
            high_fdp_rate=float(high.mean()) if len(valid) else np.nan))
    donors=pd.DataFrame(donor_rows)
    donors.to_csv(final/'diagnostic_donor_performance.csv',index=False)
    rows=[]
    rng=np.random.default_rng(2026103102)
    fields=['valid_fraction','sensitivity','false_positive_rate','alert_rate','high_fdp_rate']
    for key,g in donors.groupby(keys,dropna=False):
        row=dict(zip(keys,key),donors=len(g),
                 interval_scope='target-donor bootstrap of fixed held-out predictions; no threshold refitting')
        draws=rng.integers(0,len(g),size=(1000,len(g)))
        for field in fields:
            values=g[field].to_numpy(float)
            row[field]=float(np.nanmean(values)) if np.isfinite(values).any() else np.nan
            sample=values[draws]
            count=np.isfinite(sample).sum(1)
            boot=np.divide(np.nansum(sample,axis=1),count,out=np.full(len(sample),np.nan),where=count>0)
            finite=boot[np.isfinite(boot)]
            lo,hi=np.quantile(finite,[.025,.975]) if len(g)>=3 and len(finite) else (np.nan,np.nan)
            row[field+'_ci_low'],row[field+'_ci_high']=float(lo),float(hi)
        rows.append(row)
    pd.DataFrame(rows).to_csv(final/'diagnostic_setting_performance.csv',index=False)


def make_plots(final,summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    colors={'tmm':'#2976a3','marker_matched':'#24866c','marker_denominator_only':'#b27627',
            'nnls_linear':'#8264a7','oracle':'#555555'}
    for metric,outcome in [('fdp','null'),('power','up025')]:
        data=summary[summary.outcome.eq(outcome)&summary.method.isin(colors)&summary.support_rule.eq('allow_zero')]
        settings=list(data.groupby(['reference_source','target_source'],sort=True))
        if not settings:
            continue
        fig,axes=plt.subplots(1,len(settings),figsize=(4*len(settings),4.2),squeeze=False,sharey=True)
        for ax,((reference,target),part) in zip(axes[0],settings):
            for method,g in part.groupby('method'):
                values=g.set_index('subset_kind')[metric].reindex(['legacy','random','matched'])*100
                ax.plot(range(3),values,marker='o',color=colors[method],label=method.replace('_',' '))
            ax.set_xticks(range(3),['Predefined','Random','Matched'],rotation=20)
            ax.set_title(reference+' to '+target)
            ax.spines[['top','right']].set_visible(False)
            ax.set_ylim(0,100)
            if metric=='fdp':
                ax.axhline(5,color='#888888',linestyle=':',linewidth=.8)
        axes[0,0].set_ylabel(('Empirical FDR' if metric=='fdp' else 'Power')+' (%)')
        handles,labels=axes[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False)
        fig.subplots_adjust(top=.83,bottom=.30,wspace=.2)
        for extension in ('pdf','png'):
            fig.savefig(final/('matched_'+metric+'.'+extension),dpi=220)
        plt.close(fig)


def summarize(output,cfg,mode):
    from run_coupling import tasks
    from enhancement_experiment import folder_for,METHODS
    from support import checkpoint
    output=Path(output)
    final=output/'final'
    final.mkdir(exist_ok=True)
    plan=list(tasks(cfg))
    folders=[folder_for(output,task) for task in plan if checkpoint(folder_for(output,task),cfg['signature'])]
    simulation_complete=len(folders)==len(plan)
    require(mode not in ('all','simulation') or simulation_complete,'Cannot summarize incomplete simulation')
    if simulation_complete:
        merged={name:pd.concat([read_nonempty(f/(name+'.csv.gz')) for f in folders],ignore_index=True)
                for name in ('metrics','sensitivity','approximation','invariance','stress','coefficient_contrasts')}
        metrics=merged['metrics']
        expected=len(cfg['units'])*cfg['settings']['repetitions']*len(cfg['settings']['subset_ids'])*3*(len(METHODS)+1)
        require(len(metrics)==expected and not metrics.duplicated([*PAIR,'method']).any(),'Incomplete/duplicate metric grid')
        summary=donor_tables(metrics,final)
        direct_controls(metrics,final)
        differences=method_differences(metrics)
        diff_fields=[f+'_minus_tmm' for f in FIELDS]
        differences.groupby(DONOR,dropna=False)[diff_fields].mean().reset_index().to_csv(final/'paired_method_differences.csv',index=False)
        comparisons=subset_comparisons(metrics)
        if len(comparisons):
            keys=['unit_id','target_donor','reference_source','target_source','comparison',
                  'base_subset','subset_kind','support_rule','outcome','method']
            fields=[f+'_difference' for f in FIELDS]
            means=comparisons.groupby(keys,dropna=False)[fields].mean().reset_index()
            counts=comparisons.groupby(keys,dropna=False).agg(planned_pairs=('rep','size'),
                           common_feasible_pairs=('common_design_feasible','sum')).reset_index()
            means.merge(counts,on=keys,validate='one_to_one').to_csv(final/'paired_subset_differences.csv',index=False)
        stress=merged['stress']
        skeys=['unit_id','reference_source','target_source','target_donor','base_subset','subset_kind','support_rule']
        sr=stress.groupby(skeys,dropna=False).agg(planned_draws=('rep','size'),
            feasible_draws=('intervention_feasible','sum'),
            requested_counts=('requested_budget_total','mean'),
            unsupported_sample_classes=('unsupported_sample_class_count','sum'),
            effective_relative_amplitude=('effective_relative_amplitude','mean'),
            target_mean_selected_count=('target_mixture_mean_selected_count','mean'),
            target_selected_zero_fraction=('target_mixture_selected_zero_fraction','mean'),
            target_control_selected_zero_fraction=('target_control_selected_zero_fraction','mean'),
            added_counts_on_observed_zero=('added_counts_on_observed_zero','sum'),
            total_added_counts=('total_added_counts','sum')).reset_index()
        sr['zero_injection_fraction']=sr.added_counts_on_observed_zero/sr.total_added_counts.replace(0,np.nan)
        sr.to_csv(final/'perturbation_support.csv',index=False)
        sensitivity=merged['sensitivity']
        diagnostic_variants=[s for s in cfg['settings']['subset_ids']
                             if s.split('__')[0] in cfg['settings']['diagnostic_subsets']]
        expected_scores=(len(cfg['units'])*min(cfg['settings']['diagnostic_repetitions'],cfg['settings']['repetitions'])*
                         len(diagnostic_variants)*2*len(METHODS))
        require(len(sensitivity)==expected_scores and not sensitivity.duplicated([*PAIR,'method']).any(),
                'Missing or duplicate planned diagnostic records')
        heldout,performance=validate_diagnostic(sensitivity)
        diagnostic_details(heldout,final)
        performance.to_csv(final/'diagnostic_performance.csv',index=False)
        # Compact scores retain fold membership and allow inspection of false alerts.
        heldout.to_csv(final/'diagnostic_heldout.csv.gz',index=False,compression='gzip')
        sensitivity.groupby(['reference_source','target_source','subset_kind','support_rule','outcome','method','status'],
            dropna=False).size().rename('records').reset_index().to_csv(final/'diagnostic_status_counts.csv',index=False)
        for name in ('approximation','invariance','coefficient_contrasts'):
            frame=merged[name]
            if frame.empty:
                pd.DataFrame([dict(status='NO_RECORDS')]).to_csv(final/(name+'_summary.csv'),index=False)
                continue
            keys=[k for k in ['reference_source','target_source','subset_kind','support_rule','outcome','method','log_count_step'] if k in frame]
            numeric=[k for k in frame.select_dtypes(include='number').columns if k not in keys and k not in ('rep','effect')]
            agg=frame.groupby(keys,dropna=False)[numeric].agg(['mean','max','count'])
            agg.columns=['__'.join(c) for c in agg.columns]
            agg.reset_index().to_csv(final/(name+'_summary.csv'),index=False)
            if 'status' in frame:
                frame.groupby(keys+['status'],dropna=False).size().rename('records').reset_index().to_csv(
                    final/(name+'_status.csv'),index=False)
        make_plots(final,summary)
    clinical=output/'clinical'
    clinical_complete=bool(checkpoint(clinical,cfg['signature']))
    require(mode not in ('all','clinical') or cfg['settings']['software_fixture'] or clinical_complete,
            'Clinical phase incomplete')
    require(simulation_complete or clinical_complete,'No completed analysis module to export')
    if clinical_complete:
        for name in ('clinical_fraction_scale.csv','method_stability.csv','design_diagnostics.csv','sensitivity.csv',
                     'approximation.csv','invariance.csv','composition_predictions.csv','rna_predictions.csv',
                     'fraction_yields.csv','analysis_cohort.csv','cohort_exclusions.csv','protocol.json'):
            shutil.copy2(clinical/name,final/('clinical_'+name if not name.startswith('clinical_') else name))
    for name in ('reference_marker_features.csv','matching_pairs.csv','marker_matching_balance.csv'):
        if (output/'matching'/name).exists():
            shutil.copy2(output/'matching'/name,final/name)
    for name in ('protocol.json','runtime.json','selected_markers.json','R_versions.json'):
        if (output/name).exists():
            shutil.copy2(output/name,final/name)
    shutil.copy2(Path(__file__).with_name('RESULTS_GUIDE.md'),final/'READ_RESULTS.md')
    save_json(final/'analysis_scope.json',dict(simulation_complete=simulation_complete,
        clinical_complete=clinical_complete,software_fixture=cfg['settings']['software_fixture'],
        cross_engine=['dtangle','linear NNLS'],matched_sets=cfg['settings']['random_subsets'],
        not_performed=['MuSiC/Bisque native analyses','new cohort','clinical truth-based FDR',
                       'new drift directions','real-bulk spike-in'],
        source_input_boundary='exact verified prepared banks required; no claim of a public raw-to-bank regeneration',
        rates='simulation error rates conditional on feasible intervention; infeasibility reported separately',
        plots='generated exclusively from this run; no embedded publication numbers'))
    return dict(simulation_complete=simulation_complete,clinical_complete=clinical_complete,
                completed_batches=len(folders),planned_batches=len(plan))
