"""Fixed-cohort GSE125583 sensitivity application; no truth-based clinical claims."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import bulk_data
from common import require, save_json, digest, read_json


def compare_fraction_scales(y,x,cov,q,p,genes):
    import revision_core as C
    from extension_stats import fit_hc3
    from coupling_math import geometry_change
    frames, fits = {}, {}
    for name, fractions in [('rna_contribution',q),('nuclear_fraction',p)]:
        fit=fit_hc3(y,x,fractions,covariates=cov)
        valid=fit['basic'] & fit['good']
        pv=np.ones(len(genes))
        pv[valid]=2*stats.t.sf(abs(fit['beta'][valid]/fit['se'][valid]),fit['df'])
        frames[name]=pd.DataFrame(dict(gene_id=genes,fraction_scale=name,
            beta=np.where(valid,fit['beta'],np.nan),se=np.where(valid,fit['se'],np.nan),
            p=pv,q=C.M.bh(pv),available=valid))
        fits[name]=fit
    a,b=frames['rna_contribution'],frames['nuclear_fraction']
    valid=a.available & b.available
    sa=a.available & a.q.le(.05)
    sb=b.available & b.q.le(.05)
    row=dict(planned_genes=len(genes),rna_available=int(a.available.sum()),
        nuclear_available=int(b.available.sum()),rna_discoveries=int(sa.sum()),nuclear_discoveries=int(sb.sum()),
        gained_with_nuclear=int((sb&~sa).sum()),lost_with_nuclear=int((sa&~sb).sum()),
        threshold_crossing_fraction=float((sa^sb).mean()),
        median_absolute_beta_difference=float(abs(b.beta[valid]-a.beta[valid]).median()),
        median_nuclear_over_rna_se=float((b.se[valid]/a.se[valid]).median()),
        interpretation='fixed outcome and test family; no clinical effect truth')
    try:
        row.update(geometry_change(x,q,p,cov))
    except (ValueError,np.linalg.LinAlgError) as exc:
        row['geometry_status']='UNAVAILABLE:'+str(exc)
    return row,pd.concat(frames.values(),ignore_index=True)


def analyze(args,jobs,r,cfg):
    import revision_core as C
    import revision_driver as D
    from extension_stats import fit_hc3
    from experiment import METHODS, INVARIANT, read_predictions, evaluate_probe
    from support import checkpoint, seal
    root=args.output/'clinical'
    root.mkdir(exist_ok=True)
    if checkpoint(root,cfg['signature']):
        return
    raw,metadata,mapping,source=bulk_data.prepare(args,root/'inputs')
    meta,x,cov,excluded=bulk_data.cohort(metadata)
    require(len(meta)==288 and int((x==1).sum())==218,'Clinical sample set differs from the manuscript')
    raw=raw.loc[:,meta.sample_id]
    unique=mapping[mapping.one_to_one].set_index('ensembl')
    references=[]
    all_markers=set()
    for job in (j for j in jobs if 'target_source' not in j):
        markers=job['panels']['balanced']
        genes=job['source']['genes'][markers]
        all_markers.update(mapping.loc[mapping.ensembl.isin(genes),'entrez'].tolist())
        owner=C.E.official_markers(job).cell_type.to_numpy(int)
        selected,entrez,classes=[],[],[]
        for i,gene in enumerate(genes):
            eid=unique.loc[gene,'entrez'] if gene in unique.index else ''
            if eid:
                all_markers.add(eid)
            if eid and eid in raw.index:
                selected.append(markers[i]); entrez.append(eid); classes.append(int(owner[i]))
        require(min(classes.count(k) for k in range(1,8))>=10,'Insufficient frozen clinical markers')
        ref=job['base']
        A=ref['A'][selected]
        sums=ref['counts'][selected].sum(0)
        full_by_class=sums/A.sum(0)
        anchor=float(sums.sum()/full_by_class.sum())
        references.append(dict(job=job,name=job['source']['label'],A=A,anchor=anchor,
                               marker_rows=raw.index.get_indexer(entrez),owners=classes,
                               full_mass=float(full_by_class.sum()),full_yield=ref['pooled_yield']))
    total=raw.sum(0).to_numpy(float)
    eligible=((raw.to_numpy()/total*1e6>=1).sum(1)>=np.ceil(len(meta)/2)) & ~raw.index.isin(all_markers)
    test_rows=np.flatnonzero(eligible)
    genes=raw.index[eligible].to_numpy()
    require(len(genes)==16161,'Clinical test family differs from manuscript; do not silently refilter')
    h=cfg['settings']['derivative_step']
    steps=[0.,-h,h,-h/2,h/2,*cfg['settings']['probe_steps']]
    require(len(references)==2,'Expected LIBD and SEAAD clinical references')
    association,designs,diagnostics,changes,predictions,checks=[],[],[],[],[],[]
    scale_rows,scale_associations,rna_predictions,yields=[],[],[],[]
    y=None
    for ref in references:
        work=root/ref['name']
        work.mkdir(exist_ok=True)
        # Translate the clinical excluded gene set to reference IDs before fixing its mass anchor.
        excluded_ens=set(mapping.loc[mapping.entrez.isin(genes),'ensembl'])
        source_genes=ref['job']['source']['genes']
        reference_excluded=np.flatnonzero(np.isin(source_genes,list(excluded_ens)))
        non_anchor=float(1-ref['job']['base']['counts'][reference_excluded].sum()/ref['full_mass'])
        require(0<non_anchor<=1,'Invalid clinical reference exclusion anchor')
        C.write_bin(work/'raw.bin',raw.to_numpy(float))
        C.write_bin(work/'A.bin',ref['A'])
        C.write_bin(work/'total.bin',total)
        axes=dict(normalization_genes=len(raw),marker_rows=(ref['marker_rows']+1).tolist(),
                  marker_owner=ref['owners'],test_rows=(test_rows+1).tolist())
        entry=dict(key='clinical',n=len(meta),raw='raw.bin',A='A.bin',total='total.bin',
                   group=x.tolist(),marker_anchor=ref['anchor'],non_test_anchor=non_anchor,steps=steps)
        C.write_json(work/'jobs.json',dict(axes=axes,entries=[entry]))
        D.run_r([r['executable'],'--vanilla',str(Path(__file__).with_name('bridge.R')),str(work),r['library']],
                r['library'],work/'bridge.log',args.r_timeout)
        pp,rr=read_predictions(work,{('clinical',method,t) for method in METHODS for t in steps})
        denom=pp['clinical','tmm',0.][1]
        require(np.isfinite(denom).all() and (denom>0).all(),'Invalid common clinical outcome denominator')
        response=np.log2(raw.iloc[test_rows].to_numpy(float).T/denom[:,None]*1e6+1)
        if y is None:
            y=response
        else:
            require(np.array_equal(y,response),'Clinical outcome differs between references')
        marker=raw.iloc[ref['marker_rows']].to_numpy(float).T
        ref['marker_count']=len(ref['marker_rows'])
        ref['non_test_anchor']=non_anchor
        for method in METHODS:
            q,d,reason,scale=pp['clinical',method,0.]
            yield_=ref['full_yield']*(ref['A'].sum(0) if method=='marker_matched' else 1.)
            p=C.E.convert_available_rows(q,yield_)
            summary,scale_table=compare_fraction_scales(y,x,cov,q,p,genes)
            q_definition='marker_RNA_contribution' if method=='marker_matched' else 'full_reference_RNA_contribution'
            scale_rows.append(dict(reference=ref['name'],method=method,q_definition=q_definition,**summary))
            scale_table['reference']=ref['name']; scale_table['method']=method
            scale_associations.append(scale_table)
            qtable=pd.DataFrame(q,columns=C.M.CT)
            qtable.insert(0,'sample_id',meta.sample_id.tolist())
            qtable['reference']=ref['name']; qtable['method']=method; qtable['fraction_definition']=q_definition
            rna_predictions.append(qtable)
            yields.extend(dict(reference=ref['name'],method=method,cell_class=name,yield_value=float(value))
                          for name,value in zip(C.M.CT,yield_))
            fit=fit_hc3(y,x,p,covariates=cov)
            valid=fit['basic'] & fit['good']
            pv=np.ones(len(genes))
            pv[valid]=2*stats.t.sf(abs(fit['beta'][valid]/fit['se'][valid]),fit['df'])
            association.append(pd.DataFrame(dict(reference=ref['name'],method=method,gene_id=genes,
                beta=np.where(valid,fit['beta'],np.nan),se=np.where(valid,fit['se'],np.nan),
                p=pv,q=C.M.bh(pv),available=valid)))
            designs.append(dict(reference=ref['name'],method=method,available_genes=int(valid.sum()),
                planned_genes=len(genes),prediction_reason=reason,
                **{k:fit[k] for k in ('reason','df','composition_rank','exposure_r2','exposure_residual_ss','leverage')}))
            table=pd.DataFrame(p,columns=C.M.CT)
            table.insert(0,'sample_id',meta.sample_id.tolist())
            table['reference']=ref['name']; table['method']=method
            predictions.append(table)
            probes,denominators={},{0.:d}
            for step in steps:
                if step==0.:
                    continue
                qt,dt,_,_=pp['clinical',method,step]
                probes[step]=C.E.convert_available_rows(qt,yield_)
                denominators[step]=dt
                if method in INVARIANT:
                    finite=np.isfinite(q).all() and np.isfinite(qt).all()
                    delta=float(np.max(abs(q-qt))) if finite else np.nan
                    require(np.array_equal(np.isfinite(q),np.isfinite(qt)) and (not finite or delta<1e-10),
                            'Clinical complete-exclusion invariance failed')
                    checks.append(dict(reference=ref['name'],method=method,step=step,prediction_finite=finite,max_change=delta))
            score,approx=evaluate_probe(dict(group=x,y=y),p,probes,denominators,marker,scale,h,fit,cov)
            diagnostics.append(dict(reference=ref['name'],method=method,**score))
            changes.extend(dict(reference=ref['name'],method=method,**a) for a in approx)
        rr.to_csv(work/'prediction_manifest.csv',index=False)
    assoc=pd.concat(association,ignore_index=True)
    stability=[]
    for reference,frame in assoc.groupby('reference'):
        baseline=frame[frame.method.eq('tmm')].set_index('gene_id')
        for method,g in frame.groupby('method'):
            other=g.set_index('gene_id').reindex(baseline.index)
            a=baseline.available & baseline.q.le(.05)
            b=other.available & other.q.le(.05)
            valid=baseline.available & other.available
            stability.append(dict(reference=reference,method=method,planned_genes=len(genes),
                discoveries=int(b.sum()),gained=int((b&~a).sum()),lost=int((a&~b).sum()),
                threshold_crossings=int((a^b).sum()),threshold_crossing_fraction=float((a^b).mean()),
                median_absolute_beta_change=float(abs(baseline.beta[valid]-other.beta[valid]).median()),
                jointly_significant_opposite_sign=int((a&b&(np.sign(baseline.beta)!=np.sign(other.beta))).sum())))
    tables=dict(disease_associations=assoc,method_stability=pd.DataFrame(stability),
        design_diagnostics=pd.DataFrame(designs),sensitivity=pd.DataFrame(diagnostics),
        approximation=pd.DataFrame(changes),invariance=pd.DataFrame(checks),
        composition_predictions=pd.concat(predictions),analysis_cohort=meta,cohort_exclusions=excluded,
        clinical_fraction_scale=pd.DataFrame(scale_rows),rna_predictions=pd.concat(rna_predictions),
        fraction_yields=pd.DataFrame(yields),fraction_scale_associations=pd.concat(scale_associations))
    files=[]
    for name,table in tables.items():
        filename=name+'.csv.gz' if name in ('disease_associations','fraction_scale_associations') else name+'.csv'
        table.to_csv(root/filename,index=False,compression='gzip' if filename.endswith('.gz') else None)
        files.append(filename)
    provenance=dict(source=source,sample_ids=meta.sample_id.tolist(),test_genes=genes.tolist(),
        outcome='Fixed original all-gene TMM log2 CPM + 1; NEVER changed by probes',
        covariates=['standardized_age','source_recorded_sex'],
        perturbation='All tested non-marker counts multiplied by exp(t*AD); sensitivity only',
        exclusion_reference_anchor=[dict(reference=r['name'],anchor=r['non_test_anchor'],markers=r['marker_count']) for r in references],
        true_fdr_evaluated=False,diagnostic_validity='sensitivity, not accuracy or causal effect verification')
    save_json(root/'protocol.json',provenance)
    seal(root,[*files,'protocol.json'],cfg['signature'],samples=len(meta),genes=len(genes))
