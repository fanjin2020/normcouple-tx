"""Metadata-first bulk repair. Processed expression never becomes raw counts."""
from pathlib import Path
from types import SimpleNamespace
import json
import shutil
import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
import real_bulk as B
from bulk_access import prepare_manifest


def metadata_first(args, output, r_info):
    """Create cohort_exclusions before validating raw counts, even on a new run."""
    output.mkdir(parents=True,exist_ok=True)
    rda=B.acquire(B.BULK_NAME,args,output,B.BULK_SHA256)
    snmeta=B.acquire(B.META_NAME,args,output,B.S.PUBLIC_SHA256[B.META_NAME])
    D.run_r([r_info['executable'],'--vanilla',str(Path(B.__file__).parent/'export_bulk.R'),str(rda),str(output)],
            r_info['library'],output/'export.log',timeout=600)
    cells=pd.read_csv(snmeta)
    meta,paired=B.audit_samples(pd.read_csv(output/'bulk_metadata.csv'),cells)
    m,x,cov,names,excluded=B.design(meta)
    C.write_csv(output/'cohort_exclusions.csv',excluded)
    C.write_csv(output/'bulk_snRNA_pairing_audit.csv',paired)
    C.write_json(output/'public_inputs.json',[dict(name=p.name,sha256=C.sha(p),url=B.BASE+p.name) for p in (rda,snmeta)])
    processed=pd.read_csv(output/'bulk_processed_expression.csv',index_col=0)
    C.require(set(m.SampleID)<=set(processed.columns),'Processed sample labels do not match metadata')
    table,diag=B.association_table(processed.loc[:,m.SampleID].to_numpy().T,processed.index,x,cov,None,
                                   'deposited_processed_expression_ONLY')
    table.to_csv(output/'processed_disease_associations_DESCRIPTIVE.csv.gz',index=False,compression='gzip')
    C.write_json(output/'processed_analysis_scope.json',dict(status='COMPLETE_DESCRIPTIVE_ONLY',
        diagnosis_groups=m.Diagnosis.value_counts().to_dict(),covariates=names,diagnostics=diag,
        measurement='Deposited normExpr.reg; previous preprocessing retained',
        normalization_comparison=False,true_fdr_evaluated=False,raw_counts_reconstructed=False,
        warning='Association on deposited processed expression, not validation of a normalization intervention.'))
    return meta,cells,paired,table


def marker_descriptive(table,jobs,output):
    """Report the deposited-expression scale; no fold-change reconstruction."""
    t=table.copy()
    t['gene_id']=t.gene_id.astype(str).str.replace(r'\.\d+$','',regex=True)
    C.require(t.gene_id.is_unique,'Processed gene IDs become duplicate after version removal')
    parts=[]; counts=[]
    for job in jobs:
        if 'target_source' in job: continue
        ids=job['source']['genes'][job['panels']['balanced']]
        part=pd.DataFrame({'gene_id':ids}).merge(t,on='gene_id',how='left',validate='one_to_one',indicator=True)
        part['reference']=job['source']['label']
        part['measurement_scope']='Coefficient on deposited processed-expression scale; not a raw marker-stress dose'
        parts.append(part)
        counts.append(dict(reference=job['source']['label'],planned_markers=len(ids),
                           matched=int(part._merge.eq('both').sum())))
    if parts:
        C.write_csv(output/'marker_disease_associations_DESCRIPTIVE.csv',pd.concat(parts,ignore_index=True))
    C.write_json(output/'marker_processed_overlap.json',counts)


def run(args,output,jobs,r_info):
    application=output/'real_bulk'
    try:
        meta,cells,paired,processed=metadata_first(args,application,r_info)
        marker_descriptive(processed,jobs,application)
        options=SimpleNamespace(**vars(args))
        options.run_root=args.project_root.parent/'paper3_bulk_repro_results_v1'
        options.output=output
        # The repaired metadata stage makes this independent of old compact packs.
        manifest,status=prepare_manifest(options,output,output/'raw_count_access')
        if manifest is None:
            result=dict(status='PARTIAL_DESCRIPTIVE_ONLY_RAW_COUNTS_REQUIRED',raw_count_status=status,
                        descriptive_application='COMPLETE',normalization_comparison='NOT_RUN',
                        true_fdr_evaluated=False,simulations_need_rerun=False,
                        action='Supply authorized raw gene counts plus provenance, then rerun --only bulk into the same output.')
        else:
            m=json.loads(manifest.read_text(encoding='utf-8'))
            result=B.raw_application(options,application,jobs,meta,cells,paired,r_info,m)
            associations=pd.read_csv(application/'disease_associations.csv')
            C.require(associations.method.nunique()==11,'Incomplete raw-bulk comparison')
            C.require(associations.groupby('method').size().nunique()==1,'Different bulk test-gene families')
            result.update(true_fdr_evaluated=False,raw_count_provenance_verified=True,
                          all_planned_methods_evaluated=True,descriptive_application='COMPLETE')
        C.write_json(application/'STATUS.json',result)
        return result
    except Exception as exc:
        result=dict(status='BULK_STAGE_FAILED',error=type(exc).__name__+': '+str(exc),
                    normalization_comparison='NOT_CONFIRMED_COMPLETE',simulations_need_rerun=False)
        C.write_json(application/'STATUS.json',result)
        return result


def collect(output,final):
    dest=final/'real_bulk'; dest.mkdir(parents=True,exist_ok=True)
    excluded={'bulk_processed_expression.csv','bulk_counts_for_R.csv','bulk_raw_counts.csv'}
    for p in (output/'real_bulk').glob('*'):
        if p.is_file() and p.name not in excluded and p.suffix in ('.csv','.json','.txt','.gz'):
            shutil.copy2(p,dest/p.name)
    for p in (output/'raw_count_access').glob('*.json'):
        shutil.copy2(p,dest/p.name)
