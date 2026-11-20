"""Publication table for fixed-response dtangle/TMM dose and gene-level probes."""
from pathlib import Path

import numpy as np
import pandas as pd

from .io import read_json, sha, verify_receipt

REFS = ('LIBD', 'SEAAD')
COHORTS = {'GSE125583': (16161, 218, 70), 'GSE53697': (14745, 9, 8)}
SETS = tuple(k+'_'+r for r in ('low','middle','high') for k in ('single','five')) + ('test_family',)
TITLE = 'Realized count dose and gene-level sensitivity in two brain cohorts'
NOTE = ('GSE125583: 218 cases/70 controls; GSE53697: 9/8. Input change is imposed in cases. Count share and raw library changes are sample medians; effective library size equals raw size times the TMM factor. Each reference uses its own dtangle/TMM baseline. Own and family crossings count selected and all tested genes changing BH status (q <= 0.05). Signed beta changes apply only to singletons. LIBD fits in GSE53697 lack marker coverage (Table S18); their blank results are unavailable, not zero.')
COLUMNS = ['cohort','probe','selected_genes','selected_ids','AD_input_change_pct','AD_count_share_pct','AD_library_change_pct',
           'AD_effective_q25_pct','AD_effective_median_pct','AD_effective_q75_pct','control_library_change_pct',
           'control_effective_q25_pct','control_effective_median_pct','control_effective_q75_pct'] + [
           ref+'_'+stat for ref in REFS for stat in ('singleton_beta_change','own_BH_crossings','family_BH_crossings')]


def _index(frame, keys):
    return frame.set_index(keys,verify_integrity=True).sort_index()


def marker_roster(frame):
    keys=['cohort','reference','cell_class','gene_id']
    if frame[keys].isna().any().any() or frame.duplicated(['cohort','reference','gene_id']).any():
        raise ValueError('Ambiguous clinical marker assignments')
    if not frame.marker_probe.isin([True,False]).all():
        raise ValueError('Invalid clinical marker selection indicator')
    frame=frame.sort_values(keys)
    g=frame.groupby(keys[:3]).gene_id
    selected=g.cumcount() < (g.transform('size')+9)//10
    if not np.array_equal(selected.to_numpy(),frame.marker_probe.to_numpy()):
        raise ValueError('Clinical marker selection differs from the prespecified within-class rule')
    return frame.loc[selected,keys].reset_index(drop=True)


def build(directory, clinical_contrasts):
    directory=Path(directory)
    files=read_json(directory/'receipt.json')['files']
    required={'protocol.json','dose.csv','probe_effects.csv','selected_genes.csv','self_effects.csv','markers.csv'}
    if not required <= files.keys():
        raise ValueError('Incomplete receipt for clinical dose and single-gene evidence')
    verify_receipt(directory,files)
    expected=pd.read_csv(Path(__file__).resolve().parents[1]/'publication/input_definitions/clinical_probe_markers.csv')
    pd.testing.assert_frame_equal(marker_roster(pd.read_csv(directory/'markers.csv')),expected)
    protocol=read_json(directory/'protocol.json')
    def load(name, keys):
        f=pd.read_csv(directory/(name+'.csv'),keep_default_na=False,na_values=['','NaN','nan'],float_precision='round_trip')
        if 'method' in f: f=f.loc[f.method.eq('dtangle_tmm')]
        return _index(f,['cohort','reference']+keys)
    dose=load('dose',['probe','group']); effects=load('probe_effects',['probe'])
    detail=load('self_effects',['probe','gene_id']); selected=load('selected_genes',['selection','gene_id'])
    anchors=_index(clinical_contrasts.query("method=='dtangle_tmm' and probe=='test_counts'"),['cohort','reference'])
    result=[]; identities={}
    for cohort,(family,cases,controls) in COHORTS.items():
        inputs=[protocol['clinical_inputs'][cohort+'_'+r] for r in REFS]
        keys=('raw','cohort_sha256','outcome_sha256')
        if any(inputs[0][k]!=inputs[1][k] for k in keys):
            raise ValueError('References do not share raw samples and outcomes: '+cohort)
        identities[cohort]={k:inputs[0][k] for k in keys}
        for group in SETS:
            count=family if group=='test_family' else (1 if group.startswith('single_') else 5)
            ids=None
            if count!=family:
                ids=selected.loc[(cohort,REFS[0],group)].index
                if len(ids)!=count or not ids.equals(selected.loc[(cohort,REFS[1],group)].index):
                    raise ValueError('Selected genes differ by reference: '+cohort+' '+group)
            labels='all tested genes' if ids is None else ';'.join(ids)
            for direction,step in [('down',-.1),('up',.1)]:
                probe=group+'_'+direction
                rows=dose.loc[(cohort,REFS[0],probe)]
                pd.testing.assert_frame_equal(rows,dose.loc[(cohort,REFS[1],probe)],check_exact=False,rtol=1e-11,atol=1e-13)
                ad,ctl=rows.loc['AD'],rows.loc['control']
                if (list(rows.index)!=['AD','control'] or ad.samples!=cases or ctl.samples!=controls or
                        not rows.genes.eq(count).all() or ctl.library_relative_change_median!=0 or
                        not np.isclose(ad.library_relative_change_median,ad.tested_input_share_median*np.expm1(step),atol=1e-13,rtol=1e-10)):
                    raise ValueError('Invalid library dose or study population: '+cohort+' '+probe)
                row=[cohort,probe,count,labels,100*np.expm1(step),100*ad.tested_input_share_median,100*ad.library_relative_change_median]
                for i,g in enumerate([ad,ctl]):
                    if i: row.append(100*g.library_relative_change_median)
                    q=g[['effective_relative_change_'+s for s in ('q25','median','q75')]].to_numpy(float)
                    if not np.isfinite(q).all() or (np.diff(q)<-1e-13).any():
                        raise ValueError('Invalid effective-denominator quartiles: '+probe)
                    row.extend(100*q)
                for ref in REFS:
                    f=effects.loc[(cohort,ref,probe)]; anchor=anchors.loc[(cohort,ref)]
                    if (f.available!=anchor.available or f.tested_family!=family or not np.isclose(f.log_step,step)):
                        raise ValueError('Probe does not match its clinical specification: '+probe+' '+ref)
                    if not f.available:
                        if f.reason!='insufficient_marker_coverage' or pd.notna(f.bh_crossings) or pd.notna(f.own_beta_rms):
                            raise ValueError('Invalid unavailable clinical record: '+probe+' '+ref)
                        row.extend([np.nan]*3); continue
                    if (f.paired_tests!=family or f.own_paired_genes!=count or not 0<=f.bh_crossings<=family or
                            f.bh_crossings!=int(f.bh_crossings)):
                        raise ValueError('Incomplete paired tests: '+probe+' '+ref)
                    delta=np.nan
                    if count==family:
                        own=int(f.bh_crossings)
                        if step>0 and (anchor.family_size!=family or f.bh_crossings!=anchor.bh_crossings or
                                not np.isclose(f.own_beta_rms,anchor.beta_rms,atol=1e-12,rtol=1e-9)):
                            raise ValueError('Full-family probe differs from clinical comparison: '+ref)
                    else:
                        s=detail.loc[(cohort,ref,probe)]
                        if (not ids.equals(s.index) or not s.available.eq(True).all() or not s.bh_crossed.isin([True,False]).all() or
                                not np.isfinite(s.beta_change).all() or not np.isclose(np.sqrt(np.mean(s.beta_change**2)),f.own_beta_rms,atol=1e-14,rtol=1e-9)):
                            raise ValueError('Mismatched own-gene effects: '+probe+' '+ref)
                        own=int(s.bh_crossed.sum())
                        if count==1: delta=float(s.beta_change.iloc[0])
                    if own>f.bh_crossings: raise ValueError('Own crossings exceed family crossings')
                    row.extend([delta,own,int(f.bh_crossings)])
                result.append(row)
    return pd.DataFrame(result,columns=COLUMNS),dict(receipt_sha256=sha(directory/'receipt.json'),
        producer_framework_sha256=protocol['software']['manifest'],clinical_inputs=identities)
