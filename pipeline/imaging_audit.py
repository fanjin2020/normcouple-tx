"""Frozen dtangle specifications versus an orthogonal block-level RNA imaging assay."""
import os
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from .download import download
from .io import ROOT, atomic_json, read_json, sha, verify_receipt

UPSTREAM='LieberInstitute/Human_DLPFC_Deconvolution'
COMMIT='96131ef436fdaa9a9b2129983ac7bd6868204fc4'
FILES={
    'DLPFC_bulk_counts.txt':('processed-data/08_bulk_deconvolution/07_deconvolution_CIBERSORTx_prep/',
        '290861e3a23fd366fea2bb963c29657cd596f2481a9981bac9669f131420a8c3'),
    'sample_info.csv':('processed-data/00_data_prep/','84f27820c79aff05788974553bc3284f769ffe36dc501a7f143200f3cb5cfb10'),
    'HALO_cell_type_proportions.csv':('processed-data/03_HALO/08_explore_proportions/',
        '4d3c4bb8a310a39b1cdd92e7f4535fed5bdfb969d25cfd6361e720cf7dd6d74a')}
CLASSES=('Astro','EndoMural','Excit','Inhib','Micro','OligoOPC')
TABLE_TITLE = 'Donor-weighted errors against tissue-block imaging proportions'
TABLE_NOTE = ('Huuki-Myers et al. prefrontal cortex: 14 bulk libraries, seven matched blocks, five held-out donors. '
    'Only reference-derived transfer panels met the ten-marker minimum in every class. Technical libraries and '
    'eligible imaging sections were averaged within block. Errors were averaged within donor and then equally '
    'across donors with measurements for that class. MAE and RMSE are in fraction units; delta is MRN minus TMM. '
    'Oligodendrocytes and OPCs are combined. These measurements evaluate proportions, not disease associations.')


def paired_table(blocks, audit):
    """Summarize measured fractions on identical block/class support for both dtangle fits."""
    group=['reference','panel','cell_type']
    keys=group+['donor','block']
    methods={'dtangle_marker','dtangle_tmm'}
    fits={(p['reference'],p['panel']) for p in audit['panels'] if p['status']=='estimable'}
    if (blocks.empty or blocks[keys+['method']].isna().any().any() or
            blocks.duplicated(keys+['method']).any() or set(blocks.method)!=methods or
            set(map(tuple,blocks[['reference','panel']].drop_duplicates().to_numpy()))!=fits):
        raise ValueError('Incomplete or duplicate imaging fits')
    if (not blocks[['measured','prediction']].apply(lambda s:s.between(0,1).all()).all() or
            blocks[['libraries','sections']].isna().any().any() or
            (blocks[['libraries','sections']]<1).any().any()):
        raise ValueError('Invalid imaging observations')
    pairing=blocks.groupby(keys)
    if (not pairing.size().eq(2).all() or
            not blocks.groupby(['block','cell_type'])[['measured','sections']].nunique().eq(1).all().all() or
            not blocks.groupby('block')[['donor','libraries']].nunique().eq(1).all().all()):
        raise ValueError('Imaging methods do not share the same tissue measurements')
    support=blocks[['donor','block','cell_type']].drop_duplicates()
    if (not blocks.groupby(['reference','panel']).size().eq(2*len(support)).all() or
            set(support.cell_type)!=set(CLASSES) or
            (support.donor.nunique(),support.block.nunique(),len(support),
             blocks.drop_duplicates('block').libraries.sum()) !=
            tuple(audit[k] for k in ('matched_donors','matched_blocks','measured_block_class_pairs','matched_libraries'))):
        raise ValueError('Imaging block counts disagree with the pairing audit')
    error=blocks.prediction-blocks.measured
    donor=blocks.assign(mae=error.abs(),mse=error**2).groupby(group+['donor','method']).agg(
        mae=('mae','mean'),mse=('mse','mean'),blocks=('block','size')).reset_index()
    by=donor.groupby(group+['method']).agg(mae=('mae','mean'),mse=('mse','mean'),
        donors=('donor','size'),blocks=('blocks','sum'))
    t=by.xs('dtangle_tmm',level='method'); m=by.xs('dtangle_marker',level='method')
    if not t[['donors','blocks']].equals(m[['donors','blocks']]):
        raise ValueError('Imaging donor support differs by method')
    return t[['donors','blocks']].assign(tmm_mae=t.mae,mrn_mae=m.mae,
        mrn_minus_tmm_mae=m.mae-t.mae,tmm_rmse=np.sqrt(t.mse),mrn_rmse=np.sqrt(m.mse)).reset_index()


def publication_table(directory):
    directory=Path(directory)
    required={'imaging_blocks.csv','imaging_match.json','imaging_error.csv','protocol.json'}
    manifest=read_json(directory/'receipt.json')['files']
    if not required<=manifest.keys():raise ValueError('Missing imaging evidence receipt')
    verify_receipt(directory,{name:manifest[name] for name in required})
    audit=read_json(directory/'imaging_match.json')
    sources=read_json(directory/'protocol.json')['imaging_source_sha256']
    if sources!={name:h for name,(_,h) in FILES.items()}:
        raise ValueError('Imaging inputs differ from the pinned public study')
    table=paired_table(pd.read_csv(directory/'imaging_blocks.csv',float_precision='round_trip'),audit)
    original=pd.read_csv(directory/'imaging_error.csv',float_precision='round_trip')
    keys=['reference','panel','cell_type']
    for method,short in [('dtangle_marker','mrn'),('dtangle_tmm','tmm')]:
        reported=original.loc[original.method.eq(method)].set_index(keys,verify_integrity=True).sort_index()
        rebuilt=table.set_index(keys,verify_integrity=True).sort_index()
        if (not reported.index.equals(rebuilt.index) or not np.array_equal(reported[['donors','blocks']],rebuilt[['donors','blocks']]) or
                not np.allclose(reported[['mae','rmse']],rebuilt[[short+'_mae',short+'_rmse']],rtol=1e-11,atol=1e-14)):
            raise ValueError('Imaging aggregate does not agree with the paired block records')
    return table,dict(source_sha256=sources,matched={k:audit[k] for k in
        ('matched_libraries','matched_blocks','matched_donors','measured_block_class_pairs')})


def sources(cache,allow):
    cache=Path(cache)
    for name,(prefix,expected) in FILES.items():
        path=cache/name
        if path.exists():
            if sha(path)!=expected: raise ValueError('Imaging cache checksum differs: '+str(path))
            continue
        if not allow: raise ValueError('Missing pinned public imaging source: '+str(path))
        for mirror,base in enumerate((f'https://raw.githubusercontent.com/{UPSTREAM}/{COMMIT}/',
                                      f'https://cdn.jsdelivr.net/gh/{UPSTREAM}@{COMMIT}/')):
            try:
                tmp=cache/('mirror'+str(mirror))/name
                download(base+prefix+name,tmp,expected,attempts=2,timeout=45)
                os.replace(tmp,path)
                break
            except OSError:
                if mirror: raise
    return {n:h for n,(_,h) in FILES.items()}


def matched_inputs(cache):
    counts=pd.read_csv(cache/'DLPFC_bulk_counts.txt',sep='\t',index_col=0)
    indices=counts.index.astype(str)
    if not indices.is_unique or not indices.str.fullmatch(r'ENSG\d+').all():
        raise ValueError('Nonunique or unexpected count gene identifiers')
    meta=pd.read_csv(cache/'sample_info.csv')
    if not meta.sample_id.is_unique: raise ValueError('Ambiguous author block mapping')
    units=pd.read_csv(ROOT/'publication/input_definitions/reference_target_units.csv')
    held=set(units[(units.reference_source=='LIBD')&(units.target_source=='LIBD')].target_donor)
    if len(held)!=6: raise ValueError('Frozen held-out donor set differs')
    libs=pd.Series(counts.columns,name='library').to_frame()
    labels=libs.library.str.extract(r'^([^_]+)_(Br\d+)_(Ant|Mid|Post)_(Bulk|Cyto|Nuc)$')
    if labels.isna().any().any(): raise ValueError('Unrecognized upstream library layout')
    libs['block']=labels[1]+'_'+labels[2].str.lower()
    libs['donor']=labels[1]; libs['extraction']=labels[3]
    check=libs.merge(meta[['sample_id','subjects']],left_on='block',right_on='sample_id',how='left',validate='many_to_one')
    if (check.subjects.notna() & check.subjects.ne(check.donor)).any(): raise ValueError('Inconsistent author donor crosswalk')
    eligible=(check.extraction=='Bulk') & check.donor.isin(held) & check.subjects.notna()
    libraries=check.loc[eligible,['library','donor','block']].copy()
    image=pd.read_csv(cache/'HALO_cell_type_proportions.csv')
    if image.duplicated(['SAMPLE_ID','cell_type']).any(): raise ValueError('Duplicate imaging section and class')
    keep=image.Confidence.isin(['High','OK']) & image.cell_type.isin(CLASSES)
    section=image.loc[keep,['Sample','SAMPLE_ID','cell_type','prop']].rename(columns={'Sample':'block'})
    if not np.isfinite(section.prop).all() or not section.prop.between(0,1).all(): raise ValueError('Invalid measured proportions')
    truth=section.groupby(['block','cell_type'],as_index=False).agg(measured=('prop','mean'),sections=('SAMPLE_ID','nunique'))
    truth=truth[truth.block.isin(libraries.block)].copy()
    if set(truth.cell_type)!=set(CLASSES) or not len(libraries): raise ValueError('Incomplete matched imaging strata')
    audit=dict(input_bulk_libraries=len(libs),total_extraction=int(check.extraction.eq('Bulk').sum()),
        heldout_total_no_exact_block=int((check.extraction.eq('Bulk') & check.donor.isin(held)&check.subjects.isna()).sum()),
        eligible_bulk_libraries=len(libraries),eligible_blocks=libraries.block.nunique(),
        matched_blocks=truth.block.nunique(),matched_libraries=int(libraries.block.isin(truth.block).sum()),
        matched_donors=libraries[libraries.block.isin(truth.block)].donor.nunique(),
        measured_block_class_pairs=len(truth),expected_heldout_donors=sorted(held))
    return counts.loc[:,libraries.library], libraries, truth, audit


def score(fractions,libraries,truth):
    if fractions.duplicated(['sample_id','method']).any(): raise ValueError('Duplicated imaging predictions')
    if set(fractions.method)!={'dtangle_marker','dtangle_tmm'}: raise ValueError('Missing paired imaging methods')
    for _,f in fractions.groupby('method'):
        if set(f.sample_id)!=set(libraries.library): raise ValueError('Library identities differ')
    fractions=fractions.copy()
    p=fractions[['Astro','EndoMural','Excit','Inhib','Micro','Oligo','OPC']].to_numpy(float)
    if not np.isfinite(p).all() or (p<0).any() or not np.allclose(p.sum(1),1,rtol=0,atol=1e-8):
        raise ValueError('Invalid imaging predictions')
    fractions['OligoOPC']=fractions.Oligo+fractions.OPC
    f=fractions.melt(id_vars=['sample_id','method'],value_vars=list(CLASSES),var_name='cell_type',value_name='prediction')
    f=f.merge(libraries,left_on='sample_id',right_on='library',validate='many_to_one')
    block=f.groupby(['donor','block','cell_type','method'],as_index=False).agg(prediction=('prediction','mean'),libraries=('library','nunique'))
    block=block.merge(truth,on=['block','cell_type'],validate='many_to_one')
    if len(block)!=2*len(truth): raise ValueError('Missing common block/class support')
    block['error']=block.prediction-block.measured
    block['abs_error']=block.error.abs(); block['squared_error']=block.error**2
    donor=block.groupby(['cell_type','method','donor'],as_index=False).agg(mae=('abs_error','mean'),mse=('squared_error','mean'),blocks=('block','nunique'))
    result=donor.groupby(['cell_type','method'],as_index=False).agg(mae=('mae','mean'),mse=('mse','mean'),donors=('donor','nunique'),blocks=('blocks','sum'))
    result['rmse']=np.sqrt(result.pop('mse'))
    base=result[result.method=='dtangle_tmm'].set_index('cell_type').mae
    result['paired_mae_change_vs_tmm']=result.mae-result.cell_type.map(base)
    return result,block[['donor','block','cell_type','method','libraries','sections','measured','prediction']]


def reference_receipts(cross,protocol):
    out={}
    for ref in ('LIBD','SEAAD'):
        prefix='references/'+ref+'/'
        expected={n:protocol['reference_hashes'][prefix+n] for n in
                  ('reference.mtx','cells.csv','genes.csv','source_receipt.json')}
        folder=cross/'references'/ref
        verify_receipt(folder,expected)
        old=read_json(folder/'source_receipt.json')
        if old['files']!={k:v for k,v in expected.items() if k!='source_receipt.json'}:
            raise ValueError('Disagreement between reference seals: '+ref)
        out[ref]=expected
    return out


def transfer_panel(matrix,ids,meta,classes,universe,excluded):
    """Train on reference donors and gene names, without bulk values or imaging outcomes."""
    from scipy import sparse
    ids=pd.Index(ids); k=len(classes)
    if not ids.is_unique or len(ids)!=matrix.shape[0] or len(meta)!=matrix.shape[1]:
        raise ValueError('Reference matrix does not match its registered rows and columns')
    if not np.isfinite(matrix.data).all() or np.any(matrix.data<0):
        raise ValueError('Invalid reference counts')
    donors=sorted(meta.donor_id.unique())
    if len(donors)!=4 or set(meta.cell_type)!=set(classes) or set(donors)&set(excluded):
        raise ValueError('Need four reference-only donors and all frozen classes')
    gi=np.asarray([donors.index(d)*k+classes.index(c) for d,c in zip(meta.donor_id,meta.cell_type)])
    nc=np.bincount(gi,minlength=len(donors)*k)
    if min(nc)<1: raise ValueError('Missing donor/class cells in reference')
    h=sparse.csr_matrix((np.ones(len(gi)),(np.arange(len(gi)),gi)),shape=(len(gi),len(nc)))
    pooled=(matrix@h).toarray()
    mass=pooled.sum(0)
    if (mass<=0).any(): raise ValueError('Empty donor/class reference pseudobulk')
    fy=(mass/nc).reshape(len(donors),k).mean(0)
    per=(pooled/mass).reshape(len(ids),len(donors),k)
    avg=per.mean(1)
    lead=avg.argmax(1)
    order=np.sort(avg,axis=1)
    score=np.log2((order[:,-1]*1e6+1)/(order[:,-2]*1e6+1))
    observed=ids.isin(universe)
    rows=[]; owners=[]; available=[]
    for i in range(k):
        ix=np.flatnonzero(observed & (lead==i) & (score>0) & ((per[:,:,i]*1e6>=10).sum(axis=1)>=3))
        available.append(len(ix))
        top=ix[np.lexsort((ids[ix].to_numpy(),-score[ix]))[:25]]
        rows.extend(top); owners.extend([i+1]*len(top))
    panel=avg[rows,:]
    anchor=float(np.dot(fy,panel.sum(0))/fy.sum())
    return dict(ids=ids[rows].tolist(),owner=owners,classes=classes,A=panel,full_yield=fy,anchor=anchor,
                donors=donors,candidates=available)


def compute(cache,clinical,cross,work,rscript,libs):
    from scipy.io import mmread
    counts,libraries,truth,audit=matched_inputs(cache)
    scores=[]; pairs=[]; panel_info=[]; selected=[]
    for ref in ('LIBD','SEAAD'):
        bank=clinical/'GSE125583'/ref; plan=read_json(bank/'plan.json')
        classes=plan['classes']
        marker=np.asarray(plan['genes'])[np.asarray(plan['marker_rows'])-1]
        active=np.flatnonzero(np.isin(marker,counts.index))
        owners=np.asarray(plan['marker_owner'])[active]
        if len(set(classes))!=7 or np.any(owners<1) or np.any(owners>7): raise ValueError('Invalid original marker labels')
        source=cross/'references'/ref
        trained=transfer_panel(mmread(source/'reference.mtx').tocsr(),pd.read_csv(source/'genes.csv').gene_id,
                 pd.read_csv(source/'cells.csv'),classes,counts.index,audit['expected_heldout_donors'])
        refwork=work/ref; refwork.mkdir(parents=True)
        alternate=refwork/'reference_panel'; alternate.mkdir()
        cols=np.asarray(trained['owner'],dtype=int)
        a=trained['A']; full=trained['full_yield']
        atomic_json(alternate/'plan.json',dict(classes=classes,marker_rows=list(range(1,len(cols)+1)),
           marker_owner=cols.tolist(),anchor=trained['anchor'],full_yield=full.tolist(),marker_yield=(full*a.sum(0)).tolist()))
        a.astype('<f8').ravel(order='F').tofile(alternate/'A.bin')
        for label,src,idx,genes,own in [('original_clinical',bank,active,marker[active],owners),
                ('reference_transfer',alternate,np.arange(len(cols)),trained['ids'],cols)]:
            ns=np.bincount(own,minlength=8)[1:]
            if len(ns)!=7: raise ValueError('Unexpected panel class labels')
            supported=bool(min(ns)>=10)
            panel_info.append(dict(reference=ref,panel=label,status='estimable' if supported else 'insufficient_marker_support',
                mapped=dict(zip(classes,map(int,ns))),minimum_per_class=10,reference_donors=trained['donors'],
                eligible_reference_candidates=dict(zip(classes,map(int,trained['candidates']))) if label=='reference_transfer' else None))
            for gid,owner in zip(genes,own):
                selected.append(dict(reference=ref,panel=label,cell_class=classes[owner-1],gene_id=gid))
            if not supported: continue
            folder=refwork/label; folder.mkdir()
            counts.to_numpy(dtype='<f8').ravel(order='F').tofile(folder/'imaging.bin')
            where=counts.index.get_indexer(genes)
            if (where<0).any() or len(set(genes))!=len(genes): raise ValueError('Invalid shared marker identities')
            atomic_json(folder/'imaging_plan.json',dict(genes=len(counts),samples=list(libraries.library),
                 panel_indices=(idx+1).tolist(),count_rows=(where+1).tolist()))
            subprocess.run([rscript,str(ROOT/'pipeline/probe_engine.R'),'imaging',str(src),str(folder),*map(str,libs)],check=True)
            summary,block=score(pd.read_csv(folder/'imaging_fractions.csv'),libraries,truth)
            scores.append(summary.assign(reference=ref,panel=label)); pairs.append(block.assign(reference=ref,panel=label))
    # Unsupported panels have an explicit disposition; they do not get zero-valued error rows.
    sc=['cell_type','method','mae','donors','blocks','rmse','paired_mae_change_vs_tmm','reference','panel']
    bl=['donor','block','cell_type','method','libraries','sections','measured','prediction','reference','panel']
    return (pd.concat(scores) if scores else pd.DataFrame(columns=sc),
            pd.concat(pairs) if pairs else pd.DataFrame(columns=bl),dict(audit,panels=panel_info),
            pd.DataFrame(selected,columns=['reference','panel','cell_class','gene_id']))
