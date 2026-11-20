"""Public bulk cohorts with explicit count columns and sample identities."""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .download import download
from .io import ROOT, read_json, sha


def external_sources(cache, allow_downloads=False):
    spec = read_json(ROOT/'publication/input_definitions/external_bulk_sources.json')
    cache = Path(cache)
    for name, item in spec['files'].items():
        path = cache/name
        if path.exists():
            if sha(path) != item['sha256']:
                raise ValueError('External bulk source checksum differs: ' + name)
        elif allow_downloads:
            download(item['url'], path, item['sha256'])
        else:
            raise ValueError('Missing external bulk cache: ' + str(path))
    return spec


def soft_samples(path):
    rows, current = [], None
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        for text in stream:
            line = text.rstrip()
            if line.startswith('^SAMPLE = '):
                if current:
                    rows.append(current)
                current = dict(sample_id=line.split(' = ', 1)[1], descriptions=[])
            elif current is not None and line.startswith('!Sample_title = '):
                current['title'] = line.split(' = ', 1)[1]
            elif current is not None and line.startswith('!Sample_description = '):
                current['descriptions'].append(line.split(' = ', 1)[1])
            elif current is not None and line.startswith('!Sample_characteristics_ch1 = '):
                key, value = line.split(' = ', 1)[1].split(':', 1)
                key = key.strip().lower().replace(' ', '_')
                if key in current:
                    raise ValueError('Duplicate sample characteristic: ' + key)
                current[key] = value.strip()
    if current:
        rows.append(current)
    frame = pd.DataFrame(rows)
    if frame.empty or not frame.sample_id.is_unique:
        raise ValueError('Invalid public sample identities')
    return frame


def external_cohort(cache, allow_downloads=False):
    spec = external_sources(cache, allow_downloads)
    meta = soft_samples(Path(cache)/'GSE53697_family.soft.gz')
    frame = pd.read_csv(Path(cache)/'GSE53697_RNAseq_AD.txt.gz', sep='\t', dtype={'GeneID':str})
    if not frame.GeneID.is_unique or not frame.GeneID.str.fullmatch(r'\d+').all():
        raise ValueError('External genes require unique Entrez IDs')
    columns, records = [], []
    for row in meta.to_dict('records'):
        codes = [x for x in row['descriptions'] if re.fullmatch(r'[AC]\d+', x)]
        if len(codes) != 1:
            raise ValueError('Missing unique deposited column identifier')
        status = row['disease_status'].lower()
        group = {'control':'Control', "alzheimer's disease":'AD', 'ad':'AD',
                 "advanced alzheimer's disease":'AD',
                 "alzheimer's disease (ad)":'AD'}.get(status)
        if group is None:
            raise ValueError('Unknown external diagnosis: ' + status)
        column = codes[0] + '_raw'
        if column not in frame or (codes[0].startswith('A') != (group=='AD')):
            raise ValueError('Metadata and raw count column disagree')
        columns.append(column)
        records.append(dict(sample_id=row['sample_id'],case_id=codes[0],group=group,
                            tissue=row['tissue_subtype']))
    cohort = pd.DataFrame(records)
    batch_file = ROOT/'publication/input_definitions/GSE53697_batch.csv'
    batches = pd.read_csv(batch_file)
    if not batches.case_id.is_unique or set(batches.case_id)!=set(cohort.case_id):
        raise ValueError('External batch and expression sample identities differ')
    cohort = cohort.merge(batches,on='case_id',how='left',validate='one_to_one',sort=False)
    if not cohort.sequencing_batch.isin([1,2]).all():
        raise ValueError('Unknown external sequencing batch')
    raw = frame.set_index('GeneID')[columns].copy()
    raw.columns = cohort.sample_id
    if (len(cohort)!=17 or cohort.group.value_counts().to_dict()!={'AD':9,'Control':8}
            or not cohort.case_id.is_unique or not np.isfinite(raw.to_numpy()).all()
            or (raw.to_numpy()<0).any() or (raw.sum(0)<=0).any()):
        raise ValueError('External count or cohort validation failed')
    spec = dict(spec, fractional_entries=int((raw.to_numpy()!=np.rint(raw.to_numpy())).sum()),
                covariates=['sequencing_batch'], covariate_table_sha256=sha(batch_file),
                covariate_scope='Sequencing batch from Supplementary file 1A; age/PMI patient-to-library mapping unresolved and sex unreported',
                inference='Fixed-outcome sensitivity; no known clinical effects or FDR truth')
    return raw, cohort, spec


def gene_mapping():
    path = ROOT/'kernels/coupling/gene_id_mapping_snapshot.csv.gz'
    frame = pd.read_csv(path, dtype={'entrez':str,'ensembl':str,'symbol':str}, keep_default_na=False)
    keep = frame[frame.one_to_one].copy()
    if not keep.entrez.is_unique or not keep.ensembl.is_unique:
        raise ValueError('Mapping is not one-to-one')
    return keep.set_index('entrez')


def expression_family(raw, excluded):
    total = raw.sum(0).to_numpy(float)
    keep = (raw.to_numpy(float)/total*1e6 >= 1).sum(1) >= np.ceil(raw.shape[1]/2)
    keep &= ~raw.index.isin(excluded)
    if keep.sum()<100:
        raise ValueError('Insufficient expressed non-marker test family')
    return np.flatnonzero(keep)
