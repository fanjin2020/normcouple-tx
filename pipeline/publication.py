"""Build manuscript data from completed computations, with no archived-fit fallback."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
from .io import ROOT, atomic_json, inside, read_json, receipt, separate, sha, verify_receipt

PUB = ROOT / 'publication'
TABLE_ALIASES = {'normalization.csv': 2, 'pseudocount.csv': 3, 'dose_grid.csv': 4,
                 'common_support_summary.csv': 5, 'common_support_contrasts.csv': 6,
                 'power_components.csv': 8, 'clinical_design.csv': 9, 'clinical_stability.csv': 10,
                 'actual_span.csv': 12, 'clinical_finite.csv': 13,
                 'diagnostic_qualification.csv': 14, 'clinical_derivatives.csv': 15,
                 'cross_method.csv': 16}
EXTRA_TABLES = {
    17: dict(file='cross_precision.csv',title='Paired Monte Carlo precision',
             note='Within-dataset differences; conditional on the observed donors and references. The fixed repetition count is not chosen by the observed interval.'),
    18: dict(file='clinical_estimators.csv',title='Clinical estimator and gene-policy comparison',
             note='Fixed baseline outcomes and BH family within each cohort. Discovery counts are not false-discovery estimates.'),
    19: dict(file='clinical_estimator_contrasts.csv',title='Paired clinical sensitivity across estimators',
             note='Baseline methods are compared with dtangle/TMM; each perturbation is compared with its own unperturbed method.'),
    20: dict(file='clinical_probes.csv'),21: dict(file='imaging_comparison.csv')}
AUXILIARY = ('dose_contrasts.csv', 'power_trajectories.csv', 'power_component_donors.csv', 'power_paired_contrasts.csv',
             'common_support_donor_subsets.csv', 'matching_balance.csv',
             'clinical_coefficient_pairs.csv.gz', 'clinical_markers.csv', 'diagnostic_performance.csv',
             'clinical_linear_checks.csv', 'analytic_curves.csv')

# Logical views, in manuscript order. The catalog resolves shared physical CSVs.
FIGURE_INPUTS = {
    'Fig1': ('analytic_curves.csv',), 'Fig2': ('normalization.csv',),
    'Fig3': ('pseudocount.csv',), 'Fig4': ('dose_grid.csv', 'dose_contrasts.csv'),
    'Fig5': ('common_support_summary.csv',), 'Fig6': ('common_support_donor_subsets.csv',),
    'Fig7': ('power_trajectories.csv', 'power_components.csv', 'power_component_donors.csv'),
    'Fig8': ('clinical_coefficient_pairs.csv.gz', 'clinical_stability.csv'),
    'Fig9': ('clinical_finite.csv', 'clinical_probes.csv'),
    'FigS1': ('dose_grid.csv', 'dose_contrasts.csv'),
    'FigS2': ('dose_retention.csv', 'dose_by_class.csv'),
    'FigS3': ('common_support_contrasts.csv',), 'FigS4': ('clinical_design.csv',),
    'FigS5': ('clinical_markers.csv',), 'FigS6': ('actual_span.csv',),
    'FigS7': ('diagnostic_qualification.csv', 'clinical_derivatives.csv'),
    'FigS8': ('cross_method.csv',), 'FigS9': ('cross_precision.csv',),
    'FigS10': ('clinical_estimator_contrasts.csv', 'clinical_estimators.csv'),
    'FigS11': ('imaging_comparison.csv',),
}


def figure_ids(views, estimator_only=False):
    if estimator_only:
        return {'FigS8', 'FigS9', 'FigS10'}
    keys = set(FIGURE_INPUTS)
    if 'clinical_estimators.csv' not in views:
        keys -= {'FigS9', 'FigS10'}
    if 'imaging_comparison.csv' not in views:
        keys.remove('FigS11')
    return keys


def figure_index(views, estimator_only=False, *, partial=False):
    scope = figure_ids(views, estimator_only)
    result = {}
    for key, inputs in FIGURE_INPUTS.items():
        if key not in scope:
            continue
        if key == 'Fig9' and 'clinical_probes.csv' not in views:
            inputs = ('clinical_finite.csv',)
        missing = set(inputs).difference(views)
        if missing:
            if partial:
                continue
            raise ValueError(f"{key}: missing figure inputs: {', '.join(sorted(missing))}")
        result[key] = dict(sources=sorted({views[s]['file'] for s in inputs}),
                           views={s: views[s] for s in inputs})
    return result


def read(path):
    return pd.read_csv(path, keep_default_na=False, na_values=['', 'nan', 'NaN'], float_precision='round_trip')


def attach_paired(arms, pairs, keys, reference):
    """Add within-repetition contrast intervals without changing arm rows or marginal intervals."""
    fields = {'mean': 'paired_change', 'low': 'paired_low', 'high': 'paired_high'}
    block = pairs[keys + list(fields)].rename(columns=fields).copy()
    block['paired_reference'] = reference
    return arms.merge(block, on=keys, how='left', validate='one_to_one', sort=False)


def load_frame(data, logical):
    """Plots and tables share physical CSV files through explicit, checked views."""
    data = Path(data)
    catalog = read_json(data / 'catalog.json')['views']
    if logical not in catalog:
        raise ValueError('Unregistered publication source: ' + logical)
    view = catalog[logical]
    frame = read(inside(data, view['file']))
    if view.get('filter') == 'all_classes':
        frame = frame[frame.cell_class.astype(str).eq('all')]
    elif view.get('filter') == 'individual_classes':
        frame = frame[~frame.cell_class.astype(str).eq('all')]
        frame = frame.copy()
        frame['cell_class'] = frame.cell_class.astype(int)
    return frame.reset_index(drop=True)


def analytic_curves():
    signal = np.logspace(-3, 3, 240)
    values = {'signal_over_pseudocount': signal, 'derivative_log_denominator': -signal/(1+signal)/np.log(2)}
    for step in (.1, .25, .5):
        values['change_' + str(step)] = np.log2(signal/(1+step)+1)-np.log2(signal+1)
    return pd.DataFrame(values)


def cross_method_table(summary):
    keys = ['reference_source','target_source','effect','method','probe']
    values = summary.pivot(index=keys, columns='metric', values='donor_equal_mean')
    metadata = summary.groupby(keys)[['donors','datasets']].first()
    result = metadata.join(values[['fdp','power','coverage','availability','proportion_rmse',
                                  'max_fraction_change','weight_l2','beta_rms','bh_crossings']])
    for metric in ['fdp','beta_rms']:
        subset = summary[summary.metric.eq(metric)].set_index(keys)
        result[metric+'_mc_se'] = subset.conditional_mc_se
    for metric in ['max_fraction_change','weight_identity_error']:
        subset = summary[summary.metric.eq(metric)].set_index(keys)
        result[metric+'_maximum'] = subset.max_value
    return result.reset_index()


def clinical_weight_geometry(fractions, cohort):
    classes = ['Astro','EndoMural','Excit','Inhib','Micro','OPC','Oligo']
    x = cohort.group.eq('AD').to_numpy(float)
    cov = np.column_stack([np.ones(len(x)),cohort.age.to_numpy(float),cohort.sex.eq('M').to_numpy(float)])
    records = []
    for ref in ['LIBD','SEAAD']:
        residuals, infos = [], []
        for method in ['tmm','marker_matched']:
            f = fractions[fractions.reference.eq(ref) & fractions.method.eq(method)]
            if not f.sample_id.is_unique or set(f.sample_id) != set(cohort.sample_id):
                raise ValueError('Clinical composition and cohort identifiers differ')
            props = f.set_index('sample_id').loc[cohort.sample_id,classes].to_numpy(float)
            c = np.column_stack([cov,props])
            u,s,_ = np.linalg.svd(c,full_matrices=False)
            basis = u[:,s > s[0]*max(c.shape)*np.finfo(float).eps]
            z = x-basis@(basis.T@x)
            residuals.append(z)
            infos.append(float(z@z))
        if min(infos) <= 0:
            raise ValueError('Clinical exposure not estimable')
        distance = float(np.linalg.norm(residuals[1]/infos[1]-residuals[0]/infos[0]))
        error = abs(distance-np.linalg.norm(residuals[1]-residuals[0])/np.sqrt(np.prod(infos)))
        records.append(dict(reference=ref,method='marker_matched',coefficient_weight_l2_vs_tmm=distance,weight_identity_error=error))
    return pd.DataFrame(records)


def collect(paths):
    """Read only the declared analysis outputs. Missing stages are errors."""
    paths = {k: Path(v) for k, v in paths.items()}
    sources = {}
    def source(stage, name):
        p = inside(paths[stage], name)
        if not p.is_file():
            raise ValueError('Missing analysis output: ' + str(p))
        sources[stage + '/' + name] = sha(p)
        return read(p)
    f = {}
    norm = source('normalization', 'tables/source_summary.csv')
    norm = norm[norm.condition.eq('up05')].copy()
    names = {'dt_observed_a1_c1000000': 'Observed total', 'dt_fixed_a1_c1000000': 'Fixed total',
             'dt_tmm_a1_c1000000': 'TMM', 'dt_deseq2_a1_c1000000': 'DESeq2 poscounts'}
    n = norm[norm.method.isin(names) & norm.metric.isin(
        ['fdp','null_coverage','power','proportion_rmse','availability','mean_se','abs_leakage'])].copy()
    n['method'] = n.method.map(names)
    ref = 'dt_observed_a1_c1000000'
    paired = source('normalization', 'tables/paired_summary.csv')
    paired = paired.loc[paired.condition.eq('up05') & paired.method_b.eq(ref) &
                        paired.method_a.isin(names) & paired.metric.isin(n.metric)].copy()
    paired['method'] = paired.method_a.map(names)
    f['normalization.csv'] = attach_paired(n.rename(columns={'mc_low': 'low', 'mc_high': 'high'}),
        paired.rename(columns={'mc_low': 'low', 'mc_high': 'high'}),
        ['reference_source','target_source','method','metric'], names[ref])
    pc = norm[norm.method.str.match(r'dt_(observed|fixed)_a') & ~norm.method.str.endswith('_native') &
              ~norm.method.str.contains('a0.01_', regex=False) &
              norm.metric.isin(['fdp','null_coverage','proportion_rmse','availability'])].copy()
    pc['denominator'] = np.where(pc.method.str.contains('observed'), 'Observed total', 'Fixed total')
    pc['pseudocount'] = pc.method.str.extract(r'_a([\d.]+)_')[0].astype(float)
    pc['scale'] = pc.method.str.extract(r'_c([\d]+)')[0].astype(float)
    pc['ratio'] = pc.pseudocount / pc.scale
    f['pseudocount.csv'] = pc.rename(columns={'mc_low': 'low', 'mc_high': 'high'})
    grid = source('marker_grid', 'final_results/grid_summary.csv')
    f['dose_grid.csv'] = grid[grid.method.isin(['tmm','marker_matched','marker_denominator_only','oracle']) &
        grid.grid_id.ne('f10_a050_random') & grid.metric.isin(['fdp','power','coverage','proportion_rmse','availability'])].copy()
    contrast = source('marker_grid', 'final_results/grid_paired_vs_TMM.csv')
    f['dose_contrasts.csv'] = contrast[contrast.method.eq('marker_matched') & contrast.reference.eq('tmm') &
        contrast.grid_id.ne('f10_a050_random') & contrast.metric.isin(['fdp','power','coverage'])].copy()
    f['dose_grid.csv'] = attach_paired(f['dose_grid.csv'], f['dose_contrasts.csv'],
        ['reference_source','target_source','grid_id','outcome','method','metric'], 'tmm')
    f['power_paired_contrasts.csv'] = contrast[contrast.grid_id.eq('zero') & contrast.metric.eq('power') &
        contrast.reference_source.eq(contrast.target_source) & contrast.outcome.ne('null') &
        contrast.method.isin(['marker_matched','marker_denominator_only'])].copy()
    mechanism = source('mechanism', 'final_results/summary.csv')
    f['power_trajectories.csv'] = mechanism[mechanism.reference_source.eq(mechanism.target_source) &
        mechanism.condition.isin(['up01','up025','up05']) & mechanism.method.isin(['tmm','marker_matched']) &
        mechanism.rule.eq('full_six') & mechanism.metric.eq('power')].copy()
    f['actual_span.csv'] = mechanism[mechanism.condition.eq('up05') &
        mechanism.method.eq('mead_reference_converted') & mechanism.rule.isin(['full_six','adaptive']) &
        mechanism.metric.isin(['availability','coverage','delivered_coverage','fdp'])].copy()
    power = source('marker_grid', 'final_results/power_mechanism_gene_audit.csv')
    power = power[power.outcome.eq('up025') & power.truth.ne(0)].copy()
    power['effect_group'] = 'alternative'
    power['empirical_sd'] = np.sqrt(power.error_variance)
    keys = ['reference_source','target_source','outcome','method','effect_group']
    metrics = ['bias','leakage_mean','noise_mean','error_variance','leakage_variance','noise_variance',
        'twice_covariance','empirical_sd','rms_se','mean_residual_variance','mean_max_sample_hc3_contribution',
        'mean_t','mean_absolute_t','conditional_rejection','mean_beta','paired_beta_difference','paired_se_difference']
    donor = power.groupby(keys + ['target_donor'], as_index=False, dropna=False)[metrics].mean()
    f['power_component_donors.csv'] = donor
    f['power_components.csv'] = donor.groupby(keys, as_index=False, dropna=False)[metrics].mean()
    base = 'common_support'
    s = source(base, 'common_support/method_summary.csv')
    f['common_support_summary.csv'] = s
    f['common_support_donor_subsets.csv'] = source(base, 'common_support/donor_subset_results.csv')
    keys = ['reference_source','target_source','subset_kind','support_rule','outcome','method']
    measures = ['fdp','power','coverage','delivered_coverage','availability','proportion_rmse','bias','mean_se',
                'coefficient_error_rms','exposure_r2','exposure_residual_ss','leverage']
    p = s.set_index(keys)[measures]
    if not p.index.is_unique:
        raise ValueError('Duplicate common-support method/setting keys')
    contrasts = []
    for comparator in ['tmm','marker_denominator_only','tmm_exclude_both']:
        left = p.xs('marker_matched', level='method')
        right = p.xs(comparator, level='method')
        if set(right.index) != set(left.index):
            raise ValueError('Unpaired common-support settings')
        right = right.reindex(left.index)
        delta = (left - right).reset_index()
        delta['comparison'] = 'marker_matched_minus_' + comparator
        contrasts.append(delta)
    f['common_support_contrasts.csv'] = pd.concat(contrasts, ignore_index=True)
    d = source(base, 'common_support/dose_retention.csv')
    dose = d.groupby(['reference_source','target_source'], as_index=False).agg(
        units=('unit_id','size'), draws=('draws','sum'), zero_budget_draws=('zero_common_budget_draws','sum'),
        requested_counts=('original_counts','sum'), retained_counts=('retained_counts','sum'))
    dose['retained_fraction'] = dose.retained_counts / dose.requested_counts
    bc = source(base, 'common_support/budget_by_class.csv.gz').groupby(
        ['reference_source','target_source','cell_class'], as_index=False)[['requested_counts','retained_counts']].sum()
    bc['retained_fraction'] = bc.retained_counts / bc.requested_counts
    f['dose_table.csv'] = pd.concat([dose.assign(cell_class='all'), bc], ignore_index=True)
    f['matching_balance.csv'] = source(base, 'matching/marker_matching_balance.csv')
    for name, rel in {
        'clinical_stability.csv': 'corrections/clinical_method_stability.csv',
        'clinical_design.csv': 'corrections/clinical_design_diagnostics.csv',
        'clinical_finite.csv': 'clinical_multistep/finite_perturbations.csv',
        'clinical_derivatives.csv': 'clinical_multistep/derivative_steps.csv',
        'diagnostic_qualification.csv': 'corrections/diagnostic_qualification_counts.csv',
        'diagnostic_performance.csv': 'corrections/diagnostic_qualified_performance.csv',
        'clinical_linear_checks.csv': 'clinical_multistep/linear_approximation.csv'}.items():
        f[name] = source(base, rel)
    associations = source(base, 'corrections/clinical_disease_associations.csv.gz')
    if associations.duplicated(['reference','method','gene_id']).any():
        raise ValueError('Duplicate clinical gene/model keys')
    pairs = []
    for ref, g in associations.groupby('reference'):
        t = g[g.method.eq('tmm')].set_index('gene_id')
        m = g[g.method.eq('marker_matched')].set_index('gene_id')
        if set(t.index) != set(m.index) or not len(t):
            raise ValueError('Clinical comparisons do not share the same gene family')
        m = m.reindex(t.index)
        pairs.append(pd.DataFrame({'reference': ref, 'gene_id': t.index, 'tmm_beta': t.beta, 'tmm_se': t.se,
            'tmm_q': t.q, 'mrn_beta': m.beta, 'mrn_se': m.se, 'mrn_q': m.q}).reset_index(drop=True))
        stability = f['clinical_stability.csv']
        mask = stability.reference.eq(ref) & stability.method.eq('marker_matched')
        stability.loc[mask, 'median_se_ratio'] = float((m.se / t.se).median())
    f['clinical_coefficient_pairs.csv.gz'] = pd.concat(pairs, ignore_index=True)
    fractions = source('coupling', 'clinical/composition_predictions.csv')
    cohort = source('coupling', 'clinical/analysis_cohort.csv')
    geometry = clinical_weight_geometry(fractions, cohort)
    f['clinical_stability.csv'] = f['clinical_stability.csv'].merge(geometry,on=['reference','method'],how='left',validate='one_to_one')
    f['cross_method.csv'] = cross_method_table(source('cross_method','final/summary.csv'))
    f['cross_precision.csv'] = source('cross_method','final/paired_precision.csv')
    f['clinical_estimators.csv'] = source('clinical_estimators','final/method_summary.csv')
    f['clinical_estimator_contrasts.csv'] = source('clinical_estimators','final/paired_sensitivity.csv')
    markers = source('marker_associations', 'final_results/marker_gene_associations.csv.gz')
    markers = markers[markers.model.eq('primary') & markers.normalizer.eq('TMM_all') & markers.scope.eq('owner')]
    f['clinical_markers.csv'] = markers[['reference','cell_type','gene_id','available','log2FC',
        'qlf_p_q_planned_family','treat_p_50_q_planned_family']].copy()
    rows = []
    for (ref, cell), g in markers.groupby(['reference','cell_type']):
        valid = g[g.available.eq(True)]
        rows.append(dict(reference=ref, cell_type=cell, planned=len(g), tested=len(valid),
                         absolute_FC_at_least_1_5=int(valid.log2FC.abs().ge(np.log2(1.5)).sum()),
                         BH_discoveries=int(valid.qlf_p_q_planned_family.le(.05).sum()),
                         minimum_q=float(valid.qlf_p_q_planned_family.min()) if len(valid) else np.nan))
    f['marker_summary.csv'] = pd.DataFrame(rows)
    f['analytic_curves.csv'] = analytic_curves()
    return f, sources


def table_frames(frames):
    schema = read_json(PUB / 'table_schema.json')
    result = {}
    logical = {number: name for name, number in TABLE_ALIASES.items()}
    logical.update({7: 'dose_table.csv', 11: 'marker_summary.csv'})
    for item in schema['tables']:
        number = item['number']
        if number in EXTRA_TABLES:
            continue
        frame = pd.DataFrame(schema['resources'], columns=item['columns']) if number == 1 else frames[logical[number]]
        result[number] = frame[item['columns']].copy()
        if result[number].empty:
            raise ValueError('Empty required manuscript table: ' + str(number))
    for number,item in EXTRA_TABLES.items():
        if item['file'] in frames:
            result[number] = frames[item['file']].copy()
    if 'clinical_estimators.csv' in frames:
        cohorts = frames['clinical_estimators.csv'].groupby('cohort').samples.unique()
        if any(len(values)!=1 for values in cohorts):
            raise ValueError('Clinical sample counts differ across reference specifications')
        if 'GSE53697' in cohorts and 'GSE53697' not in set(result[1].Resource):
            row = {'Resource':'GSE53697', 'Analysis':'Independent clinical bulk sensitivity',
                   'Donors':int(cohorts['GSE53697'][0]), 'Nuclei':np.nan,
                   'Planned reference markers':np.nan,
                   'Public source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53697'}
            result[1] = pd.concat([result[1],pd.DataFrame([row])],ignore_index=True)
    return result


def write_data(destination, tables, auxiliary):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    views = {}
    for number, frame in tables.items():
        name = 'TableS%d.csv' % number
        frame.to_csv(destination / name, index=False, lineterminator='\n')
        views[name] = {'file': name}
    for name, number in TABLE_ALIASES.items():
        if number in tables:
            views[name] = {'file': 'TableS%d.csv' % number}
    for number,item in EXTRA_TABLES.items():
        if number in tables:
            views[item['file']] = {'file':'TableS%d.csv' % number}
    if 7 in tables:
        views['dose_retention.csv'] = {'file': 'TableS7.csv', 'filter': 'all_classes'}
        views['dose_by_class.csv'] = {'file': 'TableS7.csv', 'filter': 'individual_classes'}
    for name in AUXILIARY:
        if name not in auxiliary:
            continue
        compression = {'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None
        auxiliary[name].to_csv(destination / name, index=False, compression=compression, lineterminator='\n')
        views[name] = {'file': name}
    atomic_json(destination / 'catalog.json', {'schema': 1, 'views': views})


def workbook(data, output):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Border, Font, Side
        from openpyxl.utils import get_column_letter
    except ImportError as error:
        raise ValueError('Excel export needs optional openpyxl; CSV analysis and plots do not. See requirements-export.txt.') from error
    schema = read_json(Path(data).parent/'tables.json')
    book = Workbook()
    toc = book.active; toc.title = 'Guide'
    text=Font(name='Arial',size=10,color='000000')
    bold=Font(name='Arial',size=10,bold=True,color='000000')
    title=Font(name='Arial',size=13,bold=True,color='000000')
    outer=Side(style='medium',color='000000')
    inner=Side(style='thin',color='000000')
    def rules(sheet, head, bottom, width):
        for j in range(1,width+1):
            c=sheet.cell(head,j);c.border=Border(top=outer,bottom=inner);c.font=bold
            sheet.cell(bottom,j).border=Border(bottom=outer)
        sheet.sheet_view.showGridLines=False
    def banner(sheet,row,last,message,font,height):
        sheet.merge_cells(start_row=row,end_row=row,start_column=1,end_column=last)
        c=sheet.cell(row,1,message);c.font=font;c.alignment=Alignment(wrap_text=True,vertical='center')
        sheet.row_dimensions[row].height=height
    banner(toc,1,6,'Online Resource 2 | Neuroinformatics',title,25)
    banner(toc,2,6,'Normalization coupling and marker restriction in brain transcriptome inference',text,25)
    banner(toc,3,6,'Jin Fan, Chang Liu, Ying Cao, Xiankun Zhang, Chen Song\nCorrespondence: Xiankun Zhang (zhxkun@tust.edu.cn) and Chen Song (songchen@tust.edu.cn)\nSchool of Artificial Intelligence, Tianjin University of Science and Technology, Tianjin 300457, China',text,48)
    banner(toc,4,6,'Rates are fractions, unless a column explicitly states percent. Blank cells are unavailable or not applicable.',text,25)
    for j,v in enumerate(('Table','Title','Rows','Columns','Related figures','Interpretation'),1):toc.cell(6,j,v)
    toc.freeze_panes='A7'
    for j,width in enumerate((10,58,10,10,32,95),1):toc.column_dimensions[get_column_letter(j)].width=width
    for item in schema['tables']:
        number = item['number']; name = 'TableS%d' % number
        frame = load_frame(data, name + '.csv')
        if number == 7:
            # The pooled stratum is text; individual class identifiers are Excel numbers.
            frame['cell_class'] = frame.cell_class.map(lambda value: 'all' if value == 'all' else int(value))
        toc.append(['S%d'%number,item['title'],len(frame),len(frame.columns),item.get('figures'),item['note']])
        for c in toc[toc.max_row]:
            c.font=text;c.alignment=Alignment(wrap_text=True,vertical='center')
        toc.row_dimensions[toc.max_row].height=max(35,15*(2+len(item['note'])//82))
        sheet = book.create_sheet('S%d'%number)
        banner(sheet,1,min(8,len(frame.columns)),name.replace('TableS','Table S')+' | '+item['title'],title,25)
        banner(sheet,2,min(8,len(frame.columns)),item['note'],text,max(48,14*(2+len(item['note'])//75)))
        for j,v in enumerate(frame.columns,1):
            c=sheet.cell(5,j,v);c.alignment=Alignment(wrap_text=True,vertical='center')
        sheet.row_dimensions[5].height=42
        for row in frame.itertuples(index=False, name=None):
            sheet.append([None if pd.isna(v) else v for v in row])
        for row in sheet.iter_rows(min_row=6):
            for c in row:
                c.font=text
                if isinstance(c.value,float) and c.value and abs(c.value)<1e-4:c.number_format='0.000E+00'
        rules(sheet,5,sheet.max_row,len(frame.columns))
        sheet.freeze_panes = 'A6'
        for j,name in enumerate(frame.columns,1):
            sample=frame.iloc[:,j-1].dropna().astype(str)
            longest=0 if sample.empty else max(sample.str.len())
            width=min(78,max(15,len(str(name))+2,longest+2))
            sheet.column_dimensions[get_column_letter(j)].width=width
    rules(toc,6,toc.max_row,6)
    views=read_json(Path(data)/'catalog.json')['views']
    numbers={t['number'] for t in schema['tables']}
    estimators=numbers=={16,17,18,19}
    mapping=figure_index(views,estimators,partial=not estimators and 1 not in numbers)
    if mapping:
        header=toc.max_row+2
        banner(toc,header,6,'Direct inputs used to draw each figure',title,25)
        note='The list below identifies the numerical inputs to each figure.'
        if 'Fig9' in mapping and 'clinical_probes.csv' in mapping['Fig9']['views']:
            note+=' Table S20 supplies the single-gene and full-family probes in Fig. 9c,d and also reports the five-gene interventions.'
        banner(toc,header+1,6,note,text,34)
        for fig,desc in mapping.items():
            row=toc.max_row+1
            toc.cell(row,1,'Fig. '+fig[3:]).font=text
            toc.merge_cells(start_row=row,start_column=2,end_row=row,end_column=6)
            toc.cell(row,2,'; '.join(desc['sources'])).font=text
        for j in range(1,7):
            toc.cell(header,j).border=Border(top=outer)
            toc.cell(header+1,j).border=Border(bottom=inner)
            toc.cell(toc.max_row,j).border=Border(bottom=outer)
    book.save(output)


def identity_check(frame):
    col = 'exact_identity_max_error'
    values = pd.to_numeric(frame[col], errors='raise')
    mask = values.notna()
    if not mask.any() or not np.isfinite(values[mask]).all() or (values[mask] < 0).any():
        raise ValueError('Invalid clinical coefficient-identity records')
    row = frame.loc[values.idxmax()]
    return dict(max_absolute_error=float(row[col]), reference=str(row.reference),
                method=str(row.method), log_count_step=float(row.log_count_step),
                finite_records=int(mask.sum()), total_records=len(frame))


def render(data, output, formats, software_acceptance=False, estimator_only=False):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    subprocess.run([sys.executable, str(PUB / 'reproduce_figures.py'), '--data', str(data),
                    '--output', str(output), '--formats', *formats,
                    *(['--software-acceptance'] if software_acceptance else []),
                    *(['--estimators-only'] if estimator_only else [])], env=env, check=True)
    expected = figure_ids(read_json(Path(data)/'catalog.json')['views'], estimator_only)
    for extension in formats:
        actual = {p.stem for p in Path(output).glob('*.' + extension)}
        if actual != expected:
            raise ValueError('Figure numbering is incomplete for ' + extension)


def publish(destination, frames, tables, provenance, formats=('pdf',), make_workbook=False):
    destination = separate(destination)
    from .io import digest
    recipe = [Path(__file__), PUB/'table_schema.json']
    if formats:
        recipe += [*sorted(PUB.glob('figure*.py')), PUB/'reproduce_figures.py', PUB/'figure_captions.json']
    provenance = dict(provenance, publication_recipe=digest({p.name:sha(p) for p in recipe}))
    if destination.exists() and any(destination.iterdir()):
        old = read_json(destination / 'export_manifest.json')
        if old['provenance'] != provenance or old['formats'] != list(formats) or old['workbook'] != make_workbook:
            raise ValueError('Publication output belongs to a different export; choose a new output directory')
        verify_receipt(destination, old['files'])
        print('VERIFIED PUBLICATION OUTPUT', destination)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.' + destination.name + '-', dir=destination.parent) as temporary:
        stage = Path(temporary) / 'final'
        write_data(stage / 'data', tables, frames)
        schema = read_json(PUB / 'table_schema.json')
        schema['tables'] = [x for x in schema['tables'] if x['number'] in tables]
        if 13 in tables:
            check = identity_check(frames['clinical_linear_checks.csv'])
            note = next(t for t in schema['tables'] if t['number']==13)
            note['note'] += (' Maximum absolute difference between Eq. (5) and the refitted coefficient change: '
                '{max_absolute_error:.3g} ({reference}, {method}, log step {log_count_step:g}).').format(**check)
        atomic_json(stage / 'tables.json', schema)
        if formats:
            render(stage / 'data', stage / 'figures', formats,
                   software_acceptance='software_acceptance' in provenance.get('mode',''),
                   estimator_only=set(tables)=={16,17,18,19})
        if make_workbook:
            workbook(stage / 'data', stage / 'Numerical_Tables.xlsx')
        files = receipt(stage, [p.relative_to(stage).as_posix() for p in stage.rglob('*') if p.is_file()])
        atomic_json(stage / 'export_manifest.json', dict(provenance=provenance, formats=list(formats),
                    workbook=make_workbook, files=files))
        if destination.exists():
            destination.rmdir()  # Only an empty destination is allowed here.
        os.replace(stage, destination)
    print('PUBLICATION OUTPUT', destination, flush=True)
    print(str(len(tables))+' supplementary CSVs; '+str(sum(k in frames for k in AUXILIARY))+' additional figure/diagnostic datasets.' +
          (' Figures recorded in figure_sources.json.' if formats else ' No figures requested.'))
    return destination


def attach_probes(tables, provenance, evidence):
    from .clinical_probes import build
    from .imaging_audit import publication_table
    tables[20], provenance['clinical_probes'] = build(evidence,tables[19])
    tables[21], provenance['imaging_comparison'] = publication_table(evidence)


def export_run(work, destination, formats=('pdf',), make_workbook=False, software_acceptance=False, evidence=None):
    from .run import STAGES, validate_run
    work = Path(work).resolve()
    identity = validate_run(work,allow_acceptance=software_acceptance)
    destination = separate(destination, *(work / s for s in STAGES), work / '.runtime', work / '.pipeline', work / '.input-cache')
    frames, sources = collect({s: work / s for s in STAGES})
    tables=table_frames(frames)
    provenance=dict(mode='software_acceptance_not_research' if software_acceptance else 'computed',
                    run_signature=identity, analysis_sources=sources)
    if evidence is not None:
        if software_acceptance: raise ValueError('Clinical evidence requires a full run')
        separate(destination,evidence)
        attach_probes(tables,provenance,evidence)
        for cohort,hashes in provenance['clinical_probes']['clinical_inputs'].items():
            for ref in ('LIBD','SEAAD'):
                bank=work/'clinical_estimators'/cohort/ref
                for key,p in [('raw',bank/'raw.bin'),('cohort_sha256',bank.parent/'cohort.csv'),('outcome_sha256',bank/'outcomes.bin')]:
                    if sha(p)!=hashes[key]: raise ValueError('Clinical evidence has different inputs: '+cohort+' '+key)
    return publish(destination,frames,tables,provenance,formats,make_workbook)


def export_estimators(work, destination, formats=('pdf',), make_workbook=False):
    from .io import digest
    from .run import completed,framework_identity,stage_files
    work=Path(work).resolve()
    cfg=read_json(work/'pipeline.json')
    identity=cfg.pop('signature')
    if digest(cfg)!=identity or cfg['framework_source']!=framework_identity():
        raise ValueError('Estimator run code/configuration differs')
    smoke=cfg.get('profile')=='acceptance'
    for stage in ('cross_method','clinical_estimators'):
        if not completed(work,stage,identity,strict=True):
            raise ValueError('Incomplete estimator stage: '+stage)
        stage_files(work,stage,smoke)
    locations={16:work/'cross_method/final/summary.csv',
               17:work/'cross_method/final/paired_precision.csv',
               18:work/'clinical_estimators/final/method_summary.csv',
               19:work/'clinical_estimators/final/paired_sensitivity.csv'}
    tables={n:read(path) for n,path in locations.items()}
    tables[16]=cross_method_table(tables[16])
    return publish(destination,{},tables,dict(mode='software_acceptance_estimators' if smoke else 'computed_estimators',
        run_signature=identity,analysis_sources={str(p.relative_to(work)):sha(p) for p in locations.values()}),
        formats,make_workbook)


def export_published(destination, formats=('pdf',), make_workbook=False, evidence=None):
    data = PUB / 'data'
    manifest = read_json(PUB / 'published_manifest.json')
    verify_receipt(data, manifest['files'])
    tables = {int(i): load_frame(data, 'TableS%s.csv' % i) for i in manifest['table_rows']}
    frames = {name: load_frame(data, name) for name in AUXILIARY}
    provenance=dict(mode='published_data_replay',source_manifest=sha(PUB/'published_manifest.json'))
    if evidence is not None:
        separate(destination,evidence)
        archived={i:tables[i] for i in (20,21)}
        attach_probes(tables,provenance,evidence)
        if provenance['clinical_probes']['clinical_inputs'] != manifest['clinical_probe_provenance']['clinical_inputs']:
            raise ValueError('Clinical evidence has different raw samples or outcomes from the archived study')
        if provenance['imaging_comparison'] != manifest['imaging_provenance']:
            raise ValueError('Imaging evidence uses different tissue blocks or study inputs')
        for i in archived:
            pd.testing.assert_frame_equal(archived[i],tables[i],check_dtype=False,rtol=1e-10,atol=1e-13)
        provenance['mode']='published_data_with_verified_clinical_probes'
    return publish(destination,frames,tables,provenance,formats,make_workbook)
