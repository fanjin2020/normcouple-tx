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
AUXILIARY = ('dose_contrasts.csv', 'power_trajectories.csv', 'power_component_donors.csv', 'power_paired_contrasts.csv',
             'common_support_donor_subsets.csv', 'matching_balance.csv',
             'clinical_coefficient_pairs.csv.gz', 'clinical_markers.csv', 'diagnostic_performance.csv',
             'clinical_linear_checks.csv', 'analytic_curves.csv')


def read(path):
    return pd.read_csv(path, keep_default_na=False, na_values=['', 'nan', 'NaN'], float_precision='round_trip')


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
    f['normalization.csv'] = n.rename(columns={'mc_low': 'low', 'mc_high': 'high'})
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
        frame = pd.DataFrame(schema['resources'], columns=item['columns']) if number == 1 else frames[logical[number]]
        result[number] = frame[item['columns']].copy()
        if result[number].empty:
            raise ValueError('Empty required manuscript table: ' + str(number))
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
        views[name] = {'file': 'TableS%d.csv' % number}
    views['dose_retention.csv'] = {'file': 'TableS7.csv', 'filter': 'all_classes'}
    views['dose_by_class.csv'] = {'file': 'TableS7.csv', 'filter': 'individual_classes'}
    for name in AUXILIARY:
        compression = {'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None
        auxiliary[name].to_csv(destination / name, index=False, compression=compression, lineterminator='\n')
        views[name] = {'file': name}
    atomic_json(destination / 'catalog.json', {'schema': 1, 'views': views})


def workbook(data, output):
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Border, Font, Side
    except ImportError as error:
        raise ValueError('Excel export needs optional openpyxl; CSV analysis and plots do not. See requirements-export.txt.') from error
    schema = read_json(PUB / 'table_schema.json')
    book = Workbook()
    toc = book.active; toc.title = 'Contents'; toc.append(['Table', 'Title', 'Rows'])
    line = Side(style='thin', color='000000')
    for item in schema['tables']:
        number = item['number']; name = 'TableS%d' % number
        frame = load_frame(data, name + '.csv')
        toc.append([name, item['title'], len(frame)])
        sheet = book.create_sheet(name)
        sheet.append([item['title']]); sheet.append([item['note']]); sheet.append(list(frame.columns))
        for row in frame.itertuples(index=False, name=None):
            sheet.append([None if pd.isna(v) else v for v in row])
        for cell in sheet[3]:
            cell.font = Font(name='Arial', size=10, bold=True)
            cell.border = Border(top=line, bottom=line)
        for cell in sheet[sheet.max_row]:
            cell.border = Border(bottom=line)
        sheet.freeze_panes = 'A4'
        sheet.auto_filter.ref = 'A3:' + sheet.cell(sheet.max_row, sheet.max_column).coordinate
        for col in sheet.columns:
            sheet.column_dimensions[col[0].column_letter].width = min(32, max(14, len(str(col[2].value)) + 2))
        sheet['A2'].alignment = Alignment(wrap_text=True)
    book.save(output)


def render(data, output, formats):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    subprocess.run([sys.executable, str(PUB / 'reproduce_figures.py'), '--data', str(data),
                    '--output', str(output), '--formats', *formats], env=env, check=True)
    expected = {'Fig' + str(i) for i in range(1, 10)} | {'FigS' + str(i) for i in range(1, 9)}
    for extension in formats:
        actual = {p.stem for p in Path(output).glob('*.' + extension)}
        if actual != expected:
            raise ValueError('Figure numbering is incomplete for ' + extension)


def publish(destination, frames, tables, provenance, formats=('pdf',), make_workbook=False):
    destination = separate(destination)
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
        atomic_json(stage / 'tables.json', schema)
        if formats:
            render(stage / 'data', stage / 'figures', formats)
        if make_workbook:
            workbook(stage / 'data', stage / 'Numerical_Tables.xlsx')
        files = receipt(stage, [p.relative_to(stage).as_posix() for p in stage.rglob('*') if p.is_file()])
        atomic_json(stage / 'export_manifest.json', dict(provenance=provenance, formats=list(formats),
                    workbook=make_workbook, files=files))
        if destination.exists():
            destination.rmdir()  # Only an empty destination is allowed here.
        os.replace(stage, destination)
    print('PUBLICATION OUTPUT', destination, flush=True)
    print('16 supplementary CSVs; 11 additional figure/diagnostic datasets.' +
          (' 9 main and 8 supplementary figures.' if formats else ' No figures requested.'))
    return destination


def export_run(work, destination, formats=('pdf',), make_workbook=False):
    from .run import STAGES, validate_run
    work = Path(work).resolve()
    identity = validate_run(work)
    destination = separate(destination, *(work / s for s in STAGES), work / '.runtime', work / '.pipeline', work / '.input-cache')
    frames, sources = collect({s: work / s for s in STAGES})
    return publish(destination, frames, table_frames(frames),
                   dict(mode='computed', run_signature=identity, analysis_sources=sources), formats, make_workbook)


def export_published(destination, formats=('pdf',), make_workbook=False):
    data = PUB / 'data'
    manifest = read_json(PUB / 'published_manifest.json')
    verify_receipt(data, manifest['files'])
    tables = {i: load_frame(data, 'TableS%d.csv' % i) for i in range(1, 17)}
    frames = {name: load_frame(data, name) for name in AUXILIARY}
    return publish(destination, frames, tables,
                   dict(mode='published_data_replay', source_manifest=sha(PUB / 'published_manifest.json')),
                   formats, make_workbook)
