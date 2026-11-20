from pathlib import Path
import csv
import json
import subprocess
import sys
import unittest
import numpy as np
import pandas as pd
from pipeline.publication import load_frame,table_frames
from pipeline.io import read_json, verify_receipt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'publication/data'


def data(name):
    return load_frame(DATA, name)


class PublicationTests(unittest.TestCase):
    def test_intervals_are_for_paired_contrasts(self):
        n = data('TableS2.csv')
        a = n.query("reference_source == 'LIBD' and target_source == 'LIBD' and method == 'TMM' and metric == 'fdp'").iloc[0]
        self.assertEqual(a.paired_reference, 'Observed total')
        self.assertEqual(np.round(100*a[['paired_change','paired_low','paired_high']].to_numpy(float),2).tolist(), [-3.44,-3.83,-3.04])
        self.assertTrue(n.loc[n.method.eq('Observed total'),'paired_change'].isna().all())
        g = data('TableS4.csv')
        keys = ['reference_source','target_source','grid_id','outcome','method','metric']
        a = g[g.paired_reference.notna()].set_index(keys,verify_integrity=True)
        b = data('dose_contrasts.csv').set_index(keys,verify_integrity=True).loc[a.index]
        self.assertEqual(len(a),len(b))
        self.assertTrue(a.paired_reference.eq('tmm').all())
        for a_col, b_col in [('paired_change','mean'),('paired_low','low'),('paired_high','high')]:
            np.testing.assert_allclose(a[a_col], b[b_col], equal_nan=True, rtol=0, atol=0)

    def test_figure_captions_match_current_numbering(self):
        captions=json.loads((ROOT/'publication/figure_captions.json').read_text(encoding='utf-8'))
        self.assertEqual(set(captions),{f'Fig{i}' for i in range(1,10)}|{f'FigS{i}' for i in range(1,12)})
        self.assertIn('Eq. (1)',captions['Fig1'])
        self.assertTrue(all(not text.endswith('.') for text in captions.values()))

    def test_compact_table_rows_equal_source_csv(self):
        content=json.loads((ROOT/'publication/table_schema.json').read_text(encoding='utf-8'))
        self.assertEqual([x['number'] for x in content['tables']],list(range(1,22)))
        manifest=read_json(ROOT/'publication/published_manifest.json')
        verify_receipt(DATA,manifest['files'])
        for table in content['tables']:
            frame=data('TableS%d.csv'%table['number'])
            self.assertEqual(frame.columns.tolist(),table['columns'])
            self.assertEqual(len(frame),manifest['table_rows'][str(table['number'])])

    def test_clinical_summary_from_gene_results(self):
        pairs=data('clinical_coefficient_pairs.csv.gz')
        expected={'LIBD':(2752,2031,48,769),'SEAAD':(3272,1345,105,2032)}
        for ref,vals in expected.items():
            g=pairs[pairs.reference.eq(ref)]
            self.assertEqual(len(g),16161)
            a,b=g.tmm_q.le(.05),g.mrn_q.le(.05)
            self.assertEqual((int(a.sum()),int(b.sum()),int((b&~a).sum()),int((a&~b).sum())),vals)
            self.assertEqual(int((a&b&(g.tmm_beta*g.mrn_beta<0)).sum()),0)

    def test_common_support_headline_and_complete_denominators(self):
        f=data('common_support_summary.csv')
        g=f[f.reference_source.eq('LIBD')&f.target_source.eq('GSE174367')&f.outcome.eq('null')&f.support_rule.eq('positive_only')]
        expected={('legacy','tmm'):31.866666667,('legacy','marker_matched'):5.2,
                  ('random','tmm'):49.644444444,('random','marker_matched'):22.222222222}
        for (group,method),value in expected.items():
            a=g[g.subset_kind.eq(group)&g.method.eq(method)]
            self.assertEqual(len(a),1)
            self.assertAlmostEqual(a.fdp.iloc[0]*100,value,places=2)
        self.assertTrue(f.planned_records.eq(f.feasible_records).all())

    def test_finite_clinical_invariance(self):
        f=data('clinical_finite.csv')
        invariant=f[f.method.isin(['tmm_exclude_both','marker_denominator_only','marker_matched'])]
        self.assertEqual(len(invariant),18)
        self.assertTrue(invariant.threshold_crossings.eq(0).all())
        self.assertTrue(invariant.coefficient_change_rms.eq(0).all())
        counts=f[f.method.eq('tmm')&f.log_count_step.eq(.1)].set_index('reference').threshold_crossings.to_dict()
        self.assertEqual(counts,{'LIBD':3191,'SEAAD':3711})

    def test_qualified_diagnostic_has_no_predictive_discrimination(self):
        f=data('diagnostic_qualification.csv')
        self.assertTrue(f[f.method.isin(['tmm','tmm_exclude_factor'])].qualified_scores.eq(0).all())
        self.assertEqual(f.qualified_scores.sum(),10646)
        self.assertTrue(f.qualified_scores.eq(f.numerical_invariance_records).all())
        p=data('diagnostic_performance.csv')
        self.assertEqual(len(p),50)
        self.assertTrue(p.auc.eq(.5).all())
        self.assertTrue(p.alert_rate.eq(0).all())

    def test_gene_families_are_complete(self):
        definitions=ROOT/'publication/input_definitions'
        t=pd.read_csv(definitions/'complete_test_gene_families.csv')
        m=pd.read_csv(definitions/'fixed_marker_panels.csv')
        for ref in ['LIBD','SEAAD']:
            genes=set(t[t.reference_source.eq(ref)].gene_id)
            self.assertEqual(len(genes),100)
            panel=set(m[m.reference.eq(ref)].gene_id)
            self.assertFalse(genes & panel)

    def test_retained_dose(self):
        f=data('dose_retention.csv')
        self.assertEqual(f.requested_counts.sum(),30235885)
        self.assertEqual(f.retained_counts.sum(),29154250)
        self.assertEqual(f.draws.sum(),2100)
        self.assertEqual(f.zero_budget_draws.sum(),0)

    def test_figure6_all_points_fit_axis_contract(self):
        f=data('common_support_donor_subsets.csv')
        f=f[f.method.eq('marker_matched')&f.outcome.eq('null')&f.support_rule.eq('positive_only')]
        p=f.groupby(['reference_source','target_source','base_subset']).fdp.mean()
        self.assertEqual(len(p),28)
        self.assertGreater(p.max(),.45)
        self.assertLessEqual(p.max(),1.)

    def test_variance_decomposition(self):
        f=data('power_components.csv')
        np.testing.assert_allclose(f.error_variance,f.leakage_variance+f.noise_variance+f.twice_covariance,rtol=1e-10,atol=1e-12)
        f=f[f.reference_source.eq('SEAAD')&f.target_source.eq('SEAAD')].set_index('method')
        share=(f.loc['marker_matched','leakage_variance']-f.loc['tmm','leakage_variance'])/(f.loc['marker_matched','error_variance']-f.loc['tmm','error_variance'])
        self.assertAlmostEqual(share*100,91.65,places=2)
        p=data('power_paired_contrasts.csv')
        np.testing.assert_allclose(p.low,p['mean']-1.96*p.mcse,atol=1e-12)
        np.testing.assert_allclose(p.high,p['mean']+1.96*p.mcse,atol=1e-12)
        row=p[p.reference_source.eq('SEAAD')&p.outcome.eq('up025')&p.method.eq('marker_matched')].iloc[0]
        self.assertAlmostEqual(row['mean']*100,-27.4083333333,places=8)
        self.assertAlmostEqual(row['low']*100,-28.6209791781,places=8)
        self.assertAlmostEqual(row['high']*100,-26.1956874886,places=8)

    def test_gene_level_probes_and_realized_denominators(self):
        both=data('TableS20.csv').set_index(['cohort','probe'],verify_integrity=True)
        self.assertEqual(len(both),28)
        f=both.loc['GSE125583']
        self.assertEqual(len(f),14)
        singles=f[f.selected_genes.eq(1)]
        self.assertEqual(singles.selected_ids.nunique(),3)
        for ref in ('LIBD','SEAAD'):
            self.assertTrue(singles[ref+'_family_BH_crossings'].eq(0).all())
            self.assertTrue(singles[ref+'_singleton_beta_change'].abs().gt(0).all())
            self.assertTrue(f.loc[f.selected_genes.ne(1),ref+'_singleton_beta_change'].isna().all())
            self.assertTrue(f.loc[f.selected_genes.eq(5),ref+'_own_BH_crossings'].eq(0).all())
        hi=singles.loc['single_high_up']
        self.assertEqual(hi.selected_ids,'ENSG00000149485')
        self.assertLess(hi.LIBD_singleton_beta_change,0)
        self.assertGreater(hi.SEAAD_singleton_beta_change,0)
        full=f.loc['test_family_up']
        np.testing.assert_allclose([full.AD_count_share_pct,full.AD_library_change_pct,full.AD_effective_median_pct,
                                    full.control_library_change_pct,full.control_effective_median_pct],
                                   [99.25286,10.43851,9.70946,0,2.41358],atol=.00001)
        self.assertEqual(tuple(full[['LIBD_family_BH_crossings','SEAAD_family_BH_crossings']]),(3191,3711))
        self.assertEqual(tuple(f.loc['test_family_down',['LIBD_family_BH_crossings','SEAAD_family_BH_crossings']]),(2721,3819))

    def test_secondary_clinical_probes_keep_unavailable_fits_missing(self):
        f=data('TableS20.csv').query("cohort=='GSE53697'").set_index('probe',verify_integrity=True)
        self.assertEqual(len(f),14)
        self.assertTrue(f.filter(regex='^LIBD_').isna().all().all())
        self.assertTrue(f.SEAAD_family_BH_crossings.eq(0).all())
        self.assertTrue(f.loc[f.selected_genes.eq(1),'SEAAD_singleton_beta_change'].abs().gt(0).all())
        full=f.loc['test_family_up']
        self.assertEqual(full.selected_genes,14745)
        np.testing.assert_allclose([full.AD_count_share_pct,full.AD_effective_median_pct],[99.37770,10.19509],atol=.00001)

    def test_clinical_probes_require_original_measurements(self):
        import tempfile
        from pipeline.clinical_probes import build
        with tempfile.TemporaryDirectory() as d:
            (Path(d)/'receipt.json').write_text('{"files":{}}',encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'Incomplete receipt'):
                build(d,data('TableS19.csv'))

    def test_help_is_dependency_light(self):
        p=subprocess.run([sys.executable,str(ROOT/'reproduce.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)

    def test_paired_imaging_support_and_mixed_directions(self):
        f=data('TableS21.csv')
        self.assertEqual(len(f),12)
        self.assertEqual(set(f.panel),{'reference_transfer'})
        self.assertEqual(f.groupby('reference').cell_type.nunique().tolist(),[6,6])
        self.assertTrue(set(f.donors)=={3,5} and set(f.blocks)=={4,6})
        np.testing.assert_allclose(f.mrn_mae-f.tmm_mae,f.mrn_minus_tmm_mae,atol=1e-14)
        self.assertEqual((f.mrn_minus_tmm_mae<0).sum(),5)
        np.testing.assert_allclose(100*f.mrn_minus_tmm_mae.agg(['min','max']),[-.732282,.931016],atol=.00001)

    def test_imaging_rejects_missing_or_misaligned_members(self):
        from pipeline.imaging_audit import paired_table,CLASSES
        rows=[dict(reference='R',panel='p',method=m,cell_type=c,donor=d,block=b,libraries=1,
                   sections=1,measured=.25,prediction=(.25 if d=='a' else .75))
              for c in CLASSES for d,b in [('a','a1'),('a','a2'),('b','b1')]
              for m in ('dtangle_marker','dtangle_tmm')]
        blocks=pd.DataFrame(rows)
        audit=dict(panels=[dict(reference='R',panel='p',status='estimable')],matched_donors=2,
                   matched_blocks=3,matched_libraries=3,measured_block_class_pairs=len(rows)//2)
        out=paired_table(blocks,audit)
        self.assertTrue(out.tmm_mae.eq(.25).all())  # donor-equal; pooled blocks would give 1/6
        with self.assertRaises(ValueError):paired_table(blocks.iloc[1:],audit)
        blocks.loc[0,'measured']=.1
        with self.assertRaises(ValueError):paired_table(blocks,audit)

    def test_clinical_probe_roster_rule(self):
        from pipeline.clinical_probes import marker_roster
        f=pd.DataFrame(dict(cohort='C',reference='R',cell_class='A',
            gene_id=['g%02d'%i for i in range(11)],marker_probe=[True,True]+[False]*9))
        self.assertEqual(marker_roster(f.sample(frac=1,random_state=3)).gene_id.tolist(),['g00','g01'])
        f.loc[[1,2],'marker_probe']=[False,True]
        with self.assertRaises(ValueError):marker_roster(f)

    def test_workbook_has_three_rules_and_no_interior_grid(self):
        import importlib.util,tempfile
        if importlib.util.find_spec('openpyxl') is None:self.skipTest('Optional workbook dependency')
        from openpyxl import load_workbook
        from pipeline.publication import workbook
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);(p/'data').mkdir()
            (p/'tables.json').write_text(json.dumps({'tables':[dict(number=21,title='Imaging',note='Six classes')]}))
            (p/'data/catalog.json').write_text(json.dumps({'views':{'TableS21.csv':{'file':'TableS21.csv'}}}))
            (p/'data/TableS21.csv').write_text('group,n\nA,1\nB,2\n')
            workbook(p/'data',p/'book.xlsx')
            wb=load_workbook(p/'book.xlsx')
            try:
                for sh,head,last in [(wb['Guide'],6,7),(wb['S21'],5,7)]:
                    self.assertFalse(sh.sheet_view.showGridLines)
                    self.assertEqual(sh.cell(head,1).border.top.style,'medium')
                    self.assertEqual(sh.cell(head,1).border.bottom.style,'thin')
                    self.assertEqual(sh.cell(last,2).border.bottom.style,'medium')
                self.assertIsNone(wb['S21']['A6'].border.left.style)
            finally:wb.close()

    def test_tables_only_workbook_maps_all_figures_and_preserves_class_types(self):
        import importlib.util,tempfile
        if importlib.util.find_spec('openpyxl') is None:self.skipTest('Optional workbook dependency')
        from openpyxl import load_workbook
        from pipeline.publication import export_published,figure_index
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)/'out'
            export_published(out,formats=(),make_workbook=True)
            self.assertFalse((out/'figures').exists())
            index=figure_index(read_json(DATA/'catalog.json')['views'])
            wb=load_workbook(out/'Numerical_Tables.xlsx',read_only=True)
            try:
                actual=[(r[0],r[1]) for r in wb['Guide'].values if str(r[0]).startswith('Fig. ')]
                self.assertEqual(actual,[('Fig. '+k[3:],'; '.join(v['sources'])) for k,v in index.items()])
                rows=list(wb['S7'].values)
                j=rows[4].index('cell_class')
                classes=[row[j] for row in rows[5:]]
                self.assertEqual(classes.count('all'),4)
                self.assertEqual([v for v in classes if type(v) is int],list(data('dose_by_class.csv').cell_class))
            finally:wb.close()

    def test_computed_resource_catalog_includes_external_clinical_cohort(self):
        catalog=read_json(DATA/'catalog.json')
        frames={name:data(name) for name in catalog['views']}
        frames['dose_table.csv']=data('TableS7.csv')
        frames['marker_summary.csv']=data('TableS11.csv')
        before=table_frames(frames)
        frames['clinical_estimators.csv']=pd.DataFrame(dict(cohort=['GSE125583','GSE53697'],samples=[288,17]))
        after=table_frames(frames)
        self.assertEqual(len(after[1]),len(before[1]))
        self.assertEqual((after[1].Resource=='GSE53697').sum(),1)
        self.assertEqual(after[1].iloc[-1].Resource,'GSE53697')
        self.assertEqual(after[1].iloc[-1].Donors,17)


if __name__=='__main__':
    unittest.main()
