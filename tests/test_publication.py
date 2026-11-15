from pathlib import Path
import csv
import json
import subprocess
import sys
import unittest
import numpy as np
import pandas as pd
from pipeline.publication import load_frame
from pipeline.io import read_json, verify_receipt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'publication/data'


def data(name):
    return load_frame(DATA, name)


class PublicationTests(unittest.TestCase):
    def test_figure_captions_match_current_numbering(self):
        captions=json.loads((ROOT/'publication/figure_captions.json').read_text(encoding='utf-8'))
        self.assertEqual(set(captions),{f'Fig{i}' for i in range(1,10)}|{f'FigS{i}' for i in range(1,8)})
        self.assertIn('Eq. (1)',captions['Fig1'])
        self.assertTrue(all(not text.endswith('.') for text in captions.values()))

    def test_compact_table_rows_equal_source_csv(self):
        content=json.loads((ROOT/'publication/table_schema.json').read_text(encoding='utf-8'))
        self.assertEqual([x['number'] for x in content['tables']],list(range(1,16)))
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

    def test_help_is_dependency_light(self):
        p=subprocess.run([sys.executable,str(ROOT/'reproduce.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)


if __name__=='__main__':
    unittest.main()
