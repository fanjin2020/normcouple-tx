from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'publication'))
import figure_data as fd
from pipeline.io import read_json
from pipeline.publication import figure_ids, figure_index, load_frame


def data(key):
    return load_frame(ROOT / 'publication/data', key)


class FigureInputTests(unittest.TestCase):
    def test_equal_donor_random_summary_and_actual_subsets(self):
        x = fd.subset_points(data('common_support_donor_subsets.csv'))
        a = x[(x.reference_source=='LIBD') & (x.target_source=='GSE174367')].set_index('base_subset')
        self.assertEqual(len(a),9)
        np.testing.assert_allclose(a.loc[['random_001','random_002','random_003'],'fdp'],
                                   [403/750,61/750,36/750])
        np.testing.assert_allclose(a.loc[['legacy','random_mean'],'fdp'],[.052,2/9])
        np.testing.assert_allclose(a.loc[['legacy','random_mean'],'proportion_rmse'],
                                   [.1112,.1108],atol=.00005)

    def test_donor_omission_or_duplication_is_not_silently_averaged(self):
        a=data('common_support_donor_subsets.csv')
        a=a.query("method=='marker_matched' and outcome=='null' and (support_rule=='positive_only' or base_subset=='zero')")
        with self.assertRaisesRegex(ValueError,'target donors'):
            fd.subset_points(a.iloc[1:])
        with self.assertRaisesRegex(ValueError,'repeated'):
            fd.subset_points(pd.concat([a,a.iloc[[0]]]))
        with self.assertRaisesRegex(ValueError,'three random'):
            fd.subset_points(a[~a.base_subset.eq('random_002')])

    def test_fig9_rejects_inconsistent_probe_inputs(self):
        curve=data('clinical_finite.csv')
        probes=data('clinical_probes.csv')
        self.assertEqual(len(fd.fusiform_doses(curve,probes)),14)
        for key,column,value,message in [
            ('test_family_up','SEAAD_family_BH_crossings',0,'S13 and S20'),
            ('test_family_up','AD_input_change_pct',50,'S13 and S20'),
            ('test_family_down','control_effective_q25_pct',1000,'interval'),
            ('single_high_down','selected_ids','wrong-gene','identity')]:
            changed=probes.copy()
            changed.loc[changed.cohort.eq('GSE125583') & changed.probe.eq(key),column]=value
            with self.subTest(column=column), self.assertRaisesRegex(ValueError,message):
                fd.fusiform_doses(curve,changed)

    def test_figure_index_uses_physical_sources_and_order_for_each_scope(self):
        views=read_json(ROOT/'publication/data/catalog.json')['views']
        full=figure_index(views)
        order=[f'Fig{i}' for i in range(1,10)]+[f'FigS{i}' for i in range(1,12)]
        self.assertEqual(list(full),order)
        for key,sources in [('Fig9',[13,20]),('FigS2',[7]),('FigS10',[18,19]),('FigS11',[21])]:
            self.assertEqual(full[key]['sources'],[f'TableS{i}.csv' for i in sources])
        self.assertEqual(set(full['FigS2']['views']),{'dose_retention.csv','dose_by_class.csv'})
        reduced=figure_index({k:v for k,v in views.items() if k not in {'clinical_probes.csv','imaging_comparison.csv'}})
        self.assertEqual(list(reduced),order[:-1])
        self.assertEqual(reduced['Fig9']['sources'],['TableS13.csv'])
        self.assertEqual(list(figure_index(views,True)),['FigS8','FigS9','FigS10'])
        with self.assertRaisesRegex(ValueError,'Fig2: missing'):
            figure_index({k:v for k,v in views.items() if k!='normalization.csv'})

    def test_bh_membership_keeps_full_family(self):
        p=data('clinical_coefficient_pairs.csv.gz')
        s=data('clinical_stability.csv')
        truth={'LIBD':[1983,769,48,13361],'SEAAD':[1240,2032,105,12784]}
        for ref in truth:
            g=p[p.reference.eq(ref)]
            r=s[s.reference.eq(ref) & s.method.eq('marker_matched')].iloc[0].copy()
            _, counts=fd.clinical_status(g,r)
            self.assertEqual(counts.tolist(),truth[ref])
            r['lost'] += 1
            with self.assertRaisesRegex(ValueError,'disagree'):
                fd.clinical_status(g,r)
            missing=g.copy();missing.iloc[0,missing.columns.get_loc('mrn_q')]=np.nan
            with self.assertRaisesRegex(ValueError,'invalid'):
                fd.clinical_status(missing,r)

    def test_imaging_sign_and_units(self):
        f=data('imaging_comparison.csv');a=fd.imaging_rows(f)
        z=a[['LIBD','SEAAD']].to_numpy()
        self.assertEqual(int((z<0).sum()),5)
        self.assertEqual(int((z>0).sum()),7)
        self.assertAlmostEqual(z.min(),-.73,places=2)
        self.assertAlmostEqual(z.max(),.93,places=2)
        self.assertEqual(sorted(set(map(int,a.donors))),[3,5])
        bad=f.copy();bad.loc[0,'mrn_minus_tmm_mae'] *= -1
        with self.assertRaisesRegex(ValueError,'minus TMM'):
            fd.imaging_rows(bad)

    def test_figure_scopes_do_not_borrow_unavailable_imaging(self):
        full={'clinical_estimators.csv','clinical_probes.csv','imaging_comparison.csv'}
        self.assertEqual(len(figure_ids(full)),20)
        self.assertNotIn('FigS11',figure_ids(full - {'imaging_comparison.csv'}))
        self.assertEqual(figure_ids(full,True),{'FigS8','FigS9','FigS10'})

    def test_numeric_zero_is_distinct_from_unavailable(self):
        df=data('clinical_estimator_contrasts.csv')
        q=df[(df.cohort=='GSE53697') & (df.probe=='test_counts')]
        self.assertTrue(q.loc[~q.available,'crossing_fraction'].isna().all())
        self.assertTrue(q.loc[q.available,'crossing_fraction'].eq(0).all())
        self.assertTrue(q.reason.str.contains('HC3',na=False).any())
        self.assertTrue(q.reason.str.contains('MARKER_COVERAGE',na=False).any())
