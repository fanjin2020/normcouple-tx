import unittest
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
from coupling_math import coefficient,exact_change,coefficient_derivative,local_report
from report import validate_diagnostic,read_nonempty


class Mathematics(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(301)
        self.p=self.rng.dirichlet(np.ones(7)*3,90)
        self.x=self.rng.integers(0,2,90).astype(float)
        self.y=self.rng.normal(size=(90,12))+self.x[:,None]*.3
        self.d=self.rng.normal(0,.001,self.p.shape)
        self.d-=self.d.mean(1,keepdims=True)

    def test_exact_identity(self):
        actual,pred=exact_change(self.x,self.y,self.p,self.p+self.d)
        np.testing.assert_allclose(actual,pred,atol=1e-12)

    def test_identity_with_covariates(self):
        cov=self.rng.normal(size=(90,2))
        actual,pred=exact_change(self.x,self.y,self.p,self.p+self.d,cov)
        np.testing.assert_allclose(actual,pred,atol=1e-12)

    def test_analytic_derivative(self):
        h=1e-3
        actual=(coefficient(self.x,self.y,self.p+h*self.d)-coefficient(self.x,self.y,self.p-h*self.d))/(2*h)
        pred=coefficient_derivative(self.x,self.y,self.p,self.d)
        np.testing.assert_allclose(actual,pred,atol=1e-9,rtol=1e-6)

    def test_invariant_diagnostic(self):
        h=.001
        probes={t:self.p.copy() for t in (-h,h,-h/2,h/2,.01,.05,.1)}
        score,rows=local_report(self.x,self.y,self.p,probes,h,np.ones(12))
        self.assertEqual(score['score'],0.)
        self.assertTrue(all(r['exact_identity_max_error']==0 for r in rows))

    def test_local_approximation(self):
        h=.001
        probes={t:self.p+t*self.d for t in (-h,h,-h/2,h/2,.01,.05,.1)}
        score,rows=local_report(self.x,self.y,self.p,probes,h,np.ones(12))
        self.assertTrue(score['stable_nuisance_rank'])
        self.assertLess(rows[0]['linear_error_rms'],rows[-1]['linear_error_rms'])

    def test_rank_deficient_derivative_rejected(self):
        p=np.tile(np.ones(7)/7,(90,1))
        with self.assertRaisesRegex(ValueError,'rank'):
            coefficient_derivative(self.x,self.y,p,self.d)

    def test_exposure_not_estimable(self):
        with self.assertRaisesRegex(ValueError,'estimable'):
            coefficient(np.ones(90),self.y,self.p)

    def test_noisy_composition_has_no_independence_assumption(self):
        p1=self.p+self.d*self.x[:,None]
        actual,pred=exact_change(self.x,self.y,self.p,p1)
        np.testing.assert_allclose(actual,pred,atol=1e-12)


class Validation(unittest.TestCase):
    def frame(self):
        return pd.DataFrame([dict(target_source='EXTERNAL',target_donor=str(d),
            reference_source=r,method=m,outcome='null',score=(d/10 if m=='tmm' else 0),
            fdp=(.4 if d>4 else 0),rep=k) for d in range(10) for r in ['LIBD','SEAAD']
            for m in ['tmm','marker_matched'] for k in range(2)])

    def test_donor_not_reference_split(self):
        rows,performance=validate_diagnostic(self.frame())
        self.assertTrue(rows.groupby('donor_group').fold.nunique().eq(1).all())
        self.assertTrue(rows.evaluated.all())
        self.assertEqual(len(performance),10)

    def test_zero_score_is_not_calibration(self):
        rows,_=validate_diagnostic(self.frame())
        m=rows[rows.method.eq('marker_matched')]
        self.assertFalse(m.alert.any())
        self.assertTrue(m.high_realized_fdp.any())

    def test_insufficient_donors(self):
        rows,p=validate_diagnostic(self.frame().query('target_donor == "1"'))
        self.assertTrue(rows.empty)
        self.assertEqual(p.iloc[0].status,'INSUFFICIENT_DONORS')

    def test_checkpoint_corruption_fails(self):
        import sys
        sys.path.insert(0,str(Path(__file__).parent/'frozen'))
        from support import seal,checkpoint
        with tempfile.TemporaryDirectory() as name:
            folder=Path(name)
            (folder/'x.csv').write_text('a\n1\n')
            seal(folder,['x.csv'],'sig')
            self.assertTrue(checkpoint(folder,'sig'))
            (folder/'x.csv').write_text('a\n2\n')
            with self.assertRaisesRegex(ValueError,'Damaged'):
                checkpoint(folder,'sig')

    def test_null_scenario_survives_csv(self):
        with tempfile.TemporaryDirectory() as name:
            path=Path(name)/'metrics.csv.gz'
            pd.DataFrame([dict(outcome='null',score=np.nan),dict(outcome='up025',score=0.)]).to_csv(path,index=False,compression='gzip')
            table=read_nonempty(path)
            self.assertEqual(table.outcome.tolist(),['null','up025'])
            self.assertTrue(np.isnan(table.score.iloc[0]))


if __name__=='__main__':
    unittest.main()
