import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'engine')]
from mrn_tools import mrn_inputs,yield_convert,actual_span_hc3
from theory import log_input_bound,coefficient_bound,variance_decomposition
from extension_stats import fit_hc3
import revision_driver as D
import revision_core as C
import grid
from bulk_access import normalize_counts
from summarize import equal_donor


class CoreTests(unittest.TestCase):
    def setUp(self): self.rng=np.random.default_rng(14)

    def test_mrn_reference_yield(self):
        a=self.rng.uniform(.0001,.001,(30,7)); y=self.rng.uniform(1000,9000,7)
        p=self.rng.dirichlet(np.ones(7),8)
        rna=p*y; rna/=rna.sum(1,keepdims=True)
        marker_q=rna*a.sum(0); marker_q/=marker_q.sum(1,keepdims=True)
        self.assertTrue(np.allclose(yield_convert(marker_q,y*a.sum(0)),p))
        counts=self.rng.poisson(3,(8,30))
        out=mrn_inputs(counts,a,y,.02)
        self.assertTrue(np.allclose(out['reference_log2'],np.log2(20000*(a/a.sum(0)).T+1)))
        # Disjoint test counts never enter this API or change its marker input.
        self.assertTrue(np.array_equal(out['target_log2'],mrn_inputs(counts.copy(),a,y,.02)['target_log2']))

    def test_zero_library_rejected(self):
        with self.assertRaises(ValueError): mrn_inputs(np.zeros((2,30)),np.ones((30,7))*.001,np.ones(7),.01)

    def test_missing_donor_not_silently_dropped(self):
        f=pd.DataFrame({'group':['x']*4,'target_donor':['a','a','b','b'],'value':[1.,2.,np.nan,np.nan]})
        s,_=equal_donor(f,['group'])
        self.assertTrue(np.isnan(s.iloc[0]['mean']))
        self.assertEqual(s.iloc[0]['contributing_donor_mean'],1.5)
        self.assertEqual(s.iloc[0]['contributing_donors'],1)

    def test_finite_input_bound(self):
        for a in (.01,1.,1e4):
            x=self.rng.poisson(3,(40,20)); d=self.rng.uniform(10,100,40); e=d*np.exp(self.rng.uniform(-4,4,40))
            check=log_input_bound(x,d,e,pseudocount=a)
            self.assertLessEqual(check['max_violation'],1e-10)

    def test_api_matches_frozen_full_rank(self):
        x=self.rng.normal(size=80); y=self.rng.normal(size=(80,6)); p=self.rng.dirichlet(np.ones(7),80)
        new=actual_span_hc3(y,x,p); old=fit_hc3(y,x,p)
        for k in ('beta','se','exposure_residual_ss','exposure_r2'):
            self.assertTrue(np.allclose(new[k],old[k],atol=1e-12,rtol=1e-12),k)
        self.assertTrue(np.allclose(new['hc3_contributions'].sum(0),new['se']**2))

    def test_actual_span_lower_rank(self):
        x=self.rng.normal(size=80); y=self.rng.normal(size=(80,3)); p=self.rng.dirichlet(np.ones(3),80)
        p=np.column_stack([p,p[:,0],np.zeros((80,3))])
        new=actual_span_hc3(y,x,p); old=fit_hc3(y,x,p,rule='adaptive')
        self.assertTrue(new['exposure_estimable']); self.assertEqual(new['composition_rank'],2)
        self.assertTrue(np.allclose(new['beta'],old['beta']))

    def test_unidentifiable_exposure(self):
        x=np.repeat([0.,1.],40); y=self.rng.normal(size=(80,3))
        fit=actual_span_hc3(y,x,covariates=x[:,None])
        self.assertFalse(fit['exposure_estimable']); self.assertFalse(fit['available'].any())

    def test_coefficient_bound_and_variance(self):
        x=self.rng.normal(size=80); tau=np.array([0.,.2,.5]); y=x[:,None]*tau+self.rng.normal(size=(80,3))
        p=self.rng.dirichlet(np.ones(7),80); q=self.rng.dirichlet(np.ones(7),80)
        a,b=fit_hc3(y,x,p),fit_hc3(y,x,q)
        check=coefficient_bound(y,x,tau,a,b)
        self.assertTrue(check['eligible']); self.assertLess(check['max_violation'],1e-10)
        l=self.rng.normal(size=200); e=self.rng.normal(size=200)+.4*l
        v=variance_decomposition(l+e,l,e)
        self.assertLess(v['variance_identity_error'],1e-12)
        self.assertGreater(v['twice_covariance'],0)

    def test_checkpoint_integrity(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); (root/'a.txt').write_text('a')
            D.seal(root,['a.txt'],'fixed'); self.assertTrue(D.valid_checkpoint(root,'fixed'))
            (root/'a.txt').write_text('bad')
            with self.assertRaises(ValueError): D.valid_checkpoint(root,'fixed')

    def test_grid_marker_selection_and_baseline(self):
        with tempfile.TemporaryDirectory() as td:
            job=D.synthetic_jobs(Path(td))[0]; unit=C.plan_units([job],smoke=True)[0]
            self.assertEqual(len(grid.stress_grid()),11)
            selected=[set(grid.selected_rows(job,f,'gene_id')) for f in (.05,.1,.2)]
            self.assertTrue(selected[0]<=selected[1]<=selected[2])
            base=C.draw_base(job,unit,0,C.SEEDS)
            point=next(p for p in grid.stress_grid() if p['grid_id']=='f10_a050')
            new,ix=grid.stress_draw(job,base,point)
            from extension_driver import stress_additions
            old,oldix=stress_additions(job,base)
            self.assertTrue(np.array_equal(new,old)); self.assertEqual(ix.tolist(),oldix)
            data=C.intervene(job,base,C.conditions(True)[1],C.SEEDS)
            C.verify_baseline(job,unit,0,C.SEEDS,data,.5)

    def test_bulk_counts_validate_and_refuse_transformation(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            counts=pd.DataFrame({'s1':np.ones(100),'s2':np.ones(100)},index=[f'ENSG{i:011d}' for i in range(100)])
            counts.to_csv(root/'counts.csv')
            pd.DataFrame({'SampleID':['s1','s2'],'included':[True,True]}).to_csv(root/'cohort.csv',index=False)
            result=normalize_counts(root/'counts.csv',root/'cohort.csv',root/'out.csv')
            self.assertEqual(result['cohort_samples'],2)
            counts.iloc[0,0]=.5; counts.to_csv(root/'counts.csv')
            with self.assertRaises(ValueError): normalize_counts(root/'counts.csv',root/'cohort.csv',root/'out.csv')
            counts.to_csv(root/'processed.csv')
            with self.assertRaises(ValueError): normalize_counts(root/'processed.csv',root/'cohort.csv',root/'out.csv')


if __name__=='__main__': unittest.main()
