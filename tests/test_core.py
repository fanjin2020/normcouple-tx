import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mrn_tools',ROOT/'publication/mrn_tools.py')
M=importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


class CoreTests(unittest.TestCase):
    def test_mrn_definition(self):
        counts=np.array([[1,4,2],[2,1,7]],float)
        ref=np.array([[.05,.08],[.08,.02],[.01,.04]])
        yield_=np.array([30.,50.])
        result=M.mrn_inputs(counts,ref,yield_,.2)
        np.testing.assert_allclose(result['target_log2'],np.log2(2e5*counts/counts.sum(1)[:,None]+1))
        np.testing.assert_allclose(result['marker_yield'],yield_*ref.sum(0))

    def test_mrn_does_not_take_test_gene_counts(self):
        markers=np.array([[5,2],[3,4]])
        ref=np.array([[.2,.1],[.1,.2]])
        before=M.mrn_inputs(markers,ref,[20,30],.3)
        unrelated_counts=np.array([[2,9],[4,12]])
        unrelated_counts[1]*=4
        after=M.mrn_inputs(markers,ref,[20,30],.3)
        np.testing.assert_array_equal(before['target_log2'],after['target_log2'])

    def test_empty_marker_library_rejected(self):
        with self.assertRaises(ValueError):
            M.mrn_inputs([[0,0]],[[.2,.1],[.1,.2]],[20,30],.3)

    def test_yield_conversion(self):
        p=M.yield_convert([[.25,.75]],[10,30])
        np.testing.assert_allclose(p,[[.5,.5]])
        with self.assertRaises(ValueError):
            M.yield_convert([[.2,.3]],[10,30])

    def test_hc3_matches_direct_full_rank_ols(self):
        rng=np.random.default_rng(123)
        n=80
        x=np.r_[np.zeros(n//2),np.ones(n//2)]
        c=rng.normal(size=(n,2))
        y=.3*x[:,None]+c@rng.normal(size=(2,3))+rng.normal(size=(n,3))
        fit=M.actual_span_hc3(y,x,covariates=c)
        design=np.column_stack([np.ones(n),c,x])
        pinv=np.linalg.pinv(design)
        beta=pinv@y
        residual=y-design@beta
        h=np.diag(design@pinv)
        se=np.sqrt(((pinv[-1,:,None]*residual/(1-h[:,None]))**2).sum(0))
        np.testing.assert_allclose(fit['beta'],beta[-1],atol=1e-12)
        np.testing.assert_allclose(fit['se'],se,atol=1e-12)
        self.assertTrue(fit['available'].all())

    def test_duplicate_nuisance_columns_keep_same_fit(self):
        rng=np.random.default_rng(123)
        x=rng.normal(size=50)
        c=rng.normal(size=(50,2))
        y=rng.normal(size=(50,2))
        a=M.actual_span_hc3(y,x,covariates=c)
        b=M.actual_span_hc3(y,x,covariates=np.column_stack([c,c[:,0]]))
        np.testing.assert_allclose(a['beta'],b['beta'],atol=1e-12)
        np.testing.assert_allclose(a['se'],b['se'],atol=1e-12)

    def test_aliased_exposure_unavailable(self):
        x=np.arange(20,dtype=float)
        fit=M.actual_span_hc3(np.sin(x)[:,None],x,covariates=x[:,None])
        self.assertFalse(fit['exposure_estimable'])
        self.assertFalse(fit['available'].any())

    def test_input_derivative(self):
        x,d,c,a=7.,100.,1e3,1.
        h=1e-5
        observed=(np.log2(c*x/(d*np.exp(h))+a)-np.log2(c*x/(d*np.exp(-h))+a))/(2*h)
        expected=-1/np.log(2)*(c*x)/(c*x+a*d)
        self.assertAlmostEqual(observed,expected,places=9)


if __name__=='__main__':
    unittest.main()
