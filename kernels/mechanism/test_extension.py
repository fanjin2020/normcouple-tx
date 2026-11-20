"""Focused tests; artificial data are software tests, not research evidence."""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from revision_core import M, SEEDS, draw_base, intervene, conditions, plan_units
from revision_driver import synthetic_jobs
from extension_stats import fit_hc3, evaluate, mechanism_checks
from extension_driver import marker_properties, stress_additions, conditions_for
from real_bulk import audit_samples, design


class ExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parent=Path(__file__).resolve().parent/"test_work"
        parent.mkdir(exist_ok=True)
        cls.tmp=tempfile.TemporaryDirectory(dir=parent)
        cls.job=synthetic_jobs(Path(cls.tmp.name))[0]
        cls.unit=plan_units([cls.job],True)[0]
        cls.base=draw_base(cls.job,cls.unit,0,SEEDS)
        cls.data=intervene(cls.job,cls.base,conditions()[1],SEEDS)

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def test_mechanism(self):
        self.assertLess(max(mechanism_checks().values()),1e-8)

    def test_same_full_rank_fit(self):
        d=self.data;p=d["p"]
        old=M.model_fit(d["y"],d["group"],p,d["nuisance"],d["noise"],d["truth"])
        for rule in ("full_six","adaptive"):
            f=fit_hc3(d["y"],d["group"],p,rule=rule)
            np.testing.assert_allclose(f["beta"],old["beta"],atol=1e-10)
            np.testing.assert_allclose(f["se"],old["se"],atol=1e-10)
            self.assertEqual(f["df"],old["df"])

    def test_low_rank_exposure_still_estimable(self):
        rng=np.random.default_rng(44);n=80
        t=rng.uniform(.1,.7,n);p=np.column_stack([t,1-t,np.zeros((n,5))])
        x=np.tile([0,1],40);y=rng.normal(size=(n,3))+x[:,None]
        full=fit_hc3(y,x,p,rule="full_six")
        adapt=fit_hc3(y,x,p,rule="adaptive")
        self.assertFalse(full["basic"]);self.assertTrue(adapt["basic"])
        self.assertEqual(adapt["composition_rank"],1);self.assertEqual(adapt["df"],77)
        X=np.column_stack([np.ones(n),x,t]);inv=np.linalg.pinv(X)
        residual=y-X@(inv@y);h=np.sum(X*inv.T,1)
        np.testing.assert_allclose(adapt["beta"],(inv@y)[1],atol=1e-12)
        np.testing.assert_allclose(adapt["se"],np.sqrt(np.sum((inv[1,:,None]*residual/(1-h[:,None]))**2,0)),atol=1e-12)

    def test_exposure_confounding_rejected(self):
        x=np.tile([0.,1.],40);p=np.column_stack([x,1-x,np.zeros((80,5))])
        f=fit_hc3(np.ones((80,3)),x,p)
        self.assertEqual(f["reason"],"exposure_not_estimable")

    def test_unavailable_denominator(self):
        row,a,_=evaluate(self.data,np.full((80,7),np.nan))
        self.assertEqual(row["availability"],0);self.assertEqual(row["delivered_coverage"],0)
        self.assertTrue(np.isnan(row["coverage"]));self.assertTrue((a["p"]==1).all())

    def test_marker_yield_conversion(self):
        A,sm,anchor=marker_properties(self.job,self.base)
        self.assertTrue(0<anchor<1)
        mu=self.base["ref"]["mean_counts"][self.job["panels"]["balanced"]]
        np.testing.assert_allclose(mu.sum(0),sm)
        p=self.data["p"];q=M.normalize(p*sm)
        np.testing.assert_allclose(M.normalize(q/sm),p)

    def test_marker_stress_disjoint(self):
        stress,selected=stress_additions(self.job,self.base)
        self.assertEqual(stress.shape,(80,300))
        self.assertTrue((stress[self.base["group"]==0]==0).all())
        self.assertFalse(set(self.job["test"]) & set(self.job["panels"]["balanced"][selected]))
        self.assertEqual(len(conditions_for(self.job)),9)

    def test_donor_not_cell_replicates(self):
        meta=pd.DataFrame({"SampleID":["a","b"],"Case.Year":[2000,2001],"Case.Num":[1,1],"Region":["FC","FC"],
                           "Diagnosis":["AD","Control"],"Age":[80,82],"Sex":["F","M"]})
        cells=pd.DataFrame({"SampleID":["a","a","b"],"Diagnosis":["AD","AD","Control"],"Age":[80,80,82],"Sex":["F","F","M"]})
        _,pairs=audit_samples(meta,cells)
        self.assertEqual(len(pairs),2);self.assertTrue(pairs.metadata_agree.all())
        cells.loc[0,"Diagnosis"]="Control"
        with self.assertRaises(ValueError):audit_samples(meta,cells)

    def test_fwl_saved(self):
        r,a,f=evaluate(self.data,self.data["p"],"adaptive")
        np.testing.assert_allclose(a["beta"]-a["truth"],a["leakage"]+a["noise_projection"],atol=1e-10)
        self.assertLess(r["identity_error"],1e-8)
        self.assertEqual(len(f["residualized_exposure"]),80)
        self.assertIn("design_condition_number",r)
        self.assertFalse({"condition","rep","method","rule"} & r.keys())

    def test_metadata_serialization(self):
        _,selected=stress_additions(self.job,self.base)
        self.assertEqual(json.loads(json.dumps(selected)),selected)

    def test_both_disease_groups_required(self):
        n=40
        meta=pd.DataFrame(dict(SampleID=["s"+str(i) for i in range(n)],donor_key=range(n),
            Region="FC",Diagnosis="AD",Age=80,Sex="F",PMI=4,RIN=7))
        meta["Sequencing.Group"]="a"
        with self.assertRaisesRegex(ValueError,"disease groups"):design(meta)


if __name__=="__main__":unittest.main(verbosity=2)
