"""Numerical/provenance tests. R integration is `run_revision.py --synthetic`."""
import copy
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from revision_core import E, M, SEEDS, conditions, plan_units, draw_base, intervene, verify_baseline, write_json, prepare_job
from revision_driver import synthetic_jobs, valid_checkpoint, seal, summarize_genes
from revision_metrics import evaluate, pair_one, BASELINE


class RevisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        cls.root=Path(cls.temp.name)
        cls.job=synthetic_jobs(cls.root)[0]
        cls.unit=plan_units([cls.job],True)[0]
        cls.base=draw_base(cls.job,cls.unit,0,SEEDS)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_frozen_generator(self):
        for condition in conditions()[:4]:
            data=intervene(self.job,self.base,condition,SEEDS)
            checks=verify_baseline(self.job,self.unit,0,SEEDS,data,condition["effect"])
            self.assertLess(max(checks.values()),1e-10)

    def test_external_generator(self):
        job=dict(self.job,target_source=self.job["source"],alt_mask=np.arange(100)<20)
        base=draw_base(job,self.unit,0,SEEDS)
        for condition in conditions(True):
            data=intervene(job,base,condition,SEEDS)
            checks=verify_baseline(job,self.unit,0,SEEDS,data,condition["effect"])
            self.assertLess(max(checks.values()),1e-10)

    def test_balanced_intervention(self):
        condition=next(c for c in conditions() if c["name"]=="balanced05")
        data=intervene(self.job,self.base,condition,SEEDS)
        self.assertTrue((data["additions"]>=0).all())
        self.assertEqual((data["truth"]>0).sum(),10)
        self.assertEqual((data["truth"]<0).sum(),10)
        mass=data["expected_added_mass"]
        self.assertAlmostEqual(mass[data["group"]==0].mean(),mass[data["group"]==1].mean())
        np.testing.assert_allclose(data["y"],data["nuisance"]+data["group"][:,None]*data["truth"]+data["noise"])

    def test_gene_selection_and_panel_disjointness(self):
        j=self.job
        self.assertFalse(set(j["mead_genes"]) & set(j["test"]))
        self.assertFalse(set(j["panels"]["balanced"]) & set(j["test"]))
        mean=j["base"]["mean_counts"].mean(1)[j["test"]]
        self.assertLessEqual(mean[j["injections"]["low"]].max(),mean[j["injections"]["high"]].min())

    def test_truth_and_null_outcome_invariance(self):
        null=intervene(self.job,self.base,conditions()[0],SEEDS)
        for condition in conditions()[1:]:
            d=intervene(self.job,self.base,condition,SEEDS)
            mask=d["truth"]==0
            np.testing.assert_array_equal(d["y"][:,mask],null["y"][:,mask])
            np.testing.assert_allclose(d["mixture"]*d["total"][:,None],null["mixture"]*null["total"][:,None])
            np.testing.assert_array_equal(d["p"],null["p"])

    def test_nnls_scaling_invariance(self):
        data=intervene(self.job,self.base,conditions()[1],SEEDS)
        a,_,_=E.repaired_prediction(data["A"],data["mixture"],"drop_zero_weighted")
        b,_,_=E.repaired_prediction(data["A"],data["mixture"]*(data["total"]/data["original_total"])[:,None],"drop_zero_weighted")
        np.testing.assert_allclose(a,b,atol=1e-10,rtol=1e-9)

    def test_fwl_identity(self):
        for condition in conditions():
            d=intervene(self.job,self.base,condition,SEEDS)
            row,arr=evaluate(d,d["p"])
            self.assertEqual(row["complete_fit"],1)
            self.assertLess(row["identity_error"],1e-8)
            np.testing.assert_allclose(arr["beta"]-arr["truth"],arr["leakage"]+arr["noise_projection"],atol=1e-8)

    def test_unavailable_pvalues_keep_family(self):
        d=intervene(self.job,self.base,conditions()[1],SEEDS)
        row,arr=evaluate(d,np.full((80,7),np.nan),"deliberate_failure")
        self.assertEqual(row["availability"],0)
        self.assertEqual(row["delivered_coverage"],0)
        self.assertTrue(np.isnan(row["coverage"]))
        self.assertTrue((arr["p"]==1).all() and (arr["q"]==1).all())
        self.assertEqual(len(arr["p"]),100)

    def test_method_pair_keeps_failure(self):
        d=intervene(self.job,self.base,conditions()[1],SEEDS)
        good,g=evaluate(d,d["p"])
        bad,b=evaluate(d,np.full((80,7),np.nan))
        pairs=pair_one({BASELINE:good,"failed":bad},{BASELINE:g,"failed":b})
        self.assertEqual(next(p["delta"] for p in pairs if p["metric"]=="availability"),-1)
        self.assertTrue(np.isnan(next(p["delta"] for p in pairs if p["metric"]=="coverage")))

    def test_resume_integrity(self):
        path=self.root/"checkpoint";path.mkdir(exist_ok=True)
        (path/"data.csv").write_text("a,b\n1,2\n")
        seal(path,["data.csv"],"sig")
        self.assertTrue(valid_checkpoint(path,"sig"))
        with self.assertRaises(ValueError):valid_checkpoint(path,"changed")
        (path/"data.csv").write_text("damaged")
        with self.assertRaises(ValueError):valid_checkpoint(path,"sig")

    def test_condition_plan(self):
        self.assertEqual(len(conditions()),7)
        self.assertEqual([c["name"] for c in conditions(True)],["null","up05"])

    def test_draw_reproducibility(self):
        b=draw_base(self.job,self.unit,0,SEEDS)
        np.testing.assert_array_equal(b["draw"]["counts"],self.base["draw"]["counts"])
        self.assertEqual(b["allocation"],self.base["allocation"])

    def test_groupby_split_internal_and_external_seed_parity(self):
        frame = pd.DataFrame({"source": ["SYNTHETIC"], "split": np.array([0], dtype=np.int64)})
        (_, split), _ = next(iter(frame.groupby(["source", "split"], sort=True)))
        self.assertIsInstance(split, np.integer)
        for external in (False, True):
            with self.subTest(external=external):
                clean = dict(self.job)
                if external:
                    clean.update(target_source=self.job["source"], alt_mask=np.arange(100)<20)
                dirty = dict(clean, split=split)
                fixed = prepare_job(dirty, self.root, len(self.job["mead_genes"]))
                self.assertIs(type(fixed["split"]), int)
                unit = plan_units([fixed], True)[0]
                self.assertEqual(unit, plan_units([clean], True)[0])
                for rep in (0, 1):
                    expected = draw_base(clean, unit, rep, SEEDS)
                    actual = draw_base(fixed, unit, rep, SEEDS)
                    self.assertEqual(M.digest(actual["key"]), M.digest(expected["key"]))
                    self.assertEqual(actual["allocation"], expected["allocation"])
                    for name in ("counts", "mass", "cell_truth"):
                        np.testing.assert_array_equal(actual["draw"][name], expected["draw"][name])
                    for condition in conditions(True):
                        data = intervene(fixed, actual, condition, SEEDS)
                        checks = verify_baseline(fixed, unit, rep, SEEDS, data, condition["effect"])
                        self.assertLess(max(checks.values()), 1e-10)


if __name__=="__main__":
    unittest.main(verbosity=2)
