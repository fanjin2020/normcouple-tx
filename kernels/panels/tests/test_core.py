import os
from io import StringIO
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "engine")]
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
from extension_driver import stress_additions
from selection import CASES, METHODS, freeze, stress, subsets
from summarize import mc, read_metrics, validate_frame, wilson
from support import atomic_json, checkpoint, digest, file_sha, seal, task_work, tasks_for


class SelectionTests(unittest.TestCase):
    def test_fixed_stratified_sizes(self):
        genes = np.array(["g%03d" % i for i in range(300)])
        owner = np.arange(300) % 7 + 1
        first = subsets(genes, owner, "panel")
        second = subsets(genes, owner, "panel")
        self.assertEqual(len(first), 52)
        for key, indices in first.items():
            np.testing.assert_array_equal(indices, second[key])
            self.assertEqual(len(indices), len(set(indices)))
            if key == "zero":
                self.assertEqual(len(indices), 0)
            else:
                for k in range(1, 8):
                    self.assertEqual(sum(owner[indices] == k), int(np.ceil(sum(owner == k) * .1)))
        self.assertGreater(len({tuple(first[k]) for k in first if k.startswith("random")}), 45)

    def test_row_order_invariance(self):
        genes = np.array(["g%03d" % i for i in range(300)])
        owner = np.arange(300) % 7 + 1
        a = subsets(genes, owner, "fixed-panel", 2)
        b = subsets(genes[::-1], owner[::-1], "fixed-panel", 2)
        for key in a:
            self.assertEqual(set(genes[a[key]]), set(genes[::-1][b[key]]))

    def test_missing_class_rejected(self):
        with self.assertRaises(ValueError):
            subsets(np.array(["a", "b"]), np.array([1, 2]), "bad")


class CheckpointTests(unittest.TestCase):
    def test_partial_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "x.csv").write_text("a\n1\n")
            self.assertIsNone(checkpoint(root, "abc"))
            seal(root, ["x.csv"], "abc", work_items=1)
            self.assertIsNotNone(checkpoint(root, "abc"))
            with self.assertRaises(ValueError):
                checkpoint(root, "other")
            (root / "x.csv").write_text("a\n2\n")
            with self.assertRaises(ValueError):
                checkpoint(root, "abc")

    def test_task_partition(self):
        cfg = dict(settings=dict(repetitions=11, batch_size=5, subset_batch_size=5,
                                subset_ids=list(range(52))), units=[dict(unit_id="u")])
        tasks = list(tasks_for(cfg))
        self.assertEqual(len(tasks), 33)
        self.assertEqual(sum(task_work(cfg, t) for t in tasks), 11 * 52)

    def test_signature_ignores_dict_order(self):
        self.assertEqual(digest(dict(a=1,b=2)), digest(dict(b=2,a=1)))


class MetricTests(unittest.TestCase):
    def test_identical_paired_series(self):
        result = mc(np.zeros(200))
        self.assertEqual(result["mean"], 0)
        self.assertEqual(result["mc_se"], 0)
        self.assertEqual(result["finite_repetitions"], 200)
        self.assertTrue(np.isnan(mc([np.nan])["mean"]))

    def test_wilson_boundaries(self):
        lo, hi = wilson(0, 50)
        self.assertAlmostEqual(lo, 0)
        self.assertGreater(hi, 0)
        lo, hi = wilson(50, 50)
        self.assertLess(lo, 1)
        self.assertAlmostEqual(hi, 1)

    def test_complete_grid_required(self):
        cfg = dict(settings=dict(subset_ids=["zero"], repetitions=2))
        rows = []
        for case, _ in CASES:
            for method in METHODS:
                for rep in range(2):
                    rows.append(dict(subset_id="zero",outcome=case,method=method,rep=rep,
                        availability=1.,valid_tests=100,planned_tests=100,covered_tests=95,
                        delivered_coverage=.95,fdp=.1,false_discoveries=1,discoveries=10,
                        power=np.nan if case=="null" else .5))
        frame = pd.DataFrame(rows)
        frame = read_metrics(StringIO(frame.to_csv(index=False)))
        self.assertEqual(int(frame.outcome.eq("null").sum()), 8)
        self.assertEqual(len(validate_frame(frame,cfg,dict(unit_id="fixture"))),24)
        with self.assertRaises(ValueError):
            validate_frame(frame.iloc[:-1],cfg,dict(unit_id="fixture"))
        with self.assertRaises(ValueError):
            validate_frame(pd.concat([frame,frame.iloc[:1]]),cfg,dict(unit_id="fixture"))


class FrozenFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.jobs = D.synthetic_jobs(Path(cls.temp.name))
        freeze(cls.jobs, 2)
        cls.job = cls.jobs[0]
        cls.unit = C.plan_units(cls.jobs, smoke=True)[0]
        cls.base = C.draw_base(cls.job, cls.unit, 0, C.SEEDS)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_legacy_rng_and_counts_unchanged(self):
        new, _ = stress(self.job, self.base, "legacy")
        old, selected = stress_additions(self.job, self.base)
        np.testing.assert_array_equal(old, new)
        np.testing.assert_array_equal(selected, self.job["stress_subsets"]["legacy"])

    def test_zero_and_random_reproducibility(self):
        zero, _ = stress(self.job, self.base, "zero")
        self.assertEqual(zero.sum(), 0)
        a, detail = stress(self.job, self.base, "random_001")
        b, _ = stress(self.job, self.base, "random_001")
        np.testing.assert_array_equal(a, b)
        self.assertEqual(a[self.base["group"] == 0].sum(), 0)
        self.assertGreater(detail["selected_markers"], 0)

    def test_frozen_outcomes(self):
        for name, effect in CASES:
            data = C.intervene(self.job,self.base,dict(name=name,effect=effect,pattern="up",panel="frozen"),C.SEEDS)
            if name != "up025":
                C.verify_baseline(self.job,self.unit,0,C.SEEDS,data,effect)
            self.assertEqual(int((data["truth"] != 0).sum()), 0 if not effect else 20)


if __name__ == "__main__":
    unittest.main()
