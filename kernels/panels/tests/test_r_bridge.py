"""Optional real-R equivalence test against the reviewed uncached implementation."""
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "engine")]
import numpy as np
import pandas as pd
import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from selection import CASES, freeze, stress


@unittest.skipUnless(os.environ.get("PAPER3_TEST_RSCRIPT") and os.environ.get("PAPER3_TEST_RLIB"),
                     "Set PAPER3_TEST_RSCRIPT and PAPER3_TEST_RLIB for real-R equivalence checks")
class BridgeEquivalenceTests(unittest.TestCase):
    def test_cached_predictions_match_frozen_uncached_bridge(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            jobs = D.synthetic_jobs(root / "inputs")
            freeze(jobs, 1)
            job = jobs[0]
            unit = C.plan_units(jobs, smoke=True)[0]
            base = C.draw_base(job, unit, 0, C.SEEDS)
            A, _, anchor = marker_properties(job, base)
            new = root / "new"
            old = root / "old"
            new.mkdir()
            old.mkdir()
            C.write_bin(new / "raw.bin", base["draw"]["counts"][:, job["norm_genes"]].T)
            C.write_bin(new / "A.bin", A)
            C.write_bin(new / "total.bin", base["draw"]["mass"].sum(1))
            outcomes, selections, conditions = [], [], []
            for name, effect in CASES:
                data = C.intervene(job, base, dict(name=name, effect=effect, pattern="up", panel="frozen"), C.SEEDS)
                C.write_bin(new / (name + ".bin"), data["additions"].T)
                outcomes.append(dict(name=name, additions=name + ".bin"))
            for subset_id in ("zero", "legacy", "random_001"):
                added, _ = stress(job, base, subset_id)
                C.write_bin(new / (subset_id + ".bin"), added.T)
                selections.append(dict(name=subset_id, stress=subset_id + ".bin"))
                for outcome in outcomes:
                    conditions.append(dict(name=subset_id + "_" + outcome["name"],
                        additions=outcome["additions"], stress=subset_id + ".bin"))
            for path in new.iterdir():
                shutil.copy2(path, old / path.name)
            shared = dict(key="r0", n=unit["n"], raw="raw.bin", A="A.bin", total="total.bin", marker_anchor=anchor)
            C.write_json(new / "jobs.json", dict(reference_folder=job["reference_folder"],
                         entries=[dict(shared, outcomes=outcomes, subsets=selections)]))
            C.write_json(old / "jobs.json", dict(reference_folder=job["reference_folder"],
                         entries=[dict(shared, conditions=conditions)]))
            r = os.environ["PAPER3_TEST_RSCRIPT"]
            lib = os.environ["PAPER3_TEST_RLIB"]
            D.run_r([r,"--vanilla",str(ROOT / "bridge.R"),str(new),lib],lib,root / "new.log",timeout=180)
            D.run_r([r,"--vanilla",str(ROOT / "tests" / "fixtures" / "legacy_grid_bridge.R"),str(old),lib],
                    lib,root / "old.log",timeout=180)
            new_results = pd.read_csv(new / "outputs.csv", keep_default_na=False)
            old_results = pd.read_csv(old / "outputs.csv", keep_default_na=False).set_index(["key", "method"])
            self.assertEqual(len(new_results), 27)
            for row in new_results.itertuples():
                other = old_results.loc[(row.key, row.method)]
                self.assertEqual(row.reason, other.reason)
                a = np.loadtxt(new / row.file, delimiter=",")
                b = np.loadtxt(old / other.file, delimiter=",")
                np.testing.assert_allclose(a, b, atol=1e-12, rtol=1e-12, equal_nan=True)
            checks = pd.read_csv(new / "cache_checks.csv")
            self.assertTrue(checks.marker_input_identical.all())
            self.assertEqual(int(checks.fresh_repeat_checked.sum()), 2)
            checked = checks[checks.fresh_repeat_checked & checks.prediction_finite]
            self.assertTrue(checked.fresh_repeat_max_difference.le(1e-12).all())


if __name__ == "__main__":
    unittest.main()
