"""Regression checks for literal null labels, finite-case weighting and repair."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import revision_summary as S
import repair_revision as R


def fixtures(unit, repetitions=4):
    records, pairs, genes, mechanism = [], [], [], []
    for condition in [item["name"] for item in unit["conditions"]]:
        for method in sorted(S.methods_for(condition)):
            for rep in range(repetitions):
                row = dict(condition=condition, method=method, rep=rep, reason="", n_tests=100,
                           n_alternatives=0 if condition == "null" else 20, complete_fit=1)
                row.update({metric: .1 for metric in S.METRICS})
                row.update(availability=1., null_availability=1.)
                if condition == "null":
                    row["power"] = np.nan
                records.append(row)
                if condition != "null":
                    mechanism.append(dict(condition=condition, method=method, rep=rep, proportion_max_change=0.))
            for gene in range(100):
                genes.append(dict(condition=condition, method=method, gene_id=f"g{gene}",
                    truth=0. if condition == "null" or gene >= 20 else .5,
                    planned=repetitions, available=repetitions, bias=.01, empirical_sd=.1,
                    mean_se=.1, rms_se=.1, mean_leakage=.01, mean_noise_projection=0.,
                    coverage=.95, delivered_coverage=.95, availability=1.))
        keys, _ = S.expected_index(dict(conditions=[dict(name=condition)]), range(repetitions), "paired_replicates")
        for cond, a, b, metric in keys:
            denominator = "all_planned" if metric in S.ALL_PLANNED else (
                "common_complete_fits" if metric in S.CONDITIONAL else "common_available_genes")
            for rep in range(repetitions):
                pairs.append(dict(condition=cond, method_a=a, method_b=b, metric=metric,
                    denominator=denominator, rep=rep, delta=np.nan if condition == "null" and metric == "power" else 0.))
    return {"replicates": pd.DataFrame(records), "paired_replicates": pd.DataFrame(pairs),
            "per_gene_summary": pd.DataFrame(genes), "mechanism_replicates": pd.DataFrame(mechanism)}


def sealed(folder, names, signature):
    R.write_json(folder / "complete.json", dict(signature=signature, files={name: R.sha(folder / name) for name in names}))


class SummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.unit = dict(unit_id="revision_123abc", reference_source="TEST", target_source="TEST",
                         target_donor="A", conditions=[dict(name="null"), dict(name="up05")])

    def tearDown(self):
        self.temp.cleanup()

    def test_reader_preserves_null_and_na_labels_but_missing_numeric_values(self):
        path = self.root / "data.csv"
        path.write_text("condition,target_donor,value\nnull,NA,\nup05,NULL,2\n")
        frame = S.read_result_csv(path)
        self.assertEqual(frame.condition.tolist(), ["null", "up05"])
        self.assertEqual(frame.target_donor.tolist(), ["NA", "NULL"])
        self.assertTrue(np.isnan(frame.value.iloc[0]))

    def test_method_counts(self):
        self.assertEqual(len(S.methods_for("null")), 24)
        self.assertEqual(len(S.methods_for("up05")), 24)
        self.assertEqual(len(S.methods_for("up01")), 14)

    def test_complete_data_matches_original_paired_mc(self):
        matrix = np.array([[1., 2.], [4., 3.], [8., 5.], [3., 4.]])
        frame = pd.DataFrame(dict(rep=np.repeat(range(4), 2), target_donor=["A", "B"] * 4, x=matrix.ravel()))
        result = S.donor_equal_mc(frame, "x", ["A", "B"], 4)
        self.assertAlmostEqual(result["mean"], matrix.mean())
        self.assertAlmostEqual(result["mcse"], matrix.mean(axis=1).std(ddof=1) / 2)

    def test_unequal_success_rates_keep_equal_donor_weight_and_ratio_influence(self):
        matrix = np.array([[1., 8.], [2., 12.], [3., np.nan], [4., np.nan]])
        frame = pd.DataFrame(dict(rep=np.repeat(range(4), 2), target_donor=["A", "B"] * 4, x=matrix.ravel()))
        result = S.donor_equal_mc(frame, "x", ["A", "B"], 4)
        self.assertAlmostEqual(result["mean"], 6.25)
        self.assertNotAlmostEqual(result["mean"], np.nanmean(matrix, axis=1).mean())
        influence = (np.array([-1.5, -.5, .5, 1.5]) + np.array([-4., 4., 0., 0.])) / 2
        self.assertAlmostEqual(result["mcse"], influence.std(ddof=1) / 2)
        self.assertEqual(result["min_available_repetitions_per_donor"], 2)

    def test_missing_entire_donor_does_not_silently_drop_it(self):
        frame = pd.DataFrame(dict(rep=range(4), target_donor=["A"] * 4, x=[1.] * 4))
        result = S.donor_equal_mc(frame, "x", ["A", "B"], 4)
        self.assertTrue(np.isnan(result["mean"]))
        self.assertEqual(result["contributing_donors"], 1)
        self.assertEqual(result["summary_status"], "donor_without_finite_values")

    def test_one_valid_observation_has_no_mc_interval(self):
        frame = pd.DataFrame(dict(rep=[0, 1], target_donor=["A", "A"], x=[3., np.nan]))
        result = S.donor_equal_mc(frame, "x", ["A"], 2)
        self.assertEqual(result["mean"], 3.)
        self.assertTrue(np.isnan(result["mcse"]))

    def test_missing_condition_and_duplicate_reps_are_rejected(self):
        frame = fixtures(self.unit)["replicates"]
        S.validate_unit_table(frame, self.unit, list(range(4)), "replicates")
        with self.assertRaisesRegex(ValueError, "Missing/unexpected"):
            S.validate_unit_table(frame[frame.condition.ne("null")], self.unit, list(range(4)), "replicates")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            S.validate_unit_table(pd.concat([frame, frame.iloc[:1]]), self.unit, list(range(4)), "replicates")

    def test_gene_label_recovery_is_guarded_by_plan_truth_and_family(self):
        genes = fixtures(self.unit)["per_gene_summary"]
        genes.loc[genes.condition.eq("null"), "condition"] = np.nan
        fixed, count = S.recover_gene_null_label(genes, self.unit, 4)
        self.assertEqual(count, 2400)
        self.assertEqual(set(fixed.condition), {"null", "up05"})
        genes.loc[genes.condition.isna(), "truth"] = .5
        with self.assertRaisesRegex(ValueError, "nonzero-truth"):
            S.recover_gene_null_label(genes, self.unit, 4)

    def make_old_result(self):
        source = self.root / "old"
        source.mkdir()
        signature = "a" * 64
        R.write_json(source / "protocol.json", dict(version="1.0.1", signature=signature, units=[self.unit],
            settings=dict(repetitions=4, batch_size=2, software_smoke=False, synthetic=False)))
        R.write_json(source / "RUN_STATUS.json", dict(status="COMPLETE", signature=signature, units=1))
        frames = fixtures(self.unit)
        folder = source / "units" / self.unit["unit_id"]
        folder.mkdir(parents=True)
        for start in (0, 2):
            chunk = folder / f"chunk_{start:05d}"
            chunk.mkdir()
            for name in R.CHUNK_FILES:
                (chunk / name).write_text("TEST FIXTURE ONLY\n")
            for name in ("replicates", "paired_replicates", "mechanism_replicates"):
                frame = frames[name]
                S.write_csv(chunk / (name + ".csv"), frame[frame.rep.isin([start, start + 1])])
            R.write_json(chunk / "baseline_checks.json", [{"rep": start, "fixture": True}])
            R.write_json(chunk / "r_controls.json", {"fixture": True})
            (chunk / "bridge.log").write_text("")
            sealed(chunk, R.CHUNK_FILES, signature)
        for name, frame in frames.items():
            frame = frame.copy()
            frame.loc[frame.condition.eq("null"), "condition"] = np.nan
            S.write_csv(folder / (name + ".csv"), frame)
        R.write_json(folder / "unit.json", self.unit)
        sealed(folder, R.UNIT_FILES, signature)
        return source

    def test_summary_only_repair_restores_all_null_groups_and_preserves_source(self):
        source = self.make_old_result()
        before = {p.relative_to(source): R.sha(p) for p in source.rglob("*") if p.is_file()}
        output = self.root / "fixed"
        with contextlib.redirect_stdout(io.StringIO()):
            R.repair(source, output, plots=False)
        self.assertEqual(before, {p.relative_to(source): R.sha(p) for p in source.rglob("*") if p.is_file()})
        status = S.read_result_csv(output / "tables" / "method_status.csv")
        self.assertEqual(set(status.condition), {"null", "up05"})
        self.assertEqual(len(status), 48)
        self.assertTrue(status.planned_replicates.eq(4).all())
        self.assertEqual(R.read_json(output / "SUMMARY_PROVENANCE.json")["recovered_gene_labels"], 2400)
        self.assertFalse((output / "tables" / "replicates.csv").exists())
        for item in R.read_json(output / "MANIFEST.json")["files"]:
            self.assertEqual(R.sha(output / item["path"]), item["sha256"])

    def test_repair_refuses_corrupt_chunk_or_existing_output(self):
        source = self.make_old_result()
        output = self.root / "fixed"
        output.mkdir()
        with self.assertRaisesRegex(ValueError, "Output exists"):
            R.repair(source, output, plots=False)
        chunk = source / "units" / self.unit["unit_id"] / "chunk_00000" / "replicates.csv"
        chunk.write_text("corrupt")
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "Damaged source"):
            R.repair(source, self.root / "different", plots=False)
        self.assertEqual(R.read_json(self.root / "different" / "EXPORT_STATUS.json")["status"], "INCOMPLETE")


if __name__ == "__main__":
    unittest.main(verbosity=2)
