"""Conditional Monte Carlo summaries, with donor and subset denominators explicit."""
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd
from scipy import stats

from selection import CASES, METHODS
from support import atomic_json, checkpoint, read_json, task_folder, tasks_for

HERE = Path(__file__).resolve().parent
METRICS = ["fdp", "power", "availability", "coverage", "delivered_coverage", "proportion_rmse",
           "bias", "abs_leakage", "mean_se", "exposure_r2", "exposure_residual_ss", "leverage"]
COUNTS = ["planned_tests", "valid_tests", "covered_tests"]
COMPARATORS = ["tmm", "marker_denominator_only"]
NAMES = {"fdp": "fdr", "coverage": "conditional_coverage", "proportion_rmse": "composition_rmse"}


def read_metrics(path):
    # pandas' default NA tokens include the literal scenario name "null".
    return pd.read_csv(path, keep_default_na=False, na_values=["", "NaN", "nan"])


def average(values, axis):
    count = np.isfinite(values).sum(axis)
    mean = np.divide(np.nansum(values, axis=axis), count,
                     out=np.full(count.shape, np.nan, dtype=float), where=count > 0)
    return mean, count


def mc(values):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    n = len(values)
    if not n:
        return dict(mean=np.nan, mc_se=np.nan, mc_low=np.nan, mc_high=np.nan, finite_repetitions=0)
    se = values.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
    half = stats.t.ppf(.975, n - 1) * se if n > 1 else np.nan
    return dict(mean=float(values.mean()), mc_se=float(se), mc_low=float(values.mean() - half),
                mc_high=float(values.mean() + half), finite_repetitions=n)


def wilson(success, total):
    if not total:
        return np.nan, np.nan
    z = stats.norm.ppf(.975)
    p = success / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
    return float(center - half), float(center + half)


class Accumulator:
    def __init__(self, shape):
        self.total = np.zeros(shape)
        self.count = np.zeros(shape, dtype=np.int32)

    def add(self, values):
        good = np.isfinite(values)
        self.total += np.where(good, values, 0)
        self.count += good

    def mean(self):
        return np.divide(self.total, self.count, out=np.full(self.total.shape, np.nan), where=self.count > 0)


def validate_frame(frame, cfg, unit):
    s = cfg["settings"]
    keys = ["subset_id", "outcome", "method", "rep"]
    expected = pd.MultiIndex.from_product([s["subset_ids"], [n for n, _ in CASES], METHODS,
                                           range(s["repetitions"])], names=keys)
    indexed = frame.set_index(keys)
    if not indexed.index.is_unique or len(indexed) != len(expected) or len(expected.difference(indexed.index)):
        raise ValueError("Incomplete/duplicated repetition grid: " + unit["unit_id"])
    frame = indexed.reindex(expected)
    if not np.allclose(frame.availability, frame.valid_tests / frame.planned_tests):
        raise ValueError("Availability denominator mismatch")
    if not np.allclose(frame.delivered_coverage, frame.covered_tests / frame.planned_tests):
        raise ValueError("Delivered coverage denominator mismatch")
    if not np.allclose(frame.fdp, frame.false_discoveries / np.maximum(frame.discoveries, 1)):
        raise ValueError("FDP denominator mismatch")
    if frame.loc[pd.IndexSlice[:, "null", :, :], "power"].notna().any():
        raise ValueError("Power must be undefined in global-null scenarios")
    return frame


def summarize(output, cfg):
    output = Path(output)
    s = cfg["settings"]
    tasks = list(tasks_for(cfg))
    for task in tasks:
        if not checkpoint(task_folder(output, task), cfg["signature"]):
            raise ValueError("Cannot publish partial summaries; unfinished batch: " + str(task_folder(output, task)))
    by_unit = {u["unit_id"]: [] for u in cfg["units"]}
    for task in tasks:
        by_unit[task[0]["unit_id"]].append(task_folder(output, task))
    settings, pairs, unit_rows, stress_rows, reasons = {}, {}, [], [], []
    oracle_error = identity_error = 0.0
    metric_rows = fresh_checks = 0
    shape = (len(s["subset_ids"]), len(CASES), len(METHODS), s["repetitions"], len(METRICS) + len(COUNTS))
    pair_shape = (shape[0], shape[1], len(COMPARATORS), shape[3], len(METRICS))
    for number, unit in enumerate(cfg["units"], 1):
        folders = by_unit[unit["unit_id"]]
        raw = pd.concat([read_metrics(p / "replicates.csv.gz") for p in folders], ignore_index=True)
        frame = validate_frame(raw, cfg, unit)
        metric_rows += len(raw)
        finite = raw.identity_error.dropna()
        if len(finite):
            identity_error = max(identity_error, float(finite.max()))
        values = frame[METRICS + COUNTS].to_numpy().reshape(shape)
        oracle = values[:, :, METHODS.index("oracle")]
        oracle_diff = np.abs(oracle - oracle[:1])
        if np.isfinite(oracle_diff).any():
            oracle_error = max(oracle_error, float(np.nanmax(oracle_diff)))
        if oracle_error > 1e-12:
            raise ValueError("Oracle outcomes changed across stress subsets")
        key = unit["reference_source"], unit["target_source"]
        if key not in settings:
            settings[key] = Accumulator(shape)
            pairs[key] = Accumulator(pair_shape)
        settings[key].add(values)
        deltas = np.stack([values[:, :, METHODS.index("marker_matched"), :, :len(METRICS)] -
                           values[:, :, METHODS.index(comp), :, :len(METRICS)] for comp in COMPARATORS], axis=2)
        pairs[key].add(deltas)
        means, counts = average(values, 3)
        attrs = {k: unit[k] for k in ("unit_id", "reference_source", "target_source", "target_donor")}
        for si, subset_id in enumerate(s["subset_ids"]):
            for ci, (outcome, effect) in enumerate(CASES):
                for mi, method in enumerate(METHODS):
                    row = dict(attrs, subset_id=subset_id, outcome=outcome, effect=effect, method=method,
                               planned_repetitions=s["repetitions"])
                    for ki, metric in enumerate(METRICS + COUNTS):
                        name = NAMES.get(metric, metric)
                        row[name] = means[si, ci, mi, ki]
                        row[name + "_finite_repetitions"] = counts[si, ci, mi, ki]
                    unit_rows.append(row)
        diag = pd.concat([pd.read_csv(p / "stress_diagnostics.csv.gz") for p in folders], ignore_index=True)
        if diag.duplicated(["rep", "subset_id"]).any() or len(diag) != shape[0] * shape[3]:
            raise ValueError("Missing/duplicate stress diagnostic")
        grouped = diag.drop(columns="rep").groupby("subset_id", sort=False).mean().reset_index()
        for k, value in attrs.items():
            grouped[k] = value
        stress_rows.extend(grouped.to_dict("records"))
        failing = raw[raw.reason.notna() & raw.reason.ne("")]
        if len(failing):
            counted = failing.groupby(["subset_id", "outcome", "method", "reason"]).size().rename("repetitions").reset_index()
            for k, value in attrs.items():
                counted[k] = value
            reasons.extend(counted.to_dict("records"))
        for folder in folders:
            checks = pd.read_csv(folder / "cache_checks.csv.gz")
            if not checks.marker_input_identical.all():
                raise ValueError("Cache input identity failed")
            fresh_checks += int(checks.fresh_repeat_checked.sum())
        print("SUMMARY unit %d/%d" % (number, len(cfg["units"])), flush=True)
    summary_rows, paired_rows = [], []
    for key, acc in settings.items():
        n_donors = sum((u["reference_source"], u["target_source"]) == key for u in cfg["units"])
        mean = acc.mean()
        # Conditional coverage is a ratio over available donor-tests within each replicate.
        available = acc.total[..., len(METRICS) + COUNTS.index("valid_tests")]
        covered = acc.total[..., len(METRICS) + COUNTS.index("covered_tests")]
        mean[..., METRICS.index("coverage")] = np.divide(covered, available,
                out=np.full(covered.shape, np.nan), where=available > 0)
        paired = pairs[key].mean()
        for si, subset_id in enumerate(s["subset_ids"]):
            for ci, (outcome, effect) in enumerate(CASES):
                attrs = dict(reference_source=key[0], target_source=key[1], subset_id=subset_id,
                             outcome=outcome, effect=effect, target_donors=n_donors,
                             planned_repetitions=s["repetitions"])
                for mi, method in enumerate(METHODS):
                    row = dict(attrs, method=method)
                    for ki, metric in enumerate(METRICS):
                        name = NAMES.get(metric, metric)
                        for suffix, value in mc(mean[si, ci, mi, :, ki]).items():
                            row[name + "_" + suffix] = value
                        row[name + "_mean_finite_donors"] = float(acc.count[si, ci, mi, :, ki].mean())
                    summary_rows.append(row)
                for mi, comp in enumerate(COMPARATORS):
                    row = dict(attrs, comparison="marker_matched_minus_" + comp)
                    for ki, metric in enumerate(METRICS):
                        name = NAMES.get(metric, metric)
                        for suffix, value in mc(paired[si, ci, mi, :, ki]).items():
                            row[name + "_" + suffix] = value
                        row[name + "_mean_common_donors"] = float(pairs[key].count[si, ci, mi, :, ki].mean())
                    paired_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    paired = pd.DataFrame(paired_rows)
    random = paired[paired.subset_id.str.startswith("random_")].copy()
    distribution = []
    keys = ["reference_source", "target_source", "outcome", "effect", "comparison"]
    for values, g in random.groupby(keys, sort=False):
        row = dict(zip(keys, values))
        valid = g.fdr_mean.notna()
        n = int(valid.sum())
        wins = int((g.loc[valid, "fdr_mean"] < 0).sum())
        low, high = wilson(wins, n)
        row.update(random_subsets=n, fdr_improved_sets=wins, fraction_fdr_improved=wins / n if n else np.nan,
                   subset_sampling_ci_low=low, subset_sampling_ci_high=high,
                   sets_paired_mc_upper_below_zero=int((g.fdr_mc_high < 0).sum()))
        for metric in ("fdr", "power", "availability", "delivered_coverage", "composition_rmse"):
            v = g[metric + "_mean"].dropna()
            row[metric + "_delta_median"] = v.median()
            row[metric + "_delta_q10"] = v.quantile(.1)
            row[metric + "_delta_q90"] = v.quantile(.9)
        match = summary[(summary.reference_source == row["reference_source"]) &
                        (summary.target_source == row["target_source"]) & (summary.outcome == row["outcome"]) &
                        (summary.method == "marker_matched") & summary.subset_id.str.startswith("random_")]
        row["mrn_sets_fdr_point_at_most_005"] = int((match.fdr_mean <= .05).sum())
        row["mrn_sets_fdr_mc_upper_at_most_005"] = int((match.fdr_mc_high <= .05).sum())
        distribution.append(row)
    stress_table = pd.DataFrame(stress_rows)
    mass = stress_table.groupby(["reference_source", "target_source", "subset_id"], as_index=False).mean(numeric_only=True)
    associations = []
    merged = random.merge(mass, on=["reference_source", "target_source", "subset_id"], validate="many_to_one")
    for values, g in merged.groupby(keys, sort=False):
        for x in ("selected_mean_count_per_case_marker", "selected_fraction_of_case_marker_mass",
                  "realized_added_fraction_of_case_library"):
            for y in ("fdr_mean", "power_mean"):
                finite = g[[x, y]].dropna()
                rho = stats.spearmanr(finite[x], finite[y]).statistic if len(finite) >= 3 and finite[x].nunique() > 1 and finite[y].nunique() > 1 else np.nan
                associations.append(dict(zip(keys, values), feature=x, metric=y, subsets=len(finite),
                                         spearman_rho=rho, interpretation="descriptive; not causal or mass-matched"))
    with tempfile.TemporaryDirectory(prefix="summary_build_", dir=output) as temporary:
        final = Path(temporary)
        summary.to_csv(final / "subset_method_results.csv", index=False)
        paired.to_csv(final / "paired_MRN_differences.csv", index=False)
        pd.DataFrame(distribution).to_csv(final / "random_subset_benefit_distribution.csv", index=False)
        pd.DataFrame(unit_rows).to_csv(final / "donor_subset_results.csv.gz", index=False, compression="gzip")
        stress_table.to_csv(final / "stress_expression_and_mass.csv", index=False)
        pd.DataFrame(associations).to_csv(final / "identity_mass_descriptive_associations.csv", index=False)
        pd.DataFrame(reasons, columns=["unit_id", "reference_source", "target_source", "target_donor", "subset_id",
                                     "outcome", "method", "reason", "repetitions"]).to_csv(final / "failure_counts.csv", index=False)
        for name in ("protocol.json", "selected_markers.json"):
            shutil.copy2(output / name, final / name)
        shutil.copy2(HERE / "README.md", final / "ANALYSIS_README.md")
        make_plots(summary, paired, final, s["software_fixture"])
        atomic_json(final / "QC.json", dict(software_fixture=s["software_fixture"], checksum_verified_batches=len(tasks),
            metric_rows=metric_rows, random_subsets=s["random_subsets"], reference_target_units=len(cfg["units"]),
            fresh_cache_refits_checked=fresh_checks, maximum_oracle_change_across_subsets=oracle_error,
            maximum_fwl_identity_error=identity_error, complete_grid=True,
            inference_scope="Monte Carlo conditional on fixed donor banks; subset CI conditional on realized paired simulations"))
        write_conclusion(summary, pd.DataFrame(distribution), final, s["software_fixture"])
        destination = output / "final_results"
        destination.mkdir(exist_ok=True)
        for path in final.iterdir():
            os.replace(path, destination / path.name)
    return destination


def write_conclusion(summary, distribution, final, synthetic):
    lines = ["# Randomized stress-subset robustness", "",
             "SOFTWARE FIXTURE ONLY. NOT RESEARCH EVIDENCE." if synthetic else "Complete pre-fixed research grid.", "",
             "Differences are MRN minus comparator. Negative FDR difference is relative improvement, not proof of FDR control.",
             "Positive power difference is a power gain. Global-null power is undefined.",
             "Subset intervals describe sampling of perturbed gene sets conditional on the realized simulations, not patients or disease prevalence.", ""]
    for row in distribution.itertuples():
        lines.append("- %s -> %s, %s, %s: FDR improved in %d/%d random subsets; median FDR difference %.4f; median power difference %.4f; MRN FDR point estimate <=0.05 in %d/%d subsets." %
                     (row.reference_source, row.target_source, row.outcome, row.comparison,
                      row.fdr_improved_sets, row.random_subsets, row.fdr_delta_median, row.power_delta_median,
                      row.mrn_sets_fdr_point_at_most_005, row.random_subsets))
    lines += ["", "Use the donor table to inspect heterogeneity. External targets reused under different references are not additional independent donors.",
              "Expression/count-mass associations are descriptive; this experiment does not isolate a mass-matched gene-identity effect.",
              "No real bulk validation, new theoretical theorem, public code archive, or wet experiment is delivered by this extension."]
    (final / "RESULTS_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_plots(summary, paired, final, synthetic):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    settings = summary[["reference_source", "target_source"]].drop_duplicates().itertuples(index=False, name=None)
    settings = list(settings)
    fig, axes = plt.subplots(len(settings), 3, figsize=(11, 2.8 * len(settings)), squeeze=False, constrained_layout=True)
    colors = ["#32668c", "#15836a", "#b97832"]
    methods = ["tmm", "marker_matched", "marker_denominator_only"]
    for i, (ref, target) in enumerate(settings):
        for j, (outcome, effect) in enumerate(CASES):
            ax = axes[i, j]
            g = summary[(summary.reference_source == ref) & (summary.target_source == target) &
                        (summary.outcome == outcome) & summary.subset_id.str.startswith("random_")]
            for k, method in enumerate(methods):
                v = g[g.method == method].fdr_mean.to_numpy()
                x = k + np.linspace(-.15, .15, len(v))
                ax.scatter(x, v, s=13, alpha=.65, color=colors[k])
                if len(v):
                    ax.plot([k - .18, k + .18], [np.median(v)] * 2, color="black", lw=1.5)
            ax.axhline(.05, color="#993c3c", ls="--", lw=1)
            ax.set(xticks=range(3), xticklabels=["TMM", "MRN", "Denom-only"], ylabel="Estimated FDR",
                   title=f"{ref} -> {target}; effect {effect:g}", ylim=(-.02, 1.02))
    if synthetic:
        fig.suptitle("Software fixture only")
    fig.savefig(final / "random_subset_FDR.png", dpi=300)
    fig.savefig(final / "random_subset_FDR.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(len(settings), 2, figsize=(9, 2.8 * len(settings)), squeeze=False, constrained_layout=True)
    for i, (ref, target) in enumerate(settings):
        for j, (outcome, effect) in enumerate(CASES[1:]):
            ax = axes[i, j]
            g = paired[(paired.reference_source == ref) & (paired.target_source == target) &
                       (paired.outcome == outcome) & paired.subset_id.str.startswith("random_") &
                       (paired.comparison == "marker_matched_minus_tmm")]
            ax.scatter(g.fdr_mean, g.power_mean, color="#32668c", s=24, alpha=.7)
            ax.axvline(0, color="gray", lw=.8)
            ax.axhline(0, color="gray", lw=.8)
            ax.set(xlabel="FDR difference: MRN - TMM", ylabel="Power difference: MRN - TMM",
                   title=f"{ref} -> {target}; effect {effect:g}")
    if synthetic:
        fig.suptitle("Software fixture only")
    fig.savefig(final / "random_subset_tradeoff.png", dpi=300)
    fig.savefig(final / "random_subset_tradeoff.pdf")
    plt.close(fig)
