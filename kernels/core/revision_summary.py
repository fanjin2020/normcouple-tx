"""Shared, simulation-free result validation and donor-equal aggregation."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

SUMMARY_VERSION = "1.0.2"
BASELINE = "dt_observed_a1_c1000000"
ALL_PLANNED = ["availability", "fdp", "any_false", "power", "delivered_coverage",
               "null_delivered_coverage", "null_availability"]
CONDITIONAL = ["coverage", "null_coverage", "gene_error_mean", "abs_leakage", "proportion_rmse", "mean_se"]
METRICS = ALL_PLANNED + CONDITIONAL + ["group_error_contrast", "identity_error"] + [
    "error_contrast_" + ct for ct in ("Astro", "EndoMural", "Excit", "Inhib", "Micro", "OPC", "Oligo")]
TABLE_NAMES = ("replicates", "paired_replicates", "per_gene_summary", "mechanism_replicates")
WEIGHTING = "equal_target_donor_after_within_donor_available_repetitions"
MC_METHOD = "paired_repetition_ratio_influence;fixed_donor_banks"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_result_csv(path):
    # Empty numeric fields are missing; literal null/NA donor and condition labels are not.
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def write_csv(path, frame):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def methods_for(condition):
    common = {"oracle", "unadjusted", "nnls_observed", "nnls_fixed", "mead_native",
              "mead_reference_converted", "mead_lambda0_native", "mead_lambda0_reference_converted",
              "dt_tmm_a1_c1000000", "dt_deseq2_a1_c1000000",
              BASELINE + "_native", "dt_fixed_a1_c1000000_native"}
    settings = [("1", "1000000")]
    if condition in ("null", "up05"):
        settings += [("0.1", "1000000"), ("10", "1000000"), ("1", "10000"),
                     ("1", "100000000"), ("0.01", "10000")]
    return common | {f"dt_{kind}_a{a}_c{c}" for kind in ("observed", "fixed") for a, c in settings}


def pair_specs(condition):
    methods = methods_for(condition)
    pairs = {(method, BASELINE) for method in methods if method != BASELINE}
    pairs |= {(method, method.replace("dt_observed_", "dt_fixed_")) for method in methods
              if method.startswith("dt_observed_")}
    pairs |= {(method, method.replace("mead_", "mead_lambda0_", 1))
              for method in ("mead_native", "mead_reference_converted")}
    return pairs


def expected_index(unit, reps, kind):
    conditions = [item["name"] for item in unit["conditions"]]
    if kind == "paired_replicates":
        metrics = ALL_PLANNED + CONDITIONAL + ["beta_common_gene_mean", "leakage_common_gene_mean", "se_common_gene_mean"]
        keys = [(condition, a, b, metric) for condition in conditions
                for a, b in sorted(pair_specs(condition)) for metric in metrics]
        return keys, ["condition", "method_a", "method_b", "metric"]
    if kind == "mechanism_replicates":
        conditions = [condition for condition in conditions if condition != "null"]
    return [(condition, method) for condition in conditions for method in sorted(methods_for(condition))], ["condition", "method"]


def validate_unit_table(frame, unit, reps, kind, n_genes=100):
    require(not frame.condition.isna().any(), f"Missing condition labels in {unit['unit_id']}/{kind}")
    keys, columns = expected_index(unit, reps, kind)
    found = set(frame[columns].itertuples(index=False, name=None))
    require(found == set(keys), f"Missing/unexpected planned method-condition keys: {unit['unit_id']}/{kind}")
    if kind == "per_gene_summary":
        require(not frame.duplicated(columns + ["gene_id"]).any(), "Duplicate gene summary keys")
        require(frame.groupby(columns, dropna=False).size().eq(n_genes).all(), "Incomplete tested-gene family")
        require(frame.gene_id.nunique() == n_genes, "Gene sets differ across method/condition summaries")
        require(frame.planned.eq(len(reps)).all(), "Wrong planned gene repetitions")
        require(frame.loc[frame.condition.eq("null"), "truth"].eq(0).all(), "Nonzero truth in null condition")
        require(frame.available.between(0, len(reps)).all(), "Invalid available gene count")
    else:
        require(not frame.duplicated(columns + ["rep"]).any(), f"Duplicate replicate keys: {kind}")
        require(frame.rep.isin(reps).all(), "Unexpected replicate index")
        require(frame.groupby(columns, dropna=False).size().eq(len(reps)).all(), f"Incomplete repetitions: {kind}")
        if kind == "replicates":
            require(frame.n_tests.eq(n_genes).all(), "Unexpected test family size")
            require(frame.loc[frame.condition.eq("null"), "n_alternatives"].eq(0).all(), "Null condition has alternatives")
            require(frame.availability.between(0, 1).all(), "Invalid availability")
        elif kind == "paired_replicates":
            expected = np.where(frame.metric.isin(ALL_PLANNED), "all_planned",
                        np.where(frame.metric.isin(CONDITIONAL), "common_complete_fits", "common_available_genes"))
            require((frame.denominator == expected).all(), "Paired denominator labels differ from metric definition")


def recover_gene_null_label(frame, unit, repetitions):
    frame = frame.copy()
    blank = frame.condition.isna()
    if blank.any():
        require("null" in {item["name"] for item in unit["conditions"]}, "No null condition planned")
        require(not frame.condition.eq("null").any(), "Mixed blank and explicit null gene groups; inspect manually")
        missing = frame.loc[blank]
        require(missing.truth.eq(0).all(), "Cannot relabel nonzero-truth blank conditions")
        require(set(missing.method) == methods_for("null"), "Blank labels do not match the planned null method set")
        require(missing.groupby("method").size().eq(100).all(), "Blank labels do not contain full null gene families")
        frame.loc[blank, "condition"] = "null"
    validate_unit_table(frame, unit, list(range(repetitions)), "per_gene_summary")
    return frame, int(blank.sum())


def donor_equal_mc(frame, column, donors, repetitions):
    """Mean of within-donor finite-repetition means, with paired MC delta SE.

    Missing conditional fits are NOT zeros. For donor j, the influence for
    repetition r is I_rj*(x_rj-theta_j)/p_j. Averaging these influences across
    donors retains paired-repetition covariance and unequal success rates.
    A donor with no finite observations makes the full-bank estimand undefined.
    """
    require(not frame.duplicated(["rep", "target_donor"]).any(), "Duplicate donor/replicate cell")
    matrix = frame.pivot(index="rep", columns="target_donor", values=column).reindex(
        index=range(repetitions), columns=donors).to_numpy(dtype=float)
    require(not np.isinf(matrix).any(), "Infinite metric in summary: " + column)
    valid = np.isfinite(matrix)
    counts = valid.sum(axis=0)
    n_donors = len(donors)
    result = dict(mean=np.nan, mcse=np.nan, mc_low=np.nan, mc_high=np.nan,
        contributing_replicates=int(valid.any(axis=1).sum()), planned_repetitions=repetitions,
        target_donors=n_donors, contributing_donors=int((counts > 0).sum()),
        min_available_repetitions_per_donor=int(counts.min()), weighting=WEIGHTING,
        mcse_method=MC_METHOD, interval_scope="MC_only_conditional_fixed_banks",
        summary_status="no_finite_values" if not valid.any() else "donor_without_finite_values")
    if not (counts > 0).all():
        return result
    means = np.nansum(matrix, axis=0) / counts
    result.update(mean=float(means.mean()), summary_status="estimable")
    if repetitions > 1 and (counts >= 2).all():
        probability = counts / repetitions
        influence = (np.where(valid, matrix - means, 0.) / probability).mean(axis=1)
        se = float(influence.std(ddof=1) / np.sqrt(repetitions))
        result.update(mcse=se, mc_low=result["mean"] - 1.96 * se, mc_high=result["mean"] + 1.96 * se)
    else:
        result["summary_status"] = "mean_only_insufficient_repetitions_for_mcse"
    return result


def aggregate_frames(frames, units, repetitions, tables, include_raw=True):
    tables = Path(tables)
    tables.mkdir(parents=True, exist_ok=True)
    for name, frame in frames.items():
        if include_raw or name not in ("replicates", "paired_replicates"):
            write_csv(tables / (name + ".csv"), frame)
    groups = ["reference_source", "target_source", "condition", "method"]
    reps = frames["replicates"]
    plan = pd.DataFrame(units)
    donor_banks = {}
    for key, block in plan.groupby(["reference_source", "target_source"]):
        donors = block.target_donor.tolist()
        require(len(donors) == len(set(donors)), "This stage requires one unit per reference/target donor")
        donor_banks[key] = sorted(donors)
    donor = reps.groupby(groups + ["target_donor"], dropna=False)[METRICS].mean().reset_index()
    write_csv(tables / "donor_summary.csv", donor)
    summary = []
    for key, frame in reps.groupby(groups, dropna=False):
        for metric in METRICS:
            summary.append(dict(zip(groups, key), metric=metric,
                **donor_equal_mc(frame, metric, donor_banks[key[:2]], repetitions)))
    summary = pd.DataFrame(summary)
    write_csv(tables / "source_summary.csv", summary)
    pairgroups = groups[:-1] + ["method_a", "method_b", "metric", "denominator"]
    paired = []
    for key, frame in frames["paired_replicates"].groupby(pairgroups, dropna=False):
        paired.append(dict(zip(pairgroups, key), **donor_equal_mc(frame, "delta", donor_banks[key[:2]], repetitions),
            planned_pairs=len(frame), available_pairs=int(frame.delta.notna().sum())))
    write_csv(tables / "paired_summary.csv", pd.DataFrame(paired))
    failed = reps.loc[reps.availability < 1].copy()
    failed["reason"] = failed.reason.fillna("")
    fails = failed.groupby(groups + ["reason"], dropna=False).size().rename("replicates").reset_index()
    write_csv(tables / "failure_reasons.csv", fails)
    status = reps.groupby(groups, dropna=False).agg(planned_replicates=("rep", "size"),
        complete_fits=("complete_fit", "sum"), availability=("availability", "mean")).reset_index()
    status["status"] = np.where(status.availability == 0, "ALL_UNAVAILABLE", "evaluated")
    write_csv(tables / "method_status.csv", status)
    expected = {(u["reference_source"], u["target_source"], c["name"], method)
                for u in units for c in u["conditions"] for method in methods_for(c["name"])}
    for frame in (summary, status, donor):
        require(set(frame[groups].itertuples(index=False, name=None)) == expected,
                "Final report is missing planned method-condition groups")
    return summary, donor, status


def make_plots(tables, summary, donor):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keep = [BASELINE, "dt_fixed_a1_c1000000", "dt_tmm_a1_c1000000", "dt_deseq2_a1_c1000000",
            "nnls_observed", "mead_native", "mead_reference_converted", "oracle"]
    tables = Path(tables)
    selected = summary[summary.condition.eq("up05") & summary.method.isin(keep)]
    for (ref, target), block in selected.groupby(["reference_source", "target_source"]):
        fig, axes = plt.subplots(1, 3, figsize=(14, 5), layout="constrained")
        for ax, metric in zip(axes, ["fdp", "null_coverage", "availability"]):
            sub = block[block.metric.eq(metric)].set_index("method").reindex(keep)
            ax.barh(np.arange(len(keep)), sub["mean"], color="#187d87")
            ax.set_yticks(np.arange(len(keep)), keep, fontsize=8)
            ax.set_xlim(0, 1)
            ax.set_title("Mean FDP (FDR estimate)" if metric == "fdp" else metric)
        fig.suptitle(ref + " -> " + target + " (equal target-donor means)")
        fig.savefig(tables / f"summary_{ref}_{target}.pdf")
        fig.savefig(tables / f"summary_{ref}_{target}.png", dpi=160)
        plt.close(fig)
    selected = donor[donor.condition.eq("up05") & donor.method.isin(keep)]
    for (ref, target), block in selected.groupby(["reference_source", "target_source"]):
        fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
        for i, method in enumerate(keep):
            values = block[block.method.eq(method)].fdp.to_numpy()
            ax.scatter(np.full(len(values), i), values, s=20, alpha=.7)
        ax.set_xticks(range(len(keep)), keep, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Mean FDP by target donor")
        ax.set_ylim(0, 1)
        ax.set_title(ref + " -> " + target)
        fig.savefig(tables / f"donors_{ref}_{target}.pdf")
        plt.close(fig)


def write_summary_notes(path):
    Path(path).write_text(
        "# Summary version " + SUMMARY_VERSION + "\n\n"
        "Literal null condition labels are preserved. All planned method/condition/repetition keys are checked.\n"
        "Source means average within each target donor first, then give every planned donor equal weight.\n"
        "Conditional means use finite repeats within each donor; unavailable fits are not imputed as zero.\n"
        "If any donor has no finite value, the full-bank mean is undefined and is labeled accordingly.\n"
        "MC standard errors use paired-repetition ratio influences I*(x-donor_mean)/p_available,\n"
        "averaged over donors, then sample SD/sqrt(planned repetitions). They are delta-method MC\n"
        "errors conditional on the fixed donor banks, not donor-population confidence intervals.\n"
        "With complete data these reduce to the original paired-repetition mean and MC SE.\n"
        "At least two finite repetitions per donor are required for a reported MC SE.\n"
        "Paired method contrasts use the same weighting on within-repetition paired differences.\n"
        "All-planned availability/FDP/power/delivered coverage retain failure denominators.\n"
        "Mean FDP is an FDR estimate. any_false is a different probability under alternatives.\n"
        "Plots label mean FDP explicitly; older plots used any_false in that panel.\n"
        "MEAD is a point-estimator comparison, without propagation of proportion covariance to HC3.\n",
        encoding="utf-8")
