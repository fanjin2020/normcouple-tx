"""Freeze randomized STRESS subsets, never change the normalization panel."""
import numpy as np
import revision_core as C

SELECTION_SEED = 2026101407
STRESS_SEED = 2026101409
FRACTION = 0.10
AMPLITUDE = 0.5
METHODS = ["oracle", "tmm", "marker_matched", "marker_denominator_only"]
CASES = [("null", 0.0), ("up025", 0.25), ("up05", 0.5)]


def subsets(genes, owners, panel_key, number=50):
    genes, owners = np.asarray(genes), np.asarray(owners)
    C.require(len(genes) == len(set(genes.tolist())), "Duplicate marker identifiers")
    C.require(set(owners.tolist()) == set(range(1, 8)), "Need all seven marker classes")
    ordered = [np.array(sorted(np.flatnonzero(owners == k), key=lambda i: genes[i]), dtype=int)
               for k in range(1, 8)]
    sizes = [max(1, int(np.ceil(len(ix) * FRACTION))) for ix in ordered]
    out = {"zero": np.array([], dtype=int),
           "legacy": np.concatenate([ix[:n] for ix, n in zip(ordered, sizes)])}
    # Independent draws across subset IDs; no outcome-dependent replacement of a subset.
    for j in range(number):
        parts = [C.M.rng_for(SELECTION_SEED, "stress_subset", panel_key, j, k).choice(ix, n, replace=False)
                 for k, (ix, n) in enumerate(zip(ordered, sizes), 1)]
        out["random_%03d" % (j + 1)] = np.sort(np.concatenate(parts))
    return out


def freeze(jobs, number):
    manifest = []
    original = {j["source"]["label"]: j["panels"]["balanced"] for j in jobs if "target_source" not in j}
    for index, job in enumerate(jobs):
        marker = job["panels"]["balanced"]
        genes = job["source"]["genes"][marker]
        owner = C.E.official_markers(job).cell_type.to_numpy(int)
        panel_key = [job["source"]["label"], genes.tolist(), owner.tolist()]
        job["stress_subsets"] = subsets(genes, owner, panel_key, number)
        full = original[job["source"]["label"]]
        excluded = np.setdiff1d(full, marker)
        target = job.get("target_source", job["source"])["label"]
        for name, indices in job["stress_subsets"].items():
            manifest.append(dict(job_index=index, reference_source=job["source"]["label"],
                target_source=target, subset_id=name, selection_seed=SELECTION_SEED,
                requested_fraction=0.0 if name == "zero" else FRACTION,
                amplitude=0.0 if name == "zero" else AMPLITUDE,
                panel_size=len(marker), original_panel_size=len(full),
                annotation_coverage=len(marker) / len(full),
                annotation_excluded_gene_ids=job["source"]["genes"][excluded].tolist(),
                selected_count=len(indices), realized_fraction=len(indices) / len(marker),
                selected_gene_ids=genes[indices].tolist(),
                selected_class_ids=owner[indices].tolist(),
                per_class_counts=[int(np.sum(owner[indices] == k)) for k in range(1, 8)],
                selection_rule="uniform within retained annotated class; same panel for all subsets; no replacement within class"))
    return manifest


def stress(job, base, subset_id):
    marker = job["panels"]["balanced"]
    indices = job["stress_subsets"][subset_id]
    original = base["draw"]["counts"][:, marker]
    expected = np.zeros_like(original, dtype=float)
    expected[:, indices] = AMPLITUDE * original[:, indices] * base["group"][:, None]
    if subset_id == "legacy":
        rng = C.M.rng_for(2026101303, "marker_stress", base["key"])
    else:
        rng = C.M.rng_for(STRESS_SEED, "random_marker_stress", base["key"], subset_id)
    value = rng.poisson(expected)
    C.require(not np.any(value[base["group"] == 0]), "Controls received marker stress")
    C.require(not np.any(value[:, np.setdiff1d(np.arange(len(marker)), indices)]), "Unselected markers changed")
    case = base["group"] == 1
    selected_mass = original[case][:, indices].sum(1)
    marker_mass = original[case].sum(1)
    total = base["draw"]["mass"].sum(1)[case]
    details = dict(selected_markers=len(indices), panel_markers=len(marker),
        realized_fraction=len(indices) / len(marker),
        selected_mean_count_per_case_marker=float(selected_mass.mean() / max(1, len(indices))),
        selected_baseline_case_mass=float(selected_mass.mean()),
        selected_fraction_of_case_marker_mass=float(np.mean(selected_mass / marker_mass)),
        expected_added_case_mass=float(expected[case].sum(1).mean()),
        realized_added_case_mass=float(value[case].sum(1).mean()),
        realized_added_fraction_of_case_library=float(np.mean(value[case].sum(1) / total)))
    return value, details
