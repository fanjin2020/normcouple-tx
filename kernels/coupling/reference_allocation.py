"""A new additive intervention, not imputation of an undefined observed-weight rule."""
import hashlib
import numpy as np
from common import require, save_json

SEED = 2026102703


def allocate_fixed_weights(weights, selected, owners, budget, rng):
    weights = np.asarray(weights, float)
    owners, selected, budget = np.asarray(owners), np.asarray(selected), np.asarray(budget)
    require(weights.ndim == 1 and owners.shape == weights.shape, "Invalid allocation weights/classes")
    require(np.isfinite(weights).all() and (weights >= 0).all(), "Invalid fixed allocation weights")
    require(np.isin(owners, np.arange(1, 8)).all(), "Invalid marker classes")
    require(selected.ndim == 1 and np.issubdtype(selected.dtype, np.integer)
            and len(np.unique(selected)) == len(selected) and (selected >= 0).all()
            and (selected < len(weights)).all(), "Invalid selected indices")
    require(budget.ndim == 2 and budget.shape[1] == 7 and np.isfinite(budget).all()
            and (budget >= 0).all() and np.array_equal(budget, np.rint(budget))
            and (budget < np.iinfo(np.int64).max).all(), "Invalid integer budgets")
    value = np.zeros((len(budget), len(weights)), dtype=np.int64)
    for k in range(1, 8):
        ix = np.sort(selected[owners[selected] == k])
        require(len(ix) > 0, "Selected set lacks a marker class")
        w = weights[ix]
        require(w.sum() > 0, "ZERO_FROZEN_REFERENCE_SUPPORT: no allocation is defined; no uniform fallback")
        probabilities = w / w.sum()
        for sample in np.flatnonzero(budget[:, k-1] > 0):
            value[sample, ix] = rng.multinomial(int(budget[sample, k-1]), probabilities)
    actual = np.column_stack([value[:, owners == k].sum(1) for k in range(1, 8)])
    require(np.array_equal(actual, budget), "Fixed-weight class budget changed")
    return value


def frozen_weights(job, owners, policy):
    marker = np.asarray(job["panels"]["balanced"])
    require(policy in ("reference", "uniform"), "Unknown additive allocation policy")
    if policy == "uniform":
        return np.ones(len(marker))
    # The original reference bank is fixed before target draws and direct effects.
    # Do not use the target's observed counts, disease effects, or true proportions.
    A = np.asarray(job["base"]["A"], float)[marker]
    require(A.shape == (len(marker), 7) and np.isfinite(A).all() and (A >= 0).all(),
            "Invalid frozen reference profile")
    return A[np.arange(len(marker)), owners-1]


def reference_stress(job, base, subset_id):
    import selection as S
    import revision_core as C
    policy = job["budget_allocation"]
    marker = job["panels"]["balanced"]
    original = base["draw"]["counts"][:, marker]
    owners = C.E.official_markers(job).cell_type.to_numpy(int)
    legacy, _ = S.stress(job, base, "legacy")
    budget = np.column_stack([legacy[:, owners == k].sum(1) for k in range(1, 8)])
    selected = np.asarray(job["stress_subsets"][subset_id], int)
    weights = frozen_weights(job, owners, policy)
    _, detail = S.stress(job, base, subset_id)
    if subset_id == "zero":
        value = np.zeros_like(original, dtype=np.int64)
        budget = np.zeros_like(budget)
    else:
        # Legacy and random sets are BOTH reallocated by the same fixed-weight rule.
        value = allocate_fixed_weights(weights, selected, owners, budget,
            C.M.rng_for(SEED, "additive_"+policy+"_budget", base["key"], subset_id))
    case = base["group"] == 1
    require(not value[~case].any(), "Controls received marker perturbation")
    actual = np.column_stack([value[:, owners == k].sum(1) for k in range(1, 8)])
    require(np.array_equal(actual, budget), "Realized sample/class budget changed")
    observed_selected = float(original[case][:, selected].sum())
    added = int(value.sum())
    selected_legacy = job["stress_subsets"]["legacy"]
    expected = float((.5*original[case][:, selected_legacy].sum(1)).mean()) if subset_id != "zero" else 0.
    library = base["draw"]["mass"].sum(1)[case]
    wbytes = np.asarray(weights, dtype="<f8").tobytes()
    detail.update(
        expected_added_case_mass=expected,
        realized_added_case_mass=float(value[case].sum(1).mean()),
        realized_added_fraction_of_case_library=float(np.mean(value[case].sum(1)/library)),
        intervention="additive fixed-"+policy+" allocation under legacy sample-by-class integer budgets",
        allocation_policy=policy,
        allocation_weights_sha256=hashlib.sha256(wbytes).hexdigest(),
        budget_sha256=hashlib.sha256(budget.astype("<i8").tobytes()).hexdigest(),
        max_class_budget_error=int(np.max(abs(actual-budget))),
        effective_relative_amplitude=added/observed_selected if observed_selected else np.nan,
        effective_relative_amplitude_defined=bool(observed_selected),
        added_counts_on_observed_zero=int(value[original == 0].sum()),
        total_added_counts=added,
        intervention_feasible=True, infeasibility_reason="",
        unsupported_sample_class_count=0, unsupported_sample_count=0,
        unsupported_budget_total=0, requested_budget_total=int(budget.sum()),
        unsupported_sample_classes="[]")
    return value, detail


def preflight(jobs, cfg, destination, stress):
    """Test all planned marker sets over every unit without starting R fits."""
    import revision_core as C
    records = []
    repetitions = min(2, cfg["settings"]["repetitions"])
    for unit in cfg["units"]:
        job = jobs[unit["job_index"]]
        for rep in range(repetitions):
            base = C.draw_base(job, unit, rep, cfg["settings"])
            hashes = set()
            for subset in cfg["settings"]["subset_ids"]:
                value, detail = stress(job, base, subset)
                require(value is not None and detail["intervention_feasible"], "Additive allocation preflight failed")
                require(detail["max_class_budget_error"] == 0, "Preflight budget mismatch")
                if subset != "zero":
                    hashes.add(detail["budget_sha256"])
                records.append(dict(unit_id=unit["unit_id"], rep=rep, subset_id=subset,
                    feasible=True, max_class_budget_error=0,
                    added_counts_on_observed_zero=detail["added_counts_on_observed_zero"]))
            require(len(hashes) == 1, "Preflight budgets differ across marker identities")
    result = dict(status="PASS", scientific_validation=False,
                  scope="construction only; every unit and subset, first two repetitions; no R fits",
                  pairs=len(records), units=len(cfg["units"]), subsets=len(cfg["settings"]["subset_ids"]),
                  allocation_policy=cfg["settings"]["allocation_policy"], records=records)
    save_json(destination, result)
    return result
