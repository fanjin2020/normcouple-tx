"""Frozen-engine adapters and paired count interventions; no downloads."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("paper3_frozen_stage2", HERE / "vendor" / "stage2_frozen.py")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
E, M = S.E, S.M
E.GATES = ["basic"]
SEEDS = dict(paired_seed=2026101103, external_seed=2026101119, nuclei=200)
MEAD_COMMIT = "d66fa10491fb86de720a173e64b5cb26eeeefcd7"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return E.file_sha(Path(path))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    E.write_json(path, value)


def write_csv(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    E.write_csv(path, value)


def write_bin(path, matrix):
    # R reads doubles in column-major order; keep the byte order explicit.
    np.asarray(matrix, dtype="<f8").ravel(order="F").tofile(path)


def full_inputs(folders, gene_map):
    """Use all unambiguous shared genes, preserving the frozen panel and roles."""
    config = S.E.BUNDLED_INPUTS["config"]
    genes, axes, aliases = M.common_axes(folders, gene_map)
    sources = {label: M.read_source(label, folder, axes[label], genes,
               config["sources"][label]["source_id"]) for label, folder in folders.items()}
    roles = pd.DataFrame(E.BUNDLED_INPUTS["donor_roles"])
    panels = pd.DataFrame(E.BUNDLED_INPUTS["marker_panels"])
    tests = pd.DataFrame(config["tested_genes_native"])
    for frame in (roles, panels, tests):
        frame["split"] = frame["split"].astype(int)
    panels["order"] = panels["order"].astype(int)
    index, jobs = pd.Index(genes), []
    for (label, split), rows in roles.groupby(["source", "split"], sort=True):
        if split != 0:
            continue
        selected = panels[panels.source.eq(label) & panels.split.eq(split) & panels.panel.eq("balanced")]
        markers = index.get_indexer(selected.sort_values("order").gene_id)
        native = tests[tests.source.eq(label) & tests.split.eq(split)].gene_id.tolist()
        test = index.get_indexer([aliases[label].get(g) for g in native])
        require(len(markers) == 300 and len(test) == 100 and min(np.r_[markers, test]) >= 0,
                "Frozen panel/test identity could not be recovered")
        require(not set(markers) & set(test), "Frozen marker/test overlap")
        jobs.append(dict(source=sources[label], split=int(split),
                         refs=rows.loc[rows.role.eq("reference"), "donor_id"].tolist(),
                         targets=rows.loc[rows.role.eq("target"), "donor_id"].tolist(),
                         test=test, panels={"balanced": markers}))
    return sources, jobs


def prepare_job(job, output, mead_max_genes=2000):
    # pandas group keys can be np.int64; the frozen RNG hashes JSON integers.
    require(isinstance(job["split"], (int, np.integer)), "Split must be an integer")
    job["split"] = int(job["split"])
    source = job["source"]
    base = M.reference(source, job["refs"])
    markers, test = job["panels"]["balanced"], job["test"]
    target = job.get("target_source", source)
    present = target.get("present", np.ones(len(source["genes"]), bool))
    norm_genes = np.flatnonzero(present & np.isfinite(base["mean_counts"]).all(1)
                               & (base["mean_counts"].mean(1) > 0))
    # A broad reference-defined panel, not pure class markers; test genes are
    # excluded before every MEAD reference normalization and fit.
    mean = base["mean_counts"].mean(1)
    variability = base["mean_counts"].var(1) / np.maximum(mean, 1e-12)
    eligible = np.setdiff1d(norm_genes, test)
    order = np.lexsort((source["genes"][eligible], -variability[eligible]))
    mead_genes = np.sort(eligible[order[:mead_max_genes]])
    require(len(mead_genes) >= 100, "Too few disjoint genes for MEAD")
    require(set(markers) <= set(norm_genes), "Normalization universe omits a marker")
    test_scale = mean[test]
    require((test_scale > 0).all(), "Zero frozen test-gene scale")
    frozen_alt = np.flatnonzero(job.get("alt_mask", np.arange(len(test)) < 20))
    ranked = np.lexsort((source["genes"][test], test_scale))
    n_alt = min(20, len(test) // 4)
    job.update(base=base, norm_genes=norm_genes, mead_genes=mead_genes,
               injections={"frozen": frozen_alt, "low": ranked[:n_alt], "high": ranked[-n_alt:]})
    fingerprint = M.digest(dict(reference=source["label"], refs=job["refs"],
        markers=source["genes"][markers].tolist(), tests=source["genes"][test].tolist(),
        normalizers=source["genes"][norm_genes].tolist(), mead=source["genes"][mead_genes].tolist()))[:20]
    folder = Path(output) / "reference_inputs" / fingerprint
    folder.mkdir(parents=True, exist_ok=True)
    cells = np.concatenate([np.concatenate(source["masks"][d]) for d in job["refs"]])
    ref_counts = source["x"][mead_genes][:, cells].tocsc()
    if not (folder / "reference.mtx").exists():
        mmwrite(folder / "reference.mtx", ref_counts, field="integer")
    meta = source["meta"].iloc[cells].copy().reset_index(drop=True)
    meta["ref_index"] = meta.donor_id.map({d: i for i, d in enumerate(job["refs"])})
    write_csv(folder / "reference_cells.csv", meta[["cell_id", "donor_id", "cell_type", "ref_index"]])
    write_csv(folder / "mead_genes.csv", {"gene_id": source["genes"][mead_genes]})
    write_csv(folder / "normalization_genes.csv", {"gene_id": source["genes"][norm_genes]})
    selections = []
    for panel, indices in job["injections"].items():
        selections.extend(dict(panel=panel, order=i, gene_id=source["genes"][test[j]],
                               reference_mean_count=float(test_scale[j])) for i, j in enumerate(indices))
    write_csv(folder / "injection_sets.csv", selections)
    marks = E.official_markers(job)
    owner = marks.cell_type.to_numpy(dtype=int)
    norm_index = pd.Index(norm_genes)
    write_json(folder / "axes.json", dict(normalization_genes=len(norm_genes), classes=M.CT,
        marker_rows=(norm_index.get_indexer(markers) + 1).tolist(), marker_owner=owner.tolist(),
        mead_rows=(norm_index.get_indexer(mead_genes) + 1).tolist(),
        test_rows=(norm_index.get_indexer(test) + 1).tolist(),
        reference_cells=len(meta), mead_genes=len(mead_genes), mead_commit=MEAD_COMMIT))
    require((norm_index.get_indexer(test) >= 0).all(), "Normalization universe omits a test gene")
    job["reference_folder"] = str(folder)
    return job


def conditions(external=False):
    result = [dict(name="null", effect=0., pattern="up", panel="frozen"),
              dict(name="up05", effect=.5, pattern="up", panel="frozen")]
    if not external:
        result += [dict(name="up01", effect=.1, pattern="up", panel="frozen"),
                   dict(name="up025", effect=.25, pattern="up", panel="frozen"),
                   dict(name="balanced05", effect=.5, pattern="balanced", panel="frozen"),
                   dict(name="low05", effect=.5, pattern="up", panel="low"),
                   dict(name="high05", effect=.5, pattern="up", panel="high")]
    return result


def draw_base(job, unit, rep, options):
    source, target = job["source"], job.get("target_source", job["source"])
    external = "target_source" in job
    seed = options["external_seed"] if external else options["paired_seed"]
    key = ([S.SOURCE, source["label"], job["split"], unit["target_donor"], unit["n"], rep]
           if external else ["counts", source["label"], job["split"], unit["target_donor"], unit["n"], rep])
    prefix = "locked_external" if external else "innovation_v1"
    rng = M.rng_for(seed, prefix + "_design", key)
    group = np.repeat([0., 1.], unit["n"] // 2)
    rng.shuffle(group)
    weights = rng.multinomial(4, np.full(4, .25))
    ref = M.reference(source, job["refs"], weights)
    p = rng.dirichlet(M.concentration(job["base"], "reference_anchored"), unit["n"])
    p[:, 0] *= np.exp(group * unit["composition_shift"])
    p = M.normalize(p)
    draw = M.empirical_draw(target, unit["target_donor"], p, options["nuclei"],
                           M.rng_for(seed, prefix + "_counts", key))
    donor = target["donors"].index(unit["target_donor"])
    test = job["test"]
    scale = job["base"]["mean_counts"].mean(1)[test]
    mu = target["nums"][donor] / target["ncells"][donor]
    nuisance = (draw["cell_truth"] @ mu[test].T) / scale
    return dict(draw=draw, group=group, p=draw["cell_truth"], ref=ref,
                weights=weights, scale=scale, nuisance=nuisance, key=key,
                seed=seed, prefix=prefix, allocation=";".join(map(str, weights)))


def intervene(job, base, condition, options):
    test = job["test"]
    truth = np.zeros(len(test))
    chosen = job["injections"][condition["panel"]]
    truth[chosen] = condition["effect"]
    intercept = np.zeros(len(test))
    if condition["pattern"] == "balanced" and condition["effect"]:
        cut = len(chosen) // 2
        up, down = chosen[:cut], chosen[cut:]
        # Both arms receive valid Poisson additions. Equal expected added mass
        # across groups is defined using reference scales, never target outcomes.
        magnitude = condition["effect"] * base["scale"][up].sum() / base["scale"][down].sum()
        truth[down] = -magnitude
        intercept[down] = magnitude
    expected = options["nuclei"] * base["scale"][None, :] * (
        intercept[None, :] + base["group"][:, None] * truth[None, :])
    require((expected >= -1e-10).all(), "Negative Poisson expectation")
    old_style = condition["panel"] == "frozen" and condition["pattern"] == "up"
    rng = M.rng_for(base["seed"], base["prefix"] + "_spike", base["key"]) if old_style else M.rng_for(
        2026101207, "revision_spike", base["key"], condition["name"])
    additions = rng.poisson(np.maximum(expected, 0))
    y = (base["draw"]["counts"][:, test] + additions) / (options["nuclei"] * base["scale"])
    nuisance = base["nuisance"] + intercept[None, :]
    noise = y - nuisance - base["group"][:, None] * truth
    markers = job["panels"]["balanced"]
    original_total = base["draw"]["mass"].sum(1)
    total = original_total + additions.sum(1)
    return dict(y=y, group=base["group"], truth=truth, nuisance=nuisance, noise=noise,
                p=base["p"], ref=base["ref"], A=base["ref"]["A"][markers],
                mixture=base["draw"]["counts"][:, markers] / total[:, None],
                original_total=original_total, total=total, additions=additions,
                allocation=base["allocation"], intercept=intercept,
                expected_added_mass=expected.sum(1))


def verify_baseline(job, unit, rep, options, data, effect):
    case = dict(module="external" if "target_source" in job else "compare", n=unit["n"],
                target_donor=unit["target_donor"], composition_shift=unit["composition_shift"], effect=effect)
    old = S.scenario_data(job, case, rep, options)
    errors = {}
    for key in ("group", "p", "truth", "A", "mixture", "y", "nuisance", "noise"):
        difference = float(np.max(np.abs(old[key] - data[key])))
        require(np.allclose(old[key], data[key], rtol=1e-10, atol=1e-12), "Frozen generator mismatch: " + key)
        errors[key] = difference
    require(old["allocation"] == data["allocation"], "Reference draws changed")
    return errors


def plan_units(jobs, smoke=False):
    result = []
    for i, job in enumerate(jobs):
        for donor in (job["targets"][:1] if smoke else job["targets"]):
            unit = dict(job_index=i, reference_source=job["source"]["label"],
                        target_source=job.get("target_source", job["source"])["label"],
                        target_donor=donor, split=int(job["split"]), n=80, composition_shift=1.)
            unit["unit_id"] = "revision_" + M.digest(unit)[:16]
            unit["conditions"] = conditions("target_source" in job)
            result.append(unit)
    return result
