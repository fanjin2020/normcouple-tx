"""Bounded Paper 3 extension. All output and protocol identities are separate."""
import argparse
import copy
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
import platform
import shutil
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from revision_core import (S, E, M, SEEDS, require, sha, write_json, write_csv,
    write_bin, full_inputs, prepare_job, plan_units, draw_base, intervene, verify_baseline)
import revision_driver as D
from extension_stats import evaluate, mechanism_checks
from real_bulk import run_bulk, acquire, find_named

HERE = Path(__file__).resolve().parent
VERSION = "1.0.0"
JOBS, OPTIONS = [], {}
ARRAY_KEYS = ["truth", "beta", "se", "leakage", "noise_projection", "available", "covered", "p", "q"]
METRICS = ["fdp", "any_false", "power", "coverage", "delivered_coverage", "null_coverage",
           "null_delivered_coverage", "availability", "proportion_rmse", "abs_leakage", "mean_se"]
CONDITIONAL = {"coverage", "null_coverage", "abs_leakage", "mean_se"}


def source_identity():
    manifest = json.loads((HERE / "extension_manifest.json").read_text())
    for name, expected in manifest.items():
        require(sha(HERE / name) == expected, "Source checksum mismatch: " + name)
    return manifest


def choose_r(args):
    if args.r_library:
        return
    candidates = [args.project_root.parent / "paper3_revision_results_v1" / "R_library"]
    if args.old_results:
        candidates.insert(0, args.old_results / "R_library")
    for c in candidates:
        if (c / "MEAD").is_dir() and (c / "dtangle").is_dir():
            args.r_library = c
            print("REUSE R library", c, flush=True)
            break


def inputs(args, output):
    if args.input_bundle:
        return D.resolve_all_inputs(args)
    folders, gene_map, identities = E.resolve_inputs(args)
    paths = []
    roots = [args.external_cache, args.project_root, args.project_root.parent / "paper3_revision_results_v1"]
    for name, expected in ((S.META_NAME, S.PUBLIC_SHA256[S.META_NAME]), (S.COUNT_NAME, D.H5_SHA)):
        try:
            p = D.cached_file(name, roots, expected)
        except FileNotFoundError:
            p = acquire(name, args, output, expected)
        paths.append(p)
        identities.append(dict(name="external/" + name, path=str(p), sha256=sha(p), bytes=p.stat().st_size))
    return folders, gene_map, *paths, identities


def export_portable(identities, dest):
    dest.mkdir(parents=True, exist_ok=True)
    if any(dest.iterdir()) and not (dest / "INPUT_BUNDLE_OWNER.json").exists():
        raise ValueError("Refusing non-owned input export directory: " + str(dest))
    write_json(dest / "INPUT_BUNDLE_OWNER.json", {"package": "paper3_mechanism_inputs", "version": VERSION})
    entries = []
    for item in identities:
        name = "gene_identity_map.csv" if item["name"] == "gene_map" else item["name"]
        target = dest / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            require(sha(target) == item["sha256"], "Portable input conflict: " + str(target))
        else:
            temp = target.with_suffix(target.suffix + ".partial")
            shutil.copy2(item["path"], temp)
            require(sha(temp) == item["sha256"], "Portable copy checksum failed")
            os.replace(temp, target)
        entries.append(dict(relative_path=name, sha256=item["sha256"], bytes=target.stat().st_size))
    write_json(dest / "inputs_manifest.json", entries)
    (dest / "README.txt").write_text("Verified prepared inputs, separate from compact results.\n"
        "Use --input-bundle PATH on another machine. Check source licenses/consent before public deposition.\n"
        "This local export does not create a public DOI or permission to redistribute.\n", encoding="utf-8")


def conditions_for(job):
    units = plan_units([job])
    conditions = copy.deepcopy(units[0]["conditions"])
    for c in conditions:
        c["marker_stress"] = False
    for base in conditions[:2]:
        c = dict(base, name="marker_shift50_" + base["name"], marker_stress=True)
        conditions.append(c)
    return conditions


def marker_properties(job, base):
    marker = job["panels"]["balanced"]
    ref = base["ref"]
    A = ref["A"][marker]
    marker_yield = ref["pooled_yield"] * A.sum(0)
    sums = ref["counts"][marker].sum(0)
    anchor = float(sums.sum() / (sums / A.sum(0)).sum())
    require(np.isfinite(anchor) and 0 < anchor <= 1 and (marker_yield > 0).all(), "Invalid reference marker yield")
    return A, marker_yield, anchor


def stress_additions(job, base):
    marker = job["panels"]["balanced"]
    owner = E.official_markers(job).cell_type.to_numpy(int)
    selected = []
    for k in range(1, 8):
        ix = np.flatnonzero(owner == k)
        order = np.argsort(job["source"]["genes"][marker[ix]])
        selected.extend(ix[order[:max(1, int(np.ceil(len(ix) * .1)))]] )
    expectation = np.zeros((len(base["group"]), len(marker)))
    expectation[:, selected] = .5 * base["draw"]["counts"][:, marker[selected]] * base["group"][:, None]
    return M.rng_for(2026101303, "marker_stress", base["key"]).poisson(expectation), [int(i) for i in selected]


def run_batch(job, unit, reps, folder, options):
    folder.mkdir(parents=True, exist_ok=True)
    print("BATCH START", unit["unit_id"], "repetitions", reps[0], "to", reps[-1],
          "of", options["repetitions"], flush=True)
    records, arrays, props_all, z_all, index, pairs, invariant_rows = [], [], [], [], [], [], []
    start_clock = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="bridge_", dir=folder) as temporary:
        work = Path(temporary)
        entries, payload = [], {}
        for rep in reps:
            base = draw_base(job, unit, rep, options)
            A, marker_yield, anchor = marker_properties(job, base)
            stress, selected = stress_additions(job, base)
            key = f"r{rep}"
            write_bin(work / (key + "_raw.bin"), base["draw"]["counts"][:, job["norm_genes"]].T)
            write_bin(work / (key + "_A.bin"), A)
            write_bin(work / (key + "_total.bin"), base["draw"]["mass"].sum(1))
            write_bin(work / (key + "_stress.bin"), stress.T)
            conds = []
            for condition in unit["conditions"]:
                # Preserve the baseline test-gene draws even in marker stress arms.
                c = dict(condition)
                if c["marker_stress"]: c["name"] = c["name"].removeprefix("marker_shift50_")
                data = intervene(job, base, c, options)
                payload[(rep, condition["name"])] = data, marker_yield, base["ref"]["pooled_yield"]
                if rep == 0 and not condition["marker_stress"] and condition["name"] in ("null", "up05"):
                    verify_baseline(job, unit, rep, options, data, c["effect"])
                name = key + "_" + condition["name"] + ".bin"
                write_bin(work / name, data["additions"].T)
                conds.append(dict(name=condition["name"], additions=name, marker_stress=condition["marker_stress"]))
            entries.append(dict(key=key, n=unit["n"], weights=base["weights"].tolist(),
                raw=key+"_raw.bin", A=key+"_A.bin", total=key+"_total.bin", stress=key+"_stress.bin",
                marker_anchor=anchor, conditions=conds))
        write_json(work / "jobs.json", dict(reference_folder=job["reference_folder"], entries=entries,
            mead_cache=str(Path(options["output"]) / "units" / unit["unit_id"] / "mead_reference_cache")))
        info = options["r_info"]
        D.run_r([info["executable"], "--vanilla", str(HERE / "bridge_extension.R"), str(work), info["library"]],
                info["library"], folder / "bridge.log", options["r_timeout"])
        result = {}
        for row in pd.read_csv(work / "outputs.csv", keep_default_na=False).itertuples():
            p = np.loadtxt(work / row.file, delimiter=",", ndmin=2)
            require(p.shape == (unit["n"], 7), "Prediction shape mismatch")
            result[(row.key, row.arm)] = p, row.reason
        for rep in reps:
            null_preds = {}
            for condition in unit["conditions"]:
                name = condition["name"]
                data, marker_yield, full_yield = payload[(rep, name)]
                predictions = {("oracle", "full_six"): (data["p"], "")}
                for method in ("observed", "fixed", "tmm", "marker_matched", "marker_unmatched_scale", "marker_denominator_only"):
                    q, reason = result[(f"r{rep}_{name}", method)]
                    yield_ = marker_yield if method in ("marker_matched", "marker_unmatched_scale") else full_yield
                    predictions[(method, "full_six")] = E.convert_available_rows(q, yield_), reason
                if not condition["marker_stress"]:
                    mead, reason = result[(f"r{rep}", "mead_native")]
                    for method, p in (("mead_native", mead), ("mead_reference_converted", E.convert_available_rows(mead, full_yield))):
                        for rule in ("full_six", "adaptive"):
                            predictions[(method, rule)] = p, reason
                current = {}
                for (method, rule), (p, reason) in predictions.items():
                    row, arr, fit = evaluate(data, p, rule, reason)
                    attrs = dict(rep=rep, condition=name, method=method, rule=rule)
                    records.append(dict(attrs, **row)); index.append(attrs)
                    arrays.append(np.column_stack([arr[k] for k in ARRAY_KEYS]))
                    props_all.append(p); z_all.append(fit["residualized_exposure"])
                    current[(method, rule)] = row
                    if name == "null": null_preds[method, rule] = p
                    if name == "up05" and method.startswith("marker_"):
                        old = null_preds[method, rule]
                        finite = np.isfinite(p).all() and np.isfinite(old).all()
                        error = float(abs(p-old).max()) if finite else np.nan
                        invariant_rows.append(dict(rep=rep, method=method, finite=finite, max_change=error))
                        require(not finite or error < 1e-8, "Marker-only invariance failed: " + method)
                reference = current[("observed", "full_six")]
                comparisons = [((method, rule), ("observed", "full_six")) for method, rule in current if method != "oracle"]
                comparisons += [((method, "adaptive"), (method, "full_six")) for method in ("mead_native", "mead_reference_converted")
                                if (method, "adaptive") in current]
                for key_a, key_b in comparisons:
                    a, b = current[key_a], current[key_b]
                    common = a["availability"] == 1 and b["availability"] == 1
                    for metric in METRICS:
                        val = a[metric] - b[metric] if metric not in CONDITIONAL or common else np.nan
                        pairs.append(dict(rep=rep, condition=name, method=key_a[0], rule=key_a[1],
                            reference=key_b[0], reference_rule=key_b[1], metric=metric, value=val,
                            denominator="common_complete_fits" if metric in CONDITIONAL else "all_planned"))
        write_csv(folder / "replicates.csv", records)
        write_csv(folder / "pairs.csv", pairs)
        write_csv(folder / "fit_index.csv", index)
        write_csv(folder / "marker_invariance.csv", invariant_rows)
        np.save(folder / "gene_fits.npy", np.stack(arrays), allow_pickle=False)
        np.save(folder / "proportions.npy", np.stack(props_all), allow_pickle=False)
        np.save(folder / "residualized_exposure.npy", np.stack(z_all), allow_pickle=False)
        write_json(folder / "axes.json", dict(array_keys=ARRAY_KEYS, cell_types=M.CT,
            genes=job["source"]["genes"][job["test"]].tolist(), array_order="fit_index row, gene or mixture, field",
            stress_marker_rows=selected, elapsed_seconds=time.monotonic()-start_clock))
        for file in ("normalization_factors.csv", "outputs.csv", "r_controls.json"):
            shutil.copy2(work / file, folder / file)
    files = [p.name for p in folder.iterdir() if p.is_file() and p.name != "complete.json"]
    D.seal(folder, files, options["signature"])
    print("BATCH DONE", unit["unit_id"], reps[0], "seconds", round(time.monotonic()-start_clock, 1), flush=True)


def run_unit(unit):
    job, options = JOBS[unit["job_index"]], OPTIONS
    folder = Path(options["output"]) / "units" / unit["unit_id"]
    folder.mkdir(parents=True, exist_ok=True)
    chunks = []
    for start in range(0, options["repetitions"], options["batch_size"]):
        c = folder / f"chunk_{start:05d}"
        if not D.valid_checkpoint(c, options["signature"]):
            run_batch(job, unit, list(range(start, min(start+options["batch_size"], options["repetitions"]))), c, options)
        else:
            print("RESUME", unit["unit_id"], start, flush=True)
        chunks.append(c)
    write_json(folder / "unit.json", unit)
    return chunks


def summarize(output, units, options):
    records, pairs, invariants, genes = [], [], [], []
    for unit in units:
        folder = output / "units" / unit["unit_id"]
        chunks = [folder / f"chunk_{i:05d}" for i in range(0, options["repetitions"], options["batch_size"])]
        attrs = {k:unit[k] for k in ("unit_id", "reference_source", "target_source", "target_donor")}
        ix, aa = [], []
        for chunk in chunks:
            require(D.valid_checkpoint(chunk, options["signature"]), "Missing checkpoint")
            for name, target in (("replicates", records), ("pairs", pairs), ("marker_invariance", invariants)):
                f = pd.read_csv(chunk / (name+".csv"), keep_default_na=False, na_values=[""])
                require(not f.condition.isna().any() if "condition" in f else True, "Null condition parsed as missing")
                for k,v in attrs.items(): f[k]=v
                target.append(f)
            ix.append(pd.read_csv(chunk / "fit_index.csv", keep_default_na=False))
            aa.append(np.load(chunk / "gene_fits.npy", allow_pickle=False))
        ix = pd.concat(ix, ignore_index=True); a = np.concatenate(aa)
        names = json.loads((chunks[0] / "axes.json").read_text())["genes"]
        for (condition, method, rule), loc in ix.groupby(["condition", "method", "rule"]).groups.items():
            indices = np.asarray(loc)
            require(len(indices) == options["repetitions"] and ix.loc[indices].rep.nunique() == options["repetitions"], "Incomplete gene repetitions")
            block = a[indices]; truth = block[0,:,0]; avail = block[:,:,5].sum(0)
            mean = lambda v: np.divide(np.nansum(v,0),np.isfinite(v).sum(0),out=np.full(v.shape[1],np.nan),where=np.isfinite(v).sum(0)>0)
            beta = block[:,:,1]; mu = mean(beta)
            sd = np.sqrt(np.divide(np.nansum((beta-mu)**2,0),avail-1,out=np.full(len(names),np.nan),where=avail>1))
            f = pd.DataFrame(dict(gene_id=names, truth=truth, bias=mu-truth, empirical_sd=sd,
                mean_se=mean(block[:,:,2]), rms_se=np.sqrt(mean(block[:,:,2]**2)),
                availability=avail/options["repetitions"],
                coverage=np.divide(block[:,:,6].sum(0),avail,out=np.full(len(names),np.nan),where=avail>0),
                delivered_coverage=block[:,:,6].mean(0)))
            for k,v in dict(attrs, condition=condition, method=method, rule=rule).items(): f[k]=v
            genes.append(f)
    records = pd.concat(records, ignore_index=True); pairs = pd.concat(pairs, ignore_index=True)
    group = ["reference_source", "target_source", "condition", "method", "rule"]
    require(not records.duplicated(["unit_id","rep","condition","method","rule"]).any(), "Duplicate repetition records")
    actual = records.groupby(["unit_id","condition","method","rule"]).rep.nunique()
    require(actual.eq(options["repetitions"]).all(), "Unequal repetition denominators")
    final = output / "final_results"; final.mkdir(exist_ok=True)
    # Per-repetition tables can be large; the collection folder holds summaries.
    write_csv(output / "audit_tables" / "replicate_diagnostics.csv", records)
    write_csv(output / "audit_tables" / "paired_replicates.csv", pairs)
    write_csv(final / "gene_summary.csv", pd.concat(genes, ignore_index=True))
    write_csv(final / "marker_invariance.csv", pd.concat(invariants, ignore_index=True))
    donor_rows, summary = [], []
    for keys, f in records.groupby(group, dropna=False):
        for metric in METRICS:
            vals, variances = [], []
            for donor, d in f.groupby("target_donor"):
                x = pd.to_numeric(d[metric], errors="raise").dropna().to_numpy()
                value = float(x.mean()) if len(x) else np.nan
                variance = float(x.var(ddof=1)/len(x)) if len(x)>1 else np.nan
                vals.append(value); variances.append(variance)
                donor_rows.append(dict(zip(group,keys), target_donor=donor, metric=metric, mean=value, contributing_repetitions=len(x)))
            valid = np.isfinite(vals); n=int(valid.sum())
            value=float(np.nanmean(vals)) if n else np.nan
            se=float(np.sqrt(np.sum(np.asarray(variances)[valid]))/n) if n and np.isfinite(np.asarray(variances)[valid]).all() else np.nan
            summary.append(dict(zip(group,keys), metric=metric, mean=value, mcse=se, low=value-1.96*se, high=value+1.96*se,
                planned_donors=len(vals), contributing_donors=n, uncertainty="MC conditional on fixed donor banks"))
    paired_summary = []
    pg = group + ["reference", "reference_rule", "metric", "denominator"]
    for keys, f in pairs.groupby(pg, dropna=False):
        ds = f.groupby("target_donor").value.agg(["mean","var","count"])
        valid = ds["mean"].notna(); n=int(valid.sum())
        mean=float(ds.loc[valid,"mean"].mean()) if n else np.nan
        se=float(np.sqrt((ds.loc[valid,"var"]/ds.loc[valid,"count"]).sum(min_count=n))/n) if n else np.nan
        paired_summary.append(dict(zip(pg,keys), mean=mean, mcse=se, low=mean-1.96*se, high=mean+1.96*se,
            planned_donors=len(ds), contributing_donors=n, contributing_pairs=int(ds["count"].sum())))
    write_csv(final / "summary.csv", summary); write_csv(final / "donor_summary.csv", donor_rows)
    write_csv(final / "paired_summary.csv", paired_summary)
    failures=records.groupby(group+["reason"],dropna=False).size().reset_index(name="planned_fits")
    write_csv(final / "failure_reasons.csv",failures)
    from extension_figures import figures
    figures(final)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root",type=Path,default=Path("/home/research/git/repository/paper1_public_evidence"))
    p.add_argument("--output",type=Path,default=Path("/home/research/git/repository/paper3_mechanism_results_v1"))
    for name in ("libd-dir","seaad-dir","gene-map","external-cache","input-bundle","r-library","old-results","bulk-manifest"):
        p.add_argument("--"+name,type=Path)
    p.add_argument("--rscript",default="Rscript")
    p.add_argument("--workers",type=int,default=12 if os.name!="nt" else 1)
    p.add_argument("--repetitions",type=int,default=200)
    p.add_argument("--batch-size",type=int,default=10)
    p.add_argument("--r-timeout",type=int,default=14400)
    p.add_argument("--only",choices=["all","simulations","bulk","audit"],default="all")
    for flag in ("install-r-packages","offline","synthetic","pilot-only","no-input-export"):
        p.add_argument("--"+flag,action="store_true")
    return p


def main(argv=None):
    global JOBS, OPTIONS
    args=parser().parse_args(argv)
    require(args.workers>=1 and args.batch_size>=1,"Invalid worker/batch count")
    require(os.name!="nt" or args.workers==1,"Use --workers 1 on Windows")
    require(args.synthetic or args.repetitions>=200,"Research run requires >=200 repetitions; use --pilot-only for first batch")
    out=args.output.resolve()
    if out.exists() and any(out.iterdir()):
        require((out/"OWNER.json").exists(),"Output must be new or owned")
        require(json.loads((out/"OWNER.json").read_text()).get("package")=="paper3_mechanism_extension",
                "Output belongs to a different analysis; use a new directory")
    out.mkdir(parents=True,exist_ok=True); E.hold_output_lock(out)
    write_json(out/"OWNER.json",dict(package="paper3_mechanism_extension",version=VERSION))
    source=source_identity(); choose_r(args)
    print("PREFLIGHT: checking R estimators and numerical self-tests",flush=True)
    r_info=D.prepare_r(args,out)
    rda_status={"status":"NOT_REQUESTED"}
    identities=[]
    if args.synthetic:
        jobs=D.synthetic_jobs(out)
    else:
        print("INPUTS: verifying cached matrices and preparing fixed reference/target banks",flush=True)
        args.stage2_results=None
        folders,gene_map,meta,counts,identities=inputs(args,out)
        write_json(out/"input_paths.json",identities)
        if not args.no_input_export:
            export_portable(identities,out.parent/(out.name+"_prepared_inputs"))
        sources,jobs=full_inputs(folders,gene_map)
        external=S.read_external(counts,meta,next(iter(sources.values()))["genes"],out/"external_inputs")
        jobs += S.make_external_jobs(jobs,external,out)
        jobs=[prepare_job(j,out,2000) for j in jobs]
    units=plan_units(jobs,smoke=args.synthetic)
    for unit in units: unit["conditions"]=conditions_for(jobs[unit["job_index"]])
    require(args.synthetic or len(units)==42,"Expected exactly the frozen 42 units")
    settings=dict(SEEDS,repetitions=2 if args.synthetic else args.repetitions,batch_size=2 if args.synthetic else args.batch_size,
        marker_stress_seed=2026101303,stress_fraction=.1,stress_increment=.5,primary_marker="marker_matched",
        rank_relative_tolerance=1e-10,rank_absolute_tolerance=1e-12,software_fixture=args.synthetic)
    body=dict(version=VERSION,settings=settings,units=units,code_hashes=source,
        input_hashes={i["name"]:i["sha256"] for i in identities},R=r_info["versions"],
        python={p:importlib.metadata.version(p) for p in ("numpy","pandas","scipy","h5py","matplotlib")})
    signature=M.digest(body)
    if (out/"protocol.json").exists(): require(json.loads((out/"protocol.json").read_text())["signature"]==signature,"Protocol changed; use a new output directory")
    write_json(out/"protocol.json",dict(body,signature=signature))
    write_json(out/"mechanism_checks.json",mechanism_checks())
    shutil.copy2(HERE/"MECHANISM.md",out/"MECHANISM.md")
    JOBS=jobs; OPTIONS=dict(settings,signature=signature,output=str(out),r_info=r_info,r_timeout=args.r_timeout)
    write_json(out/"RUN_STATUS.json",dict(status="RUNNING",signature=signature))
    if args.only in ("all","bulk","audit") and not args.synthetic:
        print("REAL BULK: auditing public sample identities and raw-count availability",flush=True)
        try:
            rda_status=run_bulk(args,out/"real_bulk",jobs,r_info)
        except Exception as exc:
            rda_status=dict(status="BULK_FAILED",error=str(exc))
            write_json(out/"real_bulk"/"STATUS.json",rda_status)
            print("REAL BULK ATTENTION",str(exc),flush=True)
    sim_status="NOT_REQUESTED"
    if args.only in ("all","simulations"):
        print("PILOT: testing marker invariance before expanding to",len(units),"units",flush=True)
        # Numerical gate only. No outcome-based selection of a winning normalizer.
        pilot_units=[]
        for source in sorted({u["reference_source"] for u in units}):
            pilot_units.append(next(u for u in units if u["reference_source"]==source and u["reference_source"]==u["target_source"]))
        pilot_rows=[]
        try:
            for unit in pilot_units:
                folder=out/"units"/unit["unit_id"]/"chunk_00000"
                if not D.valid_checkpoint(folder,signature):
                    run_batch(jobs[unit["job_index"]],unit,list(range(min(OPTIONS["batch_size"],OPTIONS["repetitions"]))),folder,OPTIONS)
                f=pd.read_csv(folder/"marker_invariance.csv")
                require(f.finite.all() and f.max_change.lt(1e-8).all(),"Pilot marker predictions unavailable or invariant identity failed; inspect pilot, do not expand")
                pilot_rows.extend(f.assign(unit_id=unit["unit_id"]).to_dict("records"))
        except BaseException as exc:
            write_json(out/"RUN_STATUS.json",dict(status="PILOT_INCOMPLETE",error=str(exc),real_bulk=rda_status,signature=signature))
            raise
        write_json(out/"pilot_gate.json",dict(status="PASSED",criterion="finite marker estimates and paired invariance; not better FDR",checks=pilot_rows))
        if args.pilot_only:
            write_json(out/"RUN_STATUS.json",dict(status="PILOT_ONLY",real_bulk=rda_status,signature=signature))
            return 0
        try:
            if args.workers==1:
                for unit in units: run_unit(unit)
            else:
                with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context("fork")) as pool:
                    futures=[pool.submit(run_unit,u) for u in units]
                    for future in as_completed(futures): future.result()
            summarize(out,units,OPTIONS); sim_status="COMPLETE"
        except BaseException as exc:
            write_json(out/"RUN_STATUS.json",dict(status="SIMULATION_INCOMPLETE",error=str(exc),real_bulk=rda_status,signature=signature))
            raise
    if args.only in ("bulk","audit") and (out/"final_results"/"summary.csv").exists(): sim_status="PREVIOUS_SIMULATIONS_PRESENT"
    status="SOFTWARE_FIXTURE_ONLY" if args.synthetic else (
        "COMPLETE_REQUESTED_STAGES" if rda_status["status"] in ("COMPLETE_RAW_COUNT_APPLICATION","NOT_REQUESTED") else "PARTIAL_REAL_BULK_COUNTS_OR_ACCESS_NEEDED")
    receipt=dict(status=status,simulations=sim_status,real_bulk=rda_status,signature=signature,
        public_deposition="NOT_PERFORMED",rank_adaptive_analysis="Point-estimate covariates only; not joint MEAD covariance propagation")
    write_json(out/"RUN_STATUS.json",receipt)
    final=out/"final_results"; final.mkdir(exist_ok=True)
    for name in ("RUN_STATUS.json","protocol.json","mechanism_checks.json","MECHANISM.md"):
        shutil.copy2(out/name,final/name)
    if (out/"real_bulk").exists():
        dest=final/"real_bulk";dest.mkdir(exist_ok=True)
        for p in (out/"real_bulk").iterdir():
            if p.is_file() and p.suffix in (".csv",".json",".txt") and p.name not in ("bulk_processed_expression.csv","bulk_counts_for_R.csv"):
                shutil.copy2(p,dest/p.name)
        if (out/"real_bulk"/"figures").exists():
            shutil.copytree(out/"real_bulk"/"figures",dest/"figures",dirs_exist_ok=True)
    shutil.copy2(HERE/"README.md",final/"README.md")
    print(json.dumps(receipt,indent=2),flush=True)
    print("COLLECT THIS FOLDER:",final,flush=True)
    return 0 if status in ("SOFTWARE_FIXTURE_ONLY","COMPLETE_REQUESTED_STAGES") else 2


if __name__=="__main__":
    raise SystemExit(main())
