"""Paired revision experiments using cached, identity-verified research inputs."""
import argparse
import importlib.metadata
import json
import multiprocessing as mp
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from revision_core import (HERE, S, E, M, SEEDS, MEAD_COMMIT, require, sha,
    write_json, write_csv, write_bin, full_inputs, prepare_job, plan_units,
    draw_base, intervene, verify_baseline)
from revision_metrics import evaluate, pair_one, mc_summary, ALL_PLANNED, CONDITIONAL, BASELINE
from revision_summary import (read_result_csv, validate_unit_table, aggregate_frames,
    make_plots as summary_make_plots, write_summary_notes, SUMMARY_VERSION)

VERSION = "1.0.2"
H5_SHA = "6ba98a1af8772af08c8cdd5e5e63eebb0890daed48cac2e1b8961dcc67069b77"
STAGE2_SHA = "7ea5ea22bf76660f005d1c30ef40c5d0aaa8c1fffc1c0f24f3a67cd07b2f8011"
ARRAY_KEYS = ["truth","beta","se","leakage","noise_projection","available","covered","p","q"]
JOBS, OPTIONS = [], {}


def code_hashes():
    manifest = json.loads((HERE/"bundle_manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest.items():
        require(sha(HERE/name)==expected,"Source integrity failure: "+name)
    return manifest


def run_r(args, library, log, timeout=14400):
    env = dict(os.environ, R_LIBS_USER=str(library))
    log = Path(log)
    if log.exists():
        shutil.copy2(log, log.with_name(log.name + ".previous." + str(time.time_ns())))
    try:
        with log.open("w",encoding="utf-8") as stream:
            result = subprocess.run(args,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
        require(result.returncode==0,"R infrastructure failed; read "+str(log))
    except (OSError, subprocess.TimeoutExpired, ValueError):
        if log.exists():
            with log.open(encoding="utf-8",errors="replace") as stream:
                tail = "".join(deque(stream, maxlen=100))
            print("\n--- R log tail: " + str(log) + " ---\n" + tail, file=sys.stderr, flush=True)
        raise


def prepare_r(args, output):
    exe = shutil.which(str(args.rscript))
    require(exe is not None,"Rscript not found. Install R >=4.4 and its compiler toolchain, or supply --rscript.")
    library = (args.r_library or output/"R_library").resolve()
    library.mkdir(parents=True,exist_ok=True)
    if args.install_r_packages:
        print("Installing/checking official R dependencies; log:",output/"install_R.log",flush=True)
        run_r([exe,"--vanilla",str(HERE/"install_r.R"),str(library),str(HERE/"vendor")],
              library,output/"install_R.log",timeout=14400)
    version_file = output/"R_versions_current.json"
    run_r([exe,"--vanilla",str(HERE/"r_bridge.R"),"--check",str(library),str(version_file)],
          library,output/"R_check.log",timeout=180)
    versions = json.loads(version_file.read_text())
    run_r([exe,"--vanilla",str(HERE/"r_selftest.R"),str(library),str(output/"R_selftest.json")],
          library,output/"R_selftest.log",timeout=300)
    return dict(executable=exe,library=str(library),versions=versions)


def cached_file(name, roots, expected):
    seen = set()
    for root in roots:
        if root is None or not root.exists():
            continue
        if root.is_file():
            candidates = [root] if root.name==name else []
        else:
            candidates = []
            for base, dirs, names in os.walk(root):
                dirs[:] = [d for d in dirs if d not in (".git",".venv","node_modules","__pycache__","R_library")]
                if name in names:
                    candidates.append(Path(base)/name)
        for path in candidates:
            key = str(path.resolve())
            if key in seen:
                continue
            seen.add(key)
            if path.stat().st_size==S.PUBLIC_SIZE[name] and sha(path)==expected:
                print("REUSE",path,flush=True)
                return path
    raise FileNotFoundError("Missing verified cached file: "+name+". No research download was attempted; pass --external-cache.")


def resolve_all_inputs(args):
    if args.stage2_results:
        require((args.stage2_results/"protocol.json").is_file(),
                "--stage2-results must name the directory directly containing protocol.json and cases/")
    if args.input_bundle:
        b = args.input_bundle.resolve()
        require((b/"inputs_manifest.json").is_file(),"Input bundle manifest missing")
        manifest = json.loads((b/"inputs_manifest.json").read_text())
        for item in manifest:
            path = (b/item["relative_path"]).resolve()
            require(path.is_relative_to(b) and sha(path)==item["sha256"],"Damaged portable input bundle")
        args.libd_dir,args.seaad_dir,args.gene_map = b/"LIBD",b/"SEAAD",b/"gene_identity_map.csv"
        args.external_cache = b/"external"
    folders,gene_map,identities = E.resolve_inputs(args)
    roots = [args.external_cache,args.stage2_results,args.project_root,Path.cwd(),args.project_root.parent]
    meta = cached_file(S.META_NAME,roots,S.PUBLIC_SHA256[S.META_NAME])
    counts = cached_file(S.COUNT_NAME,roots,H5_SHA)
    for path in (meta,counts):
        identities.append(dict(name="external/"+path.name,path=str(path),sha256=sha(path),bytes=path.stat().st_size))
    return folders,gene_map,meta,counts,identities


def export_inputs(identities, destination):
    destination = destination.resolve()
    require(not destination.exists() or not any(destination.iterdir()),"Export destination must be new/empty")
    destination.mkdir(parents=True,exist_ok=True)
    manifest = []
    for item in identities:
        name = "gene_identity_map.csv" if item["name"]=="gene_map" else item["name"]
        dest = destination/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(item["path"],dest)
        require(sha(dest)==item["sha256"],"Export copy checksum failed")
        manifest.append(dict(relative_path=name,sha256=item["sha256"],bytes=dest.stat().st_size))
    write_json(destination/"inputs_manifest.json",manifest)
    (destination/"README.txt").write_text(
        "Verified prepared inputs for the frozen Paper 3 benchmark.\n"
        "Run --input-bundle /path/to/this/folder. No old results are required.\n"
        "Local export is not public deposition. Check original licenses/permissions before sharing.\n"
        "Public reproducibility remains incomplete until these inputs (or their full preparation workflow) are accessible.\n",encoding="utf-8")


def synthetic_jobs(output):
    """Software fixture, deliberately not research evidence."""
    rng = np.random.default_rng(1012)
    g, donors, per_class = 700, ["ref0","ref1","ref2","ref3","target0"], 30
    base = rng.gamma(2.,1.,size=(g,7))+.1
    for k in range(7):
        base[k*42:(k+1)*42,k] += 30
    blocks,meta = [],[]
    for donor in donors:
        donor_mu = base*rng.lognormal(0,.08,size=(g,7))
        for k,ct in enumerate(M.CT):
            blocks.append(rng.poisson(donor_mu[:,k,None]*rng.lognormal(0,.2,size=(1,per_class))))
            meta.extend(dict(cell_id=f"{donor}_{k}_{i}",donor_id=donor,cell_type=ct) for i in range(per_class))
    x = sparse.csc_matrix(np.column_stack(blocks))
    source = S.source_from_counts("SYNTHETIC",np.array([f"gene_{i:04d}" for i in range(g)]),x,
                                  np.asarray(x.sum(0)).ravel(),pd.DataFrame(meta))
    # Match the pandas group-key dtype used by real input loading.
    job = dict(source=source,split=np.int64(0),refs=donors[:4],targets=donors[4:],test=np.arange(300,400),
               panels={"balanced":np.arange(300)})
    return [prepare_job(job,output,500)]


def valid_checkpoint(folder, signature):
    receipt = folder/"complete.json"
    if not receipt.exists():
        return False
    old = json.loads(receipt.read_text())
    require(old["signature"]==signature,"Checkpoint belongs to another protocol")
    for name, expected in old["files"].items():
        require((folder/name).is_file() and sha(folder/name)==expected,"Damaged checkpoint: "+str(folder/name))
    return True


def seal(folder, files, signature):
    write_json(folder/"complete.json",dict(signature=signature,files={name:sha(folder/name) for name in files}))


def batch_r(job, unit, reps, options, work):
    entries, payload, baseline_checks = [], {}, []
    for rep in reps:
        base = draw_base(job,unit,rep,options)
        key = f"r{rep}"
        write_bin(work/(key+"_raw.bin"),base["draw"]["counts"][:,job["norm_genes"]].T)
        write_bin(work/(key+"_A.bin"),base["ref"]["A"][job["panels"]["balanced"]])
        write_bin(work/(key+"_total.bin"),base["draw"]["mass"].sum(1))
        conds = []
        for condition in unit["conditions"]:
            data = intervene(job,base,condition,options)
            payload[(rep,condition["name"])] = data
            file = f"{key}_{condition['name']}_add.bin"
            write_bin(work/file,data["additions"].T)
            conds.append(dict(name=condition["name"],additions=file))
            if rep==0 and condition["name"] in ("null","up05"):
                checks = verify_baseline(job,unit,rep,options,data,condition["effect"])
                baseline_checks.append(dict(rep=rep,condition=condition["name"],max_abs_error=max(checks.values()),**checks))
        entries.append(dict(key=key,n=unit["n"],weights=base["weights"].tolist(),raw=key+"_raw.bin",A=key+"_A.bin",
                            total=key+"_total.bin",conditions=conds))
    # Cache lives within the protocol-owned unit to avoid concurrent RDS writes.
    write_json(work/"jobs.json",dict(reference_folder=job["reference_folder"],entries=entries,
        mead_cache=str(Path(options["output"])/"units"/unit["unit_id"]/"mead_reference_cache")))
    info = options["r_info"]
    run_r([info["executable"],"--vanilla",str(HERE/"r_bridge.R"),str(work),info["library"]],
          info["library"],work.parent/"bridge.log",timeout=options["r_timeout"])
    index = pd.read_csv(work/"outputs.csv",keep_default_na=False)
    require(not index.duplicated(["key","arm"]).any(),"Duplicate R output keys")
    results = {}
    for row in index.itertuples(index=False):
        q = np.loadtxt(work/row.file,delimiter=",",ndmin=2)
        require(q.shape==(unit["n"],7),"R output axes mismatch")
        results[(row.key,row.arm)] = (q,row.reason,row.warnings)
    return payload,results,baseline_checks


def run_chunk(job,unit,reps,options,folder):
    folder.mkdir(parents=True,exist_ok=True)
    records,pairs,arrays,indices,mechanism = [],[],[],[],[]
    with tempfile.TemporaryDirectory(prefix="r_work_",dir=folder) as tmp:
        work = Path(tmp)
        payload,result,checks = batch_r(job,unit,reps,options,work)
        for rep in reps:
            per_condition = {}
            for condition in unit["conditions"]:
                name = condition["name"]
                data = payload[(rep,name)]
                preds = {"oracle":(data["p"],""),"unadjusted":(None,"")}
                for kind,denom in (("observed",data["total"]),("fixed",data["original_total"])):
                    mixture = data["mixture"]*data["total"][:,None]/denom[:,None]
                    q,reason,_ = E.repaired_prediction(data["A"],mixture,"drop_zero_weighted")
                    preds["nnls_"+kind] = (E.convert_available_rows(q,data["ref"]["pooled_yield"]),reason)
                for arm in ("mead_native","mead_lambda0_native"):
                    mead,reason,_ = result[(f"r{rep}",arm)]
                    preds[arm] = (mead,reason)
                    preds[arm.replace("native","reference_converted")] = (E.convert_available_rows(mead,data["ref"]["pooled_yield"]),reason)
                prefix = f"r{rep}_{name}"
                for (key,method),(q,reason,warnings) in result.items():
                    if key==prefix:
                        preds[method] = (E.convert_available_rows(q,data["ref"]["pooled_yield"]),reason)
                        if method in (BASELINE,"dt_fixed_a1_c1000000"):
                            preds[method+"_native"] = (q,reason)
                require(BASELINE in preds,"Baseline R prediction missing")
                condition_rows,condition_arrays = {},{}
                for method,(pp,reason) in preds.items():
                    row,arr = evaluate(data,pp,reason)
                    condition_rows[method],condition_arrays[method] = row,arr
                    attrs = dict(rep=rep,condition=name,method=method)
                    records.append(dict(attrs,**row))
                    indices.append(attrs)
                    arrays.append(np.column_stack([arr[k] for k in ARRAY_KEYS]))
                pairs.extend(dict(rep=rep,condition=name,**p) for p in pair_one(condition_rows,condition_arrays))
                per_condition[name] = (preds,condition_arrays)
            # Direct paired cross-condition diagnostic on genes that remain null.
            p0,a0 = per_condition["null"]
            for name,(preds,arrs) in per_condition.items():
                if name=="null":
                    continue
                alt = payload[(rep,name)]["truth"]!=0
                for method in preds.keys() & p0.keys():
                    x,y = preds[method][0],p0[method][0]
                    err = float(np.max(abs(x-y))) if x is not None and np.isfinite(x).all() and np.isfinite(y).all() else np.nan
                    common = ~alt & (arrs[method]["available"]>0) & (a0[method]["available"]>0)
                    beta_delta = arrs[method]["beta"]-a0[method]["beta"]
                    mechanism.append(dict(rep=rep,condition=name,method=method,proportion_max_change=err,
                        common_null_genes=int(common.sum()),null_beta_mean_change=float(beta_delta[common].mean()) if common.any() else np.nan,
                        null_beta_max_change=float(np.max(abs(beta_delta[common]))) if common.any() else np.nan))
                    if (method in ("oracle","nnls_observed","nnls_fixed","dt_fixed_a1_c1000000") or method.startswith("mead_")) and np.isfinite(err):
                        require(err<1e-8,"Expected scale-invariant prediction changed: "+method)
        write_csv(folder/"replicates.csv",records)
        write_csv(folder/"paired_replicates.csv",pairs)
        write_csv(folder/"mechanism_replicates.csv",mechanism)
        write_csv(folder/"fit_index.csv",indices)
        np.save(folder/"fit_arrays.npy",np.stack(arrays),allow_pickle=False)
        write_json(folder/"baseline_checks.json",checks)
        for name in ("normalization_factors.csv","r_controls.json","outputs.csv"):
            shutil.copy2(work/name,folder/name)
    seal(folder,["replicates.csv","paired_replicates.csv","mechanism_replicates.csv","fit_index.csv","fit_arrays.npy",
                 "baseline_checks.json","normalization_factors.csv","r_controls.json","outputs.csv","bridge.log"],options["signature"])


def summarize_genes(chunks,job,planned):
    all_index = pd.concat([read_result_csv(p/"fit_index.csv") for p in chunks],ignore_index=True)
    require(not all_index.condition.isna().any(), "Missing condition in fit index")
    all_arrays = np.concatenate([np.load(p/"fit_arrays.npy",allow_pickle=False) for p in chunks],axis=0)
    require(len(all_arrays)==len(all_index),"Fit array/index mismatch")
    rows=[]
    keys = {k:i for i,k in enumerate(ARRAY_KEYS)}
    for (condition,method), idx in all_index.groupby(["condition","method"]).groups.items():
        f = all_index.loc[idx]
        require(len(f)==planned and f.rep.nunique()==planned,"Incomplete per-gene repetition set")
        a = all_arrays[np.asarray(idx)]
        truth = a[0,:,keys["truth"]]
        require(np.equal(a[:,:,keys["truth"]],truth).all(),"Estimand changes across repetitions")
        beta,se = a[:,:,keys["beta"]],a[:,:,keys["se"]]
        available = a[:,:,keys["available"]].sum(0)
        def avg(values):
            count = np.isfinite(values).sum(0)
            return np.divide(np.nansum(values,0),count,out=np.full(len(count),np.nan),where=count>0)
        mean = avg(beta)
        sd = np.sqrt(np.divide(np.nansum((beta-mean)**2,axis=0),available-1,
                     out=np.full(len(truth),np.nan),where=available>1))
        covered = a[:,:,keys["covered"]].sum(0)
        table = pd.DataFrame(dict(gene_id=job["source"]["genes"][job["test"]],truth=truth,
            bias=mean-truth,empirical_sd=sd,mean_se=avg(se),rms_se=np.sqrt(avg(se**2)),
            mean_leakage=avg(a[:,:,keys["leakage"]]),mean_noise_projection=avg(a[:,:,keys["noise_projection"]]),
            coverage=np.divide(covered,available,out=np.full(len(truth),np.nan),where=available>0),
            delivered_coverage=covered/planned,availability=available/planned,available=available,planned=planned))
        table["condition"],table["method"] = condition,method
        rows.append(table)
    return pd.concat(rows,ignore_index=True)


def run_unit(unit):
    job,options = JOBS[unit["job_index"]], OPTIONS
    folder = Path(options["output"])/"units"/unit["unit_id"]
    folder.mkdir(parents=True,exist_ok=True)
    if valid_checkpoint(folder,options["signature"]):
        print("RESUME unit",unit["unit_id"],flush=True)
        return unit["unit_id"]
    chunks=[]
    for start in range(0,options["repetitions"],options["batch_size"]):
        reps = list(range(start,min(start+options["batch_size"],options["repetitions"])))
        chunk = folder/f"chunk_{start:05d}"
        if not valid_checkpoint(chunk,options["signature"]):
            print("RUN",unit["unit_id"],unit["target_donor"],start,flush=True)
            run_chunk(job,unit,reps,options,chunk)
        chunks.append(chunk)
    for name in ("replicates.csv","paired_replicates.csv","mechanism_replicates.csv"):
        frame=pd.concat([read_result_csv(c/name) for c in chunks],ignore_index=True)
        validate_unit_table(frame,unit,list(range(options["repetitions"])),Path(name).stem)
        write_csv(folder/name,frame)
    genes=summarize_genes(chunks,job,options["repetitions"])
    validate_unit_table(genes,unit,list(range(options["repetitions"])),"per_gene_summary")
    write_csv(folder/"per_gene_summary.csv",genes)
    write_json(folder/"unit.json",unit)
    seal(folder,["replicates.csv","paired_replicates.csv","mechanism_replicates.csv","per_gene_summary.csv","unit.json"],options["signature"])
    print("DONE",unit["unit_id"],flush=True)
    return unit["unit_id"]


def aggregate(output,units,signature,smoke=False):
    frames={name:[] for name in ("replicates","paired_replicates","per_gene_summary","mechanism_replicates")}
    repetitions=json.loads((output/"protocol.json").read_text())["settings"]["repetitions"]
    for unit in units:
        folder=output/"units"/unit["unit_id"]
        require(valid_checkpoint(folder,signature),"Incomplete unit: "+unit["unit_id"])
        attrs={k:unit[k] for k in ("unit_id","reference_source","target_source","target_donor")}
        for name in frames:
            f=read_result_csv(folder/(name+".csv"))
            validate_unit_table(f,unit,list(range(repetitions)),name)
            for key,val in attrs.items(): f[key]=val
            frames[name].append(f)
    frames={k:pd.concat(v,ignore_index=True) for k,v in frames.items()}
    tables=output/"tables"; tables.mkdir(exist_ok=True)
    summary,donor,status=aggregate_frames(frames,units,repetitions,tables)
    make_plots(tables,summary,donor)
    write_summary_notes(output/"SUMMARY_METHODS.md")
    flag="SOFTWARE_SMOKE_ONLY" if smoke else "COMPLETE"
    write_json(output/"RUN_STATUS.json",dict(status=flag,signature=signature,units=len(units),summary_version=SUMMARY_VERSION,
        all_unavailable_method_conditions=int(status.status.eq("ALL_UNAVAILABLE").sum()),
        limitations=["MEAD point estimates only; no propagation of its proportion covariance to expression HC3 inference",
                    "fixed denominator uses unavailable pre-intervention information",
                    "MC intervals condition on fixed donor banks; not donor-population intervals",
                    "cached inputs still require a public access/deposition solution for external reproducibility"]))
    (output/"READ_RESULTS.md").write_text("# Revision results\n\nStatus: **"+flag+"**\n\n"
        "1. Inspect `tables/method_status.csv` and `failure_reasons.csv` before comparing methods.\n"
        "2. `paired_summary.csv`: method A minus original dtangle, on paired repetitions.\n"
        "3. `per_gene_summary.csv`: bias, empirical SD, estimated SE, coverage, and availability.\n"
        "4. `mechanism_replicates.csv`: injection versus null changes on the same mixtures.\n"
        "5. `donor_summary.csv`: donor heterogeneity; do not pool two references as independent target donors.\n"
        "6. `units/*/chunk_*/normalization_factors.csv`: actual denominators used.\n\n"
        "Failed fits remain in all-planned denominators; unavailable p values are 1 before BH. "
        "Conditional coverage is not delivered coverage. MC intervals do not quantify donor-population uncertainty.\n\n"
        "MEAD uses its official measurement-error-adjusted point estimator on a reference-selected, test-disjoint panel. "
        "Both native output and optional reference-yield conversion are reported, without asserting either scale is calibrated. "
        "Expression regression retains HC3; MEAD's own proportion covariance is not propagated. "
        "No method or parameter is selected by external outcomes.\n",encoding="utf-8")


def make_plots(tables,summary,donor):
    return summary_make_plots(tables,summary,donor)


def audit_legacy(output, units, legacy, repetitions):
    protocol=json.loads((legacy/"protocol.json").read_text())
    require(protocol["script_sha256"]==STAGE2_SHA and protocol["run_signature"]==
            "8568f6d8a2941b190ffecd0363f609613dbeff5aea101ec8932857e16bbd524f",
            "Optional legacy audit requires the reviewed stage-2 run")
    plan=pd.read_csv(legacy/"planned_cases.csv")
    aliases={BASELINE:("dtangle_official","cell_reference_converted"),
             "nnls_observed":("drop_zero_weighted","cell_reference_converted"),
             "oracle":("oracle","cell_fraction"),"unadjusted":("unadjusted","none")}
    metrics=["availability","coverage","delivered_coverage","fdp","any_false","power","gene_error_mean","proportion_rmse"]
    rows=[]
    for unit in units:
        new=read_result_csv(output/"units"/unit["unit_id"]/"replicates.csv")
        for condition,effect in (("null",0.),("up05",.5)):
            selected=plan[plan.reference_source.eq(unit["reference_source"]) & plan.target_source.eq(unit["target_source"])
                & plan.target_donor.eq(unit["target_donor"]) & plan.split.eq(unit["split"]) & plan.n.eq(unit["n"])
                & plan.composition_shift.eq(1.) & plan.effect.eq(effect) & plan.module.isin(["compare","external"])]
            require(len(selected)==1,"Legacy case missing/ambiguous for "+unit["unit_id"])
            folder=legacy/"cases"/selected.case_id.iloc[0]
            receipt=json.loads((folder/"complete.json").read_text())
            require(receipt["run_signature"]==protocol["run_signature"] and sha(folder/"replicates.csv")==receipt["files"]["replicates.csv"],"Legacy checkpoint damaged")
            old=pd.read_csv(folder/"replicates.csv")
            for method,(oldmethod,scale) in aliases.items():
                a=new[new.condition.eq(condition)&new.method.eq(method)].set_index("rep").sort_index()
                b=old[old.method.eq(oldmethod)&old.scale.eq(scale)&old.gate.eq("basic")].set_index("rep").sort_index()
                use=a.index.intersection(b.index)
                require(len(use)==min(repetitions,200),"Incomplete original baseline repetition overlap")
                for metric in metrics:
                    x,y=a.loc[use,metric].to_numpy(),b.loc[use,metric].to_numpy()
                    same=np.allclose(x,y,rtol=1e-6,atol=1e-8,equal_nan=True)
                    delta=abs(x-y)
                    rows.append(dict(unit_id=unit["unit_id"],condition=condition,method=method,metric=metric,
                        paired_replicates=len(use),max_abs_error=float(np.nanmax(delta)) if np.isfinite(delta).any() else np.nan,matched=bool(same)))
    write_csv(output/"tables"/"legacy_baseline_audit.csv",rows)
    require(all(r["matched"] for r in rows),"Legacy baseline mismatch; inspect legacy_baseline_audit.csv before scientific interpretation")


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root",type=Path,default=Path("/home/research/git/repository/paper1_public_evidence"))
    p.add_argument("--output",type=Path,default=Path("paper3_revision_results_v1"))
    for flag in ("libd-dir","seaad-dir","gene-map","external-cache","stage2-results","input-bundle","export-inputs","r-library"):
        p.add_argument("--"+flag,type=Path)
    p.add_argument("--rscript",default="Rscript")
    p.add_argument("--install-r-packages",action="store_true")
    p.add_argument("--workers",type=int,default=8 if os.name!="nt" else 1)
    p.add_argument("--repetitions",type=int,default=200)
    p.add_argument("--batch-size",type=int,default=10)
    p.add_argument("--mead-max-genes",type=int,default=2000)
    p.add_argument("--r-timeout",type=int,default=14400)
    p.add_argument("--smoke",action="store_true")
    p.add_argument("--synthetic",action="store_true",help="Run artificial software fixture, never publication evidence")
    p.add_argument("--preflight",action="store_true")
    p.add_argument("--report-only",action="store_true")
    return p


def main(argv=None):
    global JOBS,OPTIONS
    args=parser().parse_args(argv)
    require(args.workers>=1 and args.batch_size>=1,"Invalid worker/batch count")
    require(os.name!="nt" or args.workers==1,"Windows requires --workers 1; server uses fork sharing")
    require(args.smoke or args.synthetic or args.repetitions>=200,"Publication run requires at least 200 repetitions")
    output=args.output.resolve()
    if output.exists() and any(output.iterdir()):
        require((output/"OWNER.json").is_file(),"Refusing to use a non-owned/nonempty output directory")
    output.mkdir(parents=True,exist_ok=True)
    E.hold_output_lock(output)
    owner=output/"OWNER.json"
    if not owner.exists(): write_json(owner,dict(package="paper3_revision",version=VERSION))
    source=code_hashes()
    require(sha(HERE/"vendor"/"stage2_frozen.py")==STAGE2_SHA,"Frozen engine identity mismatch")
    if args.report_only:
        protocol=json.loads((output/"protocol.json").read_text())
        require(protocol["code_hashes"]==source,"Report-only source differs")
        aggregate(output,protocol["units"],protocol["signature"],protocol["settings"]["software_smoke"])
        return 0
    identities=[]
    if not args.synthetic:
        folders,gene_map,meta,counts,identities=resolve_all_inputs(args)
        if args.export_inputs:
            export_inputs(identities,args.export_inputs)
            print("EXPORTED verified prepared inputs:",args.export_inputs,flush=True)
            return 0
    r_info=prepare_r(args,output)
    if args.preflight:
        write_json(output/"preflight.json",dict(inputs=identities,r_info=r_info,synthetic=args.synthetic))
        print("PREFLIGHT PASSED; no research experiments run",flush=True)
        return 0
    if args.synthetic:
        jobs=synthetic_jobs(output)
    else:
        sources,jobs=full_inputs(folders,gene_map)
        active=next(iter(sources.values()))["genes"]
        external=S.read_external(counts,meta,active,output/"external_inputs")
        jobs += S.make_external_jobs(jobs,external,output)
        jobs=[prepare_job(j,output,args.mead_max_genes) for j in jobs]
    smoke=args.smoke or args.synthetic
    units=plan_units(jobs,smoke=smoke)
    repetitions=min(args.repetitions,3) if smoke else args.repetitions
    settings=dict(SEEDS,repetitions=repetitions,batch_size=args.batch_size,n=80,composition_shift=1.,
                  mead_max_genes=args.mead_max_genes,software_smoke=smoke,
                  mead_calc_var=False,MEAD_commit=MEAD_COMMIT,gate="basic",synthetic=args.synthetic)
    pyversions={p:importlib.metadata.version(p) for p in ("numpy","pandas","scipy","h5py","matplotlib")}
    reference_hashes={str(p.relative_to(output)).replace("\\","/"):sha(p) for j in jobs
                      for p in Path(j["reference_folder"]).iterdir() if p.is_file()}
    body=dict(version=VERSION,settings=settings,units=units,code_hashes=source,
              input_hashes={i["name"]:i["sha256"] for i in identities},reference_hashes=reference_hashes,
              python_packages=pyversions,R_packages=r_info["versions"],platform=platform.platform())
    signature=M.digest(body)
    old=output/"protocol.json"
    if old.exists():
        require(json.loads(old.read_text())["signature"]==signature,"Protocol/environment changed; choose a NEW output folder")
    write_json(old,dict(body,signature=signature))
    write_json(output/"input_paths.json",identities)
    write_json(output/"RUN_STATUS.json",dict(status="RUNNING",signature=signature,units=len(units)))
    JOBS=jobs
    OPTIONS=dict(settings,signature=signature,output=str(output),r_info=r_info,r_timeout=args.r_timeout)
    try:
        if args.workers==1:
            for unit in units: run_unit(unit)
        else:
            with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context("fork")) as pool:
                futures=[pool.submit(run_unit,unit) for unit in units]
                for future in as_completed(futures): future.result()
        aggregate(output,units,signature,smoke)
        if args.stage2_results and not args.synthetic:
            audit_legacy(output,units,args.stage2_results,repetitions)
    except BaseException as exc:
        write_json(output/"RUN_STATUS.json",dict(status="INCOMPLETE",signature=signature,error=str(exc)))
        raise
    print("FINISHED",output,flush=True)
    return 0
