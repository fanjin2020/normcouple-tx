"""GSE174367 audit and optional genuine-count application; never invent counts."""
import gzip
import json
import os
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from revision_core import S, E, M, require, sha, write_json, write_csv, write_bin
from extension_stats import fit_hc3

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/"
BULK_NAME = "GSE174367_bulkRNA_processed.rda.gz"
BULK_SHA256 = "3e453a2739da3da76827ba5d40c6e4ae3e933ecdb25f705b3e8a38fa763ef227"
META_NAME = "GSE174367_snRNA-seq_cell_meta.csv.gz"
SKIP_DIRS = {".git", ".venv", "node_modules", "R_library", "__pycache__"}


def find_named(name, roots):
    found = []
    for root in roots:
        if root is None or not Path(root).exists():
            continue
        for base, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
            if name in files:
                found.append(Path(base) / name)
    return sorted(set(p.resolve() for p in found))


def download(url, path, expected=None, max_bytes=400_000_000):
    path = Path(path)
    if path.exists():
        if expected is None or sha(path) == expected:
            return path
        raise ValueError("Cached download checksum mismatch: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Paper3-research-reproducibility/1.0"})
            with urllib.request.urlopen(request, timeout=90) as src, tmp.open("wb") as dest:
                if int(src.headers.get("Content-Length", 0)) > max_bytes:
                    raise ValueError("Download exceeds declared size limit")
                total = 0
                while True:
                    block = src.read(1024 * 1024)
                    if not block:
                        break
                    total += len(block)
                    if total > max_bytes:
                        raise ValueError("Download exceeds declared size limit")
                    dest.write(block)
            if expected is not None and sha(tmp) != expected:
                raise ValueError("Downloaded checksum mismatch: " + url)
            os.replace(tmp, path)
            return path
        except (OSError, ValueError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def acquire(name, args, output, expected=None):
    candidates = find_named(name, [args.project_root, args.external_cache])
    for path in candidates:
        if expected is None or sha(path) == expected:
            print("REUSE", path, flush=True)
            return path
    if args.offline:
        raise FileNotFoundError("Missing cached public file in offline mode: " + name)
    print("DOWNLOAD", name, flush=True)
    limit = 8_000_000_000 if name.endswith(".h5") else 400_000_000
    return download(BASE + name, output / "downloads" / name, expected, max_bytes=limit)


def audit_samples(meta, cells):
    require(meta.SampleID.is_unique, "Duplicate bulk SampleID")
    require(meta[["Case.Year", "Case.Num", "Region"]].notna().all().all(), "Missing bulk donor or region identity")
    meta = meta.copy()
    meta["donor_key"] = meta["Case.Year"].astype(int).astype(str) + "-" + meta["Case.Num"].astype(int).astype(str)
    require(meta.donor_key.is_unique, "Repeated bulk donors require a repeated-measures model; refusing independent fits")
    summaries = []
    for sid, f in cells.groupby("SampleID", sort=True):
        row = dict(SampleID=sid, nuclei=len(f))
        for key in ("Diagnosis", "Age", "Sex"):
            values = f[key].dropna().unique()
            require(len(values) == 1, "Conflicting single-nucleus donor metadata: " + sid + "/" + key)
            row["sn_" + key] = values[0]
        summaries.append(row)
    paired = meta.merge(pd.DataFrame(summaries), on="SampleID", how="outer", indicator=True, validate="one_to_one")
    both = paired._merge.eq("both")
    paired["metadata_agree"] = False
    paired.loc[both, "metadata_agree"] = (
        paired.loc[both, "Diagnosis"].eq(paired.loc[both, "sn_Diagnosis"]) &
        paired.loc[both, "Sex"].eq(paired.loc[both, "sn_Sex"]) &
        pd.to_numeric(paired.loc[both, "Age"]).eq(pd.to_numeric(paired.loc[both, "sn_Age"])))
    paired["pairing_basis"] = "Exact deposited SampleID plus diagnosis/age/sex; not inferred from ordering"
    return meta, paired


def design(meta):
    fields = ["Age", "Sex", "PMI", "RIN", "Sequencing.Group", "Diagnosis"]
    keep = meta[fields].notna().all(1) & meta.Region.eq("FC") & meta.Diagnosis.isin(["AD", "Control"])
    audit = meta[["SampleID", "donor_key", "Region", *fields]].copy()
    audit["included"] = keep
    audit["exclusion"] = ["" if ok else "non-FC, unsupported diagnosis, or missing declared covariate" for ok in keep]
    m = meta.loc[keep].reset_index(drop=True)
    require(len(m) >= 25 and set(m.Diagnosis) == {"AD", "Control"}
            and m.Diagnosis.value_counts().min() >= 8, "Insufficient complete-case disease groups")
    c = pd.get_dummies(m[["Sex", "Sequencing.Group"]].astype(str), drop_first=True, dtype=float)
    for field in ("Age", "PMI", "RIN"):
        v = m[field].to_numpy(float)
        c[field] = (v - v.mean()) / v.std(ddof=1) if v.std(ddof=1) > 0 else 0
    return m, m.Diagnosis.eq("AD").to_numpy(float), c.to_numpy(float), list(c.columns), audit


def association_table(y, genes, x, cov, props, method):
    f = fit_hc3(y, x, props, covariates=cov)
    p = np.ones(y.shape[1])
    good = f["good"] & f["basic"]
    p[good] = 2 * stats.t.sf(abs(f["beta"][good] / f["se"][good]), f["df"])
    frame = pd.DataFrame(dict(gene_id=genes, method=method, beta=np.where(good, f["beta"], np.nan),
        se=np.where(good, f["se"], np.nan), p=p, q=M.bh(p), available=good))
    diag = {k: f[k] for k in ("reason", "df", "composition_rank", "nuisance_rank", "exposure_r2", "leverage", "exposure_estimable")}
    return frame, dict(method=method, n=len(x), available_genes=int(good.sum()), planned_genes=len(good), **diag)


def raw_application(args, output, jobs, meta, cells, paired, r_info, manifest):
    from revision_driver import run_r
    path = Path(manifest["counts_path"]).expanduser().resolve()
    require(manifest.get("accession") == "GSE174367" and manifest.get("measurement") == "raw_gene_counts",
            "Raw manifest must explicitly identify GSE174367 raw_gene_counts")
    require(sha(path) == manifest["sha256"], "Raw bulk file identity mismatch")
    require(bool(manifest.get("provenance")), "Raw bulk provenance is required")
    raw = pd.read_csv(path, index_col=0)
    raw.index = raw.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    require(raw.index.is_unique and raw.columns.is_unique, "Duplicate raw gene/sample identifiers")
    values = raw.to_numpy(float)
    require(np.isfinite(values).all() and (values >= 0).all() and np.allclose(values, np.rint(values), atol=1e-6, rtol=0),
            "Only genuine nonnegative integer gene counts are accepted; no exponentiated processed matrix")
    m, x, cov, cov_names, exclusions = design(meta)
    require(set(m.SampleID) <= set(raw.columns), "Raw count columns do not cover the fixed complete-case cohort")
    raw = raw.loc[:, m.SampleID]
    total = raw.sum(0).to_numpy(float)
    require((total > 0).all(), "Empty bulk libraries")
    # Fix the outcome scale and gene family before varying composition inputs.
    all_markers = set().union(*(set(j["source"]["genes"][j["panels"]["balanced"]]) for j in jobs if "target_source" not in j))
    eligible = ((raw.to_numpy() / total * 1e6 >= 1).sum(1) >= np.ceil(len(m) / 2)) & ~raw.index.isin(all_markers)
    genes = raw.index[eligible].to_numpy()
    require(len(genes) >= 100, "Too few test genes after blind expression and marker exclusion")
    # TMM outcome is fixed across deconvolution arms, not selected by association results.
    matrix_path = output / "bulk_counts_for_R.csv"
    raw.to_csv(matrix_path)
    run_r([r_info["executable"], "--vanilla", str(Path(__file__).parent / "bulk_outcome.R"),
           str(matrix_path), str(output / "outcome_denominator.csv"), r_info["library"]],
          r_info["library"], output / "outcome.log")
    den = pd.read_csv(output / "outcome_denominator.csv").denominator.to_numpy()
    require(len(den) == len(m) and np.isfinite(den).all() and (den > 0).all(), "Invalid outcome denominator")
    y = np.log2(raw.loc[genes].to_numpy().T / den[:, None] * 1e6 + 1)
    tables, diagnostics, predictions = [], [], []
    t, d = association_table(y, genes, x, cov, None, "no_composition")
    tables.append(t); diagnostics.append(d)
    for job in (j for j in jobs if "target_source" not in j):
        work = output / ("reference_" + job["source"]["label"])
        work.mkdir(exist_ok=True)
        ref = job["base"]
        markers = job["panels"]["balanced"]
        marker_names = job["source"]["genes"][markers]
        require(set(marker_names) <= set(raw.index), "Raw bulk lacks frozen markers; no silent imputation")
        A = ref["A"][markers]
        anchor = float((ref["counts"][markers].sum(0)).sum() / (ref["counts"][markers].sum(0) / A.sum(0)).sum())
        # Real bulk uses its full deposited count universe for normalizer estimation.
        # Only marker identity must match the frozen reference; unused reference
        # genes are not silently inserted as zero-count bulk genes.
        axes = json.loads((Path(job["reference_folder"]) / "axes.json").read_text())
        axes.update(normalization_genes=len(raw), marker_rows=(raw.index.get_indexer(marker_names)+1).tolist(),
                    test_rows=[1], mead_rows=[])
        local_ref = work / "reference_axes"; local_ref.mkdir(exist_ok=True)
        write_json(local_ref / "axes.json", axes)
        write_bin(work / "raw.bin", raw.to_numpy())
        write_bin(work / "A.bin", A)
        write_bin(work / "total.bin", total)
        write_bin(work / "zero.bin", np.zeros((1, len(m))))
        write_json(work / "jobs.json", dict(reference_folder=str(local_ref), omit_mead=True, real_bulk=True,
            entries=[dict(key="bulk", n=len(m), raw="raw.bin", A="A.bin", total="total.bin", marker_anchor=anchor,
                          conditions=[dict(name="observed", additions="zero.bin", marker_stress=False)])]))
        run_r([r_info["executable"], "--vanilla", str(Path(__file__).parent / "bridge_extension.R"), str(work), r_info["library"]],
              r_info["library"], work / "bridge.log")
        for row in pd.read_csv(work / "outputs.csv", keep_default_na=False).itertuples():
            q = np.loadtxt(work / row.file, delimiter=",", ndmin=2)
            yield_ = ref["pooled_yield"] * A.sum(0) if row.arm in ("marker_matched", "marker_unmatched_scale") else ref["pooled_yield"]
            p = E.convert_available_rows(q, yield_)
            method = job["source"]["label"] + "_" + row.arm
            t, d = association_table(y, genes, x, cov, p, method)
            d["estimator_reason"] = row.reason
            tables.append(t); diagnostics.append(d)
            pp = pd.DataFrame(p, columns=M.CT); pp.insert(0, "SampleID", m.SampleID); pp["method"] = method
            predictions.append(pp)
    association = pd.concat(tables, ignore_index=True)
    prop = pd.concat(predictions, ignore_index=True)
    write_csv(output / "disease_associations.csv", association)
    write_csv(output / "design_diagnostics.csv", diagnostics)
    write_csv(output / "composition_predictions.csv", prop)
    write_csv(output / "cohort_exclusions.csv", exclusions)
    write_csv(output / "analysis_cohort.csv", m)
    compare = []
    for source in ("LIBD", "SEAAD"):
        baseline = association[association.method.eq(source + "_observed")].set_index("gene_id")
        for method, f in association[association.method.str.startswith(source + "_")].groupby("method"):
            b = f.set_index("gene_id").reindex(baseline.index)
            valid = baseline.available & b.available
            a_sig, b_sig = baseline.q.le(.05) & baseline.available, b.q.le(.05) & b.available
            union = (a_sig | b_sig).sum()
            compare.append(dict(method=method, reference=source + "_observed", common_genes=int(valid.sum()),
                sign_agreement=float((np.sign(baseline.beta[valid]) == np.sign(b.beta[valid])).mean()),
                median_absolute_beta_change=float(abs(b.beta[valid] - baseline.beta[valid]).median()),
                significant_genes=int(b_sig.sum()), gained=int((b_sig & ~a_sig).sum()), lost=int((a_sig & ~b_sig).sum()),
                significance_jaccard=float((a_sig & b_sig).sum() / union) if union else np.nan))
    write_csv(output / "method_stability.csv", compare)
    # Recovered nuclei are a paired descriptive comparator, not ground-truth tissue composition.
    mapping = {"ASC": "Astro", "END": "EndoMural", "PER.END": "EndoMural", "PER/END": "EndoMural",
               "EX": "Excit", "INH": "Inhib", "MG": "Micro", "ODC": "Oligo", "OPC": "OPC"}
    cc = cells.copy(); cc["mapped"] = cc["Cell.Type"].map(mapping)
    require(cc.mapped.notna().all(), "Unmapped single-nucleus class; inspect Cell.Type before concordance")
    usable = paired.loc[paired.metadata_agree, "SampleID"]
    recovered = pd.crosstab(cc.loc[cc.SampleID.isin(usable), "SampleID"], cc.loc[cc.SampleID.isin(usable), "mapped"]).reindex(columns=M.CT, fill_value=0)
    recovered = recovered.div(recovered.sum(1), axis=0)
    recovered.to_csv(output / "recovered_nucleus_fractions.csv")
    concordance = []
    for method, f in prop.groupby("method"):
        f = f.set_index("SampleID"); ids = f.index.intersection(recovered.index)
        for ct in M.CT:
            valid = f.loc[ids, ct].notna()
            a, b = f.loc[ids[valid], ct], recovered.loc[ids[valid], ct]
            rho = stats.spearmanr(a, b).statistic if len(a) >= 5 and a.std() > 0 and b.std() > 0 else np.nan
            concordance.append(dict(method=method, cell_type=ct, n=len(a), spearman=rho,
                interpretation="Agreement with recovered nuclei only; capture/platform bias remains"))
    write_csv(output / "paired_nucleus_concordance.csv", concordance)
    from extension_figures import bulk_figures
    bulk_figures(output)
    write_json(output / "application_protocol.json", dict(outcome="fixed TMM log2-CPM+1", covariates=cov_names,
        region="FC", disease="AD minus Control", multiple_testing="BH over the same non-marker gene family per model",
        n=len(m), genes=len(genes), count_provenance=manifest, true_fdr_evaluated=False,
        fixed_denominator_comparison=False, rank_rule="actual nuisance column space"))
    return dict(status="COMPLETE_RAW_COUNT_APPLICATION", samples=len(m), genes=len(genes))


def run_bulk(args, output, jobs, r_info):
    from revision_driver import run_r
    output.mkdir(parents=True, exist_ok=True)
    rda = acquire(BULK_NAME, args, output, BULK_SHA256)
    meta_path = acquire(META_NAME, args, output, S.PUBLIC_SHA256[META_NAME])
    run_r([r_info["executable"], "--vanilla", str(Path(__file__).parent / "export_bulk.R"), str(rda), str(output)],
          r_info["library"], output / "export.log")
    meta, paired = audit_samples(pd.read_csv(output / "bulk_metadata.csv"), pd.read_csv(meta_path))
    write_csv(output / "bulk_snRNA_pairing_audit.csv", paired)
    write_json(output / "public_inputs.json", [{"path": str(p), "sha256": sha(p), "url": BASE + p.name} for p in (rda, meta_path)])
    m, x, cov, _, excluded = design(meta)
    processed = pd.read_csv(output / "bulk_processed_expression.csv", index_col=0)
    table, diag = association_table(processed.loc[:, m.SampleID].to_numpy().T, processed.index, x, cov, None, "deposited_processed_expression_only")
    write_csv(output / "processed_expression_associations_EXPLORATORY.csv", table)
    write_csv(output / "cohort_exclusions.csv", excluded)
    manifest_path = args.bulk_manifest or args.project_root / "paper3_bulk_raw_manifest.json"
    if not manifest_path.exists():
        status = dict(status="RAW_COUNTS_REQUIRED", processed_samples=len(meta), matched_sn_samples=int(paired.metadata_agree.sum()),
            reason="GEO rda contains normExpr.reg and targets, not raw counts. Normalization comparison was NOT run.",
            required_manifest=str(manifest_path), processed_only_analysis="Exploratory; original preprocessing cannot be undone",
            public_raw_resource="https://www.synapse.org/Synapse:syn22079621",
            action="Supply a provenance-verified raw gene-count CSV and rerun --only bulk; never invert the processed expression.")
    else:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        cp = Path(manifest["counts_path"])
        if not cp.is_absolute(): manifest["counts_path"] = str((manifest_path.parent / cp).resolve())
        status = raw_application(args, output, jobs, meta, pd.read_csv(meta_path), paired, r_info, manifest)
    write_json(output / "STATUS.json", status)
    return status
