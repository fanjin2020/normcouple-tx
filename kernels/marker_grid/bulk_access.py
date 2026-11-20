"""Acquire explicit raw-count inputs without inferring provenance from integers."""
import csv
import gzip
from html.parser import HTMLParser
import os
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from common import read_json, require, sha, write_json

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/"
SYNAPSE = "https://www.synapse.org/Synapse:syn22079621"
SKIP = {".git", ".venv", "node_modules", "R_library", "__pycache__", "units", "site-packages"}


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.names = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            link = dict(attrs).get("href", "")
            if link.startswith("GSE174367_") and "/" not in link:
                self.names.append(link)


def inventory_public(output, offline=False):
    record = dict(accession="GSE174367", geo_url=GEO, synapse_url=SYNAPSE,
                  raw_counts_obtained=False, restricted_access_bypassed=False)
    if offline:
        record["status"] = "OFFLINE_NOT_QUERIED"
    else:
        try:
            req = urllib.request.Request(GEO, headers={"User-Agent": "Paper3-reproducibility-audit/1.0"})
            with urllib.request.urlopen(req, timeout=25) as response:
                body = response.read(1024 * 1024 + 1)
            require(len(body) <= 1024 * 1024, "Unexpectedly large GEO listing")
            parser = Links()
            parser.feed(body.decode("utf-8"))
            record.update(status="LISTED", files=sorted(set(parser.names)))
            record["bulk_files"] = [n for n in record["files"] if "bulk" in n.lower()]
            record["note"] = ("A directory listing does not certify measurement type. "
                              "bulkRNA_processed.rda.gz contains processed expression, not the required raw counts.")
        except (OSError, ValueError) as exc:
            record.update(status="CATALOG_UNAVAILABLE", error=type(exc).__name__)
    write_json(output / "public_resource_inventory.json", record)
    return record


def scan_local(roots, excluded, max_entries=250000):
    manifests, candidates, seen = [], [], set()
    visited = 0
    for root in roots:
        if root is None or not root.is_dir():
            continue
        root = root.resolve()
        for base, dirs, names in os.walk(root, followlinks=False):
            p = Path(base)
            dirs[:] = [d for d in sorted(dirs) if d not in SKIP and not (p / d).is_symlink()
                       and not any((p / d).resolve().is_relative_to(e) for e in excluded)]
            for name in sorted(names):
                visited += 1
                if visited > max_entries:
                    return manifests, candidates, dict(status="SCAN_LIMIT_REACHED", visited=visited)
                path = p / name
                if path.is_symlink() or str(path) in seen:
                    continue
                seen.add(str(path))
                lower = name.lower()
                if lower == "paper3_bulk_raw_manifest.json" or lower.endswith(".raw_counts.json"):
                    manifests.append(path)
                context = str(path).lower()
                if ("174367" in context or "bulk" in context) and lower.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz")):
                    if ("count" in lower or "raw" in lower) and not any(t in lower for t in ("snrna", "atac", "normexpr", "processed")):
                        candidates.append(dict(path=str(path), bytes=path.stat().st_size,
                                               status="UNVERIFIED_CANDIDATE_NOT_AUTO_USED"))
    return manifests, candidates, dict(status="COMPLETE", visited=visited)


def download_explicit(url, expected, output, max_bytes):
    parts = urllib.parse.urlsplit(url)
    require(parts.scheme == "https" and parts.hostname and not parts.username and not parts.password,
            "Count download requires HTTPS without credentials in the URL")
    require(not parts.query and not parts.fragment, "Use a stable provenance URL, not a signed/secret URL")
    require(bool(re.fullmatch(r"[a-fA-F0-9]{64}", expected or "")), "Expected download SHA256 is required")
    suffix = ".tsv.gz" if parts.path.endswith(".tsv.gz") else ".csv.gz" if parts.path.endswith(".gz") else ".tsv" if parts.path.endswith((".tsv", ".txt")) else ".csv"
    destination = output / ("raw_counts_" + expected[:16] + suffix)
    if destination.exists():
        require(sha(destination) == expected.lower(), "Cached raw-count download was modified")
        return destination
    req = urllib.request.Request(url, headers={"User-Agent": "Paper3-reproducibility-audit/1.0"})
    temp = destination.with_name(destination.name + ".partial")
    with urllib.request.urlopen(req, timeout=60) as response:
        require(urllib.parse.urlsplit(response.url).scheme == "https", "Refusing non-HTTPS redirect")
        require(int(response.headers.get("Content-Length", 0)) <= max_bytes, "Raw-count file exceeds download limit")
        size = 0
        with temp.open("wb") as target:
            for block in iter(lambda: response.read(1024 * 1024), b""):
                size += len(block)
                require(size <= max_bytes, "Raw-count download exceeds limit; raw sequencing reads are not downloaded")
                target.write(block)
    require(sha(temp) == expected.lower(), "Downloaded raw counts failed SHA256 verification")
    os.replace(temp, destination)
    return destination


def normalize_counts(path, cohort_path, destination, sample_map=None):
    """Lossless shape/ID conversion only: no rounding, imputation or inversion."""
    path = Path(path)
    sep = "\t" if path.name.lower().endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter=sep))
    require(len(header) >= 2 and all(header[1:]) and len(set(header[1:])) == len(header) - 1,
            "Missing/duplicate raw-count sample IDs")
    require(not any(t in path.name.lower() for t in ("normexpr", "tpm", "fpkm", "processed")),
            "Processed-expression input is not raw counts; do not rename or invert it")
    frame = pd.read_csv(path, sep=sep, index_col=0)
    frame.index = frame.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    require(frame.index.str.match(r"^ENSG\d+$").all() and frame.index.is_unique,
            "Unique Ensembl gene IDs required after version removal")
    require(len(frame) >= 100, "Too few raw-count genes")
    if sample_map:
        mapping = pd.read_csv(sample_map, dtype=str)
        require({"raw_sample", "SampleID"} <= set(mapping), "Sample map needs raw_sample,SampleID columns")
        require(mapping[["raw_sample", "SampleID"]].notna().all().all() and mapping.raw_sample.is_unique
                and mapping.SampleID.is_unique, "Sample mapping must be explicit and one-to-one")
        require(set(mapping.raw_sample) <= set(frame.columns), "Sample map references missing raw columns")
        frame = frame.rename(columns=dict(zip(mapping.raw_sample, mapping.SampleID)))
        require(frame.columns.is_unique, "Sample renaming created duplicate IDs")
    cohort = pd.read_csv(cohort_path, dtype={"SampleID": str})
    selected = cohort.included.astype(str).str.lower().eq("true")
    ids = cohort.loc[selected, "SampleID"].tolist()
    require(ids and len(ids) == len(set(ids)) and set(ids) <= set(frame.columns),
            "Raw counts do not cover the prespecified complete-case cohort; no sample-order matching")
    values = frame.to_numpy(dtype=float)
    require(np.isfinite(values).all() and (values >= 0).all() and (values <= 2**53).all()
            and np.equal(values, np.floor(values)).all(), "Nonnegative exact integer gene counts required")
    require((frame.loc[:, ids].sum(axis=0) > 0).all(), "A complete-case sample has an empty count library")
    temp = destination.with_name(destination.name + ".partial")
    frame.to_csv(temp, index_label="gene_id")
    os.replace(temp, destination)
    return dict(genes=len(frame), sample_columns=len(frame.columns), cohort_samples=len(ids),
                transformation="Ensembl version removal and optional explicit sample-ID mapping only")


def prepare_manifest(args, results, output):
    output.mkdir(parents=True, exist_ok=True)
    manifests, candidates, scan = scan_local([args.project_root, args.run_root], [args.output.resolve()])
    write_json(output / "local_scan.json", dict(**scan, manifests=list(map(str, manifests)), candidates=candidates))
    manifest_path = args.bulk_manifest
    if manifest_path is None and not args.raw_counts and not args.counts_url:
        unique = {}
        for p in manifests:
            try:
                m = read_json(p)
                if m.get("accession") == "GSE174367" and m.get("measurement") == "raw_gene_counts":
                    unique.setdefault(m.get("sha256", "MISSING"), []).append(p)
            except (OSError, ValueError):
                continue
        if len(unique) == 1:
            manifest_path = next(iter(unique.values()))[0]
        elif len(unique) > 1:
            return None, dict(status="AMBIGUOUS_RAW_MANIFESTS", action="Pass --bulk-manifest explicitly")
    if manifest_path:
        manifest_path = manifest_path.expanduser().resolve()
        m = read_json(manifest_path)
        require(m.get("accession") == "GSE174367" and m.get("measurement") == "raw_gene_counts", "Invalid raw-count manifest")
        require(not m.get("software_fixture"), "Software fixtures are not biological raw-count evidence")
        count_path = Path(m["counts_path"]).expanduser()
        if not count_path.is_absolute():
            count_path = manifest_path.parent / count_path
        require(count_path.is_file() and sha(count_path) == m.get("sha256"), "Raw-count manifest file/hash mismatch")
        provenance = m.get("provenance")
        source_sha = m["sha256"]
        map_path = args.sample_map
        if not map_path and m.get("sample_map"):
            map_path = Path(m["sample_map"])
            if not map_path.is_absolute():
                map_path = manifest_path.parent / map_path
    elif args.raw_counts or args.counts_url:
        provenance = args.counts_provenance
        require(provenance, "--counts-provenance must describe the actual raw gene-count resource and version")
        require(not (args.raw_counts and args.counts_url), "Choose local raw counts or one explicit URL")
        if args.raw_counts:
            count_path = args.raw_counts.expanduser().resolve()
            require(count_path.is_file(), "Raw count file missing")
        else:
            require(not args.offline, "Offline mode prohibits count download")
            count_path = download_explicit(args.counts_url, args.counts_sha256, output, args.max_download_gb * 1024**3)
        source_sha = sha(count_path)
        if args.counts_sha256:
            require(source_sha == args.counts_sha256.lower(), "Raw count SHA256 mismatch")
        map_path = args.sample_map
    else:
        inventory_public(output, args.offline)
        write_json(output / "RAW_COUNTS_MANIFEST_TEMPLATE.json", dict(accession="GSE174367",
                   measurement="raw_gene_counts", counts_path="REPLACE_WITH_ACTUAL_COUNTS.csv",
                   sha256="REPLACE_WITH_SHA256", provenance="REPLACE_WITH_RESOURCE_VERSION_AND_MEASUREMENT_DESCRIPTION"))
        return None, dict(status="RAW_COUNTS_REQUIRED", candidates=len(candidates),
                          action="Obtain authorized raw gene counts and their sample IDs; supply --bulk-manifest or --raw-counts plus provenance",
                          geo_url=GEO, synapse_url=SYNAPSE, processed_expression_not_substituted=True)
    require(isinstance(provenance, str) and provenance.strip(), "Actual raw-count provenance is required")
    normalized = output / "verified_raw_gene_counts.csv"
    info = normalize_counts(count_path, results / "real_bulk" / "cohort_exclusions.csv", normalized, map_path)
    manifest = dict(accession="GSE174367", measurement="raw_gene_counts", counts_path=str(normalized.resolve()),
                    sha256=sha(normalized), provenance=provenance, original_sha256=source_sha,
                    original_path=str(count_path), **info)
    if map_path:
        manifest.update(sample_map_original=str(map_path), sample_map_sha256=sha(map_path))
    final = output / "verified_raw_manifest.json"
    write_json(final, manifest)
    return final, dict(status="RAW_COUNTS_VALIDATED", **info, sha256=manifest["sha256"])
