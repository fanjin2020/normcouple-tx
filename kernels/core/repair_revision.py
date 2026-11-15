#!/usr/bin/env python3
"""Rebuild Paper 3 summaries from sealed results; never simulate or call R."""
import argparse
import hashlib
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from revision_summary import (SUMMARY_VERSION, TABLE_NAMES, read_result_csv, validate_unit_table,
    recover_gene_null_label, aggregate_frames, make_plots, write_summary_notes, write_csv, require)

HERE = Path(__file__).resolve().parent
UNIT_FILES = ("replicates.csv", "paired_replicates.csv", "mechanism_replicates.csv", "per_gene_summary.csv", "unit.json")
CHUNK_FILES = ("replicates.csv", "paired_replicates.csv", "mechanism_replicates.csv", "fit_index.csv",
    "fit_arrays.npy", "baseline_checks.json", "normalization_factors.csv", "r_controls.json", "outputs.csv", "bridge.log")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")


def checked_path(root, relative):
    path = root / relative
    require(path.is_file() and path.resolve().is_relative_to(root), "Missing/unsafe input: " + str(path))
    return path


def checkpoint(root, relative, signature, expected_names, verify_names):
    receipt = read_json(checked_path(root, relative / "complete.json"))
    require(receipt.get("signature") == signature, "Checkpoint signature mismatch: " + str(relative))
    files = receipt.get("files", {})
    require(set(expected_names).issubset(files), "Incomplete receipt: " + str(relative))
    for name, expected in files.items():
        require(isinstance(name, str) and Path(name).name == name and "\\" not in name,
                "Unsafe checkpoint filename")
        require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected), "Invalid SHA256")
        path = checked_path(root, relative / name)
        if name in verify_names:
            require(sha(path) == expected, "Damaged source file: " + str(path))
    return receipt


def repair(source, output, allow_smoke=False, include_replicates=False, plots=True):
    source, output = source.expanduser().resolve(), output.expanduser().resolve()
    require(not source.is_relative_to(output) and not output.is_relative_to(source), "Source/output must not overlap")
    require(not output.exists(), "Output exists; select a NEW directory (originals will not be overwritten)")
    protocol = read_json(checked_path(source, Path("protocol.json")))
    original_status = read_json(checked_path(source, Path("RUN_STATUS.json")))
    smoke = original_status.get("status") == "SOFTWARE_SMOKE_ONLY"
    require(original_status.get("status") == "COMPLETE" or (allow_smoke and smoke),
            "Source is not a completed run; use the full paper3_revision_results_v2 directory, not old v1 or compact export")
    require(smoke or not (protocol["settings"].get("synthetic") or protocol["settings"].get("software_smoke")),
            "Synthetic data must not be labeled as publication evidence")
    signature = protocol["signature"]
    require(original_status.get("signature") == signature, "Source protocol/status mismatch")
    units = protocol["units"]
    require(len(units) > 0 and len(units) == original_status.get("units"), "Invalid unit count")
    require(len({unit["unit_id"] for unit in units}) == len(units), "Duplicate unit IDs")
    require(all(re.fullmatch(r"revision_[0-9a-f]+", unit["unit_id"]) for unit in units), "Unsafe unit ID")
    repetitions = protocol["settings"]["repetitions"]
    batch_size = protocol["settings"]["batch_size"]
    require(type(repetitions) is int and repetitions > 0 and type(batch_size) is int and batch_size > 0,
            "Invalid repetition/batch counts")
    require((source / "units").is_dir(), "Full server results with units/ are required; compact export is insufficient")
    output.mkdir(parents=True)
    audit = output / "audit"
    audit.mkdir()
    write_json(output / "EXPORT_STATUS.json", dict(status="REPAIRING", analysis_rerun=False))
    frames = {name: [] for name in TABLE_NAMES}
    audit_rows, baselines, controls = [], [], []
    try:
        with (audit / "source_checkpoint_receipts.jsonl").open("w", encoding="utf-8") as receipts, \
                (audit / "bridge_logs.txt").open("w", encoding="utf-8") as logs:
            for count, unit in enumerate(units, 1):
                print(f"READ {count}/{len(units)}: {unit['unit_id']}", flush=True)
                relative = Path("units") / unit["unit_id"]
                complete = checkpoint(source, relative, signature, UNIT_FILES, UNIT_FILES)
                require(read_json(source / relative / "unit.json") == unit, "Unit metadata differs from source protocol")
                receipts.write(json.dumps(dict(path=relative.as_posix(), receipt=complete)) + "\n")
                pieces = {name: [] for name in TABLE_NAMES if name != "per_gene_summary"}
                chunks = 0
                for start in range(0, repetitions, batch_size):
                    chunk = relative / f"chunk_{start:05d}"
                    verified = {name + ".csv" for name in pieces} | {"fit_index.csv", "baseline_checks.json", "r_controls.json", "bridge.log"}
                    complete_chunk = checkpoint(source, chunk, signature, CHUNK_FILES, verified)
                    receipts.write(json.dumps(dict(path=chunk.as_posix(), receipt=complete_chunk)) + "\n")
                    rep_ids = list(range(start, min(start + batch_size, repetitions)))
                    for name in pieces:
                        frame = read_result_csv(source / chunk / (name + ".csv"))
                        validate_unit_table(frame, unit, rep_ids, name)
                        pieces[name].append(frame)
                    attrs = dict(unit_id=unit["unit_id"], chunk_start=start)
                    baselines.extend(dict(attrs, **row) for row in read_json(source / chunk / "baseline_checks.json"))
                    controls.append(dict(attrs, **read_json(source / chunk / "r_controls.json")))
                    log = source / chunk / "bridge.log"
                    if log.stat().st_size:
                        logs.write("\n### " + chunk.as_posix() + "\n")
                        with log.open(encoding="utf-8", errors="replace") as handle:
                            for line in handle:
                                logs.write(line)
                    chunks += 1
                merged = {name: pd.concat(blocks, ignore_index=True) for name, blocks in pieces.items()}
                # Reuse existing per-gene numbers only after explicit null-plan/truth validation.
                genes, relabeled = recover_gene_null_label(read_result_csv(source / relative / "per_gene_summary.csv"), unit, repetitions)
                merged["per_gene_summary"] = genes
                for name, frame in merged.items():
                    validate_unit_table(frame, unit, list(range(repetitions)), name)
                    for key in ("unit_id", "reference_source", "target_source", "target_donor"):
                        frame[key] = unit[key]
                    frames[name].append(frame)
                audit_rows.append(dict(unit_id=unit["unit_id"], batches=chunks,
                    replicate_rows=len(merged["replicates"]), paired_rows=len(merged["paired_replicates"]),
                    null_replicate_rows=int(merged["replicates"].condition.eq("null").sum()),
                    recovered_gene_labels=relabeled, unit_and_consumed_chunk_hashes="verified",
                    fit_array_hashes="not_recomputed;receipt_and_existence_checked"))
        frames = {name: pd.concat(blocks, ignore_index=True) for name, blocks in frames.items()}
        print("SUMMARIZE: all planned conditions, including null; equal target-donor weighting", flush=True)
        summary, donor, status = aggregate_frames(frames, units, repetitions, output / "tables", include_replicates)
        if plots:
            make_plots(output / "tables", summary, donor)
        write_csv(audit / "unit_repair_audit.csv", pd.DataFrame(audit_rows))
        write_csv(audit / "baseline_checks.csv", pd.DataFrame(baselines))
        write_csv(audit / "r_controls.csv", pd.DataFrame(controls))
        write_summary_notes(output / "SUMMARY_METHODS.md")
        for name in ("protocol.json", "OWNER.json", "input_paths.json", "R_versions_current.json", "R_selftest.json", "R_selftest.log", "R_check.log"):
            if (source / name).is_file():
                shutil.copyfile(checked_path(source, Path(name)), output / name)
        shutil.copyfile(source / "RUN_STATUS.json", output / "SOURCE_RUN_STATUS.json")
        refroot = source / "reference_inputs"
        if refroot.is_dir():
            for path in refroot.rglob("*"):
                if path.is_file() and path.name in ("axes.json", "injection_sets.csv", "mead_genes.csv", "normalization_genes.csv"):
                    relative = path.relative_to(source)
                    destination = output / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(checked_path(source, relative), destination)
        write_json(output / "SUMMARY_PROVENANCE.json", dict(
            summary_version=SUMMARY_VERSION, source_directory=str(source), source_signature=signature,
            source_protocol_sha256=sha(source / "protocol.json"), analysis_rerun=False,
            source_version=protocol.get("version"), summary_only=True,
            summary_code_hashes={name: sha(HERE / name) for name in ("repair_revision.py", "revision_summary.py")},
            repetitions=repetitions, units=len(units), batches=sum(row["batches"] for row in audit_rows),
            recovered_gene_labels=sum(row["recovered_gene_labels"] for row in audit_rows),
            null_replicate_rows=sum(row["null_replicate_rows"] for row in audit_rows),
            source_protocol_preserved="protocol.json describes ORIGINAL computation, not a rerun with new code",
            weighting="equal target donors after within-donor available repetition means",
            mcse="paired-repetition ratio influence; fixed-bank conditional Monte Carlo only",
            fit_array_hashes_recomputed=False, source_files_modified=False,
            created_at_utc=datetime.now(timezone.utc).isoformat()))
        require(read_json(source / "RUN_STATUS.json") == original_status and read_json(source / "protocol.json") == protocol,
                "Source changed while repair was running")
        final_status = dict(original_status, summary_version=SUMMARY_VERSION,
            summary_only_repair=True, all_planned_conditions_verified=True,
            all_unavailable_method_conditions=int(status.status.eq("ALL_UNAVAILABLE").sum()))
        write_json(output / "RUN_STATUS.json", final_status)
        (output / "README_EXPORT.md").write_text(
            "# Corrected compact Paper 3 results\n\n"
            "SUMMARY_PROVENANCE.json records a summary-only repair, NOT a new simulation.\n"
            "The copied protocol.json and SOURCE_RUN_STATUS.json describe the original computation.\n"
            "Literal null labels and all planned condition/method/repetition sets have been checked.\n"
            "Gene numbers were retained; only verified blank null labels were restored.\n"
            "Source and paired means/MC intervals were recalculated with the documented equal-donor estimand.\n"
            "See SUMMARY_METHODS.md for missingness and MC formulas; no donor-population CI is implied.\n"
            "Inspect method_status.csv and failure_reasons.csv before interpreting comparisons.\n"
            "The existing source summaries and checkpoints have not been overwritten.\n"
            "No R calls, new sampling, dependency installation, or data downloads occur.\n"
            "Large aggregate replicate tables are omitted unless --include-replicates was supplied.\n"
            "Full fit arrays, input counts, R libraries, and caches remain on the server.\n"
            "This is a review package, not a complete public reproducibility deposit.\n"
            + ("\nSOFTWARE SMOKE DATA ONLY: not publication evidence.\n" if smoke else ""), encoding="utf-8")
        write_json(output / "EXPORT_STATUS.json", dict(status="EXPORTED_SOFTWARE_SMOKE_ONLY" if smoke else "EXPORTED_COMPLETE",
            source_status=original_status["status"], summary_version=SUMMARY_VERSION, summary_only_repair=True,
            units=len(units), batches=sum(row["batches"] for row in audit_rows), analysis_rerun=False,
            include_replicates=include_replicates, signature=signature))
        files = [dict(path=path.relative_to(output).as_posix(), bytes=path.stat().st_size, sha256=sha(path))
                 for path in sorted(output.rglob("*")) if path.is_file()]
        write_json(output / "MANIFEST.json", dict(files=files, summary_only_repair=True))
        print("REPAIR FINISHED:", output, flush=True)
        print("Size: %.2f MiB; transfer this folder, no ZIP." % (sum(item["bytes"] for item in files) / 1024**2), flush=True)
        return output
    except BaseException as exc:
        write_json(output / "EXPORT_STATUS.json", dict(status="INCOMPLETE", error=str(exc), analysis_rerun=False))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("/home/research/git/repository/paper3_revision_results_v2"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--allow-smoke", action="store_true", help="Testing only")
    parser.add_argument("--include-replicates", action="store_true", help="Larger output; normally unnecessary")
    parser.add_argument("--no-plots", action="store_true", help="Skip plot creation")
    args = parser.parse_args(argv)
    output = args.output or args.source.parent / ("paper3_final_results_fixed_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f"))
    try:
        repair(args.source, output, args.allow_smoke, args.include_replicates, not args.no_plots)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print("REPAIR STOPPED:", exc, file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
