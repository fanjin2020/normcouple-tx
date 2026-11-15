# Paper 3 bounded strengthening v1

This is an additional analysis, not a replacement for reviewed manuscript
results. Default: 36 Linux CPU workers, one BLAS/OpenMP thread per worker,
200 repetitions, the same 42 reference-target units and cached count banks.
No new dataset, no GPU dependency, no ZIP output and no automatic raw FASTQ
download. Do not put its output in an old mechanism/revision result directory.

## Server run

Place the single `paper3_strengthening_oneclick.py` in the existing project and run:

```bash
python paper3_strengthening_oneclick.py --setup \
  --project-root /home/research/git/repository/paper1_public_evidence \
  --output /home/research/git/repository/paper3_strengthening_results_v1 \
  --workers 36
```

The launcher reuses a working pinned Python environment. Existing R libraries
are reused, not reinstalled. The grid needs dtangle 2.0.9, edgeR and jsonlite.
The raw bulk comparison additionally uses the already installed dependencies
of the reviewed bridge. Missing dependencies fail with a log, never trigger
an unbounded upgrade. Use `--r-library` and `--rscript` for a different location.

`--pilot-only` runs one batch per reference source first. Remove that flag to
continue. A rerun checks source/protocol/input hashes and each completed batch's
file hashes before skipping it. An incomplete batch is recalculated in full.
Changing only `--workers` does not change the analysis identity. Stop the prior
run and wait for all children to exit before resuming; do not delete run.lock.

Progress, including a simulation-only estimate from observed batch times:

```bash
python paper3_strengthening_oneclick.py --status \
  --output /home/research/git/repository/paper3_strengthening_results_v1 --workers 36
```

The first run has 42 * 40 = 1,680 batches at the default batch size 5. This is
new work, not a promise to reuse an old central-condition result as if it were
the new grid. The original central intervention and baseline RNG draws are
preserved and tested. Speed depends on CPU, fork memory sharing and storage;
256 GB RAM alone does not determine completion time. Start at 36 workers.

## Prespecified experiments

- Fractions 5%, 10%, 20% of each marker class, amplitudes 0.25, 0.5, 1.0.
- The zero-stress control and one independently prespecified random marker set
  at the central 10%/0.5 point. Class counts use ceil with at least one marker;
  actual counts and gene IDs are recorded, not described as exactly 5/10/20%.
- Null and up05 outcomes at every grid point, same mixtures and outcome noise.
- Internal zero-stress up025 control for the documented power-loss setting.
- Oracle, observed-total, TMM, MRN matched-reference/yield, denominator-only.
- Full-six regression rule retained for this comparison. The separate reusable
  actual-span function is verified against the frozen implementation. This
  stage does not rerun MEAD or claim joint MEAD uncertainty propagation.
- FDR (mean replicate FDP), coverage, delivered coverage, power, composition
  RMSE, availability, exact invariance, finite input/coefficient bounds and
  gene-level power decomposition. No outcome-based choice of a winning method.

## Bulk repair and the data boundary

The old pipeline could require cohort_exclusions.csv before generating it.
This version first reuses/downloads the hash-verified public processed RDA and
single-nucleus metadata, reconstructs the cohort and verifies explicit IDs.
It runs a clearly labeled descriptive association on deposited processed
expression, including a fixed-marker overlap report. That is NOT a raw-count
normalization comparison or an estimate of true FDR.

It then searches the specified project/previous bulk run for provenance
manifests, audits public GEO file listings, and accepts either an existing
verified manifest or explicit authorized raw gene counts:

```bash
python paper3_strengthening_oneclick.py --setup --only bulk \
  --output /home/research/git/repository/paper3_strengthening_results_v1 \
  --raw-counts /actual/path/GSE174367_bulk_raw_counts.csv \
  --counts-provenance "Actual resource, version and raw gene-count measurement" \
  --workers 36
```

Optional `--sample-map` requires two unique columns `raw_sample,SampleID`.
`--counts-url` additionally requires `--counts-sha256` and provenance; HTTPS
downloads are bounded and checksum verified. Provenance is a source statement,
not a property inferred from integer values. Gene IDs must be unambiguous
Ensembl IDs. Files with transformed expression are rejected, never inverted.

The deposited GEO bulk RDA contains normExpr.reg/targets, not raw gene counts.
The linked Synapse resource may require authorized access. This code cannot
grant access, certify unlabelled local tables, or manufacture counts. With no
valid raw input it records PARTIAL_DESCRIPTIVE_ONLY_RAW_COUNTS_REQUIRED and
continues the grid. Exit 2 after a finished grid means inspect the bulk status;
it does not mean the completed grid was lost. Supply counts later and use
`--only bulk`: the grid is not rerun.

## Files to collect

Return only `paper3_strengthening_results_v1/final_results/`. It contains
summary/donor/paired tables, the grid figure (PDF/PNG), gene-level variance
components, invariance audit, theory-check summary, protocol and bulk status.
Large per-repetition audits stay in the server output as CSV.gz; this is gzip
compression of a data table, not a ZIP bundle. No raw counts or R libraries are
copied into final_results. The submitted manuscript should be updated with
these results only after completion and scientific review.

## Reusable entry points

`mrn_tools.py` depends only on NumPy. It exports `mrn_inputs`, `yield_convert`
and `actual_span_hc3`. Supply explicit fixed marker order, reference profiles,
reference yields and anchor. `actual_span_hc3` reports estimability and HC3
availability, not a guarantee of unbiased regression. See function docstrings
and `example_reusable.py`. `THEORY.md` gives conditions and derivations.

`python -m unittest discover -s tests -v` checks algebra, reference matching,
API agreement, checkpoint integrity and input validation. `--synthetic` is a
two-repetition software fixture and is prominently labeled NON-research.
No public repository or DOI is claimed to exist.
