# Paper 3: empirical marker-stress anchor

This independent incremental analysis uses GSE174367 raw snRNA UMI data already
used by Paper 3. It does not rerun or modify the strengthening grid, manuscript,
AI declaration, or any previous output directory. No new dataset is added.

```bash
python paper3_marker_anchor_oneclick.py --setup \
  --project-root /home/research/git/repository/paper1_public_evidence \
  --output /home/research/git/repository/paper3_marker_anchor_results_v1 \
  --workers 36
```

The launcher reuses the pinned Python environment and the existing successful
R library. It reads cached public H5/metadata first; only missing exact files
are downloaded from NCBI. The H5 is approximately 274 MB, not a new FASTQ run.
Use `--offline` to disallow all downloads. `--external-cache`, `--counts`,
`--metadata`, `--rscript` and `--r-library` allow explicit existing paths.
R dependencies (edgeR, limma, statmod, jsonlite) must already be installed;
the script does not modify an R library used by another running analysis.

There are at most 28 independent cell-type/model/normalizer fits, so 36 is a
worker ceiling, not a requirement to create idle processes. BLAS uses one
thread per process. The H5 is aggregated once; it is not reread per fit.
Do not start this alongside a saturated 36-worker simulation without checking
available CPU and RAM. Actual duration depends on cache/network and R runtime;
it is not another 42-unit, hundreds-of-repetitions simulation.

Check without downloading the large count matrix:

```bash
python paper3_marker_anchor_oneclick.py --setup --preflight
```

Status:

```bash
python paper3_marker_anchor_oneclick.py --setup \
  --output /home/research/git/repository/paper3_marker_anchor_results_v1 --status
```

Rerun the same command after interruption. Only checkpointed, checksum-verified
complete fits are reused; incomplete/corrupt fits are recomputed. The output
lock prevents concurrent writers. Changed code, declared protocol or R package
versions require a new output directory, so different analyses cannot be mixed.

## Return only `final_results/`

- `marker_gene_associations.csv.gz`: every planned fixed marker, estimate,
  threshold tests, adjusted p values, missingness and the explicit family.
- `marker_dose_anchor_summary.csv`: separate up/down dose fractions, separate
  point estimates and statistically supported counts, all denominators.
- `donor_marker_libraries.csv`: donor-level captured-marker library sizes.
- `marker_library_associations.csv`: adjusted associations and HC3 intervals.
- `library_leave_one_donor_out.csv`: influence sensitivity, not new replicates.
- `design_diagnostics.csv`, `cohort_exclusions.csv`: low-cell/missing-covariate
  exclusions and unavailable cell-type designs, not silently dropped.
- `normalization_factors.csv`, `sessionInfo.txt`: sample offsets and the exact
  R environment used, plus frozen-panel source provenance.
- Three PDF/SVG/PNG figure sets, editable source data and figure legends.
- Protocol, public-input hashes, frozen markers, R versions and final status.

Raw H5 and intermediate per-fit matrices are intentionally not in this compact
folder. Keep the source and input cache on the server for reruns. A complete
status can still contain declared gaps in rare cell types; the final status
lists unavailable fits. It never means raw-bulk validation was completed.

Read METHODS.md before turning any result into a manuscript claim. In
particular, do not claim that finding disease-associated marker changes proves
MRN superiority, causal bias, or the same perturbation as a simulated spike-in.
