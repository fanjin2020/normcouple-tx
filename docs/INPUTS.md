# Scientific inputs and environment

The original analyses and reference reconstruction use eleven exact objects:

```text
paper3_inputs/
  gene_identity_map.csv
  LIBD/counts.mtx
  LIBD/cells.csv
  LIBD/genes.csv
  SEAAD/counts.mtx
  SEAAD/cells.csv
  SEAAD/genes.csv
  external/GSE174367_snRNA-seq_cell_meta.csv.gz
  external/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5
  bulk/GSE125583_RAW.tar
  bulk/GSE125583_family.soft.gz
```

`publication/input_definitions/computational_specification.json` pins prepared banks and external single-nucleus objects; `public_input_provenance.json` pins the bulk objects. Use `check-inputs` before fitting. An existing separate bulk cache can be supplied with `run --bulk-cache PATH`. The workflow creates compatibility receipts for verified bulk files in its own work directory; it does not require old download receipts or edit the original cache.

The new clinical stage additionally needs `GSE53697_RNAseq_AD.txt.gz` and
`GSE53697_family.soft.gz` in `--external-bulk-cache PATH`. Their URLs and hashes are
in `external_bulk_sources.json`. Only missing objects are downloaded when
`--allow-downloads` is explicit. The curated batch crosswalk is shipped in
`GSE53697_batch.csv` with its primary supplementary-source identity.

## Sources

| Resource | Public source | Exact rerun requirement |
| --- | --- | --- |
| LIBD | https://doi.org/10.1186/s13059-025-03552-3 | Prepared count bank, selected cell/gene axes and identity map |
| SEA-AD | https://brain-map.org/consortia/sea-ad/our-data | Prepared count bank with the study's selected nuclei and gene axes |
| GSE174367 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367 | Deposited metadata and filtered H5 matrix |
| GSE125583 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583 | Deposited raw-count tar and SOFT metadata |
| GSE53697 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53697 | Deposited unnormalized read-abundance table and SOFT metadata; fractional raw values retained |

`rebuild-inputs` reconstructs the exact LIBD/SEA-AD matrices and identity map from the files in `publication/input_definitions/raw_reference_sources.json`. All seven reconstructed files were checked against the manuscript's SHA256 identities. The recipe uses the LIBD author's HDF5-backed SingleCellExperiment, author broad labels and donor field, preserves distinct Ensembl rows, and namespaces duplicate barcodes by donor. SEA-AD uses the ten fixed DFC donor objects dated 2026-06-22 and their raw `layers/UMIs`, with the recorded subclass mapping and per-donor sampling seed. It does not invert normalized expression.

```bash
python reproduce.py rebuild-inputs --source-dir /data/public_raw \
  --output ../paper3_inputs --r-library /data/R_library --allow-downloads
python reproduce.py prepare-inputs --output ../paper3_inputs --allow-downloads
```

Place existing raw files, with their recorded basenames, in `--source-dir` to avoid downloads. Both commands are offline unless `--allow-downloads` is present. LIBD reconstruction needs HDF5Array, SingleCellExperiment and Matrix. Allow approximately 35 GB of free working space for cached raw files, extraction and reconstructed banks. The reconstruction test reused previously downloaded author files; it did not redownload their complete contents.

```bash
python reproduce.py prepare-inputs --output ../paper3_inputs \
  --cache-root /data/prepared_inputs --cache-root /data/GSE125583 --allow-downloads
```

The importer checks known input/cache layouts. It reuses exact hashes and downloads missing objects with pinned URLs only when allowed. Interrupted downloads retain resumable partial files. A final filename is committed only after SHA256 verification. Changed remote data are rejected rather than repinned.

An author-supplied JSON file can assign real download URLs to missing objects:

```json
{"files": {"LIBD/counts.mtx": {"sha256": "the exact digest from computational_specification.json", "url": "the actual deposited HTTPS object URL"}}}
```

This illustrates the structure, not a working address. `--sources FILE` accepts only registered filenames and their existing expected hashes.

## Offline input archive

```bash
python reproduce.py pack-inputs --inputs ../paper3_inputs --output ../paper3_inputs.tar
python reproduce.py prepare-inputs --archive ../paper3_inputs.tar --output ../imported_inputs
```

The archive includes the eleven original objects and has a checksum/source manifest.
Transfer the separately verified GSE53697 cache alongside it for an offline clinical
run. It is a local transfer artifact, not an automatic public deposit. Confirm
redistribution permissions before publishing it. Raw counts are excluded from the
GitHub/source ZIP.

## Software

Python dependencies for computation are pinned in `requirements.txt`: NumPy 2.2.6, pandas 2.2.3, SciPy 1.15.3, h5py 3.13.0 and Matplotlib 3.10.1. Python 3.12 is the tested local interpreter. The study R identity includes R 4.5.3, dtangle 2.0.9, edgeR 4.8.2 and jsonlite 2.0.0. Stage-specific bridges additionally check DESeq2, MEAD and dependencies using preserved protocol specifications.

The runner requires an existing R library passed with `--r-library` and an `Rscript` executable (optionally `--rscript PATH`). It does not install packages. Environment mismatches stop before fitting. Select the matching environment explicitly rather than replacing the user's current environment.

The cross-estimator comparison uses MuSiC 1.0.0 at commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c`, SingleCellExperiment 1.32.0, Matrix 1.7-6 and e1071 1.7.17 for linear SVR. Its reference moments are memoized without changing the official `music_prop` function body. Independent calls to the uncached official API validate the cached implementation. `run --music-library PATH` adds an existing library for these stages; standalone `cross-method` and `clinical` accept repeated `--r-library` options.

Primary clinical analysis uses 288 samples (218 AD, 70 controls), a fixed
16,161-gene family and TMM outcome expression across composition-covariate arms.
GSE53697 has 17 samples (9 AD, 8 controls), a cohort-specific expression family and
batch adjustment. Marker coverage and the unmapped demographic covariates are
documented in `docs/ESTIMATORS.md`. Source rates are fractions. Keep the literal
outcome label `null` as text when reading CSVs.
