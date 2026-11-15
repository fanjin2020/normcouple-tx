# Scientific inputs and environment

The full workflow needs eleven exact objects:

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

## Sources

| Resource | Public source | Exact rerun requirement |
| --- | --- | --- |
| LIBD | https://doi.org/10.1186/s13059-025-03552-3 | Prepared count bank, selected cell/gene axes and identity map |
| SEA-AD | https://brain-map.org/consortia/sea-ad/our-data | Prepared count bank with the study's selected nuclei and gene axes |
| GSE174367 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367 | Deposited metadata and filtered H5 matrix |
| GSE125583 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583 | Deposited raw-count tar and SOFT metadata |

The exact LIBD/SEA-AD matrices and identity map have no assigned download URLs in this release. Their landing pages are not substitutes. A verified cache closes the input requirement for its owner; public reproducibility additionally requires deposition of these exact objects with appropriate permissions, or a tested preparation recipe that reproduces their identities.

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

The archive includes eleven exact objects and has a checksum/source manifest. It is a local transfer artifact, not an automatic public deposit. Confirm redistribution permissions before publishing it. Raw counts are excluded from the GitHub/source ZIP.

## Software

Python dependencies for computation are pinned in `requirements.txt`: NumPy 2.2.6, pandas 2.2.3, SciPy 1.15.3, h5py 3.13.0 and Matplotlib 3.10.1. Python 3.12 is the tested local interpreter. The study R identity includes R 4.5.3, dtangle 2.0.9, edgeR 4.8.2 and jsonlite 2.0.0. Stage-specific bridges additionally check DESeq2, MEAD and dependencies using preserved protocol specifications.

The runner requires an existing R library passed with `--r-library` and an `Rscript` executable (optionally `--rscript PATH`). It does not install packages. Environment mismatches stop before fitting. Select the matching environment explicitly rather than replacing the user's current environment.

Clinical analysis uses 288 samples (218 AD, 70 controls), a fixed 16,161-gene family and TMM outcome expression across composition-covariate arms. Source rates are fractions. Keep the literal outcome label `null` as text when reading CSVs.
