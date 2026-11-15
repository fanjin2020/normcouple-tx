# Validation scope

Local validation was performed on Windows on November 15, 2026, using Python 3.12 and the pinned NumPy 2.2.6, pandas 2.2.3, SciPy 1.15.3, h5py 3.13.0 and Matplotlib 3.10.1 environment.

## Checks performed

- **48 framework/publication tests** passed: MRN transformations, actual-span HC3, gene families, headline numerical checks, source-layout hashes, dependency ordering, missing-stage rejection, fixture rejection, changed-file resume checks, safe paths, input-cache reuse, interrupted imports and downloads, checksum failures, and writer-lock status.
- **17 common-support kernel tests** passed in the materialized runtime layout. These exercise retained budgets, zero-count support, batch resume, multistep clinical probes, derivative qualification and output integrity using small software fixtures.
- All seven scientific entry points accepted their configured option names in help/CLI compatibility checks.
- The new aggregation code was run on existing server-result summaries. All fifteen supplementary tables and all eleven additional datasets matched the archived publication sources, using exact labels and numerical tolerances of relative 1e-10 and absolute 1e-12. This validates the calculation-to-publication transformations without rerunning the underlying fits.
- All nine main and seven supplementary figures were rendered from that newly aggregated data, with text-bound checks. Their source map names physical CSV files and explicit views; table/figure copies are not duplicated.
- The unified archived-data export produced 26 numerical files and 16 figures. Its optional workbook had 15 data sheets plus Contents; 37,842 data cells matched the shared CSV sources.
- Scientific kernels retain their original bytes and internal manifest identities after materialization. Identical physical source files are stored once in the upload package.

Run the public test suite with:

```bash
python -m unittest discover -s tests -v
python reproduce.py verify
python reproduce.py export --published --output ../paper3_publication_check
```

The complete source/output receipt chain is also checked by `verify --work RUN` and by `export --work RUN`. These commands reject missing stages, unsealed outputs, changed files and a different framework identity. A success status alone is not sufficient.

## Scope of evidence

No complete Linux/R-based refit or 36-worker performance run was performed on this local Windows computer. Software fixtures and reconstruction from existing server summaries are not new biological analyses. Full end-to-end verification requires the exact prepared banks, the specified R environment and execution of the seven-stage run on the server.

The exact LIBD/SEA-AD prepared matrices and gene map remain an input-access requirement for new users. Cache import, checksum-pinned URL manifests and offline input packaging are implemented; this release does not claim that missing objects have already been publicly deposited.
