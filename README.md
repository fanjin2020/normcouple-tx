# NormCouple-Tx

**Normalization-coupling diagnostics for transcriptomics.** Version 2.3.0.

For the CBC figures and the finite-dose, held-out anchor experiment, use
[`bash run_cbc.sh`](docs/CBC.md). It reuses verified inputs and does not rerun
the original benchmarks. The archived 2.2.2 study remains available through
[its Zenodo DOI](https://doi.org/10.5281/zenodo.24394454).

Run the six evidence checks from a completed estimator analysis with `bash run_evidence.sh`.
See [evidence protocol and inputs](docs/EVIDENCE.md).

NormCouple-Tx connects fixed-response analyses, compact numerical outputs and publication figures. The CBC workflow adds paired library-preparation doses, a donor-isolated comparison of fixed anchors and the original rules, and a clinical feasibility check.

Public repository: [fanjin2020/normcouple-tx](https://github.com/fanjin2020/normcouple-tx).

Authors: Jin Fan, Chang Liu, Ying Cao, Xiankun Zhang, Chen Song. Jin Fan and Chang Liu share first authorship. Correspondence: zhxkun@tust.edu.cn; songchen@tust.edu.cn.

## One entry point

`reproduce.py` connects verified inputs, dependency-ordered computation, compact tables and all manuscript figures. Full scientific fitting targets Linux. Publication export and software tests also run on Windows. No command installs packages automatically.

The archived manuscript workflows (`run`, `cross-method`, `clinical`, and
`evidence`) retain the environment in
[requirements-archived.txt](requirements-archived.txt). The CBC enhancement
requires Python 3.12.13 and the original server's core numerical versions in
[requirements.txt](requirements.txt). Run
`python reproduce.py check-enhancement-env` to verify that environment without
loading data. R software requirements for the archived workflows are in
[docs/INPUTS.md](docs/INPUTS.md). `openpyxl` is optional and used only with `--workbook`.

```bash
python reproduce.py verify
python -m unittest discover -s tests -v
```

### Reuse downloaded inputs

Paths below are examples; set them to your actual data, run and R-library directories. Multiple `--cache-root` arguments are supported. Existing caches are read, not edited.

```bash
python reproduce.py prepare-inputs --output ../paper3_inputs \
  --cache-root /data/paper3_prepared_inputs --cache-root /data/public_bulk_cache
python reproduce.py check-inputs --inputs ../paper3_inputs
```

Downloads require `--allow-downloads`. Exact LIBD/SEA-AD banks and the gene-identity map can also be rebuilt from pinned author-distributed files:

```bash
python reproduce.py rebuild-inputs --source-dir /data/public_raw \
  --output ../paper3_inputs --r-library /data/R_library --allow-downloads
python reproduce.py prepare-inputs --output ../paper3_inputs --allow-downloads
```

The reconstruction was tested against all seven prepared-bank hashes. The source manifest fixes the LIBD archive, ten SEA-AD donor files, subclass mapping and sampling seed. The second command obtains missing GEO inputs. See [docs/INPUTS.md](docs/INPUTS.md).

### Compute and export

```bash
python reproduce.py run --inputs ../paper3_inputs --work ../paper3_run \
  --r-library /data/R_library --music-library /data/music_library \
  --external-bulk-cache /data/external_bulk --allow-downloads --workers 36
```

The default `all` scope executes nine stages: normalization, mechanism, marker-grid, random-panel, coupling, common-support, marker-association, cross-estimator simulation and clinical-estimator analysis. Prerequisites are produced in the same run. The final export reads those newly computed outputs; it never falls back to archived article values. The cross-estimator stages require pinned MuSiC and e1071 installations; an additional existing library can be supplied with `--music-library`. Existing verified external data are reused without downloading.

Use `--stage estimators` to run only the new simulation and clinical analyses, with no rerun of the original seven experiments. The simulation fixes 200 repetitions per unit and effect before fitting. The clinical stage analyzes GSE125583 and GSE53697 with six estimator/gene-policy combinations and two fixed-response probes. See [docs/ESTIMATORS.md](docs/ESTIMATORS.md).

`--profile acceptance` runs the complete graph with real inputs and two simulation repetitions. Its figures are marked **software acceptance, not research results**; ordinary publication export rejects these runs. This profile tests input-to-figure connectivity, not statistical calibration.

```bash
python reproduce.py status --work ../paper3_run
python reproduce.py run --inputs ../paper3_inputs --work ../paper3_run \
  --r-library /data/R_library --workers 36 --resume
```

`run_server.sh` forwards the same arguments to `reproduce.py run`. Set `PAPER3_PYTHON` to select an existing interpreter. Each numerical worker uses one BLAS thread. `--dry-run` prints the ordered commands without fitting. Protocols, code identities and completion receipts protect resume; incomplete batches are delegated to the checkpoint-aware scientific kernels.

### Redraw the published results without refitting

```bash
python reproduce.py export --published --formats pdf png --workbook --output ../paper3_submission_v222
```

This command uses the bundled numerical data to draw the nine main figures and eleven supplementary figures as vector PDFs and 300 dpi PNGs, and writes the Guide and Tables S1-S21 to an Excel workbook. The PDFs are for standalone figures; the PNGs match the illustrations embedded in the submission documents when the rendering environment and fonts match. To generate the workbook and CSVs without drawing any figures, use `python reproduce.py tables --published --workbook --output ../paper3_tables`. Both commands use the same figure-to-data index.

This replay is marked `published_data_replay` in its manifest. To export a completed statistical run instead:

```bash
python reproduce.py export --work ../paper3_run --output ../paper3_export
```

## Final outputs

Article-data replay uses the tables and figure data bundled with this version:

- `data/TableS1.csv` through `TableS21.csv`: 21 supplementary tables.
- Eleven additional CSV datasets (one gzip-compressed) for figures and numerical checks: **32 numerical files** in total.
- `data/catalog.json`: explicit table/figure views of the shared CSV files.
- `figures/Fig1.pdf` through `Fig9.pdf` and `FigS1.pdf` through `FigS11.pdf`: **nine main and eleven supplementary figures**.
- `figures/figure_sources.json`: actual CSV files, views, dimensions and caption for each figure.
- `tables.json`: table titles, notes and columns; numerical rows remain in the CSVs.
- `export_manifest.json`: data origin and output checksums.

A full benchmark export produces S1–S19. Add `--evidence ../paper3_review/evidence` (or the `results` directory of a completed `evidence` run) to rebuild S20 from signed coefficient changes and measured doses and S21 from paired tissue-block imaging measurements. Using this switch with `tables --published` checks both against the article snapshot. The 72 clinical probe assignments are in `publication/input_definitions/clinical_probe_markers.csv`. Try `python reproduce.py example --output ../normcouple_example` for a small offline audit.

PDF is the default to avoid redundant image formats. Use `--formats pdf tiff` (or `png`, `eps`) when needed. `--workbook` adds `Numerical_Tables.xlsx` and requires `requirements-export.txt`. Working matrices, logs and checkpoint batches stay outside the final directory and are not packaged for GitHub. [docs/ANALYSES.md](docs/ANALYSES.md) lists the computation/output mapping.

## Source organization

For user-supplied data, the [`sensitivity` workflow](docs/fixed_response_audit.md) compares composition-adjusted coefficients with fixed outcomes and a fixed gene family. It reports coefficient weights, standard errors, estimability and BH status changes. `cross-method` and `clinical` expose the same pipeline modules for focused execution. Linear SVR uses LIBSVM, not NNLS; it is not presented as the official CIBERSORT implementation.

- `pipeline/`: input preparation, verification, execution, aggregation, export and source packaging.
- `kernels/`: one canonical copy of each scientific source file. `layout.json` maps files to the import layout created under `WORK/.runtime`; hashes preserve algorithms and frozen protocol identities. These are internal modules, not separate user workflows.
- `publication/`: plotting functions, shared numerical sources, captions, table schema and fixed study input definitions.
- `tests/`: numerical, input, checkpoint, download and export checks.
- `docs/`: input requirements, analysis/output mapping, figure style and validation scope.

Build the upload archive with `python reproduce.py package --output ../normcouple-tx.zip`. It contains registered package files only: no raw matrices, checkpoints, virtual environments or author-side intermediate reports. A development version is not an archived manuscript release.

## Citation and license

Use `CITATION.cff` and identify the software version and commit. The DOI above identifies version 2.2.2, not the 2.3.0 analysis. The supplied GPL-3.0-or-later license is retained; third-party licenses and data-access conditions remain applicable. See `LICENSE` and `THIRD_PARTY.md`.
