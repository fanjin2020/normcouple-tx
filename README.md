# NormCouple-Tx

**Normalization-coupling diagnostics for transcriptomics.** Version 2.1.0.

NormCouple-Tx connects input verification, statistical analysis, compact numerical outputs and manuscript figures for **Normalization coupling and marker restriction in brain transcriptome inference**, prepared for *Neuroinformatics*.

Public repository: [fanjin2020/normcouple-tx](https://github.com/fanjin2020/normcouple-tx).

Authors: Jin Fan, Chang Liu, Ying Cao, Xiankun Zhang, Chen Song. Correspondence: zhxkun@tust.edu.cn; songchen@tust.edu.cn.

## One entry point

`reproduce.py` connects verified inputs, dependency-ordered computation, compact tables and all manuscript figures. Full scientific fitting targets Linux. Publication export and software tests also run on Windows. No command installs packages automatically.

Use an existing Python environment matching `requirements.txt` for statistical fits. R software requirements are documented in [docs/INPUTS.md](docs/INPUTS.md). `openpyxl` is optional and used only with `--workbook`.

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
  --r-library /data/R_library --workers 36
```

The default `all` scope executes normalization, mechanism, marker-grid, random-panel, coupling, common-support, marker-association and cross-estimator analyses. Prerequisites are produced in the same run. The final export reads those newly computed outputs; it never falls back to archived article values. The cross-estimator stage requires the pinned MuSiC installation; an additional existing library can be supplied with `--music-library`.

```bash
python reproduce.py status --work ../paper3_run
python reproduce.py run --inputs ../paper3_inputs --work ../paper3_run \
  --r-library /data/R_library --workers 36 --resume
```

`run_server.sh` forwards the same arguments to `reproduce.py run`. Set `PAPER3_PYTHON` to select an existing interpreter. Each numerical worker uses one BLAS thread. `--dry-run` prints the ordered commands without fitting. Protocols, code identities and completion receipts protect resume; incomplete batches are delegated to the checkpoint-aware scientific kernels.

### Redraw the published results without refitting

```bash
python reproduce.py export --published --output ../paper3_published
```

This explicit replay mode uses bundled numerical data and records `published_data_replay` in its manifest. It is different from a new statistical run. To export a completed run to another location:

```bash
python reproduce.py export --work ../paper3_run --output ../paper3_export
```

## Final outputs

The default final directory contains:

- `data/TableS1.csv` through `TableS16.csv`: sixteen supplementary tables.
- Eleven additional CSV datasets (one gzip-compressed) needed for figure panels and numerical checks: **27 numerical files** in total, with no duplicate plot copies of table data.
- `data/catalog.json`: explicit table/figure views of the shared CSV files.
- `figures/Fig1.pdf` through `Fig9.pdf` and `FigS1.pdf` through `FigS8.pdf`: **nine main and eight supplementary figures**.
- `figures/figure_sources.json`: actual CSV files, views, dimensions and caption for each figure.
- `tables.json`: table titles, notes and columns; numerical rows remain in the CSVs.
- `export_manifest.json`: data origin and output checksums.

PDF is the default to avoid redundant image formats. Use `--formats pdf tiff` (or `png`, `eps`) when needed. `--workbook` adds `Numerical_Tables.xlsx` and requires `requirements-export.txt`. Working matrices, logs and checkpoint batches stay outside the final directory and are not packaged for GitHub. [docs/ANALYSES.md](docs/ANALYSES.md) lists the computation/output mapping.

## Source organization

For user-supplied data, the [`sensitivity` workflow](docs/fixed_response_audit.md) compares composition-adjusted coefficients with fixed outcomes and a fixed gene family. It reports coefficient weights, standard errors, estimability and BH status changes. `cross-method` runs the manuscript's paired dtangle/MuSiC comparison separately when only that experiment is needed.

- `pipeline/`: input preparation, verification, execution, aggregation, export and source packaging.
- `kernels/`: one canonical copy of each scientific source file. `layout.json` maps files to the import layout created under `WORK/.runtime`; hashes preserve algorithms and frozen protocol identities. These are internal modules, not separate user workflows.
- `publication/`: plotting functions, shared numerical sources, captions, table schema and fixed study input definitions.
- `tests/`: numerical, input, checkpoint, download and export checks.
- `docs/`: input requirements, analysis/output mapping, figure style and validation scope.

Build the upload archive with `python reproduce.py package --output ../paper3-normalization-coupling_v2.1.0.zip`. It contains registered package files only: no raw matrices, checkpoints, virtual environments or author-side intermediate reports.

## Citation and license

Use `CITATION.cff` and identify the software version and commit. A software DOI has not been assigned. The supplied GPL-3.0-or-later license is retained; third-party licenses and data-access conditions remain applicable. See `LICENSE` and `THIRD_PARTY.md`.
