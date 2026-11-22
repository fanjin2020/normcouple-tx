# CBC reproduction

The CBC analysis requires Python 3.12.13 and the core numerical libraries
in [requirements.txt](../requirements.txt). NumPy, pandas and SciPy versions
are taken from the original server simulation and transfer run records.
h5py and Matplotlib are pinned to the server's data I/O and figure export
versions; the original simulation record did not specify them.
The archive-only lock is [requirements-archived.txt](../requirements-archived.txt).

Run from the unpacked NormCouple-Tx source directory:

```bash
python reproduce.py check-enhancement-env
bash run_cbc.sh --pseudobulk ./../pseudobulk.npz \
  --imaging-cache ./../cbc_cache/libd --bulk-cache ./../cbc_cache/geo \
  --output ./../cbc_results --workers 36 --allow-downloads
```

Use your actual paths. Set `PAPER3_PYTHON` to choose the Python interpreter.
The equivalent Windows command starts with `python reproduce.py cbc-upgrade`.
The run is Python-only; neither command installs packages or updates Git.
The version check runs before the input or output directories are touched.
No original simulations are rerun. Reuse a directory with the same command
and `--resume`; each completed batch is verified before it is reused.

## Inputs

- The GSE174367 **reference and target source** is the existing donor-by-class
  `prepared/pseudobulk.npz` from the fixed-panel empirical analysis. Its
  SHA-256 is
  `5ea4c5ea4214fc64d43ee6f6fca3ff2bd38c203f536ba307dc8568ababdad16e`.
  It contains 18 donor by seven class profiles constructed from the 61,472
  author-annotated nuclei. The underlying inputs are the GEO matrix and cell
  annotation pinned in `publication/input_definitions`; the aggregation code is
  `kernels/markers/anchor_core.py`. Generic pseudo-bulk files are not substituted.
- The **observed dose** uses the previously pinned Huuki-Myers et al. LIBD bulk
  counts, `sample_info.csv`, imaging mapping and author `library_type.csv`.
  Existing correctly hashed files are reused. The upstream library labels are
  at the same commit as the imaging audit. Use `--library-types` if the label
  CSV resides separately.
- The **clinical mapping and support** audit uses the pinned GSE125583 raw TAR
  and the original 288 samples, original 16,161 outcome genes and frozen
  one-to-one mapping. It does not fit a clinical effect with the proposed
  anchor set.

Adding `--allow-downloads` fills only the missing, pinned public LIBD and GEO
files. It does not download or replace the donor-profile NPZ. Supply the exact
NPZ from the previous analysis. `--smoke` is a separately marked acceptance
profile, not evidence. The full run uses six reference donors, 12 held-out
donors and 60 paired draws per donor and composition direction.

## Final results

- `final/data/`: six CSVs for the observed library pairs and dose, held-out
  performance and paired contrasts, clinical count support, and the fixed
  anchor list and weights.
- `final/figures/FigS12.pdf`: observed-dose and held-out evidence. The
  figure is exported whether the candidate passes or fails the decision rule.
- `final/article/`: redraws the unchanged archived numerical inputs for the
  nine main and 11 existing supplementary figures in the CBC figure style.
  Use `--skip-replay` to omit this optional export, or `--formats pdf png`
  to include PNG copies for composition in a word processor.
- `final/summary.json` and `final/receipt.json`: design, completion,
  predeclared comparisons, numerical checks and file hashes. Working batches
  are not publication or upload files.

The TMM denominator follows the edgeR weighted-trimmed-mean definition; the
marker mean-contrast calculation is checked against dtangle on a fixed
numerical fixture. Remaining genes contribute to full library size; only the
fixed computational universe enters TMM factor estimation. The benchmark
knows the simulated RNA fractions and direct effects. The observed
paired-preparation study supplies technical doses, not disease-effect truth.

The absolute FDR check is separate from the paired FDR, coverage and power
contrasts. Poor anchor support is reported as an applicability limit; do not
remove outcome genes post hoc to make the anchor panel applicable.
