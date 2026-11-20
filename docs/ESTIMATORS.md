# Cross-estimator and clinical workflow

These analyses use the same `reproduce.py` entry point as the complete framework.

```bash
python reproduce.py run --stage estimators \
  --inputs /data/paper3_inputs --bulk-cache /data/GSE125583 \
  --external-bulk-cache /data/GSE53697 --allow-downloads \
  --work /results/paper3_estimators --workers 4 \
  --rscript /path/to/Rscript --r-library /data/R_library \
  --music-library /data/music_library
```

The `estimators` scope does not recompute the original seven experiments. Use the
same command with `--resume` after interruption. It verifies inputs and sealed
batches before reuse. No package installation or Git operation is performed.

## Simulation

The four estimators are dtangle with all-gene TMM, dtangle with marker-restricted
normalization, official MuSiC with fixed markers, and linear nu-support-vector
regression. Outcomes, exposure, samples, reference donors and the 100-gene BH
family are fixed within each baseline/probe comparison. Unavailable tests remain
in the family with p=1.

There are 42 reference-target units, two effects and 200 repetitions per effect:
16,800 datasets. The repetition count is fixed before fitting. Conditional
Monte Carlo intervals use within-dataset contrasts and equal donor weights.
They are t-based Monte Carlo approximations, not finite-sample bounds on FDR.
The FDR-change half-width target of 0.02 is reported; failing that target is not
relabelled as successful precision or used for optional stopping. These intervals
are conditional on the observed donors, not confidence intervals for a donor
population.

The SVR comparator uses e1071 1.7.17 / LIBSVM, a linear kernel, cost=1 and nu in
{0.25,0.5,0.75}. Nu minimizes reconstruction error without disease labels.
Negative coefficients are clipped and positive coefficients closed to sum one,
then converted from marker-RNA contribution to fractions using reference yields.
The optimization is support-vector regression, not NNLS. This transparent
comparator is not the official CIBERSORT implementation.

## Clinical analysis

Each cohort has a fixed baseline all-gene TMM log2(CPM+1) outcome matrix and a
fixed expressed non-marker test family. The implementation computes
log2(CPM + 1), not log2(CPM) + 1. It checks identical outcomes across both
reference banks. GSE125583 retains 288 samples and its 16,161-gene family.

The three MuSiC policies use fixed markers, all compatible genes, or genes outside
the test family. Each policy is checked against uncached official `music_prop`.
The other methods are the two dtangle specifications and linear SVR. A tested-count
probe and a class-stratified marker-count probe multiply selected AD sample counts
by exp(0.1), while the response matrix remains unchanged.

GSE53697 supplies 17 dorsolateral prefrontal cortex samples from the Mount Sinai
Brain Bank (9 AD, 8 controls). GSE125583 originated from Banner Sun Health, as
documented in its primary publication. The deposited GSE53697 `_raw` columns
contain fractional read abundances; the workflow preserves these values and does
not substitute the `_rpkm` columns. Different accessions alone are not the basis
of the cohort-source distinction.

The primary cohort adjusts for age and sex. The external cohort adjusts for RNA-seq
batch, transcribed from Supplementary file 1A (DOI 10.7554/eLife.10421.035) using its
explicit RNA-seq sample numbers. That workbook also reports age and postmortem
interval under separate within-group patient numbers, but does not provide an
explicit patient-to-library crosswalk; these fields are not inferred or joined.
Sex is not reported. This limits disease-association interpretation of the external
cohort; the diagnostic comparison keeps this same conditioning set throughout.

Marker-based specifications require at least ten mapped genes in every cell class.
GSE53697 retains only four astrocyte and eight excitatory markers from the LIBD
panel, so its four marker-based LIBD specifications are unavailable. They retain
explicit availability records, with missing discovery counts rather than zeros.
The LIBD broad/non-test MuSiC policies and all six SEA-AD specifications remain
eligible. This coverage check uses gene identities, not disease associations.
After fitting, the external SEA-AD non-test MuSiC specification has leverage one,
so its HC3 inference is unavailable for all three probes. These entries retain
missing discovery and crossing counts rather than zero findings.

Clinical outputs describe coefficient sensitivity and BH threshold crossings.
There is no clinical effect truth, so these quantities are not called clinical
FDR, false positives, or validated disease biomarkers. The external cohort is
small; its design rank, residual degrees of freedom and maximum leverage are
reported explicitly.

For report-layer reanalysis, `clinical --estimates-from COMPLETED_CLINICAL_RUN`
can reuse sealed fractions and baseline outcomes. It verifies the scientific
protocol, estimator/input code, raw matrices, marker reference and analysis plan
before reuse. HC3 and BH are recomputed. An altered estimator, input or plan is
rejected; estimator reuse is recorded in the new validation receipt.

## Compact output

`cross_method/final/summary.csv` and `paired_precision.csv` contain simulation
summaries. `clinical_estimators/final/method_summary.csv` and
`paired_sensitivity.csv` contain clinical summaries. Completion receipts and
official-API parity results are recorded separately. Detailed per-replicate and
per-gene arrays remain in the work directory, not the publication export.

The complete `run --stage all` graph maps these outputs to shared publication
CSVs (Tables S16-S19) and Figures S8-S10. `export --published` replays only the
archived manuscript snapshot; it is never used to fill missing new computations.

Sources: [GSE53697](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53697),
[Scheckel et al.](https://doi.org/10.7554/eLife.10421),
[Srinivasan et al.](https://doi.org/10.1016/j.celrep.2020.107843),
[Newman et al.](https://doi.org/10.1038/nmeth.3337).
