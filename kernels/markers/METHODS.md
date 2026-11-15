# Fixed-marker empirical anchor: locked analysis definitions

## Question and scope

Measure the direction and size of disease-associated changes in the existing
LIBD and SEA-AD reference marker panels, within deposited GSE174367 broad cell
types. This is an observational empirical anchor for a sensitivity analysis,
not an independent validation of MRN accuracy or a substitute for raw-bulk
normalization comparisons. The same resource was used in method development;
it is not described as an unused external cohort.

## Inputs and units

Use the public author-filtered Cell Ranger aggregate integer UMI matrix and
author nucleus metadata, with release SHA256 checks. Exact Barcode matches are
required. Sum all annotated nuclei within SampleID and broad cell type, keeping
whole-matrix gene totals. SampleID is the deposited donor-level analysis unit.
Do not restrict a donor to having all seven classes. Do not relabel these
single-nucleus aggregates as genuine bulk tissue libraries or raw FASTQs.

The published annotation is retained without reclustering. PER.END is mapped
to EndoMural; effects in this combined class can include subtype shifts.
Unannotated H5 columns are excluded and counted. There is no cell-level test.

Primary cell-type designs require at least 20 annotated nuclei per donor,
at least four AD and four control donors, a full-rank declared design and at
least four residual degrees of freedom. A failed design is reported, not
rescued by a data-dependent lower threshold. Missing cell-type samples do not
remove that donor from other cell types. Metadata has 18 donors (11 AD, 7
controls); actual per-type fitted counts can be lower.

Primary covariates are age, sex and batch. Extended sensitivity adds PMI and
RIN and uses complete cases. Age is modeled numerically as deposited, including
values recorded as 90; no unknown exact ages are invented. Diagnosis, metadata
agreement, group sizes, rank, exposure information and leverage are audited.
Tangle/plaque stages are not additional adjustment variables because they
define/track the disease contrast. Residual confounding remains possible.

## Gene-level anchor

The two balanced split-0 panels (300 genes per reference) are frozen before
looking at disease outcomes. Native Ensembl version suffixes are removed;
duplicated IDs stop analysis. Absent genes are not inserted as observed zeros.
The complete planned panel, annotation presence and expression-filter status
are retained in the final gene table.

Use edgeR filterByExpr with the declared design and its default expression
thresholds. Retain original whole-library UMI totals after filtering. Primary
TMM factors use all expressed genes. Sensitivity TMM factors use expressed
genes outside the union of BOTH frozen panels, with the same full-library
sizes and expression filter. Then estimate dispersion and fit robust
quasi-likelihood negative-binomial GLMs with edgeR. The installed R/package
versions are recorded and locked on first run; a version change requires a
new output directory. The locally tested combination is edgeR 4.4.2, limma
3.62.2 and statmod 1.5.1 on R 4.4.2.

Export edgeR log2FC, QL p values and glmTreat threshold p values for
log2(1.25), log2(1.5), and log2(2). These correspond to multiplicative increases
of 25%, 50%, and 100%, not log2FC values of 0.25, 0.5, and 1.0. For downregulated
genes, the same log magnitudes are reciprocal ratios (0.8, 2/3, 0.5), not
symmetric percentage decreases. Counts of point estimates beyond a threshold
are descriptive and kept separate from threshold-test support. Estimates of
relative expression depend on the normalization assumptions and do not measure
absolute transcript counts per cell.

Primary scope: each marker is tested within its frozen owner cell type.
Secondary scope: the entire marker panel within each broad cell type.
For each model/normalizer/scope, BH is applied jointly to both panels and all
seven cell types, retaining untestable planned entries as p=1. Thresholds are
separate reported families, not pooled as a single 5% guarantee across all
analyses. Untestable entries retain missing raw p values. Report planned,
annotation-present and tested denominators separately. Low power or absence
of detected effects does not establish marker invariance. No gene-binomial
confidence interval is used because genes are correlated.

## Library endpoints and sensitivity

For both panel scopes, export donor marker UMI, total UMI, annotated nucleus
count, marker fraction, marker UMI/nucleus and total UMI/nucleus. Fit declared
donor designs to log2((marker UMI+0.5)/(total UMI+1)), log2((marker UMI+0.5)/nuclei)
and log2((total UMI+0.5)/nuclei). Use HC3 standard errors and residual-df t
intervals. The +0.5 is only a declared stabilization for these descriptive
endpoints; gene-count GLMs do not add it. Export leave-one-donor-out fits for
each library endpoint and report nonestimable fits explicitly. These sensitivity
fits assess influential donors; they are not independent validation cohorts.
Library tests use BH within each model/endpoint over references, scopes and
available types. Absolute captured UMI/nucleus is capture-sensitive, not an
independent measure of biological RNA yield.

## Interpretation limits

Real disease changes can be bidirectional, correlated across genes, and include
within-class subtype composition shifts. They are not the same intervention as
adding Poisson counts to selected markers in an otherwise controlled mixture.
Empirical overlap with a simulated dose supports examining that sensitivity
range, not a claim that its FDR benefit will occur in a real cohort. This script
does not select the simulation grid after seeing associations, construct
processed bulk counts, or calculate a true real-data FDR.

## Sources

- GSE174367: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367
- Morabito et al.: https://doi.org/10.1038/s41588-021-00894-z
- edgeR methods: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
- Access checked: 2026-10-14.
