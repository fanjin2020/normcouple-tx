# Revision protocol, fixed before this new run

This is a retrospective hypothesis-driven extension, not a preregistration or an untouched external validation. No improvement is required for the run to be complete.

## Structure

The unchanged stage-2 engine, donor roles, balanced markers, tested genes, cell-class order, count sampling and seeds are retained. Every unit's first null and 0.5 repetition is independently regenerated with that engine and compared elementwise.

Internal split0 comprises two references and six targets each. External GSE174367 uses these two references and15 targets, with annotation-only overlap filtering. These are42 units, not42 independent people. Each mixture has200 sampled nuclei; n=80 mixtures, composition shift1, and200 repetitions by default. Four reference draws come from four donors. Target exposure is simulated, not clinical diagnosis. No new resource or external outcome is used to select parameters.

## Interventions

Response is count divided by200 times the original pooled-reference test-gene mean. Counts stay nonnegative integers. Null/up conditions use Poisson additions and the old random-seed scheme.

Internal conditions are null, effects0.1/0.25/0.5, balanced signed effects, and low/high reference-expression injection sets. External conditions are only null and original0.5. Low/high sets each contain20 of the frozen100 test genes, selected by pooled-reference expression alone. This tests sensitivity inside the selected family, not genome-wide representativeness.

Balanced effects add counts to ten genes in group1 and ten other genes in group0. The latter have negative between-group coefficients. Their magnitude is scaled using reference means so expected added mass is equal across groups. This is a two-group addition design, not molecular downregulation or thinning. Realized totals fluctuate. The known group-independent addition enters the nuisance term so the FWL identity remains exact.

## Normalization and dtangle

Normalization uses all annotation-present shared genes with positive original reference mean, not merely the400 selected genes. Full native library totals are retained before gene restriction. TMM/DESeq2 receive no altered-gene labels. Gene-universe lists are exported.

Official dtangle2.0.9 receives fixed class markers, fixed reference treatment, and package-default RNA-seq gamma. The fixed denominator is the unavailable pre-intervention full native total: a mechanism diagnostic, not a deployable remedy.

Null/0.5 arms cross observed/fixed denominators with (a,c) in `log2(c*u+a)`: (0.1,1e6), (1,1e6), (10,1e6), (1,1e4), (1,1e8), (0.01,1e4). Reference and target use the same a,c. Zeros remain zero before the transform; a>0 makes them finite, with no hidden epsilon. Artificial positive-input tests check pure-log sample-scaling invariance and joint-scale equivalence (0.01,1e4)=(1,1e6). Empirical zeros never enter pure log.

The ordinary (1,1e6) setting also uses edgeR TMM effective library sizes and DESeq2 `estimateSizeFactorsForMatrix(type='poscounts')`. DESeq2 factors are geometrically centered, then multiplied by the geometric mean observed native library total. TMM uses observed native library sizes and official factors. The frozen reference remains on its original full-library scale. These are specified deployable target-normalization comparisons, not exhaustive native preprocessing workflows. Strong composition shifts can violate their assumptions; improvement is not presumed.

Repaired weighted NNLS, oracle composition and unadjusted regression remain controls. Most dtangle arms use the manuscript's reference-yield conversion; original observed/fixed arms additionally retain unconverted outputs to reveal scale dependence.

## MEAD

Unmodified MEAD1.0.2 is locked by commit. Up to2000 eligible genes are ranked by across-class reference mean-count variance divided by mean (gene-ID tie break), after excluding every tested gene. This reference-only panel is not claimed to be the default/recommended native workflow.

The four reference draws retain their true donor origin in metadata; duplicated donors receive draw IDs for the outer resampling design only. Selected-gene zero-library cells are excluded and each draw-by-class must retain support. Exported official `MEAD_getX` computes reference moments/weights; `MEAD_est` computes measurement-error-adjusted point predictions. No upstream numerical code is modified.

Fixed settings: `calc_var=FALSE`, `groups=NULL`, no QP, no centering, default truncation. Native gene filtering is not invoked; this deviation and the frozen panel are explicit. Class order is checked. Native and reference-yield-converted outputs are both evaluated against realized nucleus fractions; neither scale is presumed calibrated. The selected-gene RNA yield may differ from full-library yield. An additional `lambda=0` ablation calls the SAME official estimator with the SAME panel/moments/weights but removes the reference measurement-error subtraction; it is a diagnostic ablation, not a new competing method. Its paired comparison with lambda=1 isolates that component without confounding the gene panel.

Since tested genes are excluded, injecting them leaves MEAD inputs unchanged. Cross-condition invariance is a design control, not independent evidence of general robustness. Degeneracy, zero prediction mass, numerical failures or downstream rank loss remain failed trials. MEAD proportion covariance is NOT propagated to the downstream expression regression; the question concerns estimator choice with the same HC3 analysis, not all of MEAD's theoretical guarantees.

## Metrics and provenance

All methods share mixtures, responses, noise and test family within each condition. HC3, t reference distributions, basic availability gate and BH0.05 are unchanged. Unavailable p-values are1 before full-family BH. Conditional coverage, delivered coverage and availability are separate; conditional contrasts state their denominator. Each fit retains beta, SE, leakage, noise projection, p/q and availability; FWL is verified.

Paired contrasts include each method minus original dtangle and observed minus fixed at the SAME a,c. As of summary version 1.0.2, source means first average finite repetitions within each target donor, then give every planned target donor equal weight. This deliberately replaces the earlier repetition-first conditional summary, whose donor set could vary with failures. Full-data/all-planned means are unchanged. Conditional failures are not imputed as zero. A donor with no finite observations makes the full-bank conditional mean undefined; contributor counts and reasons are reported. Donor-specific tables/plots are retained.

For donor j, let theta_j be the finite-repetition mean and p_j the fraction of finite repetitions among R planned repetitions. The paired-repetition influence is h_r = mean_j[I_rj*(x_rj-theta_j)/p_j]. The conditional Monte Carlo SE is sd(h_r)/sqrt(R), and the normal MC interval is mean_j(theta_j) +/- 1.96*SE. This ratio-estimator delta method retains within-repetition covariance across donors and handles unequal success fractions without changing the equal-donor estimand. It reduces to the original repetition-mean SE when all values are finite. At least two finite repetitions per donor are required for a reported SE. The same construction is applied to paired method differences using their specified common-fit or all-planned denominator. These intervals are conditional on fixed donor banks, not donor-population intervals.

Summary version 1.0.2 preserves the literal null condition when reading CSV and rejects missing condition/method/repetition keys. The summary-only repair reads sealed chunk tables and retains original computation provenance in a separate new folder; it does not rerun simulations or overwrite original checkpoint signatures. Blank legacy gene-summary labels are restored only after matching the planned null methods, full tested-gene family, and all-zero truth. Large fit arrays are not rehashed by the repair, and the compact report is not a full reproducibility archive.

SHA-256 identities, source manifest, package versions, gene selections, conditions and seeds are saved. Code/parameter/environment changes refuse in-place resume. Completed chunks/units have checksums. No old submission or research input is overwritten. Exporting verified prepared inputs enables portable execution; licensing and actual public deposition remain separate author actions.
