# Normalization coupling and downstream adjustment

## Proposition 1 Finite change and sensitivity

For sample i and marker g, let x_ig be a nonnegative raw count, D_i > 0 the
normalization denominator, c > 0 the multiplier, and a >= 0 the pseudocount.
Write Z_ig(D_i) = log2(c x_ig / D_i + a), on its finite domain. If an intervention
outside the fixed marker panel changes D_i to D_i + delta_i > 0 but leaves x_ig
unchanged, then exactly

Z_ig(D_i + delta_i) - Z_ig(D_i)
= log2[D_i / (D_i + delta_i)]
  + log2[(c x_ig + a (D_i + delta_i)) / (c x_ig + a D_i)].

Differentiating with respect to log D_i gives

dZ_ig / dlog(D_i) = -[c x_ig / (c x_ig + a D_i)] / log(2).

Thus marker/test disjointness alone does not isolate the marker input when the
denominator includes test genes. With a > 0, the shift generally differs between
markers. With a = 0 and strictly positive counts, it is a common sample-wise
shift. Prediction invariance then requires the estimator to be invariant to that
shift; it is not a property of every deconvolution algorithm. Official dtangle
is checked numerically on strictly positive input. Zero counts are not silently
replaced to claim the pure-log identity.

Scaling a and c by the same positive constant adds a common log offset; it does
not change a/c. Shrinking a/c approaches the common-row-shift limit for positive
counts. This controls one route of denominator dependence, not all error in the
estimated covariates. It does not imply unbiased coefficients or calibrated FDR.

## Proposition 2 A matched marker denominator

Fix a marker set M before examining target outcomes. Let reference mean counts
per nucleus be mu_gk and s_Mk = sum_(g in M) mu_gk. Define A_Mgk = mu_gk/s_Mk,
and normalize target marker counts by D_Mi = sum_(g in M) x_ig. Under an exact
linear mixture with matched reference/target recovery, the normalized marker
mixture has weights q_Mik proportional to p_ik s_Mk. Therefore the corresponding
recovered-nucleus conversion is p_ik = (q_Mik/s_Mk)/sum_j(q_Mij/s_Mj).

Using the full-library yield in this conversion would change its meaning. The
new code uses marker yields for both marker-normalized arms. When intervention
changes only non-marker counts, both x_M and D_M are unchanged, so the candidate
input is exactly invariant. This requires the same fixed panel, positive marker
totals, fixed reference fitting and no other jointly estimated normalization.
Disease effects in markers, panel reselection, unequal capture, batch effects,
reference mismatch and zeros can break the conditions or prediction accuracy.

The primary exploratory arm uses c_M = 10^6 r_ref, a = 1, where r_ref is the
pooled reference marker-count fraction. This aligns the typical scale without
target labels. A second matched-reference arm uses c_M = 10^6. A denominator-only
control uses D_M/r_ref with the original full-library reference and yield; that
control is NOT presented as a coherently scaled new estimator. All three remain
visible even when their performance is poor. No method is selected by pilot FDR.

## Proposition 3 Propagation to the expression coefficient

For y = beta x + u + epsilon, let C(P) contain an intercept and the nuisance
column space induced by the estimated proportions P (and any fixed covariates).
Let z(P) = [I - projection_C(P)]x. If z'z > 0, then

beta_hat(P) - beta = z(P)'u/[z(P)'z(P)] + z(P)'epsilon/[z(P)'z(P)].

Changes in marker input can change P, its nuisance span, z, and either projection.
The identity is checked against each estimable fit. A change in P need not change
its column space, and a smaller denominator response need not reduce either bias
or total rejection error. This proposition does not assign an optimal pseudocount.

## Reporting and rank sensitivity

The full-six rule reproduces the existing prespecified reporting requirement.
The adaptive rule spans the actually observed nuisance columns with SVD; rank is
determined at max(1e-12, largest singular value times 1e-10). It reports the
exposure coefficient only when its residual sum of squares exceeds the declared
tolerance, residual df >= 2, max leverage < 1 - 1e-9, and HC3 SE is finite and
positive. Both rules use the same MEAD point predictions and outcomes.

Delivered coverage is P(interval available AND covers the true effect); conditional
coverage uses only available intervals. Unavailable tests receive p=1 in the fixed
100-gene BH family. Exposure R2 is an alignment/collinearity diagnostic, not a
new estimator of true bias. Real disease data have no effect truth: the application
reports coefficients, agreement, sensitivity and availability, not true FDR.

References: Frisch and Waugh (1933), doi:10.2307/1907330; MacKinnon and White
(1985), doi:10.1016/0304-4076(85)90158-7; Hunt et al. (2019),
doi:10.1093/bioinformatics/bty926. Formulas here are algebraic derivations for this
specified workflow, not a claim of a new general normalization theorem.
