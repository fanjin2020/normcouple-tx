# Conditional bounds and power diagnostics

These statements are algebraic results with explicit conditions. They are not a
proof of nominal false discovery rate and do not require that MRN outperform TMM.

## 1. A finite normalization-input bound

For a fixed marker count x >= 0, scale c > 0, pseudocount a > 0 and two positive
denominators D0 and D1, define Z(D) = log2(c*x/D + a). Then

```
|Z(D1)-Z(D0)| <= [c*x / (c*x + a*min(D0,D1))] * |ln(D1/D0)| / ln(2).
```

Proof: dZ/dln(D) = -[c*x/(c*x+a*D)]/ln(2). Its magnitude on the log-denominator
interval is bounded by its value at the smaller denominator. Integrate the
derivative across that interval. The bound also covers x=0, where both sides
are zero. It does not cover a changed marker count: compare paired outcome
conditions at the SAME marker-stress level.

For fixed marker set M, raw marker counts and reference, MRN has
D_M = sum_{g in M} x_g. Perturbations outside M leave D_M and every prediction
input unchanged. The matched reference profile A_gk/s_k and yield q_k*s_k
must be retained, where s_k=sum_{g in M} A_gk and every s_k>0. The scaling
anchor is fixed from reference counts. Zero target marker libraries are not
silently replaced by pseudo-observations. Marker changes, panel changes or
reference changes invalidate this particular invariance argument.

## 2. A conditional downstream coefficient bound

Let P0 and P1 be orthogonal projectors onto the actual nuisance column spaces,
both containing the intercept. Let z_j=(I-P_j)x, q_j=z_j'z_j>0 and w_j=z_j/q_j.
For the same response y and exposure x, w_j'x=1. For any fixed coefficient tau,

```
beta1-beta0 = (w1-w0)'(y-tau*x)
|beta1-beta0| <= ||w1-w0||_2 ||y-tau*x||_2
||w1-w0||_2 <= ||z1-z0||_2/q1 + ||z0||_2 |q1-q0|/(q0*q1).
```

The first inequality is Cauchy-Schwarz. For the second write
w1-w0=(z1-z0)/q1 + z0*(1/q1-1/q0) and apply the triangle inequality.
The reverse orientation gives another valid bound; the implementation takes
the smaller of the two. Also z1-z0=(P0-P1)x, so a predictor-to-projector
stability bound, IF separately established, can be inserted here.

There is no unconditional bridge from the input bound to a coefficient bound:
it requires stability of the deconvolution map and the nuisance column space,
plus residualized-exposure information bounded away from zero. Near rank
changes such a bound can be uninformative. The code does not invent a global
Lipschitz constant or equate a small input change with calibrated inference.
Simulation truth tau is used only in this diagnostic; it is not assumed known
in a real cohort. In real data, a bound with tau=0 remains valid but can be loose.

## 3. Power cost without a presupposed explanation

The paired generator supplies y=tau*x+u+e. Since w'x=1,

```
beta-tau = w'u + w'e = leakage + noise_projection
Var(beta-tau) = Var(leakage) + Var(noise_projection)
                + 2 Cov(leakage, noise_projection).
```

The sample variance identity is checked separately for every gene, outcome,
target donor and method with at least two common available repetitions. The
same repetitions enter all three methods. Missing fits are counted and reported.
These components can be correlated; neither component alone is the variance.

For HC3, the exact variance calculation is sum_i [w_i*r_i/(1-h_i)]^2.
The reusable actual-span function returns these contributions. The grid stores
gene-level residual variance, largest sample contribution, standard error,
coefficient, signed/absolute t, BH result and both FWL projections for the
unstressed null, up05, and internal up025 controls. Design R-squared, z'z and
leverage are in paired summary tables. No claim of reduced z'z is hard-coded.

Observed correlations with marker library size are descriptive, not a causal
mediation analysis. MRN changes the denominator and reference definition, not
the physically sequenced reads. A lower marker sublibrary does not by itself
prove why a coefficient or its standard error changed.

## 4. Scope of uncertainty

Grid intervals use within-repetition paired differences, average first within
target donor and then equally across donors, and report Monte Carlo uncertainty
conditional on the fixed donor banks. They do not treat repeated mixtures as
additional biological donors. Conditional metrics report contributing donors
and repetitions; FDR/power/availability retain their all-planned denominators.
If a planned donor has no finite conditional endpoint, the full-bank mean and
interval are undefined. A separate contributing_donor_mean records the mean
over the remaining donors without silently relabeling it as the full-bank result.
