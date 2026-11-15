# Coefficient propagation under a fixed response

Let x be the exposure, y a response (or a matrix of gene responses), and W the
nuisance design including an intercept and composition covariates. Define
P = W W^+, M = I - P, z = Mx. If z'z > 0,

\[
\hat\beta(W) = (x'My)/(x'Mx).
\]

For two estimable nuisance spaces and the SAME x and y, subtraction gives the
finite identity

\[
\hat\beta_1-\hat\beta_0 =
\frac{x'(M_1-M_0)(y-x\hat\beta_0)}{x'M_1x}.
\]

Indeed x'M_0(y-x beta_0)=0. Adding and subtracting that zero in
x'M_1(y-x beta_0) proves the identity. It remains an identity across changes
in nuisance rank, provided both exposure coefficients are estimable. It is a
change of fitted coefficients, not a population-bias or causal-effect theorem.

When W has full column rank locally and changes differentiably,

\[
dP = M(dW)W^+ + (W^+)'(dW)'M,\qquad dM=-dP,
\]

\[
d\hat\beta = \frac{x'(dM)(y-x\hat\beta)}{z'z}.
\]

These formulas show why a global norm of composition error is insufficient:
the response residual, exposure direction, nuisance space and remaining
exposure information all enter the coefficient change. They do not imply
that any exposure-correlated composition difference is erroneous.

For fixed marker count m and scale c, with pseudocount a=1,

\[
Z=\log_2(cm/D+1),\quad
\frac{dZ}{d\log D}=-\frac{cm}{(cm+D)\log 2}.
\]

The experiment varies test counts by exp(t*x), leaving marker counts and y
unchanged. The input derivative is checked analytically; the composition
direction along this path is obtained by central differences of the installed
dtangle pipeline (including the fixed yield conversion). The projection
derivative then predicts the coefficient direction. This does not provide a
closed-form Jacobian of dtangle, and the finite-difference step is halved to
detect numerical/non-smooth sensitivity. TMM trimming/reference changes and
deconvolution active-set changes may violate local smoothness.

If y also changed, x'M(dy)/(x'Mx) would enter the differential; omitting it
would confound predictor construction with response changes. That is why the
paired sensitivity experiment holds the response fixed.

Removing tested counts from every predictor-construction path yields functional
invariance with respect to those counts. It does not give independence between
marker and test genes, remove shared biology, or establish FDR control. Those
are distinct statistical properties and are not asserted by this package.
