# Paper 3 review enhancements

This isolated program adds reference-matched marker stress, support-rule controls,
linear NNLS comparisons, coefficient propagation, donor-held-out diagnostics and
clinical RNA-contribution versus nuclear-fraction comparisons. It leaves previous
analyses and manuscript files unchanged. It installs no Python or R packages and
does not fetch code from GitHub.

## One-command research run

Upload `paper3_review_enhancement_oneclick.py` to the server project folder.

```bash
conda activate paper1
python /home/research/git/repository/paper1_public_evidence/paper3_review_enhancement_oneclick.py \
  --project-root /home/research/git/repository/paper1_public_evidence \
  --random-results /home/research/git/repository/paper3_random_marker_results_v1 \
  --output /home/research/git/repository/p3_review_v1 --workers 36
```

The default `--phase confirm` runs 42 reference-target units, 200 repetitions,
three outcome conditions, and eight method arms including oracle. There are three
frozen random stress subsets and three reference-feature-matched stress subsets,
in addition to the predefined subset and no-stress control. These are subsets of
the fixed marker panel, not independently discovered complete marker panels.

Each nonzero subset is evaluated with both allow-zero and positive-only allocation.
The no-stress control is computed once. This gives 15 intervention variants,
8,400 checkpoint batches and 3,024,000 planned method-outcome records. Infeasible
positive-only interventions remain in the plan with an explicit status. They are
not invented, fitted as another intervention, or reported as zero false discoveries.

For a smaller operational pilot, add `--phase pilot` and choose a different output.
It uses 10 repetitions and two random/matched subsets: 420 batches and 110,880
planned method-outcome records. The pilot cannot support final calibration claims.
No favorable-result gate chooses the confirmation settings. Changing phase requires
a new output. Changing only worker count permits resume.

## Inputs and environments

The original random-marker protocol supplies the frozen generator, input identities
and package versions. Its previous fitted results are not imported. Input matrices
and downloads may be reused after identity checks; all new fits are computed here.

Exact LIBD/SEA-AD prepared count banks, their gene map, and GSE174367 source objects
must already be available. Discovery uses the previous protocol's input paths and
the portable prepared-input directory. Overrides are `--input-bundle`, `--libd-dir`,
`--seaad-dir`, `--gene-map` and `--external-cache`. This program does not invent a
public address for unpublished prepared matrices or claim a new raw-to-bank recipe.

Missing GSE125583 public bulk objects may be downloaded. Existing objects must match
pinned source hashes and receipts. `--offline` disables those downloads as well.
Package downloads are always disabled. Python pins are numpy 2.2.6, pandas 2.2.3,
scipy 1.15.3, h5py 3.13.0 and matplotlib 3.10.1. Existing environments are checked
by imports. R requires dtangle 2.0.9, edgeR and jsonlite with the original protocol's
recorded versions. Use `--python`, `--rscript` and repeatable `--r-library` as needed.

Each worker has one BLAS/OpenMP thread. Linux fork shares prepared count banks.
The program does not silently change the requested worker count. The first complete
batch runs before launching the worker pool to catch adapter failures early.

## Scientific definitions

- Matching uses only the reference donors: log1p own-class count abundance,
  own-class detection, and own-class fraction of mean expression across classes.
  Matching stays within cell class, excludes the predefined genes, and preserves
  class quotas. Assignment uses standardized squared distance plus fixed-seed
  uniform jitter in [0, 0.1]. Weak matches and duplicated matched sets are visible
  in the exported gene assignments; there is no outcome-dependent resampling.
- Predefined, random and matched subsets receive identical sample-by-class integer
  budgets from the frozen predefined relative-stress draw. Allocation probabilities
  use fixed training-reference weights. Positive-only allocation renormalizes them
  among selected genes with positive baseline counts. An empty eligible set with
  positive budget makes that intervention infeasible, without a uniform fallback.
- The five dtangle arms are all-gene TMM, factor-only test exclusion, complete
  factor-and-library test exclusion, marker denominator only, and marker matched.
  Oracle uses known sampled nuclear proportions. The two NNLS arms use linear
  nonnegative least squares with full versus marker-normalized reference profiles
  and corresponding yield conversions. They use linear inputs, not dtangle log inputs.
- Outcomes, exposures, reference resampling and test families are fixed within
  paired comparisons. Non-marker sensitivity probes modify predictor construction
  only. Expected input invariance is checked for complete exclusion and marker
  restriction, and scale invariance is checked for linear NNLS. This does not claim
  that every deconvolution algorithm shares the same behavior.
- Local derivatives are conditional on stable nuisance rank. Step disagreement and
  finite-step approximation error are reported. Exact coefficient identities are
  software checks, not independently predicted biological effects.
- The diagnostic score is the 90th percentile of absolute local coefficient change
  divided by baseline HC3 SE. Five target-donor folds keep repeated references
  together. A training-only 90th-percentile threshold (floor 1e-6) is applied to
  held-out donors. The label is realized FDP > 0.10, not population FDR.
- Clinical analysis retains 288 samples and 16,161 test genes, one fixed all-gene
  TMM log2-CPM response, age/sex covariates and a fixed BH family. q and p adjustments
  are fitted separately. Marker-matched q is marker RNA contribution; full-profile
  q has a different reference interpretation. No clinical FDR truth is asserted.

MuSiC/Bisque native analyses, another drift direction, a new cohort and real-bulk
spike-in are not part of this bounded run. These are named as not performed in the
result scope file. Optional recommendations are not reported as completed modules.

## Progress, resume and delivery

```bash
python paper3_review_enhancement_oneclick.py --status --output /home/research/git/repository/p3_review_v1
```

Add `--verify-status` for complete batch SHA256 checks. Resume with the original
command. A completed batch is verified before reuse; an incomplete batch is rebuilt.
Never delete a lock file to allow two writers. Use `--preflight-only` to check
environments, inputs and the fixed protocol without fitting. Use `--summarize-only`
to regenerate summaries from complete batches. `--mode simulation` and
`--mode clinical` request individual modules and do not imply completion of both.

The finished run generates `OUTPUT/final` and a short adjacent archive such as
`/home/research/git/repository/p3r_0123456789.zip`. Bring back only this ZIP for
review. It contains compact numeric summaries, matching assignments, clinical
fractions and diagnostics, current-run plots, provenance and checksums. Raw matrices
and per-batch fitted records remain outside it. Full records remain in the research
output for tracing any summarized value.

```bash
python paper3_review_enhancement_oneclick.py --export-only --output /home/research/git/repository/p3_review_v1
```

`--synthetic` is a generated software fixture, explicitly labeled as such. Unit
tests can run without R; tests requiring actual dtangle/edgeR must be run with
`P3_TEST_RSCRIPT` and `P3_TEST_R_LIBRARY`. A passing Python test suite is not a full
R research rerun.
