# Paper 3: randomized perturbed-marker subsets

This is one bounded robustness experiment. It does not replace the reviewed
manuscript, modify a GitHub repository, or download new research data.

## Fixed design

- Reuse the frozen LIBD/SEA-AD references, GSE174367 targets, test genes, count
  banks, mixtures, reference resampling and random seeds: 42 reference-target
  units, 27 distinct targets, n=80 mixtures, 200 repetitions by default.
- Randomize WHICH MARKERS RECEIVE STRESS, not the normalization/deconvolution
  marker panel. Draw 50 subsets independently, uniformly within each retained
  annotated cell-class panel. Sampling is without replacement within a subset;
  subsets may overlap. Do not select a subset using observed results.
- Request 10% per class, rounded UP separately in each of seven classes. The
  realized fraction is recorded, not described falsely as exactly 10%.
- Freeze subsets within each reference/annotation panel across donors, repeats,
  outcomes and methods. Different annotation panels may have different gene
  sets. Missing external annotations are excluded before subset sampling;
  they are not treated as observed zero-expression genes. Record exclusions.
- Apply Poisson additions with mean 0.5 times selected baseline marker counts
  in exposed mixtures only. Compare direct test-gene effects 0, 0.25 and 0.5.
  Underlying mixtures/outcomes are identical across subsets and methods.
- Include zero-stress and the previous lexicographic/class-stratified stress
  subset as controls. Preserve its original stress RNG exactly.
- Methods: oracle composition, TMM, marker-matched MRN, denominator-only control.
  Regression uses the frozen full-six HC3 rule, BH threshold 0.05, and the same
  cell-yield conversion as the previous experiment. Failure is an outcome.

## Run on the existing Linux server

Place the single delivered paper3_random_marker_oneclick.py in the project:

```bash
python paper3_random_marker_oneclick.py --setup --offline \
  --project-root /home/research/git/repository/paper1_public_evidence \
  --output /home/research/git/repository/paper3_random_marker_results_v1 \
  --workers 36
```

Python environments and the existing R library are reused. Required versions:
NumPy 2.2.6, pandas 2.2.3, SciPy 1.15.3, h5py 3.13.0, matplotlib 3.10.1;
dtangle 2.0.9, edgeR and jsonlite. R versions are captured and must match on
resume. If a Python environment is absent, removing --offline allows dependency
installation only. Research-data downloads are never implemented in this entry.
R packages are not installed automatically. Use --r-library if auto-detection
does not find the old R_library, and --rscript for a nonstandard R executable.

Input overrides: --input-bundle, or --libd-dir/--seaad-dir/--gene-map and
--external-cache. All input identities remain SHA256-checked. A prepared bundle
at ../paper3_mechanism_results_v1_prepared_inputs is detected automatically.

```bash
python paper3_random_marker_oneclick.py --status \
  --output /home/research/git/repository/paper3_random_marker_results_v1
```

Status needs only standard Python, no R or numerical packages. It counts
completed receipts and file sizes; add --verify for all content hashes. ETA
uses measured completed work, excludes preparation/summary time, and is not a
promise. The initial small numerical pilot reports measured cost. --pilot-only
stops after that pilot; omit it on the next invocation to continue.

## Resume and parallelism

Linux defaults to 36 fork workers; each numerical/R worker uses one BLAS/OpenMP
thread. Counts are shared read-only by fork. GPU is not used. Memory and disk
throughput still limit speed. Windows supports one worker for software checks.
Default chunks contain five repetitions by five subsets. Completed chunks have
checksummed receipts committed LAST. A partly written chunk has no completion
receipt and is rerun from its deterministic seeds. A corrupted completed chunk
raises an error rather than being silently reused or dropped.

Stop with Ctrl+C and allow owned child processes to exit. Restart the SAME
command to resume. --workers may change after stopping; code, inputs, seeds,
batch layout, scientific settings and package versions may not change within
an output directory. Never delete run.lock while a task is alive. Old research
outputs are rejected as destinations. --summarize-only rebuilds compact results
from complete, verified checkpoints without loading count banks or invoking R.

## Compact deliverables

Transfer only OUTPUT/final_results/, not the complete working directory.
No ZIP is produced. Small CSV.gz tables use ordinary gzip for individual tables.

- RESULTS_SUMMARY.md: data-derived summary, no preset favorable conclusion.
- random_subset_benefit_distribution.csv: fraction of randomized subsets with
  lower estimated FDR, subset-sampling Wilson intervals, effect/power tradeoffs,
  and counts of subsets with MRN FDR point estimate or MC upper bound <=0.05.
- subset_method_results.csv: per-reference/target/outcome/subset metrics and
  conditional Monte Carlo t intervals, plus finite-repetition/donor counts.
- paired_MRN_differences.csv: paired MRN minus TMM/denominator-only differences.
- donor_subset_results.csv.gz: retain donor heterogeneity instead of only means.
- stress_expression_and_mass.csv and identity_mass_descriptive_associations.csv:
  expression and added-count-mass diagnostics; descriptive correlations only.
- failure_counts.csv, selected_markers.json, protocol.json, QC.json, hashes.
- Two data plots, PNG/PDF: subset FDR distribution and FDR/power tradeoffs.

## Statistical interpretation

FDR is the mean per-repetition FDP, not a pooled false-discovery ratio. Failed
tests remain in planned denominators with p=1 and no rejection. Power uses all
planned alternatives; it is undefined at global null. Delivered coverage uses
all planned tests; conditional coverage uses available tests. Setting-level
conditional coverage is computed over available donor-tests within each repeat.
Paired conditional metrics use common finite donor values; their reported
counts show selection. Unconditional availability/FDP/power/delivered coverage
include all planned units. RMSE/SE/bias can be undefined, never zero-filled.

MC intervals average fixed target donors equally within each reference-target
setting and use independent simulation repetitions, NOT millions of genes or
subsets as independent patients. The same external targets under different
references are reported separately. Do not pool them as new donors.

Subset-sampling Wilson intervals condition on the realized paired simulations
and fixed donor banks. They quantify uniform subset selection, not uncertainty
over people or disease prevalence, and do not include finite-MC classification
error. Report paired MC intervals alongside them. Point improvement is not FDR
control. Even an MC upper bound below 0.05 is conditional simulation evidence,
not a general guarantee. No multiplicity-adjusted claim is made from the number
of individual MC intervals below zero.

Raw expression/count mass can vary among equal-size subsets. Diagnostics help
interpret this, but do not isolate a count-mass-matched causal identity effect.
MRN predictions are cached only when input identity is explicitly checked;
fresh refits verify cache equivalence. Cache hits are NOT independent new
numerical invariance replications. Frozen control generator identities are
checked separately.

This extension does not deliver new real-bulk validation, a new theorem, a new
normalization method, a public DOI, or evidence of how often markers change in
real disease. It cannot ensure journal acceptance.

## Software checks

```bash
python -m unittest discover -s tests -v
python run.py --synthetic --subsets 2 --repetitions 2 --workers 1 \
  --batch-size 1 --subset-batch-size 2 --output /tmp/paper3_random_marker_fixture
```

Fixture outputs are permanently marked SOFTWARE_FIXTURE_ONLY and cannot resume
into a research run. tests/test_r_bridge.py additionally compares the cached R
implementation against the frozen uncached bridge when PAPER3_TEST_RSCRIPT and
PAPER3_TEST_RLIB are supplied. The original numerical engine is vendored without
scientific changes; source_manifest.json includes all bundled code and tests.
