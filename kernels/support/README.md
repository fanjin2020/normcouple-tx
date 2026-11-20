# Paper 3: bounded follow-up of bdd118b137

This addon reads the completed full server run and writes a separate output. It does not edit the parent run, rerun its 8,400 main batches, change manuscript numbers, install packages, or download data. It embeds the exact parent source and verifies its identity.

## Run

Upload `paper3_next_oneclick.py` to `/home/research/git/repository/paper1_public_evidence/`.

```bash
conda activate paper1
python /home/research/git/repository/paper1_public_evidence/paper3_next_oneclick.py --workers 36
```

The script locates the full completed bdd118b137 run among immediate server/project child directories. If more than one is present, specify `--previous-results /path/to/the/full/server/run`. A compact `p3_review` export is not sufficient. The default new output is `/home/research/git/repository/p3_next_v1`.

The same command resumes this addon. Only sealed, checksum-verified addon work is reused. A changed protocol or corrupt receipt stops the run rather than replacing results. Do not delete lock files. Workers can be changed without changing the protocol.

```bash
python /home/research/git/repository/paper1_public_evidence/paper3_next_oneclick.py --status
python /home/research/git/repository/paper1_public_evidence/paper3_next_oneclick.py --status --verify-status
```

For corrected reports alone use `--reports-only`; later run the normal command to perform the new computations. `--preflight-only` verifies the parent, environments and prepared banks without performing the new clinical fits or simulations. `--export-only` packages a previously generated review directory. These commands must use the same `--output` if it differs from the default.

## 1. Completed-result corrections, without new fits

- Retain the original diagnostic assessment, including its non-smooth scores. Add a separate reliability-qualified assessment requiring an estimable baseline, stable nuisance rank and derivative agreement at most 0.1.
- Set qualified scores at or below the existing 1e-6 decision floor to zero for both ranking and alerts. Report this numerical-invariance group explicitly. Do not interpret noise-ranked NNLS AUC as predictive discrimination.
- Relabel non-finite or unavailable baseline scores. Report the number of excluded and qualified scores rather than implying a missing score is a negative result.
- Report planned versus actually contributing donors, subset pairs and repetitions for each metric. Keep the original equal-subset/equal-target-donor weighting.
- Export the existing 161,610 disease-association rows and 323,220 fraction-scale association rows from the verified clinical receipt. Independently check the frozen gene family, BH calculations, availability and discovery-set summary.

The reliability filter is a post-hoc audit of the same predictions. It is not a newly validated FDR alarm, and it does not erase the original diagnostic result.

## 2. Clinical multistep check

Reuse the same verified public raw-count cache, 288 donors, age/sex adjustment, 16,161 test genes and two frozen references. Verify persistent matrices against public counts, then check newly computed baseline fractions and beta/SE/p/q against sealed parent results before interpreting perturbations.

Evaluate all derivative steps 0.01, 0.001, 0.0001 and 0.00001, with their half steps, and finite positive perturbations 0.01, 0.05 and 0.10. The response remains the original all-gene TMM log2 CPM+1. Report all step sizes, non-smoothness, linear approximation error, finite coefficient changes and threshold crossings. No best step is selected from the observed results.

The clinical output measures sensitivity, not true clinical FDR, causal effects or correctness of any changed gene list.

## 3. Common-positive-support equal-budget experiment

Default design: 42 reference-target units x 50 new repetitions (IDs 10000-10049), the seven frozen stressed marker sets, two support rules, a zero-stress arm, two outcomes (null and effect 0.25), and the existing eight inference arms. This is 2,100 new batches and 504,000 method-outcome records. The old 200 repetitions are not rerun.

For each unperturbed draw, derive the same legacy sample-by-class count budget. Retain a sample/class budget only when EVERY one of the seven frozen marker sets has a positive observed count with a positive reference weight in that sample/class. Apply that single reduced budget to EVERY stressed set under BOTH support rules. Marker identities are not redrawn and outcomes are not inspected to choose support.

The allow-zero arm can allocate within any reference-supported selected marker. The positive-only arm allocates only within selected markers already observed positive. Reference weights, donor draws, outcomes and seeds are paired. All draws, including a zero retained budget, remain in the primary analysis. A secondary positive-dose stratum is explicitly labeled.

**This changes the intervention dose and target condition.** It is not a repair of infeasible original interventions. Read `dose_retention.csv` and the class-level budget table alongside FDR and power. The old infeasibility results remain valid and are preserved. If the common dose is usually zero or tiny, this experiment cannot establish broad positive-support performance.

Fifty repetitions per unit constitute a bounded follow-up. Within-unit Monte Carlo uncertainty is reported; data-dependent subset draws are averaged inside each repetition before computing its SE. Donor reuse across reference settings does not add independent biological replicates. These results do not establish a universal 5% FDR guarantee or an automatic model-selection rule.

## Outputs and review

The script prints one short archive path, `p3n_<10-character-signature>.zip`. Transfer that archive only. It contains:

- `corrections/`: original diagnostic performance, qualified diagnostic records, effective sample counts, clinical gene-level tables and checks.
- `clinical_multistep/`: baseline identity checks, all derivative steps, finite-perturbation summaries and provenance.
- `common_support/`: method summaries, paired comparisons, Monte Carlo uncertainty, exact dose retention, and compact follow-up metrics.
- Protocol, source snapshot, input/source hashes, runtime versions and per-file integrity manifest.

Clinical binary inputs and per-batch working files remain on the server. The archive is compact relative to the full run, but retains the gene-level clinical tables needed for verification. It is not a manuscript supplement as-is.

## Environments and safety

Use the recorded existing Python environment with numpy 2.2.6, pandas 2.2.3, scipy 1.15.3, h5py 3.13.0 and matplotlib 3.10.1. Use the parent R version and package versions exactly. The launcher checks these; it never installs or downgrades packages. Explicit overrides: `--python`, `--rscript`, repeatable `--r-library`, `--cache` and the prepared-bank input options shown by `--help`.

The parent lock is held shared/read-only; the new output has an exclusive writer lock. Each R child runs in its own process group and is stopped with its owning worker on interruption. Native linear-algebra thread counts are one per worker. The initial batch runs alone to detect infrastructure problems before starting the pool.

No duration is promised before new batches have finished. `--status` estimates simulation-only remaining time from at least three completed addon batches and suppresses the estimate when no writer is active.

## Not performed

No new biological cohort, native MuSiC/Bisque comparison, additional cell-type shift directions, clinical truth-based FDR evaluation, optimal denominator estimator, automatic model selection, raw-to-prepared-bank regeneration, manuscript update or GitHub upload. The existing verified LIBD/SEAAD prepared banks are still required.

## Testing boundary

`test_followup.py` covers budget invariants, zero-support cases, reliability gating, numerical ties, effective sample counts, clinical export checks, baseline drift detection, receipts, resume and an end-to-end Python bridge fixture. The bridge fixture is labeled software-only and is not dtangle/TMM validation. Native R execution must be confirmed in the recorded server environment; the script stops on version or baseline disagreement.
