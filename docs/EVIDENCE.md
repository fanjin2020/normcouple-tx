# Six evidence checks

The server entry is `bash run_evidence.sh`. It runs `reproduce.py evidence` with
the existing scientific Python and R libraries and does not install packages or rerun
the 16,800-dataset simulation. One R process at a time bounds memory use; a wider
simulation worker pool does not speed these matrix fits. An interrupted invocation
can be repeated: finished cohort/reference pairs are reused only after checksum
verification. A changed input or framework requires a different output directory.

The default search root is `$HOME/git/repository`. The sealed inputs are
`p3_complete_20261117/clinical_final` and `p3_complete_20261117/cross_full`;
the destination is `p3_evidence_20261119`. Override with `--clinical`, `--cross`,
`--output`, `--rscript`, and repeatable `--r-library`. `PAPER3_PYTHON` overrides
the analysis Python. The standalone server handoff embeds the verified source
package and three public imaging-study files outside it; no SSH, Git or install
step is involved. In a source checkout, downloads are limited to these three files.
Use `--imaging-cache DIR` and `PAPER3_OFFLINE=1` for a prepopulated offline cache.

1. Check the manuscript's four newest tables against sealed estimator outputs.
   Rebuild S16 from the original metric summaries. Regenerate S20 and S21 from the clinical probe and imaging records, check both against the archived article tables, and render 20 figures from the verified numerical inputs.
2. Read the actual coefficient-identity errors in the clinical multistep records.
   The S13 note reports the maximum finite value and corresponding method, step
   and reference; no clinical biological truth is assumed.
3. Read the frozen simulated-probe design and clinical marker lists. For every
   perturbation, report the input fraction and realized library/effective TMM
   denominator changes separately by diagnosis (median and interquartile range).
   Counts of zero remain zero in this clinical multiplicative intervention.
4. Run a deterministic small example. For separate use:
   `python reproduce.py example --output ../example`. It writes three inputs
   for `reproduce.py sensitivity`; permuted IDs align, missing IDs fail and a
   different basis for the same nuisance space gives unchanged coefficients.
5. Choose 10th, 50th and 90th abundance-rank test genes without phenotype or
   significance information. Apply exp(+0.1) and exp(-0.1) to each singleton, its
   five-gene neighborhood, and the full tested family in AD samples. Refit actual
   dtangle/TMM and marker-restricted dtangle, keeping outcomes, covariates and BH
   family fixed. Baseline, all-family-up and marker-up predictions are checked
   against the original predictions. Report own and other genes separately;
   unavailable specifications do not generate zero-valued discoveries.
6. Compare those two frozen specifications with RNA imaging in held-out LIBD
   donors. The upstream input is the authors' released *filtered* bulk count
   matrix, not a newly sequenced cohort. Use total-tissue extracts only. Sample
   and block identifiers must match the author crosswalk exactly after normalizing
   region case; ambiguous suffixes are not guessed. The imaging comparison keeps
   High/OK sections, uses the measured fraction of all nuclei without renormalizing,
   sums predicted Oligo and OPC, and does not impute unmeasured classes. Marker
   yields are recomputed on the shared panel, keeping the frozen reference anchor.
   Average technical libraries/sections within tissue block, then average block
   errors within donor; the reported RMSE is the square root of the donor-equal MSE.
   Effective numbers of blocks and donors are reported per cell class. The released
   R count object and its text export cover the same 21,745 genes. Do not weaken a
   ten-marker minimum because a frozen clinical panel is poorly covered. Instead,
   report that panel's mapped counts and estimability. A separately identified
   reference-transfer panel is evaluated on the same blocks for both normalizations:
   train only on the four archived reference donors, using common gene IDs, at least
   10 reference CPM in the owning class in at least three donors, and a positive
   log ratio of the strongest vs second class mean CPM (add 1 to each). Keep up to
   25 highest-scoring markers per class, break ties by gene ID. Fit both scales only
   with at least 10 markers per class. Neither measured fractions nor target bulk
   expression are used to select genes. This separate experiment cannot validate
   the original clinical panel when that panel is not estimable.

Imaging source: [Huuki-Myers et al., Genome Biology](https://doi.org/10.1186/s13059-025-03552-3).
Counts, block crosswalk and imaging fractions are pinned to
[the original repository](https://github.com/LieberInstitute/Human_DLPFC_Deconvolution/tree/96131ef436fdaa9a9b2129983ac7bd6868204fc4);
the authors' count export reads the counts assay, despite its `GeneSym` column
containing Ensembl identifiers. This experiment measures proportion error, not
clinical FDR. Keep block-level measurements together when reanalyzing.

To regenerate the article tables without refitting, run `python reproduce.py tables --published --evidence ../paper3_review/evidence --output ../paper3_tables`. S20 verifies the clinical samples, probe assignments and group-specific denominator series. S21 recomputes six cell-class errors for each reference from paired block records and checks both the donor support and the archived error summaries.

Return **paper3_review.zip**: measured doses, compact own/other effect summaries,
selected genes with baseline mean CPM, selected-gene effects, marker definitions and per-class quotas, imaging errors and matched blocks,
the matching and protocol records, selected manuscript CSVs and figure/table
mapping. Full replay PDFs are in `manuscript_replay/figures`. Working matrices
are never copied into the review archive. Completion is printed only after all
checks and a ZIP integrity test; partial stages are not scientific conclusions.
