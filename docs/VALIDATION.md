# Validation scope

## Publication replay (2.2.2)

The submission contains nine main figures and eleven supplementary figures. The full replay and the tables-only path use the same figure-source index. Check the Guide, including the two logical views of S7, Fig. 9 (S13 and S20), Fig. S10 (S18 and S19) and Fig. S11 (S21). Individual S7 cell-class IDs are numbers in Excel; the pooled class is the text `all`. The renderer rejects incomplete Fig. 6 subset sets or conflicting Fig. 9 dose, family or threshold counts. These checks concern publication outputs and do not repeat the statistical fits.

## Numerical evidence in the archived tables

The figure mapping covers nine main and eleven supplementary captions. S20 is generated from the measured dose, selected-gene and paired refit records; its full-family positive controls must match S19. S21 requires matching tissue blocks and measurements in both fits; losses are averaged over blocks within donor before averaging over donors. The publication manifest retains the producing framework checksum separately from the packaging checksum. Replaying these data does not refit the clinical models.

## Previous development checks

The 2.2.0.dev0 code has passed 78 Python tests on Windows and Linux. The new R numerical
test recovered a synthetic mixture with RMSE 5.88e-5, checked scaling invariance
and fraction closure, and rejected zero-library inputs. The four-method pilot
completed 16 datasets across all four reference-target settings. In the first
360 method-probe records shared with the preceding experiment, FDR, power,
coverage, proportion RMSE, standard errors and design diagnostics agreed exactly.

The six-policy clinical analysis completed for the 288-sample GSE125583 and
17-sample GSE53697 cohorts with two reference banks. Of 72 planned method-probe
combinations, 12 lack sufficient marker coverage and three have non-estimable
HC3 variance. The remaining 57 support inference. An independent full-design SVD
implementation recalculated 935,676 gene-model records, including the three
non-estimable fits. All BH decisions agreed; maximum differences were 7.68e-13
for coefficients, 1.08e-12 for standard errors, 2.37e-12 for p values and 2.56e-12
for q values. The original clinical baseline discovery counts were reproduced.
All 22 checked cached/uncached MuSiC sample predictions agreed exactly. Estimator
outputs were reused after exact scientific-protocol, input and plan verification
for the final downstream reanalysis; this reuse is recorded in the receipts.

The 200-repetition experiment and real-input nine-stage acceptance run are still
running. Their final results and all-stage export have not been signed off.
Figures S8-S9 passed pilot-data layout checks; Figure S10 was rendered from the
independently checked clinical results. Acceptance uses fewer repetitions
and is explicitly labelled as software testing, not manuscript evidence. It does
not replace the full scientific runs. The archived manuscript data remain unchanged.

## Previously validated 2.1.0 results

Version 2.1.0 was checked on November 17, 2026. Scientific calculations ran on Linux with the pinned Python environment and R packages. Document and workbook checks ran on Windows. The remainder of this section describes that version, not a completed 2.2.0 run.

## Scientific checks

- The cross-estimator experiment completed all 210 batches: 1,680 datasets, 15,120 method-probe records and 10,080 paired contrasts. It uses 42 reference-target units and 20 repetitions per unit and effect. Its receipt records input, source and environment identities.
- Cached MuSiC reference calculations were compared with 192 uncached official API calls. The maximum prediction difference was zero. MuSiC 1.0.0 is pinned to commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c`.
- Independent aggregation reproduced all cross-estimator means and conditional Monte Carlo standard errors within 2e-15. The coefficient-weight identity had maximum absolute error 1.34e-15. All planned regressions were estimable.
- A final-launcher pilot reproduced 144 method-probe records and 96 contrast records exactly. A subsequent resume reused all four verified batches without recomputation.
- The user-data sensitivity command reproduced 161,610 gene-model results for 288 clinical samples and 16,161 genes. Availability and BH decisions agreed exactly. Maximum absolute differences were 8.52e-14 for coefficients, 1.71e-15 for standard errors, 1.63e-12 for p values and 1.90e-12 for q values.
- `rebuild-inputs` reconstructed both sparse count matrices, four cell/feature annotation files and the gene-identity map from cached author-distributed public source objects. All seven output SHA256 hashes matched the prepared-input definitions exactly. This was a reconstruction test, not a repeated download of all eleven source objects.

## Framework and publication checks

- All **65 framework/publication tests** passed on Linux. Coverage includes transformations, actual-span HC3, fraction closure, sample and family alignment, cross-estimator contracts, dependency ordering, input reconstruction, path safety, receipt verification and interrupted-run recovery.
- The preceding release's 17 common-support kernel tests also passed in its materialized runtime. The 2.1.0 release retained those kernel files unchanged. Development 2.2.0 adds real-input acceptance scopes and current-source parent verification; their runtime checks are separate.
- Aggregation from existing scientific outputs and the new cross-estimator output reproduced all 27 publication datasets: 16 supplementary tables and 11 additional figure datasets. Old table values were retained; Table S10 gained two geometry fields and Table S16 contains the new experiment.
- Published-result export generated all nine main and eight supplementary figures. Every figure has an explicit source-data mapping, with no duplicated table/plot data copies.
- The submission workbook was checked against the shared CSV sources across 39,499 header/data cells. Its original 38,839 occupied cells were preserved apart from four intended author/catalog metadata updates. The workbook has Guide and Tables S1-S16.

Run the public test suite with:

```bash
python -m unittest discover -s tests -v
python reproduce.py verify
python reproduce.py export --published --output ../paper3_publication_check
```

`verify --work RUN` and `export --work RUN` validate the source/output receipt chain and reject missing stages, unsealed outputs, changed files and a different framework identity. See `docs/ANALYSES.md` for scientific execution and `docs/INPUTS.md` for raw-source reconstruction or reuse of existing caches.

## Interpretation and execution scope

The new experiment supports a boundary between tested-gene denominator feedback and marker-count sensitivity. It does not establish universal FDR control for MuSiC or marker-restricted normalization. Its four observed reference donors are fixed; Monte Carlo standard errors describe repeated mixtures conditional on those references and target banks.

The original seven large scientific analyses were not refitted during this validation. Their outputs were reaggregated and checked, while the new cross-estimator analysis, input reconstruction and clinical sensitivity workflow were executed on the server. A single fresh invocation of the entire eight-stage framework from newly downloaded raw data was not performed. The code does not claim that an additional public archive of the prepared count banks exists; their verified reconstruction recipe and pinned public source manifest are provided instead.
