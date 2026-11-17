# Validation scope

Version 2.1.0 was checked on November 17, 2026. Scientific calculations ran on Linux with the pinned Python environment and R packages. Document and workbook checks ran on Windows.

## Scientific checks

- The cross-estimator experiment completed all 210 batches: 1,680 datasets, 15,120 method-probe records and 10,080 paired contrasts. It uses 42 reference-target units and 20 repetitions per unit and effect. Its receipt records input, source and environment identities.
- Cached MuSiC reference calculations were compared with 192 uncached official API calls. The maximum prediction difference was zero. MuSiC 1.0.0 is pinned to commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c`.
- Independent aggregation reproduced all cross-estimator means and conditional Monte Carlo standard errors within 2e-15. The coefficient-weight identity had maximum absolute error 1.34e-15. All planned regressions were estimable.
- A final-launcher pilot reproduced 144 method-probe records and 96 contrast records exactly. A subsequent resume reused all four verified batches without recomputation.
- The user-data sensitivity command reproduced 161,610 gene-model results for 288 clinical samples and 16,161 genes. Availability and BH decisions agreed exactly. Maximum absolute differences were 8.52e-14 for coefficients, 1.71e-15 for standard errors, 1.63e-12 for p values and 1.90e-12 for q values.
- `rebuild-inputs` reconstructed both sparse count matrices, four cell/feature annotation files and the gene-identity map from cached author-distributed public source objects. All seven output SHA256 hashes matched the prepared-input definitions exactly. This was a reconstruction test, not a repeated download of all eleven source objects.

## Framework and publication checks

- All **65 framework/publication tests** passed on Linux. Coverage includes transformations, actual-span HC3, fraction closure, sample and family alignment, cross-estimator contracts, dependency ordering, input reconstruction, path safety, receipt verification and interrupted-run recovery.
- The preceding release's 17 common-support kernel tests also passed in its materialized runtime. Those frozen kernel files remain unchanged.
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
