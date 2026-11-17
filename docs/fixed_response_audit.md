# Fixed-response composition audit

The `sensitivity` command compares association results across composition estimates while holding the sample set, outcomes, exposure, other covariates, and BH family fixed. It uses the existing OLS/HC3 implementation in the manuscript framework.

```bash
python reproduce.py sensitivity \
  --outcomes outcomes.csv \
  --design design.csv \
  --fractions fractions.csv \
  --baseline tmm \
  --exposure diagnosis \
  --covariates age sex \
  --output audit_results
```

Input tables:

- `outcomes.csv`: one `sample_id` column and one numeric column per gene, already on the intended regression scale. The command does not normalize or transform outcomes.
- `design.csv`: `sample_id`, the numeric exposure, and numeric covariates. Encode categorical covariates before invoking the command.
- `fractions.csv`: `sample_id`, `method`, and identically named cell-class columns for every method. Fractions must be finite and nonnegative. Row sums within 1e-6 of one are accepted and normalized to exactly one before fitting, so decimal rounding cannot add a spurious composition dimension. Each method must contain exactly the same sample IDs as the outcome table.

Outputs are `method_summary.csv`, `design_comparison.csv`, `gene_statistics.csv.gz`, and one small provenance record. Gene-level q values use the full supplied family, with unavailable tests assigned p=1. The actual nuisance column space determines estimability. Reports include rank, exposure information, availability, coefficient changes, standard-error ratios, weight-vector distance, and BH status changes. These are conditional sensitivity summaries; the command does not estimate clinical FDR or certify a preferred normalizer.

To isolate normalization feedback, generate alternative fraction matrices by modifying tested counts only in the input used for composition estimation. Keep the outcome table unchanged. Compare factor-only exclusion with complete exclusion before attributing a change to library totals. Marker perturbations address a separate source of sensitivity.

## Official second-estimator benchmark

```bash
python reproduce.py cross-method \
  --inputs /path/to/verified_inputs \
  --work /path/to/new_cross_method_run \
  --rscript /path/to/Rscript \
  --r-library /path/to/existing_R_library \
  --workers 4 --repetitions 20 --batch-size 4
```

Use `--smoke --repetitions 1` in a separate directory first. For an interrupted run, repeat the identical command with `--resume`. Only receipt-verified batches are reused. Use one batch size throughout a run.

The benchmark uses raw counts for the official MuSiC `music_prop` implementation at commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c`, with a fixed disjoint marker panel and native cell-fraction output. MuSiC uses weighted NNLS; it is a different reference-weighting and input-processing workflow, not a non-NNLS estimator. Four observed reference donors are held fixed and equally weighted. The two internal and two external settings use the original target donors, frozen marker/test definitions, and empirical pseudobulk generator. Effects are zero and 0.5. The same outcome is used under baseline, tested-count, and marker-count probes.

New benchmark results are stored under `WORK/final`: one compact summary and two compressed replicate tables. Monte Carlo uncertainty is conditional on the observed donor/reference collection; between-donor ranges are separate. The main `export --published` workflow replays the archived manuscript tables and does not silently replace them with a new benchmark run.

The official fitting function is unchanged. Reference moments from `music_basis` are cached once per fixed reference to avoid repeating the same calculation. The first batch in each setting is checked against the uncached official API on separately fitted samples; the check includes all three probes and both outcome effects. This also checks that fitting the samples together does not alter the selected-sample results.
