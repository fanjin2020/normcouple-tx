# Computation and output mapping

All commands use `reproduce.py`. Scientific fits run in `WORK`; the compact publication directory defaults to `WORK/final`.

| Stage | Scientific scope | Final numerical outputs |
| --- | --- | --- |
| normalization | Fixed/observed totals, TMM, DESeq2, pseudocount/scale grid | Table S2, S3 |
| mechanism | Paired effect/noise accounting and estimability | Table S12; power trajectories |
| marker_grid | Relative marker-count perturbations and unperturbed power | Table S4, S8; dose contrasts and donor power components |
| random_panels | Sealed paired controls for randomized panels | Inputs to the next two stages |
| coupling | Full confirmation parent, count-level clinical arrays | Inputs to common support and corrected clinical comparisons |
| common_support | Equal retained budgets, both supports, multistep clinical probes, derivative qualification | Table S5-S7, S9-S10, S13-S15; clinical and matching data |
| marker_associations | Donor-level marker/disease comparisons | Table S11; marker data |
| cross_method | Fixed-response dtangle/MuSiC comparison with raw counts and fixed reference donors | Table S16 |

Table S1 describes the fixed resources. Figure 1's curves are evaluated analytically and saved as a source CSV. `--stage common_support` automatically includes `random_panels` and `coupling`; `--stage all` includes all eight stages. Publication export requires complete receipts for all stages.

## Figure mapping

| Figure | Physical CSV source(s) |
| --- | --- |
| Fig1 | analytic_curves.csv |
| Fig2 | TableS2.csv |
| Fig3 | TableS3.csv |
| Fig4 | TableS4.csv, dose_contrasts.csv |
| Fig5 | TableS5.csv |
| Fig6 | common_support_donor_subsets.csv |
| Fig7 | TableS8.csv, power_trajectories.csv, power_component_donors.csv |
| Fig8 | TableS10.csv, clinical_coefficient_pairs.csv.gz |
| Fig9 | TableS13.csv |
| FigS1 | TableS4.csv, dose_contrasts.csv |
| FigS2 | TableS7.csv, viewed as pooled and class-specific rows |
| FigS3 | TableS6.csv |
| FigS4 | TableS9.csv |
| FigS5 | clinical_markers.csv |
| FigS6 | TableS12.csv |
| FigS7 | TableS14.csv, TableS15.csv |
| FigS8 | TableS16.csv |

Four other auxiliary CSVs retain evidence needed for the text and checks: `power_paired_contrasts.csv`, `matching_balance.csv`, `diagnostic_performance.csv`, and `clinical_linear_checks.csv`. These are distinct analyses, not extra copies of table data.

## Completion and resume

The wrapper records its source identity, exact input hashes and software versions. A stage completes only after a success state and required files are present; software fixtures cannot satisfy this condition. Later stages consume hash-verified controls and clinical arrays. Resume does not accept changed files as completed results. The kernels validate their own batch protocols before reusing work.

Aggregation reads the completed run, enforces paired keys and common gene families, and writes the paper's table schemas. A missing stage or source file stops export. `--published` is the explicit route for archived reconstruction.

Canonical kernels are byte-preserved. Their relative layouts are materialized under `WORK/.runtime` because source manifests and checkpoint signatures refer to those paths. Runtime layouts and checkpoints are working files, not extra public code packages.
