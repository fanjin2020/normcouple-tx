# Reading the compact results

1. `STATUS.json` and `analysis_scope.json`: completion and scope; pilots, fixtures
   and clinical-only runs are distinct from a complete confirmation.
2. `method_summary.csv`: equal-donor, equal-subset means. Simulation performance is
   conditional on a feasible intervention. Planned/feasible records are explicit.
   Rates are fractions, not percentages. Planned unavailable gene tests retain the
   frozen p=1 before BH rule within feasible interventions.
3. `paired_method_differences.csv`: each method minus TMM within paired replicates.
   `paired_controls.csv` adds direct MRN minus denominator-only and oracle contrasts,
   with per-subset Monte Carlo standard errors from paired repetitions.
4. `paired_subset_differences.csv`: subset versus predefined, matched versus random,
   and positive-only versus allow-zero comparisons. Differences use common feasible
   pairs; planned and common-feasible counts accompany the estimates.
5. `perturbation_support.csv`: infeasibility, relative amplitudes and the fraction
   of added counts placed on baseline zeros. Do not interpret subset identity
   without these operating characteristics and matching quality.
6. `reference_marker_features.csv`, `matching_pairs.csv`,
   `marker_matching_balance.csv`: actual reference-only selections and balance.
   Matching is approximate; a weak match is not relabeled as exact balance.
7. `donor_subset_results.csv`, `donor_results.csv`,
   `target_donor_leave_one_out.csv`: biological grouping and sensitivity. This is
   leave-one-target-donor aggregation, not re-training a reference without a donor.
   `within_unit_monte_carlo.csv` gives repetition-based means, standard errors and
   usable counts for each fixed donor/subset. Donor IQRs in `method_summary.csv`
   describe biological variation, not Monte Carlo uncertainty.
8. `approximation_summary.csv`, `coefficient_contrasts_summary.csv`,
   `invariance_summary.csv` and status tables: finite-step error, adjustment-space
   changes and controls. Numerical identities are not biological validation.
9. `diagnostic_performance.csv`, `diagnostic_heldout.csv.gz`: held-out target donors,
   false alerts, missed realized-FDP events and valid-score coverage. No FDR control
   guarantee or universally safe threshold follows from this score.
   `diagnostic_donor_performance.csv` and `diagnostic_setting_performance.csv`
   separate settings and donors. Their bootstrap intervals resample target donors
   with fixed held-out predictions; they do not refit training thresholds.
   `diagnostic_status_counts.csv` includes planned unavailable scores.
10. `clinical_fraction_scale.csv`: q versus p comparison with one outcome/family.
    Clinical association counts indicate stability, not known accuracy.
11. `clinical_method_stability.csv` and other `clinical_*` summaries: public bulk
    sensitivity, sample identities, q/p predictions and conversion yields.
12. `matched_fdp.pdf/png`, `matched_power.pdf/png`: plots made from this run's
    summaries only. Lines connect predefined/random/matched mean values, not
    trajectories or confidence intervals. They show allow-zero allocation;
    positive-only results and feasibility are in the numeric tables.

The review ZIP is the compact export, not the complete per-replicate archive.
Raw q/p association tables and batch records remain in the output directory.
The original manuscript's figures are not silently overwritten or re-used here.
