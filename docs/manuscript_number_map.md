# Where each number in the manuscript comes from

Every number in the paper comes from a script below. Each one is printed to `OUT_DIR/results/<script>.txt` and, where it says so, saved to a CSV. Paths are relative to `OUT_DIR/results`. Per-participant scores behind every analysis are in `OUT_DIR/derived`. **Test statistics** (t or F with df) for every p value in the Results and the supplementary tables are in `16_test_statistics.txt`, `drug_effect_statistics.csv` and `supplement_statistics.csv`; each line prints the p value next to its statistic so the two can be checked together.

## Methods

| Statement | Script | Output |
|---|---|---|
| 25 participants at baseline (13/12); participants with post-dose data | 11 | `11_sample_and_clinical.txt` |
| Post-dose testing mean 6.7 days (SD 1.1) after dosing | 01, 20 | `20_task_validity.txt`, "Interval" |
| Baseline testing a median of 1 day before dosing (1–8); post-dose testing a median of 7 days after (4–9); sessions 8 days apart (5–13) | 01, 20 | `derived/session_dates.csv`; `20_task_validity.txt`, "Interval" (also checks the testing log against task-file dates) |
| Days from dosing as a covariate (sensitivity) | 10, 16 | `10_drug_effects_and_coupling.txt`, "Sensitivity"; `16_test_statistics.txt` |
| Drinking arm × window model (three windows) | 11 | `11_sample_and_clinical.txt` |
| PRL: 112 responded test trials; ΔAIC and BIC for perseveration | 02, 20 | `20_task_validity.txt`, PRL |
| Rule switching: feedback valid on 80% of trials | 03, 20 | `20_task_validity.txt` (feedback validity) |
| Data quality: 6 demand selection baseline sessions at chance (accuracy .48–.54); 3 baseline and 1 post-dose rule-switching sessions at chance (.45–.53); 2 inverted sessions; 2 MSIT position responders; MSIT mapping threshold 70% | 03, 04, 05, 20 | `20_task_validity.txt`; `derived/demand_selection_scores.csv`; `derived/msit_scores.csv` |
| Construct validation: 23 cross-task pairs, 12 complete cases | 21 | `21_construct_validation.txt` |

## Results

| Section / statement | Script | Output |
|---|---|---|
| Table 2; Methods "not statistically different" (all p≥.18, Welch t / Fisher's exact) | 11 | `table2_baseline_characteristics.csv` (column p), `11_sample_and_clinical.txt` |
| Task validity paragraphs, Table 3 | 20 | `20_task_validity.txt`, `table3_task_validity.csv` |
| Construct validation paragraph, Figure S1, Table S6 | 21 | `21_construct_validation.txt`, `tableS6_*.csv`, `figures/FigureS1_construct_validation.png` |
| Drinking fell (.62 → .30, .39 → .12); arm × window; PHQ-9 means by arm and ANCOVA; PHQ × drinking | 11 (+14 for BFs, 16 for t) | `11_sample_and_clinical.txt`, `14_bayes_factors.txt`, `16_test_statistics.txt` |
| Drug effects across the battery (smallest q; MSIT and value-learning d; three effects p<.05; raw means; baseline–change r) | 10 | `10_drug_effects_and_coupling.txt`, `all_measures_results.csv` |
| BF01 range for null effects; BF10 for each effect | 14 | `14_bayes_factors.txt`, `bf_arm_effects.csv` |
| Feedback-driven vs cued shifting: 7 of 11, sign test; within-arm t tests; composites (incl. baseline-adjusted composite); switch- and repeat-trial RT | 10 | `10_drug_effects_and_coupling.txt`, `switch_repeat_effects.csv`, `stability_composite.csv`, `stability_broad_composite.csv` |
| Per-level mixed models cited in the cued-switching paragraph (Table S9) | 17 | `17_supplementary_analyses.txt`, `tableS9_demand_level_models.csv` |
| Willingness to exert control (all numbers) | 13 (+14 for BFs) | `13_demand_avoidance.txt`, `tableS5_demand_avoidance.csv` |
| Change toward stability and drinking: smallest q, 15 of 17, couplings, composite (and baseline-adjusted composite), post-dose perseveration × heavy drinking | 10 (+14 for BFs, 16 for t) | `10_drug_effects_and_coupling.txt`, `16_test_statistics.txt` |
| Working-memory span predicts drinking reduction (70 tests, q, robustness, arm-specific r_s) | 12 (+14 for BFs, 16 for t) | `12_prediction_and_wmc.txt`, `prognostic_tests.csv`, `tableS4_wmc_robustness.csv` |
| Baseline control profile and drinking (105 tests; perseveration and MSIT associations) | 12 | `12_prediction_and_wmc.txt`, `tableS3_baseline_drinking.csv` |
| Depressive symptoms and cognitive change: PHQ couplings (pooled and arm-adjusted, smallest q), 11 of 17 sign test; 35 baseline predictors of PHQ change (smallest q); operation span → PHQ change (t, BF01) | 10, 12, 14, 16 | `10_drug_effects_and_coupling.txt`, "Coupling with PHQ change"; `12_prediction_and_wmc.txt`, `prognostic_tests_phq.csv`; `14_bayes_factors.txt`; `16_test_statistics.txt` |

## Tables and figures

| Item | Script | Output |
|---|---|---|
| Table 1 | — | Classification (`meta` in `config.R`) |
| Table 2 | 11 | `table2_baseline_characteristics.csv` (arm comparison in column p) |
| Table 3 | 20 | `table3_task_validity.csv` |
| Table S1 | 10, 14, 16 | `all_measures_results.csv` (columns nP–rbc, q_arm), `bf_arm_effects.csv`, `supplement_statistics.csv` (anc_t, pairP_t, pairS_t) |
| Table S2 | 10, 16 | `all_measures_results.csv` (cpl_*, cplH_*, phq_*, q_cpl), `supplement_statistics.csv` (cpl_adj_t, phq_t) |
| Table S3 | 10, 12 | `tableS3_baseline_drinking.csv` |
| Table S4 | 12, 16 | `tableS4_wmc_robustness.csv`; t statistics in `16_test_statistics.txt`, "Prediction" |
| Table S5 | 13, 14, 16 | `tableS5_demand_avoidance.csv`, `tableS7_bayes_factors.csv`; t statistics in `16_test_statistics.txt`, "Demand avoidance" |
| Table S6 | 21 | `tableS6_construct_pairs.csv`, `tableS6_reliability.csv` |
| Table S7 | 10, 11, 13, 14 | `tableS7_bayes_factors.csv`; Cohen's d for the composites in `10_...txt`; ANCOVA d for drinking and PHQ in `11_...txt`; Bonett–Wright CIs in `14_...txt` |
| Table S8 | 17 | `tableS8_ospan_postdose.csv`, `17_supplementary_analyses.txt` |
| Table S9 | 17 | `tableS9_demand_level_models.csv`, `17_supplementary_analyses.txt` |
| Figure 1 | 30 (export only) | schematic drawn in `assets/Figure1_framework_and_findings.svg`; exported to `figures/Figure1_...png/.pdf` |
| Figures 2–7 | 30 | `figures/Figure2_...` to `figures/Figure7_...`, each as .pdf (vector), .tiff and .png (600 dpi), drawn at 7.0–7.2 in wide; `figures/Figure3_data.csv` |
| Figure S1 | 21 | `figures/FigureS1_construct_validation.png` |

## Sensitivity analysis

| Statement | Script | Output |
|---|---|---|
| Rule-switching results with the 4 sessions at chance kept (8 of 11 feedback-driven measures; Table S1 note) and with Wen et al.'s 65% criterion (9 sessions excluded; Results, Supplementary Methods) | 15 | `15_sensitivity_rule_switching.txt` |
| Discussion: change-score correlations between tasks (rule-switching cost × cued switch cost, × MSIT interference) | 10 | `10_drug_effects_and_coupling.txt`, last section |
