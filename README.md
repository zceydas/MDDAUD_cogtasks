# Psilocybin, cognitive stability and drinking reduction in AUD with depression: analysis code

Analysis code for:

> Sayalı, C., Weisman, E., & Barrett, F. S. *Psilocybin, cognitive stability and drinking reduction in co-occurring alcohol use disorder and depression: a pilot secondary analysis of a randomized, placebo-controlled trial.* medRxiv, 2026 (preprint; DOI to be added).

The scripts take the raw task files and produce every statistic, table and figure in the article and its supplement.

## Data

Participant data are not included in this repository. The data are available from the corresponding author upon reasonable request, subject to the conditions of the trial's ethics approval. The code expects the data folder laid out as follows:

```
DATA_DIR/
  III_Data/III_Data/sub-<ID>/<pre|pre dose|post|post dose>/...   raw task files, one folder per participant and session
  PRL_Analysis/dose_allocation_merged.csv                         arm allocation and screening variables
  PRL_Analysis/allocation_corrections.csv                         corrected allocation rows (see docs/data_handling.md)
  2004_tlfb_20260601.csv                                          Timeline Followback, one row per participant and window
  2015 Cog Testing Dates.xlsx                                     dates of baseline testing, drug administration and post-dose testing
  Integrative_Analysis/PHQ_scores.csv                             PHQ9 at the baseline and post-dose cognitive visits
```

## Running

```bash
pip install -r requirements.txt
Rscript install_R_packages.R
DATA_DIR=/path/to/data OUT_DIR=/path/to/outputs ./run_all.sh
```

`DATA_DIR` defaults to `./data` and `OUT_DIR` to `./outputs`; both are excluded from version control. The whole pipeline takes a few minutes. Results are written to `OUT_DIR/results` (CSV tables plus one text log per script listing every number quoted in the paper) and figures to `OUT_DIR/figures`. Model fits and permutation tests use fixed seeds, so repeated runs on the same machine give identical output. Across platforms, two things can differ slightly: the Bayes factors, which are estimated by Monte Carlo sampling (100,000 iterations; about 2–3%), and the PRL learning-rate and inverse-temperature estimates of a few sessions whose likelihood is flat in one direction (see `02_score_prl.py`), which can move with the NumPy/SciPy version. Perseveration and all other measures are unaffected.

## Pipeline

| Step | Script | What it does | Manuscript |
|---|---|---|---|
| Scoring | `01_allocation.py` | Arm allocation and screening variables, with the documented corrections; testing dates and intervals | Methods, Table 2 |
| | `02_score_prl.py` | Probabilistic reward learning: choice accuracy and Rescorla–Wagner fits (with and without perseveration) | Methods, Tables 1, 3, S1 |
| | `03_score_rule_switching.py` | Probabilistic rule switching: accuracy, pre- and post-switch accuracy and switch cost at every rule change; binomial test against chance | Tables 3, S1 |
| | `04_score_demand_selection.py` | Parametric Demand Selection Task: accuracy, RT, task-switch costs, choice probabilities, demand avoidance; learning-phase chance test | Tables 3, S1, S5 |
| | `05_score_msit.py` | MSIT: per-block key mapping, accuracy, RT, interference effects; excludes position responders | Tables 3, S1 |
| | `06_score_ospan.py` | Operation span: partial-credit score from trial-level files; recall by set size | Tables 2, 3, S4 |
| | `07_assemble_measures.py` | Combines the 34 retest measures; applies the demand selection chance exclusion | — |
| Statistics | `10_drug_effects_and_coupling.R` | ANCOVA drug effects, couplings with drinking and PHQ change, re-signing on the stability–flexibility continuum, sign tests, stability composites | Results; Tables S1, S2 |
| | `11_sample_and_clinical.R` | Baseline characteristics and their comparison between arms; drinking and PHQ outcomes | Table 2; Results |
| | `12_prediction_and_wmc.R` | 70 prognostic tests; operation span robustness checks; baseline cognition × baseline drinking | Tables S3, S4 |
| | `13_demand_avoidance.R` | Demand avoidance (willingness to exert control) | Table S5 |
| | `14_bayes_factors.R` | Bayes factors for all drug effects and main claims | Table S1 (BF10), Table S7 |
| | `15_sensitivity_rule_switching.R` | Rule-switching analyses with the sessions at chance kept, and with Wen et al.'s 65% accuracy criterion | Sensitivity analysis |
| | `16_test_statistics.R` | Test statistics (t, F with df) for every p value quoted in the Results and the supplementary tables | Results; Table 3; Tables S1, S2 |
| | `17_supplementary_analyses.R` | Post-dose operation span; per-level demand selection mixed models | Tables S8, S9 |
| Validity | `20_task_validity.py` | Baseline task effects, PRL model comparison, retest interval, data checks | Table 3; Methods |
| | `21_construct_validation.py` | Baseline construct validation (cross-task correlations, permutation test, test–retest) | Figure S1, Table S6 |
| Figures | `30_figures.R` | Figures 2–7 (PDF, TIFF and PNG at print size); export of Figure 1 | Figures 1–7 |

Figure 1 is a schematic: `assets/Figure1_framework_and_findings.svg`. Shared settings are in `code/config.py` and `code/config.R`, which also hold the session-level decisions (exclusions and corrections), each with a comment that says why it was made.

`docs/manuscript_number_map.md` lists, for each section of the paper, the script and output that produce each number. `docs/data_handling.md` lists every exclusion and correction made in scoring.

## License

MIT; see `LICENSE`.
