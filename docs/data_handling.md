# Data handling: exclusions and corrections

Every participant- or session-level decision made in scoring is listed here, with the script that applies it. Each one is also described in the docstring of that script and in the Supplementary Methods of the article.

## Sample

- **Allocation table.** The source table had two rows with the same participant number and no row for another participant. The rows were assigned using self-reported age and sex in the task files and baseline heavy-drinking days in the Timeline Followback, and the assignment was confirmed against the trial record. Both rows are in the psilocybin arm. The trial record also corrected one participant's antidepressant use. The corrected rows are kept with the data (`allocation_corrections.csv`) and applied by `01_allocation.py`.
- **Participant from another study.** One ID from a different study series is removed (`config.py`, `EXCLUDED_SUBJECTS`).

## Task scoring

| Task | Decision | Script |
|---|---|---|
| Probabilistic reward learning | Missed responses (coded 999) and trials without a chosen cue are removed before computing accuracy and before model fitting. Learning accuracy uses the chosen-cue accuracy column; the task's `Accuracy` column records the outcome. | `02_score_prl.py` |
| Probabilistic rule switching | Pre- and post-switch accuracy are computed around every change of the relevant rule (3 in the learning phase, 5 in the transfer phase). Sessions whose overall accuracy does not exceed chance (one-sided binomial test, p ≥ .05) are excluded from all rule-switching measures; `15_sensitivity_rule_switching.R` repeats the results with these sessions kept and with a 65% accuracy criterion. | `03_score_rule_switching.py`, `07_assemble_measures.py` |
| Parametric Demand Selection Task | Missed responses and presses of non-task keys count as errors and are excluded from RT. Missed decision-phase choices, which the task re-presents, are dropped. Two post-dose sessions with a reversed decision-phase key mapping (responses inverted on 95–96% of answered trials) are scored with the mapping reversed (`config.py`, `ES_RESPONSE_INVERSIONS`). Baseline sessions whose learning-phase accuracy does not exceed chance are excluded from all measures of the task. | `04_score_demand_selection.py`, `07_assemble_measures.py` |
| MSIT | The key-to-digit mapping is inferred per block (rule matching ≥ 70% of responses); one block with a mid-session change of mapping is not scored. Sessions in which interference trials were answered by the target's position (98–99% of trials in two post-dose sessions) are excluded (`config.py`, `MSIT_POSITION_RESPONDERS`). | `05_score_msit.py` |
| Operation span | Partial-credit unit score from the trial-level files. One baseline score whose trial-level file was not saved is taken from the task's own summary log, which uses the same score; one baseline file is empty and has no score. | `06_score_ospan.py` |

## Files outside the expected folders

A few task files were saved in another task's folder or directly in the session folder. Every scoring script therefore searches the whole session folder and resolves duplicate copies explicitly (see each script's `find_file()`, `find_files()` or `pick()`).

## Numerical ties

All measures are rounded to 10 decimals (`07_assemble_measures.py`, and drinking proportions in `config.R`) so that scores that are equal in exact arithmetic receive tied ranks in the Spearman analyses.
