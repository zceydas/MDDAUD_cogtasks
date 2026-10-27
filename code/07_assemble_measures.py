"""
07_assemble_measures.py — combine the per-task scores into the files used by the R analyses.

Output
  derived/all_measures_long.csv   one row per participant x session, the 34 retest measures
                                  (PRL 7, rule switching 9, demand selection 12, MSIT 6)
  derived/demand_avoidance.csv    DemandAvoid and its splits, with the chance-level flag
  derived/effort_switch_repeat.csv learning-phase switch-trial and repeat-trial median RT (seconds),
                                  used for the arm effect on switch-trial RT

Exclusion applied here
  Demand selection: baseline sessions whose learning phase was not above chance
  (ES_chance_pre) are set to missing for every demand selection measure.
  Rule switching: sessions whose overall accuracy did not exceed chance (one-sided binomial
  p >= .05 over 240 trials; 3 baseline, 1 post-dose) are set to missing for every rule-switching
  measure, because rule-dependent indices are not interpretable when the rules were not followed.
  15_sensitivity_rule_switching.R repeats the main rule-switching results with these sessions kept.
"""
import os
import numpy as np, pandas as pd
from config import DERIVED, EXCLUDED_SUBJECTS

EXCLUDE_RS_CHANCE = True   # primary analysis (see docstring)

P = pd.read_csv(os.path.join(DERIVED, "prl_scores.csv"))
P = P[["SubjectID", "Session", "PRL_RewardAcc", "PRL_PunishAcc", "PRL_alphaGain", "PRL_alphaLoss",
       "PRL_betaGain", "PRL_betaLoss", "PRL_Perseveration"]]

C = pd.read_csv(os.path.join(DERIVED, "rule_switching_scores.csv"))
cf_cols = ["CF_Accuracy", "CF_SwitchCost_learning", "CF_SwitchCost_transfer", "CF_PreSwitchAcc_learning",
           "CF_PostSwitchAcc_learning", "CF_PreSwitchAcc_transfer", "CF_PostSwitchAcc_transfer",
           "CF_RT_learning", "CF_RT_transfer"]
if EXCLUDE_RS_CHANCE:
    C.loc[C.binom_p >= .05, cf_cols] = np.nan
C = C[["SubjectID", "Session"] + cf_cols]

E = pd.read_csv(os.path.join(DERIVED, "demand_selection_scores.csv"))
es_cols = ["ES_Learn_Acc", "ES_Learn_RT", "ES_Learn_SwitchCost_Acc", "ES_Learn_SwitchCost_RT", "ES_Dec_Acc",
           "ES_Dec_RT", "ES_Dec_SwitchCost_Acc", "ES_Dec_SwitchCost_RT", "ES_ChoiceProb_L1", "ES_ChoiceProb_L4",
           "ES_MeanChosenLevel", "ES_SelectionRT"]
chance_pre = (E.Session == "pre") & E.ES_chance_pre
E.loc[chance_pre, es_cols + ["ES_Learn_SwitchRT", "ES_Learn_RepeatRT"]] = np.nan

DA = E[["SubjectID", "Session", "DemandAvoid", "Avoid_d1", "Avoid_d2", "Avoid_d3", "n_choices", "ES_chance_pre"]].copy()
DA["chance_pre"] = DA.ES_chance_pre
DA = DA[DA.DemandAvoid.notna()]          # 20040165 saved no decision-phase selection file
# Raw values are kept here (chance-level baseline sessions included, flagged by chance_pre) because
# the avoidance analyses report results both with and without them; 13_demand_avoidance.R applies the flag.
DA.to_csv(os.path.join(DERIVED, "demand_avoidance.csv"), index=False)
SR = E[["SubjectID", "Session", "ES_Learn_SwitchRT", "ES_Learn_RepeatRT", "ES_chance_pre"]].copy()   # chance-level baseline already masked above
SR.columns = ["SubjectID", "Session", "SwitchRT", "RepeatRT", "ES_chance_pre"]
SR.to_csv(os.path.join(DERIVED, "effort_switch_repeat.csv"), index=False)

M = pd.read_csv(os.path.join(DERIVED, "msit_scores.csv"))
M = M[["SubjectID", "Session", "MSIT_ControlAcc", "MSIT_InterferenceAcc", "MSIT_ControlRT", "MSIT_InterferenceRT",
       "MSIT_InterferenceEffect", "MSIT_AccuracyDrop"]]

A = P
for X in (C, E[["SubjectID", "Session"] + es_cols], M):
    A = A.merge(X, on=["SubjectID", "Session"], how="outer")
A = A[~A.SubjectID.isin(EXCLUDED_SUBJECTS)].sort_values(["SubjectID", "Session"])
# Round to 10 decimals so that scores that are equal in exact arithmetic (e.g. 44/48 reached by two different
# computations) are stored as equal and receive tied ranks in the Spearman analyses.
num = [c for c in A.columns if c not in ("SubjectID", "Session")]
A[num] = A[num].round(10)
A.to_csv(os.path.join(DERIVED, "all_measures_long.csv"), index=False)
print("all_measures_long:", A.shape, A.groupby("Session").SubjectID.count().to_dict())
