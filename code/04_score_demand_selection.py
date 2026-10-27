"""
04_score_demand_selection.py — Parametric Demand Selection Task (Sayalı & Badre, 2019).

Task
  Cued parity/magnitude judgments; demand = probability of a task switch within an 8-trial block,
  one card deck per level (level 1-4 = .13, .38, .63, .88).
  Learning phase  (Learning_Session<N>_*.xlsx, 320 trials): every level performed.
  Decision phase  (Decision_Selection_*.xlsx: 72 free choices between two decks, no incentive;
                   Decision_Execution_*.xlsx: performance of the chosen deck).

Output 1: derived/demand_selection_scores.csv (one row per participant x session)
  ES_Learn_* / ES_Dec_*   Acc (proportion correct; missed responses count as errors),
                          RT (mean over the four levels of the per-level median RT, missed and
                          non-task key presses excluded),
                          SwitchCost_Acc (repeat minus switch accuracy), SwitchCost_RT (switch minus
                          repeat median RT). Switch/repeat is defined within a block from trial 2 on.
  ES_Learn_SwitchRT, ES_Learn_RepeatRT   median RT on switch and repeat trials (learning phase)
  ES_Learn_RT_L1..L4, ES_Learn_Acc_L1..L4  median RT and accuracy at each demand level (L1 vs L4: task validity, Table 3;
                          all four: per-level models in 17_supplementary_analyses.R)
  ES_ChoiceProb_L1..L4    P(deck chosen | deck offered); ES_MeanChosenLevel; ES_SelectionRT (median)
  DemandAvoid             proportion of choices of the LOWER-demand deck in each pair (0.5 =
                          indifference); Avoid_d1..d3 by difference in level between the two decks
  ES_Learn_binomP         one-sided binomial test of learning-phase accuracy against 0.5
  ES_chance_pre           True if the BASELINE learning phase was not above chance (p >= .05).
                          Deck-demand associations cannot have been learned in such sessions, so
                          all demand selection measures from those baseline sessions are set to
                          missing in 07_assemble_measures.py.

Data-quality handling
  - CorrectAnswer == 999 marks a missed execution trial (counted as an error, RT excluded).
  - Response == 9 with a normal CorrectAnswer (39 trials across all files, RT 0.01-1.19 s, Accuracy 0)
    is a press of a non-task key: counted as an error, RT excluded, because it is not the latency of a
    task response.
  - ChosenLevel == 999 marks a missed decision; the task re-presents the trial, so misses are dropped.
  - Two post-dose sessions had their decision-phase key mapping reversed (config.ES_RESPONSE_INVERSIONS);
    accuracy is flipped for answered trials in those sessions.
  - Duplicate files ('... (1).xlsx', '... - Copy', 'part2') are ignored.
"""
import os, glob
import numpy as np, pandas as pd
from scipy.stats import binomtest
from config import RAW, DERIVED, EXCLUDED_SUBJECTS, ES_RESPONSE_INVERSIONS, session_of


def switch_repeat(T):
    acc = (T.Accuracy == 1).values
    rt = T.RT.where((T.CorrectAnswer != 999) & (T.Response != 9)).values   # missed or non-task key: no valid RT
    is_switch = np.zeros(len(T), bool); defined = np.zeros(len(T), bool)
    for b in T.BlockNo.unique():
        idx = np.where(T.BlockNo.values == b)[0]
        tt = T.TrialType.values[idx]
        for k in range(1, len(idx)):            # trial 1 of a block has no predecessor
            defined[idx[k]] = True
            is_switch[idx[k]] = tt[k] != tt[k - 1]
    sa, ra = acc[defined & is_switch].mean(), acc[defined & ~is_switch].mean()
    sr, rr = np.nanmedian(rt[defined & is_switch]), np.nanmedian(rt[defined & ~is_switch])
    lv = [np.nanmedian(rt[T.EffortLevel.values == l]) for l in range(1, 5)]
    la = [acc[T.EffortLevel.values == l].mean() for l in range(1, 5)]
    return dict(Acc=acc.mean(), RT=np.mean(lv), SwitchCost_Acc=ra - sa, SwitchCost_RT=sr - rr,
                SwitchRT=sr, RepeatRT=rr, **{f"RT_L{l}": lv[l - 1] for l in range(1, 5)}, **{f"Acc_L{l}": la[l - 1] for l in range(1, 5)})


rows = []
for d in sorted(glob.glob(os.path.join(RAW, "sub-*", "*"))):
    sess = session_of(d)
    if sess is None or not os.path.isdir(d):
        continue
    sid = int(d.split("sub-")[1].split(os.sep)[0])
    if sid in EXCLUDED_SUBJECTS:
        continue
    n = 1 if sess == "pre" else 2

    def pick(pat):
        return [x for x in glob.glob(os.path.join(d, "**", pat), recursive=True)
                if "(1)" not in os.path.basename(x) and "Copy" not in x and "part2" not in x]

    L = pick(f"Learning_Session{n}_subject*.xlsx")
    E = pick(f"Decision_Execution_Session{n}_subject*.xlsx")
    S = pick(f"Decision_Selection_Session{n}_subject*.xlsx")
    r = dict(SubjectID=sid, Session=sess)
    if L:
        T = pd.read_excel(L[0])
        m = switch_repeat(T)
        r.update({"ES_Learn_" + k: v for k, v in m.items()})
        r["ES_Learn_binomP"] = binomtest(int((T.Accuracy == 1).sum()), len(T), 0.5, alternative="greater").pvalue
    if E:
        T = pd.read_excel(E[0])
        if (str(sid), sess) in ES_RESPONSE_INVERSIONS:
            answered = T.Response != 9
            T.loc[answered, "Accuracy"] = 1 - T.loc[answered, "Accuracy"]
        m = switch_repeat(T)
        r.update({"ES_Dec_" + k: m[k] for k in ("Acc", "RT", "SwitchCost_Acc", "SwitchCost_RT")})
    if S:
        T = pd.read_excel(S[0])
        T = T[T.ChosenLevel != 999]
        for l in range(1, 5):
            offered = (T.LeftSide_Effort == l) | (T.RightSide_Effort == l)
            r[f"ES_ChoiceProb_L{l}"] = ((T.ChosenLevel == l) & offered).sum() / offered.sum()
        r["ES_MeanChosenLevel"] = T.ChosenLevel.mean()
        r["ES_SelectionRT"] = T.RT.median()
        lo = np.minimum(T.LeftSide_Effort, T.RightSide_Effort)
        hi = np.maximum(T.LeftSide_Effort, T.RightSide_Effort)
        r["DemandAvoid"] = (T.ChosenLevel == lo).mean()
        for k in (1, 2, 3):
            sel = (hi - lo) == k
            r[f"Avoid_d{k}"] = (T.ChosenLevel[sel] == lo[sel]).mean()
        r["n_choices"] = len(T)
    if len(r) > 2:
        rows.append(r)

D = pd.DataFrame(rows).sort_values(["SubjectID", "Session"])
chance = set(D[(D.Session == "pre") & (D.ES_Learn_binomP >= .05)].SubjectID)
D["ES_chance_pre"] = D.SubjectID.isin(chance)
D.to_csv(os.path.join(DERIVED, "demand_selection_scores.csv"), index=False)
print(f"demand selection: {len(D)} sessions; baseline sessions at chance: {sorted(chance)}")
