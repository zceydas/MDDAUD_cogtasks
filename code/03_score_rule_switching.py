"""
03_score_rule_switching.py — probabilistic rule switching (cognitive flexibility task).

Task (low-volatility condition of Wen, Geddert, Madlon-Kay & Egner, 2023)
  Participants sort images by one of two stimulus dimensions. Feedback is valid on 80% of trials
  (column reward_validity). The relevant dimension changes without warning: every 30 trials in the
  120-trial learning phase (3 changes) and every 20 trials in the 120-trial transfer phase, which
  uses two new dimensions (5 changes). Practice trials (practice == 'false' marks real trials) are
  removed. Trials without a response (response_time == 0) count as incorrect and have no RT.

Input : data_subj<ID>_session <N>_..._CognitiveFlexibility_<date>.csv
Output: derived/rule_switching_scores.csv, one row per participant x session:
  CF_Accuracy                         overall accuracy (240 trials)
  Acc_learning, Acc_transfer          accuracy per phase (task validity, Table 3)
  CF_PreSwitchAcc_<phase>             accuracy on the 5 trials before each rule change
  CF_PostSwitchAcc_<phase>            accuracy on the 5 trials after each rule change
  CF_SwitchCost_<phase>               pre-switch minus post-switch accuracy (averaged over changes)
  CF_RT_<phase>                       median RT (ms)
  binom_p                             one-sided binomial test of overall accuracy against 0.5
  feedback_validity                   proportion of trials with valid feedback (0.80 in every session)
  (phase = learning: trials 1-120; transfer: trials 121-240; changes are counted within a phase)

Rule changes
  Every change of the relevant rule is located directly as a change in 'sequence' within a phase, so
  that all 3 (learning) and 5 (transfer) changes enter the switch measures. (The rule index
  alternates, e.g. 0,1,0,1, so the changes cannot be found from the unique rule values.)
"""
import os, re, glob
import numpy as np, pandas as pd
from scipy.stats import binomtest
from config import RAW, DERIVED, EXCLUDED_SUBJECTS

K = 5  # trials before / after each rule change


def find_file(sub_dir, sid, sess):
    """Return the task file for one session. Files are normally in '<session>/Cognitive Flexibility/',
    but some were saved elsewhere in the session folder (20040130 baseline: the session folder itself;
    20040141 baseline: inside the PRL data folder), so the whole session folder is searched, with the
    'Cognitive Flexibility' folder preferred. Duplicates saved by the browser as '... (1).csv' are
    ignored when an unsuffixed file exists."""
    for d in glob.glob(os.path.join(sub_dir, sess + "*")):
        if os.path.basename(d) not in (sess, sess + " dose"):
            continue
        cands = [f for f in glob.glob(os.path.join(d, "**", "*.csv"), recursive=True)
                 if re.search(rf"data_subj{sid}_session.*CognitiveFlexibility", os.path.basename(f))
                 or re.search(rf"{sid}_.*cog_flex", os.path.basename(f))]
        if not cands:
            continue
        pref = [f for f in cands if "Cognitive Flexibility" in f] or cands
        plain = [f for f in pref if not re.search(r"\(\d+\)", os.path.basename(f))]
        return sorted(plain if len(plain) == 1 else pref)[0]
    return None


def switch_metrics(acc, seq, k=K):
    """Mean accuracy in the k trials before vs after every change of the relevant rule."""
    bnd = [i for i in range(1, len(seq)) if seq[i] != seq[i - 1]]
    pre = np.mean([acc[i - k:i].mean() for i in bnd])
    post = np.mean([acc[i:i + k].mean() for i in bnd])
    return pre - post, pre, post


rows = []
for sub_dir in sorted(glob.glob(os.path.join(RAW, "sub-*"))):
    sid = re.search(r"sub-(\d+)", sub_dir).group(1)
    if int(sid) in EXCLUDED_SUBJECTS:
        continue
    for sess in ("pre", "post"):
        f = find_file(sub_dir, sid, sess)
        if f is None:
            continue
        T = pd.read_csv(f)
        X = T[T.practice.astype(str).str.lower() == "false"].reset_index(drop=True)
        acc = X.response_acc.astype(str).str.lower().eq("true").values.astype(float)
        seq = X.sequence.values
        n = len(X); half = n // 2
        rt = pd.to_numeric(X.response_time, errors="coerce").replace(0, np.nan)
        scl, prl, pol = switch_metrics(acc[:half], seq[:half])
        sct, prt, pot = switch_metrics(acc[half:], seq[half:])
        k = int(acc.sum())
        rows.append(dict(SubjectID=int(sid), Session=sess, source_file=os.path.relpath(f, RAW), n_trials=n,
                         feedback_validity=float((X.reward_validity == 1).mean()),
                         CF_Accuracy=acc.mean(), Acc_learning=acc[:half].mean(), Acc_transfer=acc[half:].mean(),
                         CF_SwitchCost_learning=scl, CF_SwitchCost_transfer=sct,
                         CF_PreSwitchAcc_learning=prl, CF_PostSwitchAcc_learning=pol,
                         CF_PreSwitchAcc_transfer=prt, CF_PostSwitchAcc_transfer=pot,
                         CF_RT_learning=np.nanmedian(rt[:half]), CF_RT_transfer=np.nanmedian(rt[half:]),
                         binom_p=binomtest(k, n, 0.5, alternative="greater").pvalue))

D = pd.DataFrame(rows).sort_values(["SubjectID", "Session"])
D.to_csv(os.path.join(DERIVED, "rule_switching_scores.csv"), index=False)
print(f"rule switching: {len(D)} sessions; at chance (binomial p >= .05): {int((D.binom_p >= .05).sum())}")
