"""
02_score_prl.py — Probabilistic reward learning (PRL): choice accuracy and reinforcement-learning fits.

Task (Pessiglione et al., 2006; design as in Tichelaar et al., 2023)
  Each run presents one new cue pair on Gain trials (win 10 or nothing) and one on Loss trials
  (lose 10 or nothing), interleaved; in each pair one cue gives the better outcome on 75% of trials.
  Run 1 is training and is not analysed. Runs 2-3 are the test runs (56 trials each).

Input : raw Results{1,2}_subject<ID>_<date>.xlsx files (one per session)
Output: derived/prl_scores.csv, one row per participant x session:
  PRL_RewardAcc, PRL_PunishAcc   proportion of choices of the better cue (ChosenStimAccuracy) on
                                 Gain / Loss trials of the test runs
  PRL_alphaGain, PRL_alphaLoss   Rescorla-Wagner learning rates
  PRL_betaGain,  PRL_betaLoss    softmax inverse temperatures
  PRL_Perseveration              choice-repetition ("stickiness") parameter
  nll_stick, nll_nostick, n_trials  for the model comparison reported in the Methods

Model
  V(chosen) <- V(chosen) + alpha * (R - V(chosen)), R = Earning / 10 (Gain: 0/+1, Loss: -1/0)
  P(chosen) = 1 / (1 + exp(-(beta * (V(chosen) - V(other)) + s))), where s = +c if the same cue
  was chosen on the previous trial of that cue pair, -c otherwise (s = 0 on the first trial).
  Values start at +0.5 (Gain pairs) or -0.5 (Loss pairs) and are reset for each run x trial type.
  Bounds: alpha in [0,1], beta in [0,10], c in [-3,3]. Maximum likelihood, best of 50 random
  starts (L-BFGS-B; seed 1 for every fit; fewer starts did not always reach the maximum). In a few
  sessions the likelihood is flat in one direction (a learning rate at 0 leaves the matching inverse
  temperature undetermined), so learning-rate and inverse-temperature estimates for those sessions
  can differ slightly between NumPy/SciPy versions. Perseveration estimates are stable (to 0.003
  across start sets).

Scoring decisions
  1. Missed trials are coded 999 in RT and Accuracy (not left empty). They, and trials without a
     chosen cue, are removed before computing accuracy and before model fitting.
  2. The 'Accuracy' column records whether the trial's outcome was favourable, not whether the
     better cue was chosen. Learning accuracy is therefore taken from 'ChosenStimAccuracy'.
"""
import os, re, glob
import numpy as np, pandas as pd
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
from config import RAW, DERIVED, EXCLUDED_SUBJECTS

LB = [0, 0, 0, 0, -3]
UB = [1, 1, 10, 10, 3]


def find_files():
    """One Results file per participant x session. Results1 = baseline, Results2 = post-dose.
    Files are usually inside '<session>/Probabilistic Reward Learning/...', but for 20040130 (pre)
    the file sits directly in the session folder; when a session has more than one copy
    (20040185 pre: a duplicate in 'Working Memory'), the copy in the reward-learning folder is used."""
    found = {}
    for f in glob.glob(os.path.join(RAW, "sub-*", "**", "Results*_subject*.xlsx"), recursive=True):
        if os.path.basename(f).startswith("~$"):
            continue
        sid = int(re.search(r"sub-(\d+)", f).group(1))
        sess = {"1": "pre", "2": "post"}[re.search(r"Results(\d)_subject", os.path.basename(f)).group(1)]
        key = (sid, sess)
        in_task_folder = bool(re.search(r"eward", f, re.I))
        if key not in found or (in_task_folder and not re.search(r"eward", found[key], re.I)):
            found[key] = f
    return found


def blocks(T):
    """Test-run trials, missed trials removed, split into blocks = (run, trial type) in chronological order."""
    T = T[T.RunNo > 1]
    T = T[(T.RT != 999) & (T.Accuracy != 999) & T.ChosenStimulus.isin(["A", "B"])].copy()
    T["ch"] = np.where(T.ChosenStimulus == "B", 1, 0)
    T = T.sort_values(["RunNo", "TrialNo"])
    out = [(tt == "Gain", g.Earning.values / 10.0, g.ch.values) for (r, tt), g in T.groupby(["RunNo", "TrialType"], sort=True)]
    return out, len(T)


def nll(x, B, stick=True):
    """Negative log-likelihood of the observed choices."""
    aG, aL, bG, bL = x[:4]
    c = x[4] if stick else 0
    ll = 0.0
    for gain, R, ch in B:
        a, b = (aG, bG) if gain else (aL, bL)
        ev = np.array([0.5, 0.5]) if gain else np.array([-0.5, -0.5])
        for t in range(len(R)):
            k = ch[t]
            s = 0 if t == 0 else (c if ch[t] == ch[t - 1] else -c)
            p = 1 / (1 + np.exp(-(b * (ev[k] - ev[1 - k]) + s)))
            p = min(max(p, 1e-10), 1 - 1e-10)
            ll += np.log(p)
            ev[k] += a * (R[t] - ev[k])
    return -ll


def fit(B, stick=True, n=50, seed=1):
    """Best of n random starts; the same seed for every fit so the starts are reproducible."""
    rng = np.random.default_rng(seed)
    lb, ub = (LB, UB) if stick else (LB[:4], UB[:4])
    best = None
    for _ in range(n):
        x0 = [.1 + rng.random() * .8, .1 + rng.random() * .8, 1 + rng.random() * 8, 1 + rng.random() * 8]
        if stick:
            x0.append(-2.5 + rng.random() * 5)
        r = minimize(nll, x0, args=(B, stick), bounds=list(zip(lb, ub)), method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    return best


def score(item):
    (sid, sess), f = item
    T = pd.read_excel(f)
    B, n = blocks(T)
    r = fit(B, stick=True)
    r0 = fit(B, stick=False)
    X = T[(T.RunNo > 1) & (T.RT != 999) & (T.Accuracy != 999) & T.ChosenStimulus.isin(["A", "B"])]   # trials as in blocks()
    return dict(SubjectID=sid, Session=sess, source_file=os.path.relpath(f, RAW),
                PRL_RewardAcc=(X[X.TrialType == "Gain"].ChosenStimAccuracy == 1).mean(),
                PRL_PunishAcc=(X[X.TrialType == "Loss"].ChosenStimAccuracy == 1).mean(),
                PRL_alphaGain=r.x[0], PRL_alphaLoss=r.x[1], PRL_betaGain=r.x[2], PRL_betaLoss=r.x[3],
                PRL_Perseveration=r.x[4], nll_stick=r.fun, nll_nostick=r0.fun, n_trials=n,
                n_missed=int(((T.RunNo > 1) & (T.RT == 999)).sum()))


if __name__ == "__main__":
    files = {k: v for k, v in find_files().items() if k[0] not in EXCLUDED_SUBJECTS}
    with ProcessPoolExecutor() as ex:
        rows = list(ex.map(score, sorted(files.items())))
    D = pd.DataFrame(rows).sort_values(["SubjectID", "Session"])
    D.to_csv(os.path.join(DERIVED, "prl_scores.csv"), index=False)
    print(f"PRL: {len(D)} sessions scored")
