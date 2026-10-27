"""
05_score_msit.py — Multi-Source Interference Task (Bush & Shin, 2006).

Task
  Three digits are shown; participants report the identity of the digit that differs from the
  other two (the target). Control blocks: target matches its position, flankers are 0.
  Interference blocks: target position and identity conflict and flankers are other digits.
  Eight alternating blocks, 192 trials (PsychoPy output, one csv per session).

Key mapping
  The response keys were recorded as typed key codes, and the mapping between key code and target
  digit is not stored in the file. The mapping is therefore inferred per block by testing which rule
  reproduces >= 70% of responses:
    VALUE_NOSHIFT  key == target digit
    VALUE_SHIFT    key == target digit + 1
    POSITION_SHIFT key == target position + 2   (participant answered by position, not identity)
  Blocks below 70% for every rule are unresolved and their trials are not scored.

Threshold
  The threshold only has to establish WHICH mapping was in use, not that the participant did well;
  a stricter cutoff discards blocks because the participant made errors, which biases accuracy
  upward. Across all blocks, the best rule that predicts different keys from the winning rule never
  matches more than 50% of responses, so 70% identifies the mapping unambiguously (a stricter 80%
  would leave two blocks with a clear winner, 78% and 75%, unscored because of genuine errors).
  One block (59% vs 41%, a mid-session change of key mapping) remains unresolved.

Position responding
  Answering interference trials by the target's position means the prepotent response was not
  overridden, so these answers are not scored as correct: sessions in which an interference block is
  best explained by position responding are excluded (flag position_responder). In this dataset these
  are two post-dose sessions with 98-99% of interference trials answered by position
  (config.MSIT_POSITION_RESPONDERS). Files are looked up in '<session>/MSIT/' and, failing that, in
  the session folder itself.

Output: derived/msit_scores.csv (one row per participant x session)
  MSIT_ControlAcc, MSIT_InterferenceAcc        proportion correct (omissions count as errors)
  MSIT_ControlRT, MSIT_InterferenceRT          mean RT of correct trials (ms)
  MSIT_InterferenceEffect                      interference minus control RT (ms)
  MSIT_AccuracyDrop                            control minus interference accuracy
"""
import os, glob, re
import numpy as np, pandas as pd
from config import RAW, DERIVED, EXCLUDED_SUBJECTS, MSIT_POSITION_RESPONDERS, session_of

THRESH = 0.70   # see 'Threshold' in the docstring


def find_file(sess_dir):
    fs = sorted(glob.glob(os.path.join(sess_dir, "MSIT", "*.csv")))
    if not fs:
        fs = sorted(glob.glob(os.path.join(sess_dir, "*MSIT*.csv")))
    return fs[0] if fs else None


def score(f):
    T = pd.read_csv(f, encoding="utf-8-sig")
    T = T[T["Trials.thisIndex"].notna()].reset_index(drop=True)
    stim = T.number_stim.astype(int).map(lambda x: f"{x:03d}")
    val = np.array([[int(c) for c in s if s.count(c) == 1][0] for s in stim])            # target identity
    pos = np.array([[i for i, c in enumerate(s) if s.count(c) == 1][0] for s in stim])   # target position (0-2)
    cond = np.where(T.conds_file == "interference_stims.csv", "Interference", "Control")
    key = pd.to_numeric(T["Choice_resp.keys"], errors="coerce").values
    rt = pd.to_numeric(T["Choice_resp.rt"], errors="coerce").values
    block = T["Blocks.thisN"].values

    correct_key = np.full(len(T), np.nan)
    position_block = False
    for b in np.unique(block[~pd.isna(block)]):
        idx = block == b
        valid = idx & ~np.isnan(key)
        if valid.sum() == 0:
            continue
        rates = {"VALUE_NOSHIFT": np.mean(val[valid] == key[valid]),
                 "VALUE_SHIFT": np.mean(val[valid] + 1 == key[valid]),
                 "POSITION_SHIFT": np.mean(pos[valid] + 2 == key[valid])}
        best = max(rates, key=rates.get)          # ties resolved in the order above
        if rates[best] < THRESH:
            continue                               # unresolved block: not scored
        if best == "POSITION_SHIFT":
            position_block = True
        correct_key[idx] = {"VALUE_NOSHIFT": val, "VALUE_SHIFT": val + 1, "POSITION_SHIFT": pos + 2}[best][idx]

    scored = ~np.isnan(correct_key)
    corr = np.where(np.isnan(key), 0.0, (key == correct_key).astype(float))   # omission = error
    out = {}
    for c, lab in (("Control", "Control"), ("Interference", "Interference")):
        m = scored & (cond == c)
        out[f"MSIT_{lab}Acc"] = corr[m].mean()
        out[f"MSIT_{lab}RT"] = np.nanmean(rt[m & (corr == 1)]) * 1000
    out["MSIT_InterferenceEffect"] = out["MSIT_InterferenceRT"] - out["MSIT_ControlRT"]
    out["MSIT_AccuracyDrop"] = out["MSIT_ControlAcc"] - out["MSIT_InterferenceAcc"]
    out["n_unscored"] = int((~scored).sum())
    out["position_responder"] = position_block
    # proportion of interference trials answered by position (diagnostic)
    vi = (cond == "Interference") & ~np.isnan(key)
    out["interference_position_rate"] = float(np.mean(key[vi] == pos[vi] + 2))
    return out


rows = []
for sess_dir in sorted(glob.glob(os.path.join(RAW, "sub-*", "*"))):
    sess = session_of(sess_dir)
    if sess is None or not os.path.isdir(sess_dir):
        continue
    sid = int(re.search(r"sub-(\d+)", sess_dir).group(1))
    if sid in EXCLUDED_SUBJECTS:
        continue
    f = find_file(sess_dir)
    if f:
        rows.append(dict(SubjectID=sid, Session=sess, source_file=os.path.relpath(f, RAW), **score(f)))

D = pd.DataFrame(rows).sort_values(["SubjectID", "Session"])
flagged = set(zip(D.SubjectID[D.position_responder], D.Session[D.position_responder]))
assert flagged == MSIT_POSITION_RESPONDERS, f"position responders changed: {flagged}"
D.loc[D.position_responder, [c for c in D.columns if c.startswith("MSIT_")]] = np.nan
D.to_csv(os.path.join(DERIVED, "msit_scores.csv"), index=False)
print(f"MSIT: {len(D)} sessions; excluded (position responding): {sorted(flagged)}")
