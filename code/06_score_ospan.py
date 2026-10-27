"""
06_score_ospan.py — automated operation span (Unsworth et al., 2005).

Score
  Partial-credit unit score: for each of the 15 test sets, the number of letters recalled in the
  correct serial position divided by set size, summed over sets (range 0-15). Processing (math)
  accuracy = percentage of the 75 test math problems answered correctly; timeouts count as errors.
  Practice trials (BLOCKNUM not a positive number) are excluded.

Input
  ospan<ID>session<N>.csv trial-level files (RECALLED = 1 for a letter recalled in position,
  0 otherwise, -1 for non-recall rows; MATHACC, TIMEOUT per math problem). Files are found anywhere
  inside the session folder, because two baseline files were saved inside the PRL data folder.
  20040185 (baseline): no trial-level file was saved. Its score is taken from the task's own
  summary log (Working Memory/Results.csv), which uses the same partial-credit score; the log row
  for 20040185 is unambiguous. The log stores processing accuracy rounded to whole percent.
  20040180 (baseline): the trial-level file is empty and the log does not contain this participant;
  no score.

Output
  derived/ospan_scores.csv           one row per session with a trial-level file
  derived/ospan_baseline.csv         SubjectID, OSpan, OSpan_math (baseline; used in all analyses)
  derived/ospan_setsize.csv          baseline recall proportion by set size (Table 3)
"""
import os, re, glob
import numpy as np, pandas as pd
from config import RAW, DERIVED, EXCLUDED_SUBJECTS, session_of

LOG_ROW = re.compile(r"(\d{7,9}),(\d+\.\d{4}),(\d+\.\d{4}),(\d+\.\d{4})")


def score_file(f):
    T = pd.read_csv(f)
    if len(T) == 0:
        return None
    b = pd.to_numeric(T.BLOCKNUM, errors="coerce")
    X = T[b > 0].copy(); X["B"] = b[b > 0]
    if X.B.nunique() < 15:
        return None
    pcu, by_size = 0.0, {}
    for _, g in X.groupby("B"):
        r = g.RECALLED[g.RECALLED >= 0]
        size = int(g.BLOCKSIZE.iloc[0])
        pcu += r.sum() / size
        by_size.setdefault(size, []).append(r.sum() / size)
    math = 100 * (X[X.MATHACC >= 0].MATHACC == 1).mean()
    return dict(OSpan=round(pcu, 4), OSpan_math=round(math, 1), n_sets=int(X.B.nunique()),
                **{f"recall_size{k}": np.mean(v) for k, v in sorted(by_size.items())})


def from_log(sid, sess_dir):
    """Last row for this participant in the task's appended summary log (Results.csv).
    The log has no line breaks between appended runs, so rows are recovered by pattern."""
    f = glob.glob(os.path.join(sess_dir, "Working Memory", "Results.csv"))
    if not f:
        return None
    body = open(f[0], encoding="utf-8", errors="replace").read().split("\n", 1)[-1].replace("\r", "").replace("\n", "")
    own = [r for r in LOG_ROW.findall(body) if r[0][-8:] == str(sid)]
    if not own:
        return None
    return dict(OSpan=float(own[-1][1]), OSpan_math=float(own[-1][2]), source="Results.csv log")


rows = []
for sess_dir in sorted(glob.glob(os.path.join(RAW, "sub-*", "*"))):
    sess = session_of(sess_dir)
    if sess is None or not os.path.isdir(sess_dir):
        continue
    sid = int(re.search(r"sub-(\d+)", sess_dir).group(1))
    if sid in EXCLUDED_SUBJECTS:
        continue
    files = sorted(glob.glob(os.path.join(sess_dir, "**", "ospan*.csv"), recursive=True))
    res = None
    for f in files:
        res = score_file(f)
        if res:
            res["source"] = os.path.relpath(f, RAW)
            break
    if res is None and sess == "pre":
        res = from_log(sid, sess_dir)
    if res:
        rows.append(dict(SubjectID=sid, Session=sess, **res))

D = pd.DataFrame(rows).sort_values(["SubjectID", "Session"])
D.to_csv(os.path.join(DERIVED, "ospan_scores.csv"), index=False)
B = D[D.Session == "pre"]
B[["SubjectID", "OSpan", "OSpan_math"]].to_csv(os.path.join(DERIVED, "ospan_baseline.csv"), index=False)
B[["SubjectID"] + [c for c in D.columns if c.startswith("recall_size")]].dropna().to_csv(
    os.path.join(DERIVED, "ospan_setsize.csv"), index=False)
print(f"operation span: {len(B)} baseline scores ({int((B.source == 'Results.csv log').sum())} from the task log)")
