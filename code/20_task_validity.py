"""
20_task_validity.py — baseline task performance (Results, "Sample, task validity and clinical course";
Table 3), the PRL model comparison (Methods), the post-dose retest interval (Methods) and data checks.

All task statistics use the baseline session of the 25 randomized participants. Rule switching and demand
selection use the baseline sessions that were above chance, as in all analyses of those tasks.

Statistics
  one-sample t-tests against chance (.5) or zero; paired t-tests between conditions
  PRL model fit: likelihood-ratio test of the fitted model against random choice (chi-square, 5 df);
                 AIC and BIC of the model with vs without perseveration (paired t-test over all sessions)
  Operation span: per-participant linear slope of recall proportion on set size, t-test against zero

Output: results/table3_task_validity.csv, results/20_task_validity.txt
"""
import os, re, glob
from datetime import datetime
import numpy as np, pandas as pd
from scipy import stats
from config import RAW, DERIVED, RESULTS, EXCLUDED_SUBJECTS, session_of

LOG = open(os.path.join(RESULTS, "20_task_validity.txt"), "w")
def note(s=""):
    print(s); LOG.write(s + "\n")

A = pd.read_csv(os.path.join(DERIVED, "allocation.csv"))
A["SubjectID"] = A.StudyID.str.replace(r"\D", "", regex=True).astype(int)
ids = set(A[A.alloc.isin(["Placebo", "Psilocybin"])].SubjectID)
L = pd.read_csv(os.path.join(DERIVED, "all_measures_long.csv")); L = L[L.SubjectID.isin(ids)]
B = L[L.Session == "pre"]
T3 = []
def row(task, expected, result, n, comparison=""):
    T3.append(dict(task_measure=task, expected_effect=expected, baseline_result=result, n=n, published_comparison=comparison))

def t1(x, mu, lab):
    x = pd.Series(x).dropna(); t = stats.ttest_1samp(x, mu)
    note(f"  {lab}: n={len(x)} mean={x.mean():.3f} SD={x.std():.3f} t({len(x)-1})={t.statistic:.2f} p={t.pvalue:.2g}")
    return x.mean(), x.std(), t.pvalue, len(x)

def tpair(a, b, lab):
    x = pd.concat([a, b], axis=1).dropna(); t = stats.ttest_rel(x.iloc[:, 0], x.iloc[:, 1])
    note(f"  {lab}: n={len(x)} {x.iloc[:, 0].mean():.3f} vs {x.iloc[:, 1].mean():.3f} t({len(x)-1})={t.statistic:.2f} p={t.pvalue:.2g}")
    return t.pvalue

# ---- PRL ------------------------------------------------------------------------------------------
note("== Probabilistic reward learning (baseline)")
g = t1(B.PRL_RewardAcc, .5, "gain accuracy vs .5"); l = t1(B.PRL_PunishAcc, .5, "loss accuracy vs .5")
pgl = tpair(B.PRL_RewardAcc, B.PRL_PunishAcc, "gain vs loss accuracy")
for v in ["PRL_alphaGain", "PRL_alphaLoss", "PRL_betaGain", "PRL_betaLoss", "PRL_Perseveration"]:
    s = B[v].dropna()
    note(f"  {v}: median {s.median():.2f}, IQR {s.quantile(.25):.2f} to {s.quantile(.75):.2f}"
         + (f"; at upper bound (10): {(s >= 9.999).sum()} of {len(s)}" if "beta" in v else ""))
P = pd.read_csv(os.path.join(DERIVED, "prl_scores.csv")); P = P[P.SubjectID.isin(ids)]
Pb = P[P.Session == "pre"]
lrt = 2 * (Pb.n_trials * np.log(2) - Pb.nll_stick)               # fitted model vs random choice (p = .5 on every trial)
p_lrt = stats.chi2.sf(lrt, 5)
lik = np.exp(-Pb.nll_stick / Pb.n_trials)                        # geometric-mean per-trial choice likelihood
note(f"  model better than random choice (LRT, 5 df, p<.05): {(p_lrt < .05).sum()} of {len(Pb)}; "
     f"per-trial choice likelihood mean {lik.mean():.2f} (SD {lik.std():.2f})")
# model comparison: perseveration vs no perseveration, all sessions (baseline and post-dose)
aic_s = 2 * 5 + 2 * P.nll_stick; aic_n = 2 * 4 + 2 * P.nll_nostick
bic_s = 5 * np.log(P.n_trials) + 2 * P.nll_stick; bic_n = 4 * np.log(P.n_trials) + 2 * P.nll_nostick
ta = stats.ttest_1samp(aic_s - aic_n, 0); tb = stats.ttest_1samp(bic_s - bic_n, 0)
note(f"  adding perseveration: mean dAIC={np.mean(aic_s - aic_n):.2f}, t({len(P)-1})={ta.statistic:.2f}, p={ta.pvalue:.3f}; "
     f"mean dBIC={np.mean(bic_s - bic_n):.2f}, p={tb.pvalue:.2f} (n={len(P)} sessions)")
note(f"  test trials per session: median {int(P.n_trials.median())} (range {P.n_trials.min()}-{P.n_trials.max()})")
row("PRL: accuracy", "Above-chance learning from gains and losses",
    f"Gain {g[0]:.2f} ({g[1]:.2f}); loss {l[0]:.2f} ({l[1]:.2f}); both vs .5, p={max(g[2], l[2]):.1g}; gain > loss, p={pgl:.3f}", g[3],
    "Reward .73, punishment .74; no valence difference (Palminteri et al., 2015)")
row("PRL: model fit", "Choices predicted by the learning model",
    f"Better than random choice in {(p_lrt < .05).sum()} of {len(Pb)}; per-trial likelihood {lik.mean():.2f} ({lik.std():.2f})", len(Pb))

# check of the programmed contingency: the better cue gives the better outcome on 75% of trials
files = P[P.Session == "pre"].source_file
rows = []
for f in files:
    T = pd.read_excel(os.path.join(RAW, f)); T = T[(T.RunNo > 1) & (T.RT != 999) & T.ChosenStimulus.isin(["A", "B"])]
    best = T[T.ChosenStimAccuracy == 1]
    rows.append(((best.TrialType == "Gain") & (best.Earning > 0) | (best.TrialType == "Loss") & (best.Earning == 0)).mean())
note(f"  contingency check: better cue gave the better outcome on {np.mean(rows):.3f} of trials (range {min(rows):.2f}-{max(rows):.2f})")

# ---- rule switching ---------------------------------------------------------------------------------
note("\n== Probabilistic rule switching (baseline sessions above chance)")
C = pd.read_csv(os.path.join(DERIVED, "rule_switching_scores.csv")); C = C[C.SubjectID.isin(ids)]
note(f"  feedback validity: {C.feedback_validity.min():.3f}-{C.feedback_validity.max():.3f} across {len(C)} sessions")
Call = C[C.Session == "pre"]
ch = C[C.binom_p >= .05]
note(f"  sessions not above chance (binomial p>=.05, 240 trials): baseline {(ch.Session == 'pre').sum()} of {len(Call)}, "
     f"post {(ch.Session == 'post').sum()} of {(C.Session == 'post').sum()}; accuracy {ch.CF_Accuracy.min():.2f}-{ch.CF_Accuracy.max():.2f}; "
     f"{sorted(zip(ch.SubjectID, ch.Session))}  [excluded from all analyses]")
note(f"  baseline sessions below 65% overall accuracy (Wen et al., 2023 criterion): {(Call.CF_Accuracy < .65).sum()} of {len(Call)} "
     f"(of which {((Call.CF_Accuracy < .65) & (Call.binom_p < .05)).sum()} above chance)")
Cb = Call[Call.binom_p < .05]                                     # the sessions analysed
t1(Cb.CF_Accuracy, .5, "overall accuracy vs .5"); la = t1(Cb.Acc_learning, .5, "learning-phase accuracy vs .5")
ta_ = t1(Cb.Acc_transfer, .5, "transfer-phase accuracy vs .5"); tpair(Cb.Acc_learning, Cb.Acc_transfer, "learning vs transfer accuracy")
scl = t1(Cb.CF_SwitchCost_learning, 0, "switch cost, learning"); sct = t1(Cb.CF_SwitchCost_transfer, 0, "switch cost, transfer")
m = Cb[["CF_PreSwitchAcc_learning", "CF_PostSwitchAcc_learning", "CF_PreSwitchAcc_transfer", "CF_PostSwitchAcc_transfer"]].mean()
note(f"  pre -> post switch accuracy: learning {m.iloc[0]:.2f} -> {m.iloc[1]:.2f}; transfer {m.iloc[2]:.2f} -> {m.iloc[3]:.2f}")
row("Rule switching: switch cost", "Accuracy drops after an unannounced rule change",
    f"Learning {m.iloc[0]:.2f} -> {m.iloc[1]:.2f}; transfer {m.iloc[2]:.2f} -> {m.iloc[3]:.2f}; both p={max(scl[2], sct[2]):.1g}", len(Cb),
    f"Low-volatility condition of Wen et al. (2023); exclusion below 65% accuracy ({(Call.CF_Accuracy < .65).sum()} of {len(Call)} here)")

# ---- demand selection ------------------------------------------------------------------------------
note("\n== Parametric Demand Selection Task (baseline sessions above chance)")
E = pd.read_csv(os.path.join(DERIVED, "demand_selection_scores.csv")); E = E[E.SubjectID.isin(ids) & (E.Session == "pre")]
chance = E[E.ES_chance_pre.astype(str) == "True"]
note(f"  baseline sessions {len(E)}; not above chance {len(chance)} (accuracy {chance.ES_Learn_Acc.min():.2f}-{chance.ES_Learn_Acc.max():.2f})")
Ev = E[E.ES_chance_pre.astype(str) != "True"]
t1(Ev.ES_Learn_Acc, .5, "learning-phase accuracy"); note(f"  decision-phase accuracy {Ev.ES_Dec_Acc.mean():.2f}")
prt = tpair(Ev.ES_Learn_RT_L4, Ev.ES_Learn_RT_L1, "RT, demand level 4 vs 1 (s)")
pac = tpair(Ev.ES_Learn_Acc_L4, Ev.ES_Learn_Acc_L1, "accuracy, demand level 4 vs 1")
sc = [t1(Ev[v], 0, v) for v in ["ES_Learn_SwitchCost_RT", "ES_Dec_SwitchCost_RT", "ES_Learn_SwitchCost_Acc", "ES_Dec_SwitchCost_Acc"]]
for v in ["ES_ChoiceProb_L1", "ES_ChoiceProb_L2", "ES_ChoiceProb_L3", "ES_ChoiceProb_L4", "ES_MeanChosenLevel", "ES_SelectionRT"]:
    note(f"  {v}: {Ev[v].mean():.2f} (n={Ev[v].notna().sum()})")
D = pd.read_csv(os.path.join(DERIVED, "demand_avoidance.csv"))
D = D[D.SubjectID.isin(ids) & (D.Session == "pre") & (D.chance_pre.astype(str) != "True")]
da = t1(D.DemandAvoid, .5, "demand avoidance vs .5")
note(f"  by difference in level: {D.Avoid_d1.mean():.2f}, {D.Avoid_d2.mean():.2f}, {D.Avoid_d3.mean():.2f}")
pd31 = tpair(D.Avoid_d3, D.Avoid_d1, "avoidance, 3 vs 1 level difference")
row("Demand selection: demand effect", "Slower, less accurate responses as switch probability rises",
    f"RT {Ev.ES_Learn_RT_L1.mean():.2f} -> {Ev.ES_Learn_RT_L4.mean():.2f} s, p={prt:.1g}; accuracy {Ev.ES_Learn_Acc_L1.mean():.2f} -> {Ev.ES_Learn_Acc_L4.mean():.2f}, p={pac:.3f}",
    len(Ev), "Linear increase in RT and errors (Sayali & Badre, 2019, 2021)")
row("Demand selection: task-switch cost", "Switch trials slower and less accurate than repeat trials",
    f"RT {sc[0][0]:.2f} s (learning), {sc[1][0]:.2f} s (decision); accuracy {100*sc[2][0]:.1f} and {100*sc[3][0]:.1f} points; all p<={max(s[2] for s in sc):.1g}",
    f"{min(s[3] for s in sc)}-{max(s[3] for s in sc)}", "Signature effect of cued task switching (Monsell, 2003)")
row("Demand selection: demand avoidance", "Preference for less demanding decks",
    f"{da[0]:.2f} ({da[1]:.2f}), vs .5 p={da[2]:.2f}; no increase with demand difference (p={pd31:.2f})", da[3],
    "About half avoiders (Sayali & Badre, 2019); 67% avoidance (Sayali & Badre, 2021)")

# ---- MSIT -------------------------------------------------------------------------------------------
note("\n== MSIT (baseline)")
ie = t1(B.MSIT_InterferenceEffect, 0, "interference effect, RT (ms)"); ia = t1(B.MSIT_AccuracyDrop, 0, "interference effect, accuracy")
note("  " + ", ".join(f"{v} {B[v].mean():.3f}" for v in ["MSIT_ControlAcc", "MSIT_InterferenceAcc", "MSIT_ControlRT", "MSIT_InterferenceRT"]))
row("MSIT: interference effect", "Slower, less accurate interference trials",
    f"RT {ie[0]:.0f} ({ie[1]:.0f}) ms, p={ie[2]:.1g}; accuracy {100*ia[0]:.1f} points, p={ia[2]:.1g}", ie[3], "312 (61) ms in healthy adults (Bush et al., 2003)")

# ---- operation span ----------------------------------------------------------------------------------
note("\n== Operation span (baseline)")
O = pd.read_csv(os.path.join(DERIVED, "ospan_baseline.csv")); O = O[O.SubjectID.isin(ids)]
note(f"  partial-credit score {O.OSpan.mean():.1f} (SD {O.OSpan.std():.1f}), n={len(O)}; {100*O.OSpan.mean()/15:.0f}% of letters")
note(f"  math accuracy {O.OSpan_math.mean():.0f}% (SD {O.OSpan_math.std():.0f}); >=85%: {(O.OSpan_math >= 85).sum()} of {len(O)}")
S = pd.read_csv(os.path.join(DERIVED, "ospan_setsize.csv")); S = S[S.SubjectID.isin(ids)]
cols = [c for c in S.columns if c.startswith("recall_size")]; sizes = np.array([int(c[-1]) for c in cols])
slopes = S[cols].apply(lambda r: np.polyfit(sizes, r.values, 1)[0], axis=1); ts = stats.ttest_1samp(slopes, 0)
note(f"  recall by set size: " + ", ".join(f"{k}: {S[c].mean():.2f}" for k, c in zip(sizes, cols))
     + f"; slope {slopes.mean():.3f} per letter, t({len(S)-1})={ts.statistic:.2f}, p={ts.pvalue:.1g}, n={len(S)}")
row("Operation span: recall", "Recall declines with set size",
    f"Score {O.OSpan.mean():.1f} ({O.OSpan.std():.1f}) of 15; {100*S[cols[0]].mean():.0f}% at set size 3 to {100*S[cols[-1]].mean():.0f}% at set size 7, p={ts.pvalue:.1g}",
    f"{len(O)} ({len(S)})", "Partial-credit score 76% of letters (Redick et al., 2012; different scoring)")
row("Operation span: processing", "Math accuracy >=85%", f"{O.OSpan_math.mean():.0f}% ({O.OSpan_math.std():.0f}); {(O.OSpan_math >= 85).sum()} of {len(O)} meet criterion",
    len(O), "Criterion of Unsworth et al. (2005)")
pd.DataFrame(T3).to_csv(os.path.join(RESULTS, "table3_task_validity.csv"), index=False)

# ---- retest interval ---------------------------------------------------------------------------------
# Session dates are taken from the names of the task files used in the analyses, which carry the date the
# file was written: the PRL file (Results1_subject<ID>_<DD-Mon-YYYY>.xlsx), or, if there is none, the
# rule-switching file (..._CognitiveFlexibility_<YYYY-MM-DD>.csv). Other files in a session folder are not
# used because a few carry dates of aborted runs.
note("\n== Interval between baseline and post-dose sessions")
def file_date(f):
    m_ = re.search(r"(\d{2}-[A-Z][a-z]{2}-\d{4})", os.path.basename(f))
    if m_: return datetime.strptime(m_.group(1), "%d-%b-%Y").date()
    m_ = re.search(r"(\d{4}-\d{2}-\d{2})", os.path.basename(f))
    return datetime.strptime(m_.group(1), "%Y-%m-%d").date() if m_ else None
dates = {}
for X in (C, P):                                     # PRL dates overwrite rule-switching dates
    for _, r in X.iterrows():
        d_ = file_date(r.source_file)
        if d_: dates[(r.SubjectID, r.Session)] = d_
post_ids = set(L[L.Session == "post"].SubjectID)
iv = pd.Series({s: (dates[(s, "post")] - dates[(s, "pre")]).days for s in post_ids if (s, "pre") in dates and (s, "post") in dates})
note(f"  n={len(iv)} of {len(post_ids)} with post-dose data; median {iv.median():.0f} days (range {iv.min()}-{iv.max()})")
# Intervals relative to drug administration come from the study's testing log (01_allocation.py). The
# log's testing dates are checked against the task-file dates.
SD = pd.read_csv(os.path.join(DERIVED, "session_dates.csv"), parse_dates=["pre_date", "post_date"])
mism = [(r.SubjectID, s_) for r in SD.itertuples() for s_, col in (("pre", "pre_date"), ("post", "post_date"))
        if (r.SubjectID, s_) in dates and dates[(r.SubjectID, s_)] != getattr(r, col).date()]
note(f"  testing log vs task-file dates: {len(SD)} participants, {len(mism)} mismatches {mism if mism else ''}")
for col, lab in (("pre_to_dose", "baseline testing to drug administration"), ("dose_to_post", "drug administration to post-dose testing")):
    note(f"  {lab}: median {SD[col].median():.0f} days (range {SD[col].min()}-{SD[col].max()}; mean {SD[col].mean():.2f}, SD {SD[col].std():.2f})")
LOG.close()
