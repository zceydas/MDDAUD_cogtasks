"""
01_allocation.py — arm allocation and screening variables, with two documented corrections.

Input : PRL_Analysis/dose_allocation_merged.csv   allocation table (arm, sex, age, antidepressant use,
                                                  QIDS-SR at screening, baseline heavy-drinking days)
        PRL_Analysis/allocation_corrections.csv   corrected rows (participant-level values, kept with the data)
Output: derived/allocation.csv  (StudyID, alloc, sex, age, antidep, QIDSR, TLFB, Notes)

Also writes derived/session_dates.csv from the study's testing log (2015 Cog Testing Dates.xlsx):
SubjectID, dates of baseline testing, drug administration and post-dose testing, and the intervals
pre_to_dose, dose_to_post and pre_to_post in days (used in Methods, Participants, and in the
interval sensitivity analysis in 10_drug_effects_and_coupling.R).

Corrections
1. The source allocation table contained two rows with the same participant number and no row for
   another participant. The two rows were assigned to participants using self-reported age and sex in
   the task files and baseline heavy-drinking days in the Timeline Followback, and the assignment was
   confirmed against the trial record. Both rows are in the psilocybin arm, so allocation is not
   affected; the trial record also showed that one of the two participants was taking an antidepressant
   that the table had not recorded. The corrected rows are read from allocation_corrections.csv.
2. 2309-0021 is from a different study series and is removed.
"""
import pandas as pd
from config import ALLOC_RAW, ALLOC_FIX, DATES_RAW, DERIVED
import os

d = pd.read_csv(ALLOC_RAW, dtype=str)

fix = pd.read_csv(ALLOC_FIX, dtype=str)
cols = ["alloc", "sex", "age", "antidep", "QIDSR", "TLFB", "Notes"]
for _, r in fix.iterrows():
    assert (d.StudyID == r.StudyID).sum() == 1, f"expected exactly one row for {r.StudyID}"
    d.loc[d.StudyID == r.StudyID, cols] = r[cols].values

d = d[d.StudyID != "2309-0021"]
d.to_csv(os.path.join(DERIVED, "allocation.csv"), index=False)
print(f"allocation: {len(d)} rows; arms: {d.alloc.value_counts().to_dict()}")

# ---- session dates (testing log) ----
T = pd.read_excel(DATES_RAW).dropna(subset=["Study-ID"])
T = pd.DataFrame({"SubjectID": T["Study-ID"].str.replace("-", "").astype(int),
                  "pre_date": T["Pre Cog Date"].dt.date, "dose_date": T["Dosing Date"].dt.date,
                  "post_date": T["Post Cog Date"].dt.date})
T["pre_to_dose"] = [(b - a).days for a, b in zip(T.pre_date, T.dose_date)]
T["dose_to_post"] = [(b - a).days for a, b in zip(T.dose_date, T.post_date)]
T["pre_to_post"] = [(b - a).days for a, b in zip(T.pre_date, T.post_date)]
assert (T.pre_to_dose > 0).all() and (T.dose_to_post > 0).all(), "testing log: session order violated"
T.to_csv(os.path.join(DERIVED, "session_dates.csv"), index=False)
print(f"session dates: {len(T)} participants; dosing to post-dose testing median {T.dose_to_post.median():.0f} days "
      f"(range {T.dose_to_post.min()}-{T.dose_to_post.max()})")
