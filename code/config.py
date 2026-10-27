"""
config.py — shared paths and participant-level decisions for the Python scoring scripts.

All paths are taken from environment variables so the code runs unchanged on any machine:

    DATA_DIR  folder that holds the raw data (default: <repository>/data, which is not tracked)
                III_Data/III_Data/sub-*/...                 raw task files, one folder per participant
                PRL_Analysis/dose_allocation_merged.csv     arm allocation and screening variables
                PRL_Analysis/allocation_corrections.csv     corrected allocation rows (01_allocation.py)
                2004_tlfb_20260601.csv                      Timeline Followback summaries
                2015 Cog Testing Dates.xlsx                 dates of baseline testing, drug administration and post-dose testing
                Integrative_Analysis/PHQ_scores.csv         PHQ at the baseline and post-dose visits
    OUT_DIR   where every derived file, table and figure is written (default: <repository>/outputs)

The raw data are not part of the public repository (participant data, IRB restrictions).
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.expanduser(os.environ.get("DATA_DIR", os.path.join(REPO, "data")))
OUT_DIR = os.path.abspath(os.path.expanduser(os.environ.get("OUT_DIR", os.path.join(REPO, "outputs"))))

RAW = os.path.join(DATA_DIR, "III_Data", "III_Data")            # one folder per participant: sub-<ID>/<pre|post>[ dose]/...
ALLOC_RAW = os.path.join(DATA_DIR, "PRL_Analysis", "dose_allocation_merged.csv")
ALLOC_FIX = os.path.join(DATA_DIR, "PRL_Analysis", "allocation_corrections.csv")
DATES_RAW = os.path.join(DATA_DIR, "2015 Cog Testing Dates.xlsx")

DERIVED = os.path.join(OUT_DIR, "derived")    # per-session scores (one row per participant x session)
RESULTS = os.path.join(OUT_DIR, "results")    # statistics reported in the manuscript
FIGURES = os.path.join(OUT_DIR, "figures")
for d in (DERIVED, RESULTS, FIGURES):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# Participant-level decisions (documented in the manuscript, Methods, Data quality)
# ---------------------------------------------------------------------------
# 23090021 belongs to a different study series and is not part of this trial.
EXCLUDED_SUBJECTS = {23090021}

# Demand selection task: two post-dose sessions in which decision-phase responses were
# inverted on 95-96% of answered trials (reversed key mapping). Accuracy is flipped for
# answered trials in these sessions only. Identified from trial-level data.
ES_RESPONSE_INVERSIONS = {("20040210", "post"), ("20040313", "post")}

# MSIT: post-dose sessions in which 98-99% of interference trials were answered by the
# target's POSITION rather than its identity. These sessions are excluded (the prepotent
# response was not overridden, so the interference effect is not interpretable).
# 05_score_msit.py detects this automatically; the set below is what that detection
# returns and is asserted there so that any change in the raw data is noticed.
MSIT_POSITION_RESPONDERS = {(20040125, "post"), (20040147, "post")}


def session_of(folder_name):
    """Map a session folder name ('pre', 'pre dose', 'post', 'post dose') to 'pre'/'post'."""
    b = os.path.basename(folder_name).strip().lower()
    if b.startswith("pre"):
        return "pre"
    if b.startswith("post"):
        return "post"
    return None
