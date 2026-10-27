# =============================================================================
# config.R: shared inputs and helper functions for the R analyses (scripts 10-17 and 30)
#
# Every R script starts with  source("config.R")  and is run from the code/ folder.
#
# Paths come from environment variables, as in config.py:
#   DATA_DIR  raw data folder (default <repository>/data, which is not tracked). Read here:
#               2004_tlfb_20260601.csv               Timeline Followback summary per visit
#               Integrative_Analysis/PHQ_scores.csv   PHQ9 at the baseline (v6) and post-dose (v9) visits
#   OUT_DIR   output folder (default <repository>/outputs). Read from OUT_DIR/derived (written by scripts 01-07):
#               allocation.csv, all_measures_long.csv, ospan_baseline.csv,
#               demand_avoidance.csv, effort_switch_repeat.csv, session_dates.csv, ospan_scores.csv,
#               demand_selection_scores.csv, rule_switching_scores.csv
#             Results are written to OUT_DIR/results, figures to OUT_DIR/figures.
#
# Objects created
#   dose   SubjectID, Group (Placebo / Psilocybin)          every randomized participant in the allocation table
#   alloc  full allocation table (age, sex, antidepressant use, QIDS-SR at screening)
#   tl     Timeline Followback, one row per participant x visit (0 = baseline 90 d, 1 = 30 d after dosing, 2 = 90 d after dosing)
#   dr     drinking outcomes, one row per participant (wide)
#   phq    PHQ9 at the baseline (PHQ_v6) and post-dose (PHQ_v9) cognitive visits, PHQ_change = v9 - v6
#   L      all_measures_long.csv, one row per participant x session (pre / post), 34 retest measures
#   osp    baseline operation span (OSpan, OSpan_math)
#   meta   the 34 retest measures with task, label and stability-flexibility pole (Table 1)
# Functions
#   wide(v)              one row per participant: pre, post, change for measure v, joined to arm, drinking, PHQ
#   rankAssoc(y,x,d,cv)  Spearman rho plus rank-based (partial) regression of y on x adjusting for covariates
#   anc(w)               ANCOVA drug effect (post ~ pre + arm), d in pooled-SD units, within-arm paired tests
#   spear(x,y)           Spearman correlation on complete pairs (asymptotic p)
#   bw_ci(r,n)           Bonett-Wright 95% CI for a Spearman correlation
#   z(x)                 z-score
#   note(...)            prints a line and appends it to the script's log in OUT_DIR/results
# =============================================================================

suppressPackageStartupMessages({ library(stats) })

DATA_DIR <- path.expand(Sys.getenv("DATA_DIR", file.path("..", "data")))           # scripts run from code/
OUT_DIR  <- normalizePath(path.expand(Sys.getenv("OUT_DIR", file.path("..", "outputs"))), mustWork = FALSE)
DERIVED  <- file.path(OUT_DIR, "derived")
RESULTS  <- file.path(OUT_DIR, "results")
FIGURES  <- file.path(OUT_DIR, "figures")
for (d in c(RESULTS, FIGURES)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- logging: each script sets LOG before calling note() ---------------------
LOG <- NULL
start_log <- function(name) { LOG <<- file.path(RESULTS, name); cat("", file = LOG) }
note <- function(...) { s <- sprintf(...); cat(s, "\n"); if (!is.null(LOG)) cat(s, "\n", file = LOG, append = TRUE) }

# ---- arm allocation (01_allocation.py) ---------------------------------------
alloc <- read.csv(file.path(DERIVED, "allocation.csv"))
alloc$SubjectID <- as.numeric(gsub("\\D", "", alloc$StudyID))      # "2004-0114" -> 20040114
dose <- alloc[alloc$alloc %in% c("Placebo", "Psilocybin"), c("SubjectID", "alloc")]
names(dose)[2] <- "Group"

# ---- drinking: Timeline Followback --------------------------------------------
# Windows differ in length (90 d at baseline, 30 d and 90 d after dosing), so drinking is expressed
# as proportions of days: drinking days, heavy-drinking days, and drinks per drinking day.
# round(..., 10) removes floating-point noise so that ties are ranked identically on every machine.
tl <- read.csv(file.path(DATA_DIR, "2004_tlfb_20260601.csv"), fileEncoding = "UTF-8-BOM",
               colClasses = c(volid = "character"))
tl <- tl[!is.na(tl$visit) & tl$volid != "", ]
tl$SubjectID <- 20040000 + as.numeric(tl$volid)                     # volid "0114" -> 20040114
tl$pct_drink <- round(tl$drinking_days_total / tl$days_total, 10)
tl$pct_heavy <- round(tl$drinking_days_heavy / tl$days_total, 10)
tl$dpdd      <- ifelse(tl$drinking_days_total > 0, tl$mean_drinks_day, NA)   # undefined without drinking days

dr <- Reduce(function(a, b) merge(a, b, by = "SubjectID", all = TRUE), lapply(0:2, function(v) {
  x <- tl[tl$visit == v, c("SubjectID", "pct_drink", "pct_heavy", "dpdd")]
  names(x)[-1] <- paste0(names(x)[-1], "_v", v); x }))
for (dv in c("pct_drink", "pct_heavy")) {
  dr[[paste0(dv, "_chg1")]] <- round(dr[[paste0(dv, "_v1")]] - dr[[paste0(dv, "_v0")]], 10)  # primary: 30 d after dosing - baseline
  dr[[paste0(dv, "_chg2")]] <- round(dr[[paste0(dv, "_v2")]] - dr[[paste0(dv, "_v0")]], 10)  # 90 d after dosing - baseline
}

# ---- depressive symptoms ----------------------------------------------------
phq <- read.csv(file.path(DATA_DIR, "Integrative_Analysis", "PHQ_scores.csv"))
phq$SubjectID <- as.numeric(gsub("\\D", "", phq$StudyID))
phq$StudyID <- NULL

# ---- cognitive measures (07_assemble_measures.py, 06_score_ospan.py) -----------
L   <- read.csv(file.path(DERIVED, "all_measures_long.csv"))
osp <- read.csv(file.path(DERIVED, "ospan_baseline.csv"))
COHORT <- sort(unique(L$SubjectID[L$Session == "pre"]))            # 25 participants with a baseline session

# Demand avoidance with the chance-level exclusion applied (baseline sessions whose learning phase
# was not above chance). raw = TRUE keeps those sessions (sensitivity analyses).
demand_avoidance <- function(raw = FALSE) {
  DA <- read.csv(file.path(DERIVED, "demand_avoidance.csv"))
  if (!raw) DA$DemandAvoid[DA$Session == "pre" & DA$chance_pre %in% c("True", TRUE)] <- NA
  DA[, c("SubjectID", "Session", "DemandAvoid")]
}
add_measure <- function(L, X) merge(L, X, by = c("SubjectID", "Session"), all = TRUE)

# ---- helpers ------------------------------------------------------------------
# wide(): pre, post and change (post - pre) for one measure, one row per participant, with arm,
# drinking outcomes and PHQ attached. Participants without an arm (none in the final sample) keep Group = NA.
wide <- function(v, data = L) {
  pre  <- data[data$Session == "pre",  c("SubjectID", v)]; names(pre)[2]  <- "pre"
  post <- data[data$Session == "post", c("SubjectID", v)]; names(post)[2] <- "post"
  w <- merge(pre, post, by = "SubjectID", all = TRUE)
  w$change <- round(w$post - w$pre, 10)
  w <- merge(w, dose, by = "SubjectID", all.x = TRUE)
  w <- merge(w, dr,   by = "SubjectID", all.x = TRUE)
  merge(w, phq, by = "SubjectID", all.x = TRUE)
}

# rankAssoc(): association between x and y on complete cases.
#   rho, rho_p  Spearman correlation (unadjusted)
#   beta, p     rank-based regression: scale(rank(y)) ~ scale(rank(x)) + covariates (covariates entered raw;
#               a factor such as Group enters as a dummy). With no covariates beta equals rho.
rankAssoc <- function(y, x, d, covs = character(0)) {
  dd <- d[complete.cases(d[, c(y, x, covs)]), ]
  if (nrow(dd) < 6) return(c(n = nrow(dd), beta = NA, p = NA, rho = NA, rho_p = NA))
  dd$ry <- scale(rank(dd[[y]])); dd$rx <- scale(rank(dd[[x]]))
  co <- summary(lm(as.formula(paste(c("ry ~ rx", covs), collapse = " + ")), dd))$coefficients
  ct <- suppressWarnings(cor.test(dd[[x]], dd[[y]], method = "spearman"))
  c(n = nrow(dd), beta = co["rx", 1], p = co["rx", 4], rho = unname(ct$estimate), rho_p = ct$p.value)
}

# anc(): drug effect at post-dose testing. ANCOVA post ~ pre + Group on participants with both sessions.
#   d       adjusted arm difference (psilocybin - placebo) divided by the SD of all pre and post scores pooled
#   lo, hi  95% CI of d (CI of the coefficient divided by the same SD)
#   pP, pS  within-arm paired t-test of post vs pre
#   rbc     correlation between baseline and change (regression-to-the-mean diagnostic)
anc <- function(w) {
  d <- w[complete.cases(w[, c("pre", "post", "Group")]), ]
  d$Group <- factor(d$Group, levels = c("Placebo", "Psilocybin"))
  s <- sd(c(d$pre, d$post))
  m <- lm(post ~ pre + Group, d); co <- summary(m)$coefficients; ci <- confint(m)["GroupPsilocybin", ]
  wp <- function(g) { x <- d[d$Group == g, ]
    if (nrow(x) > 2 && sd(x$post - x$pre) > 0) t.test(x$post, x$pre, paired = TRUE)$p.value else NA }
  c(nP = sum(d$Group == "Placebo"), nS = sum(d$Group == "Psilocybin"),
    d = co["GroupPsilocybin", 1] / s, lo = ci[[1]] / s, hi = ci[[2]] / s, p = co["GroupPsilocybin", 4],
    preP = mean(d$pre[d$Group == "Placebo"]),    postP = mean(d$post[d$Group == "Placebo"]),
    preS = mean(d$pre[d$Group == "Psilocybin"]), postS = mean(d$post[d$Group == "Psilocybin"]),
    pP = wp("Placebo"), pS = wp("Psilocybin"), rbc = cor(d$pre, d$post - d$pre))
}

# Bonett & Wright (2000) approximate 95% CI for a Spearman correlation.
bw_ci <- function(r, n) { se <- sqrt((1 + r^2 / 2) / (n - 3)); tanh(atanh(r) + c(-1, 1) * 1.96 * se) }

spear <- function(x, y) { ok <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok)) }

# ---- the 34 retest measures (Table 1) -------------------------------------------
# pole: 1 = higher score means more stable, -1 = more flexible, 0 = not assigned
meta <- read.csv(text = 'var,task,label,pole
PRL_RewardAcc,PRL,Reward-learning accuracy,0
PRL_PunishAcc,PRL,Punishment-learning accuracy,0
PRL_alphaGain,PRL,Learning rate (gain),-1
PRL_alphaLoss,PRL,Learning rate (loss),-1
PRL_betaGain,PRL,Inverse temperature (gain),1
PRL_betaLoss,PRL,Inverse temperature (loss),1
PRL_Perseveration,PRL,Perseveration,1
CF_Accuracy,Rule switching,Overall accuracy,0
CF_SwitchCost_learning,Rule switching,Switch cost (learning),1
CF_SwitchCost_transfer,Rule switching,Switch cost (transfer),1
CF_PreSwitchAcc_learning,Rule switching,Pre-switch accuracy (learning),1
CF_PostSwitchAcc_learning,Rule switching,Post-switch accuracy (learning),-1
CF_PreSwitchAcc_transfer,Rule switching,Pre-switch accuracy (transfer),1
CF_PostSwitchAcc_transfer,Rule switching,Post-switch accuracy (transfer),-1
CF_RT_learning,Rule switching,Median RT (learning),0
CF_RT_transfer,Rule switching,Median RT (transfer),0
ES_Learn_Acc,Demand selection,Learning-phase accuracy,0
ES_Learn_RT,Demand selection,Learning-phase RT,0
ES_Learn_SwitchCost_Acc,Demand selection,Task-switch cost (accuracy; learning),1
ES_Learn_SwitchCost_RT,Demand selection,Task-switch cost (RT; learning),1
ES_Dec_Acc,Demand selection,Decision-phase accuracy,0
ES_Dec_RT,Demand selection,Decision-phase RT,0
ES_Dec_SwitchCost_Acc,Demand selection,Task-switch cost (accuracy; decision),1
ES_Dec_SwitchCost_RT,Demand selection,Task-switch cost (RT; decision),1
ES_ChoiceProb_L1,Demand selection,P(choose least demanding),0
ES_ChoiceProb_L4,Demand selection,P(choose most demanding),0
ES_MeanChosenLevel,Demand selection,Mean chosen demand level,0
ES_SelectionRT,Demand selection,Selection RT,0
MSIT_ControlAcc,MSIT,Control accuracy,0
MSIT_InterferenceAcc,MSIT,Interference accuracy,0
MSIT_ControlRT,MSIT,Control RT,0
MSIT_InterferenceRT,MSIT,Interference RT,0
MSIT_InterferenceEffect,MSIT,Interference effect (RT),-1
MSIT_AccuracyDrop,MSIT,Interference effect (accuracy),-1')

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
