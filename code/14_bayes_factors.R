# =============================================================================
# 14_bayes_factors.R
#
# Manuscript: Methods "Effect sizes and Bayes factors"; every BF10 / BF01 in the Results; BF10 column of
#             Table S1; Table S7.
#
# Bayes factors (BayesFactor package, default priors):
#   bfA()    drug effect: lmBF(post ~ pre + Group) / lmBF(post ~ pre)              (Rouder et al., 2012)
#   bfR()    correlation: correlationBF on ranks (Spearman), stretched-beta prior     (Ly et al., 2016)
#   bfRP()   rank regression: model with vs without the predictor of interest; continuous variables ranked
#   bfInt()  interaction with arm: model with vs without predictor x arm
#   ttestBF  one-sample (vs .5) and two-sample (composites) t tests
# BF10 = evidence for an effect; BF01 = 1/BF10. 1/3-3 inconclusive, 3-10 moderate, >10 strong.
# lmBF integrates over the effect-size parameters by Monte Carlo sampling. With the default 10,000 iterations,
# BFs near 3 varied by 5-10% between seeds, so 100,000 iterations are used and a seed is set; with this setting,
# runs with other seeds differ by about 2-3%.
#
# Outputs: results/bf_arm_effects.csv (Table S1 BF column), results/tableS7_bayes_factors.csv,
#          results/14_bayes_factors.txt
# =============================================================================
source("config.R")
suppressPackageStartupMessages(library(BayesFactor))
start_log("14_bayes_factors.txt")
set.seed(2026); options(BFprogress = FALSE)
ITER <- 1e5   # Monte Carlo iterations per lmBF model (default 10,000 gives ~5-10% error on BFs near 3)

bfA <- function(w) { d <- w[complete.cases(w[, c("pre", "post", "Group")]), ]; d$Group <- factor(d$Group)
  extractBF(lmBF(post ~ pre + Group, d, iterations = ITER) / lmBF(post ~ pre, d, iterations = ITER))$bf }
bfR <- function(x, y) { ok <- complete.cases(x, y); extractBF(correlationBF(rank(x[ok]), rank(y[ok])))$bf }
bfRP <- function(y, x, d, covs) { dd <- d[complete.cases(d[, c(y, x, covs)]), ]; dd$ry <- rank(dd[[y]]); dd$rx <- rank(dd[[x]])
  for (cv in covs) if (is.numeric(dd[[cv]])) dd[[cv]] <- rank(dd[[cv]]) else dd[[cv]] <- factor(dd[[cv]])
  f1 <- as.formula(paste(c("ry ~ rx", covs), collapse = "+")); f0 <- as.formula(paste("ry ~", paste(covs, collapse = "+")))
  c(n = nrow(dd), bf = extractBF(lmBF(f1, dd, iterations = ITER) / lmBF(f0, dd, iterations = ITER))$bf) }
bfInt <- function(y, x, d, covs) { dd <- d[complete.cases(d[, c(y, x, covs, "Group")]), ]; dd$ry <- rank(dd[[y]]); dd$rx <- rank(dd[[x]])
  dd$Group <- factor(dd$Group); for (cv in covs) dd[[cv]] <- rank(dd[[cv]])
  f1 <- as.formula(paste(c("ry ~ rx + Group + rx:Group", covs), collapse = "+")); f0 <- as.formula(paste(c("ry ~ rx + Group", covs), collapse = "+"))
  extractBF(lmBF(f1, dd, iterations = ITER) / lmBF(f0, dd, iterations = ITER))$bf }
tBF2 <- function(v, D) extractBF(ttestBF(formula = as.formula(paste(v, "~ Group")), data = transform(D, Group = factor(Group))))$bf

out <- list()
add <- function(section, claim, effect, n, bf10) {
  out[[length(out) + 1]] <<- data.frame(section = section, claim = claim, effect = effect, n = n, BF10 = bf10, BF01 = 1 / bf10,
    evidence = ifelse(bf10 > 10, "Strong for effect", ifelse(bf10 > 3, "Moderate for effect", ifelse(bf10 < 1/10, "Strong against",
               ifelse(bf10 < 1/3, "Moderate against", "Inconclusive")))))
  note("%-62s %-24s n=%-6s BF10=%-7.3g BF01=%.3g", claim, effect, n, bf10, 1 / bf10) }
ci_r <- function(r, n) { ci <- bw_ci(r, n); sprintf("%+.2f [%+.2f, %+.2f]", r, ci[1], ci[2]) }
ci_d <- function(a) sprintf("%+.2f [%.2f, %.2f]", a["d"], a["lo"], a["hi"])

# ---- 1. drug effects on all 34 measures (Table S1) ---------------------------------------------------
R <- read.csv(file.path(RESULTS, "all_measures_results.csv"))
R$BF10 <- sapply(R$var, function(v) bfA(wide(v))); R$BF01 <- 1 / R$BF10
write.csv(R[, c("var", "label", "nP", "nS", "d", "lo", "hi", "p", "q_arm", "BF10", "BF01")], file.path(RESULTS, "bf_arm_effects.csv"), row.names = FALSE)
note("== 34 drug effects: BF01 range %.2f-%.2f; BF10>3: %d; BF01>3: %d; inconclusive: %d",
     min(R$BF01), max(R$BF01), sum(R$BF10 > 3), sum(R$BF01 > 3), sum(R$BF10 <= 3 & R$BF01 <= 3))
sec <- "Drug effects at retest (ANCOVA)"
for (v in c("MSIT_InterferenceEffect", "MSIT_AccuracyDrop", "PRL_RewardAcc", "PRL_PunishAcc", "CF_SwitchCost_transfer", "CF_PostSwitchAcc_transfer",
            "PRL_Perseveration", "ES_Learn_RT", "ES_Learn_SwitchCost_RT", "ES_Learn_Acc", "ES_ChoiceProb_L1")) {
  r <- R[R$var == v, ]; add(sec, r$label, sprintf("%+.2f [%.2f, %.2f]", r$d, r$lo, r$hi), paste0(r$nP, "/", r$nS), r$BF10) }

# switch-trial and repeat-trial RT (chance-level baseline sessions masked in 07)
E <- read.csv(file.path(DERIVED, "effort_switch_repeat.csv")); LE <- add_measure(L, E[, c("SubjectID", "Session", "SwitchRT", "RepeatRT")])
for (v in c("SwitchRT", "RepeatRT")) { w <- wide(v, LE); a <- anc(w)
  add(sec, paste("Learning-phase", ifelse(v == "SwitchRT", "switch-trial", "repeat-trial"), "RT"), ci_d(a), paste0(a["nP"], "/", a["nS"]), bfA(w)) }

# demand avoidance
LA <- add_measure(L, demand_avoidance()); w <- wide("DemandAvoid", LA); a <- anc(w)
add(sec, "Demand avoidance", ci_d(a), paste0(a["nP"], "/", a["nS"]), bfA(w))

# composites (Welch comparison of change; Cohen's d printed by script 10)
S  <- read.csv(file.path(RESULTS, "stability_composite.csv"));       Sg <- S[!is.na(S$Group), ]
Z  <- read.csv(file.path(RESULTS, "stability_broad_composite.csv")); Zg <- Z[!is.na(Z$Group), ]
add(sec, "Stability composite (narrow)", "Welch t", nrow(Sg), tBF2("stab", Sg))
add(sec, "Stability composite (broad)", "Welch t", nrow(Zg), tBF2("broad", Zg))

# clinical outcomes
ph <- merge(phq, dose, by = "SubjectID"); ph$Group <- factor(ph$Group); ph <- ph[complete.cases(ph[, c("PHQ_v6", "PHQ_v9")]), ]
add(sec, "PHQ", "ANCOVA", nrow(ph), extractBF(lmBF(PHQ_v9 ~ PHQ_v6 + Group, ph, iterations = ITER) / lmBF(PHQ_v9 ~ PHQ_v6, ph, iterations = ITER))$bf)
dd <- merge(dr, dose, by = "SubjectID"); dd <- dd[dd$SubjectID %in% COHORT, ]; dd$Group <- factor(dd$Group)
for (v in c("pct_drink", "pct_heavy")) { q <- dd[complete.cases(dd[, paste0(v, c("_v0", "_v1"))]), ]
  f1 <- as.formula(paste0(v, "_v1 ~ ", v, "_v0 + Group")); f0 <- as.formula(paste0(v, "_v1 ~ ", v, "_v0"))
  add(sec, paste(v, "at 30 days"), "ANCOVA", nrow(q), extractBF(lmBF(f1, q, iterations = ITER) / lmBF(f0, q, iterations = ITER))$bf) }

# ---- 2. baseline associations -----------------------------------------------------------------------
sec <- "Baseline associations"
b <- merge(w, osp, by = "SubjectID", all.x = TRUE); x <- b$pre[!is.na(b$pre)]
add(sec, "Demand avoidance vs indifference (.5)", sprintf("mean %.2f (SD %.2f)", mean(x), sd(x)), length(x), extractBF(ttestBF(x, mu = .5))$bf)
s <- spear(b$pre, b$OSpan);  add(sec, "Demand avoidance x operation span", ci_r(s[1], s[3]), s[3], bfR(b$pre, b$OSpan))
s <- spear(b$pre, b$PHQ_v6); add(sec, "Demand avoidance x PHQ, same visit", ci_r(s[1], s[3]), s[3], bfR(b$pre, b$PHQ_v6))
s <- spear(b$post, b$PHQ_v9); add(sec, "Retest demand avoidance x retest PHQ", ci_r(s[1], s[3]), s[3], bfR(b$post, b$PHQ_v9))
pre <- LA[LA$Session == "pre", ]
s <- spear(pre$PRL_Perseveration, pre$CF_SwitchCost_transfer)
add(sec, "Perseveration x rule-switching cost (transfer)", ci_r(s[1], s[3]), s[3], bfR(pre$PRL_Perseveration, pre$CF_SwitchCost_transfer))
for (v in c("ES_Learn_Acc", "ES_Learn_RT", "ES_Learn_SwitchCost_RT", "ES_Learn_SwitchCost_Acc", "ES_Dec_Acc", "ES_Dec_RT", "ES_Dec_SwitchCost_RT", "ES_Dec_SwitchCost_Acc")) {
  s <- spear(pre$DemandAvoid, pre[[v]]); add(sec, paste("Demand avoidance x", v), sprintf("%+.2f", s[1]), s[3], bfR(pre$DemandAvoid, pre[[v]])) }

# ---- 3. cognitive change and clinical change --------------------------------------------------------
sec <- "Cognitive change and clinical change"
s <- spear(S$stab, S$pct_drink_chg1);  add(sec, "Stability composite change x drinking-day change", ci_r(s[1], s[3]), s[3], bfR(S$stab, S$pct_drink_chg1))
r <- bfRP("pct_drink_chg1", "stab", Sg, "Group"); add(sec, "  same, arm-adjusted", "rank regression", r["n"], r["bf"])
Sc <- S[!is.na(S$change.p) & !is.na(S$change.s), ]; s <- spear(Sc$stab, Sc$pct_drink_chg1)
add(sec, "  same, participants with both measures", ci_r(s[1], s[3]), s[3], bfR(Sc$stab, Sc$pct_drink_chg1))
s <- spear(Z$broad, Z$pct_drink_chg1); add(sec, "Broad composite change x drinking-day change", ci_r(s[1], s[3]), s[3], bfR(Z$broad, Z$pct_drink_chg1))
for (v in c("PRL_Perseveration", "PRL_RewardAcc")) { w2 <- wide(v); s <- spear(w2$change, w2$pct_drink_chg1)
  add(sec, paste("Change in", v, "x drinking-day change"), ci_r(s[1], s[3]), s[3], bfR(w2$change, w2$pct_drink_chg1)) }
s <- spear(S$stab, S$PHQ_change); add(sec, "Stability composite change x PHQ change", ci_r(s[1], s[3]), s[3], bfR(S$stab, S$PHQ_change))
wg <- w[!is.na(w$Group), ]
s <- spear(wg$change, wg$pct_drink_chg1); add(sec, "Demand avoidance change x drinking-day change", ci_r(s[1], s[3]), s[3], bfR(wg$change, wg$pct_drink_chg1))
s <- spear(wg$change, wg$PHQ_change);    add(sec, "Demand avoidance change x PHQ change", ci_r(s[1], s[3]), s[3], bfR(wg$change, wg$PHQ_change))
y <- merge(phq, dr, by = "SubjectID"); s <- spear(y$PHQ_change, y$pct_drink_chg1)
add(sec, "PHQ change x drinking-day change", ci_r(s[1], s[3]), s[3], bfR(y$PHQ_change, y$pct_drink_chg1))

# ---- 4. prediction from baseline --------------------------------------------------------------------
sec <- "Prediction from baseline (rank regression, adjusted for arm and baseline outcome)"
for (y in c("pct_drink_chg1", "pct_heavy_chg1", "PHQ_change")) { cv <- if (y == "PHQ_change") "PHQ_v6" else sub("chg1", "v0", y)
  r <- bfRP(y, "pre", wg, c(cv, "Group")); add(sec, paste("Demand avoidance ->", y), "main effect", r["n"], r["bf"])
  add(sec, paste("Demand avoidance x arm ->", y), "interaction", r["n"], bfInt(y, "pre", wg, cv)) }
w3 <- merge(merge(osp, dose, by = "SubjectID"), dr, by = "SubjectID"); w3$pre <- w3$OSpan
for (y in c("pct_heavy", "pct_drink")) { r <- bfRP(paste0(y, "_chg1"), "pre", w3, c("Group", paste0(y, "_v0")))
  add(sec, paste("Operation span ->", y, "change"), "main effect", r["n"], r["bf"])
  add(sec, paste("Operation span x arm ->", y, "change"), "interaction", r["n"], bfInt(paste0(y, "_chg1"), "pre", w3, paste0(y, "_v0"))) }
w4 <- merge(merge(osp, dose, by = "SubjectID"), phq, by = "SubjectID"); w4$pre <- w4$OSpan
r <- bfRP("PHQ_change", "pre", w4, c("Group", "PHQ_v6")); add(sec, "Operation span -> PHQ change", "main effect", r["n"], r["bf"])

O <- do.call(rbind, out); rownames(O) <- NULL
write.csv(O, file.path(RESULTS, "tableS7_bayes_factors.csv"), row.names = FALSE)
note("\nBF10 per measure (Table S1):")
for (i in seq_len(nrow(R))) note("  %-38s d=%+.2f p=%.3f BF10=%.2f BF01=%.2f", R$label[i], R$d[i], R$p[i], R$BF10[i], R$BF01[i])
