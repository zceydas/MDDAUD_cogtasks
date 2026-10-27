# =============================================================================
# 15_sensitivity_rule_switching.R
#
# Sensitivity analysis (Results, rule switching; Supplementary Methods).
# In the primary analysis, rule-switching sessions whose accuracy did not exceed chance (binomial p >= .05
# over 240 trials; 3 baseline, 1 post-dose) are set to missing (07_assemble_measures.py). This script
# restores those sessions from derived/rule_switching_scores.csv and repeats the rule-switching drug
# effects, their couplings with drinking change, the feedback-driven sign tests and the narrow stability
# composite, so that the effect of the exclusion can be reported.
#
# Output: results/15_sensitivity_rule_switching.txt
# =============================================================================
source("config.R")
start_log("15_sensitivity_rule_switching.txt")

C <- read.csv(file.path(DERIVED, "rule_switching_scores.csv"))
chance <- C[C$binom_p >= .05, ]
note("sessions at chance (excluded in the primary analysis): %s",
     paste(sprintf("%d %s (accuracy %.2f)", chance$SubjectID, chance$Session, chance$CF_Accuracy), collapse = "; "))
cf <- grep("^CF_", names(L), value = TRUE)
Lk <- L[, setdiff(names(L), cf)]                                    # all measures except rule switching
Lk <- merge(Lk, C[, c("SubjectID", "Session", cf)], by = c("SubjectID", "Session"), all = TRUE)
Lk[, cf] <- round(Lk[, cf], 10)

note("\n== Rule-switching measures with the chance-level sessions kept")
res <- do.call(rbind, lapply(cf, function(v) { w <- wide(v, Lk); a <- anc(w); c0 <- rankAssoc("pct_drink_chg1", "change", w)
  note("%-26s drug d=%+.2f [%.2f, %.2f] p=%.3f n=%d/%d | coupling rho=%+.2f p=%.3f n=%d", v, a["d"], a["lo"], a["hi"], a["p"], a["nP"], a["nS"], c0[4], c0[5], c0[1])
  data.frame(var = v, d = a[["d"]], cpl_rho = c0[["rho"]]) }))
m <- merge(res, meta, by = "var"); m <- m[m$pole != 0, ]
R <- read.csv(file.path(RESULTS, "all_measures_results.csv")); prl <- R[R$task == "PRL" & R$pole != 0, ]
fb_d <- c(prl$d * prl$pole, m$d * m$pole); fb_c <- c(-prl$cpl_rho * prl$pole, -m$cpl_rho * m$pole)
note("feedback-driven measures toward flexibility: %d of %d (sign test p=%.3f); stability gain with drinking reduction: %d of %d",
     sum(fb_d < 0), length(fb_d), binom.test(sum(fb_d < 0), length(fb_d))$p.value, sum(fb_c > 0), length(fb_c))

W <- merge(wide("PRL_Perseveration", Lk)[, c("SubjectID", "change", "Group", "pct_drink_chg1")],
           wide("CF_SwitchCost_transfer", Lk)[, c("SubjectID", "change")], by = "SubjectID", all = TRUE, suffixes = c(".p", ".s"))
W$stab <- round(rowMeans(cbind(z(W$change.p), z(W$change.s)), na.rm = TRUE), 10); W <- W[!is.nan(W$stab), ]; Wg <- W[!is.na(W$Group), ]
tt <- t.test(stab ~ Group, Wg); r0 <- rankAssoc("pct_drink_chg1", "stab", W); r1 <- rankAssoc("pct_drink_chg1", "stab", Wg, "Group")
note("narrow composite: Welch t(%.1f)=%.2f p=%.3f; x drinking-day change rho=%+.2f p=%.3f n=%d, arm-adjusted p=%.3f",
     tt$parameter, tt$statistic, tt$p.value, r0[4], r0[5], r0[1], r1[3])

# ---- Wen et al. (2023) criterion: exclude sessions below 65% overall accuracy ----------------------------
# The primary analysis excludes sessions not above chance (binomial). Wen et al. excluded participants below
# 65% overall accuracy; this repeats the rule-switching results with that stricter session-level criterion.
note("\n== Rule-switching measures with sessions below 65%% overall accuracy excluded (Wen et al., 2023)")
low <- C[C$CF_Accuracy < .65, ]
note("sessions below 65%%: %s", paste(sprintf("%d %s (%.2f)", low$SubjectID, low$Session, low$CF_Accuracy), collapse = "; "))
L65 <- Lk; L65[paste(L65$SubjectID, L65$Session) %in% paste(low$SubjectID, low$Session), cf] <- NA
res65 <- do.call(rbind, lapply(cf, function(v) { w <- wide(v, L65); a <- anc(w); c0 <- rankAssoc("pct_drink_chg1", "change", w)
  note("%-26s drug d=%+.2f [%.2f, %.2f] p=%.3f n=%d/%d | coupling rho=%+.2f p=%.3f n=%d", v, a["d"], a["lo"], a["hi"], a["p"], a["nP"], a["nS"], c0[4], c0[5], c0[1])
  data.frame(var = v, d = a[["d"]], cpl_rho = c0[["rho"]]) }))
m <- merge(res65, meta, by = "var"); m <- m[m$pole != 0, ]
fb_d <- c(prl$d * prl$pole, m$d * m$pole); fb_c <- c(-prl$cpl_rho * prl$pole, -m$cpl_rho * m$pole)
note("feedback-driven measures toward flexibility: %d of %d (sign test p=%.3f); stability gain with drinking reduction: %d of %d",
     sum(fb_d < 0), length(fb_d), binom.test(sum(fb_d < 0), length(fb_d))$p.value, sum(fb_c > 0), length(fb_c))
W <- merge(wide("PRL_Perseveration", L65)[, c("SubjectID", "change", "Group", "pct_drink_chg1")],
           wide("CF_SwitchCost_transfer", L65)[, c("SubjectID", "change")], by = "SubjectID", all = TRUE, suffixes = c(".p", ".s"))
W$stab <- round(rowMeans(cbind(z(W$change.p), z(W$change.s)), na.rm = TRUE), 10); W <- W[!is.nan(W$stab), ]; Wg <- W[!is.na(W$Group), ]
tt <- t.test(stab ~ Group, Wg); r0 <- rankAssoc("pct_drink_chg1", "stab", W); r1 <- rankAssoc("pct_drink_chg1", "stab", Wg, "Group")
note("narrow composite: Welch t(%.1f)=%.2f p=%.3f; x drinking-day change rho=%+.2f p=%.3f n=%d, arm-adjusted p=%.3f",
     tt$parameter, tt$statistic, tt$p.value, r0[4], r0[5], r0[1], r1[3])
