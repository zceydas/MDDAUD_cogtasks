# =============================================================================
# 12_prediction_and_wmc.R
#
# Manuscript: Results "Lower working-memory span predicted larger drinking reduction" (Figure 7, Table S4)
#             and "Baseline control profile and drinking pattern" (Table S3).
#
#   Prognostic tests   each baseline measure (34 retest measures + operation span = 35) predicting change in
#                      drinking days and in heavy-drinking days (30 d after dosing - baseline): rank-based
#                      regression adjusted for arm and baseline drinking; 70 tests, BH-FDR.
#   Table S4           robustness of the operation-span effect.
#   Depression         the same 35 baseline measures predicting PHQ change (post-dose visit - baseline visit), adjusted
#                      for arm and baseline PHQ; 35 tests, BH-FDR (results/prognostic_tests_phq.csv).
#   Table S3           baseline cognition x baseline drinking (drinking days, heavy-drinking days, drinks per
#                      drinking day), Spearman; BH-FDR across all tests.
#
# Outputs: results/prognostic_tests.csv, results/tableS3_baseline_drinking.csv, results/tableS4_wmc_robustness.csv,
#          results/12_prediction_and_wmc.txt
# =============================================================================
source("config.R")
start_log("12_prediction_and_wmc.txt")

# ---- 70 prognostic tests ------------------------------------------------------------
P <- do.call(rbind, lapply(meta$var, function(v) { wg <- wide(v); wg <- wg[!is.na(wg$Group), ]
  do.call(rbind, lapply(c("pct_drink", "pct_heavy"), function(y) {
    r <- rankAssoc(paste0(y, "_chg1"), "pre", wg, c("Group", paste0(y, "_v0")))
    data.frame(measure = v, outcome = y, n = r[["n"]], beta = r[["beta"]], p = r[["p"]]) })) }))
w <- merge(merge(osp, dose, by = "SubjectID"), dr, by = "SubjectID"); w$pre <- w$OSpan
for (y in c("pct_drink", "pct_heavy")) { r <- rankAssoc(paste0(y, "_chg1"), "pre", w, c("Group", paste0(y, "_v0")))
  P <- rbind(P, data.frame(measure = "OSpan", outcome = y, n = r[["n"]], beta = r[["beta"]], p = r[["p"]])) }
P$q <- p.adjust(P$p, "BH"); P <- P[order(P$p), ]
write.csv(P, file.path(RESULTS, "prognostic_tests.csv"), row.names = FALSE)
note("== %d prognostic tests (BH-FDR); strongest:", nrow(P))
for (i in 1:4) with(P[i, ], note("  %-24s %-9s beta=%+.2f p=%.4f q=%.3f n=%d", measure, outcome, beta, p, q, n))

# ---- Table S4: robustness of the operation-span effect ---------------------------------------
S4 <- list(); add <- function(check, y, n, txt) S4[[length(S4) + 1]] <<- data.frame(check = check, outcome = y, n = n, result = txt)
note("\n== Operation span -> drinking change (Table S4)")
for (y in c("pct_heavy", "pct_drink")) {
  yy <- paste0(y, "_chg1"); b <- paste0(y, "_v0")
  r <- rankAssoc(yy, "pre", w, c("Group", b)); add("Primary (rank-based; arm + baseline adjusted)", y, r[1], sprintf("beta=%+.2f, p=%.4f", r[2], r[3]))
  add("FDR across 70 prognostic tests", y, r[1], sprintf("q=%.3f", P$q[P$measure == "OSpan" & P$outcome == y]))
  r2 <- rankAssoc(paste0(y, "_chg2"), "pre", w, c("Group", b)); add("90-day window", y, r2[1], sprintf("p=%.4f", r2[3]))
  r3 <- rankAssoc(yy, "pre", w[w$SubjectID != 20040130, ], c("Group", b)); add("Excluding baseline-only participant (20040130)", y, r3[1], sprintf("beta=%+.2f, p=%.4f", r3[2], r3[3]))
  ps <- sapply(seq_len(nrow(w)), function(i) rankAssoc(yy, "pre", w[-i, ], c("Group", b))[3]); add("Leave-one-out (largest p)", y, nrow(w) - 1, sprintf("p=%.3f", max(ps)))
  pm <- summary(lm(as.formula(paste(yy, "~ OSpan + Group +", b)), w))$coefficients["OSpan", 4]; add("Pearson (linear; arm + baseline adjusted)", y, nrow(w), sprintf("p=%.4f", pm))
  wm <- w[w$OSpan_math >= 85, ]; r4 <- rankAssoc(yy, "pre", wm, c("Group", b)); add("Processing accuracy >=85% only", y, r4[1], sprintf("beta=%+.2f, p=%.3f", r4[2], r4[3]))
  for (g in c("Placebo", "Psilocybin")) { x <- w[w$Group == g, ]; s <- spear(x$OSpan, x[[yy]]); add(paste(g, "arm only (Spearman)"), y, s[3], sprintf("rho=%.2f, p=%.3f", s[1], s[2])) }
  dd <- w[complete.cases(w[, c(yy, "pre", "Group")]), ]; dd$ry <- rank(dd[[yy]]); dd$rx <- rank(dd$pre); dd$rb <- rank(dd[[b]])
  add("Arm x OSpan interaction", y, nrow(dd), sprintf("p=%.2f", summary(lm(ry ~ rx * Group + rb, dd))$coefficients["rx:GroupPsilocybin", 4]))
  # covariate models: linear, drinking change ~ OSpan + arm + baseline drinking + covariate
  wc <- merge(w, alloc[, c("SubjectID", "age", "sex", "QIDSR", "antidep")], by = "SubjectID")
  for (cv in c("age", "sex", "QIDSR", "antidep")) { m <- lm(as.formula(paste(yy, "~ OSpan + Group +", b, "+", cv)), wc)
    add(paste("+", cv, "(linear)"), y, nobs(m), sprintf("p=%.4f", summary(m)$coefficients["OSpan", 4])) }
}
S4 <- do.call(rbind, S4); write.csv(S4, file.path(RESULTS, "tableS4_wmc_robustness.csv"), row.names = FALSE)
for (i in seq_len(nrow(S4))) with(S4[i, ], note("  %-48s %-9s n=%-3d %s", check, outcome, n, result))
wc <- merge(w, alloc[, c("SubjectID", "age", "QIDSR")], by = "SubjectID")
note("OSpan x age r=%.2f; x QIDS-SR r=%.2f (Spearman)", cor(wc$OSpan, wc$age, method = "spearman"), cor(wc$OSpan, wc$QIDSR, method = "spearman"))
note("OSpan x baseline drinking days rho=%+.2f", spear(w$OSpan, w$pct_drink_v0)[1])

# ---- prediction of PHQ change (parallel to the drinking tests) --------------------------------------------
PQ <- do.call(rbind, lapply(meta$var, function(v) { wg <- wide(v); wg <- wg[!is.na(wg$Group), ]
  r <- rankAssoc("PHQ_change", "pre", wg, c("Group", "PHQ_v6")); data.frame(measure = v, n = r[["n"]], beta = r[["beta"]], p = r[["p"]]) }))
wq <- merge(merge(osp, dose, by = "SubjectID"), phq, by = "SubjectID"); wq$pre <- wq$OSpan
r <- rankAssoc("PHQ_change", "pre", wq, c("Group", "PHQ_v6")); PQ <- rbind(PQ, data.frame(measure = "OSpan", n = r[["n"]], beta = r[["beta"]], p = r[["p"]]))
PQ$q <- p.adjust(PQ$p, "BH"); PQ <- PQ[order(PQ$p), ]
write.csv(PQ, file.path(RESULTS, "prognostic_tests_phq.csv"), row.names = FALSE)
note("\n== %d baseline measures predicting PHQ change (arm + baseline PHQ adjusted; BH-FDR): smallest q=%.2f", nrow(PQ), min(PQ$q))
for (i in 1:3) with(PQ[i, ], note("  %-24s beta=%+.2f p=%.4f q=%.2f n=%d", measure, beta, p, q, n))
with(PQ[PQ$measure == "OSpan", ], note("  operation span -> PHQ change: beta=%+.2f p=%.2f q=%.2f n=%d", beta, p, q, n))
for (v in c("CF_SwitchCost_transfer", "PRL_RewardAcc")) { x <- merge(wide(v), osp, by = "SubjectID"); s <- spear(x$OSpan, x$change)
  note("OSpan x change in %s: rho=%+.2f p=%.3f n=%d", v, s[1], s[2], s[3]) }

# ---- Table S3: baseline cognition x baseline drinking ------------------------------------------
# All 35 baseline measures (34 retest measures + operation span) x 3 drinking variables = 105 tests.
R <- read.csv(file.path(RESULTS, "all_measures_results.csv"))
S3 <- R[, c("var", "label", "base_n", "base_drink_rho", "base_drink_p", "base_heavy_rho", "base_heavy_p", "base_dpdd_rho", "base_dpdd_p")]
o <- sapply(c("pct_drink_v0", "pct_heavy_v0", "dpdd_v0"), function(y) spear(w$OSpan, w[[y]]))
S3 <- rbind(S3, data.frame(var = "OSpan", label = "Operation span", base_n = o[3, 1], base_drink_rho = o[1, 1], base_drink_p = o[2, 1],
                           base_heavy_rho = o[1, 2], base_heavy_p = o[2, 2], base_dpdd_rho = o[1, 3], base_dpdd_p = o[2, 3]))
q <- p.adjust(c(S3$base_drink_p, S3$base_heavy_p, S3$base_dpdd_p), "BH")
write.csv(S3, file.path(RESULTS, "tableS3_baseline_drinking.csv"), row.names = FALSE)
note("\n== Baseline cognition x baseline drinking: %d tests, smallest q=%.2f (34 retest measures only: %d tests, smallest q=%.2f)",
     length(q), min(q), 3 * nrow(R), min(p.adjust(c(R$base_drink_p, R$base_heavy_p, R$base_dpdd_p), "BH")))
for (i in which(S3$base_drink_p < .05 | S3$base_heavy_p < .05 | S3$base_dpdd_p < .05)) with(S3[i, ],
  note("  %-30s drink rho=%+.2f p=%.3f | heavy rho=%+.2f p=%.3f | drinks/day rho=%+.2f p=%.3f n=%d", label, base_drink_rho, base_drink_p, base_heavy_rho, base_heavy_p, base_dpdd_rho, base_dpdd_p, base_n))
pre <- L[L$Session == "pre", ]; s <- spear(pre$PRL_Perseveration, pre$MSIT_InterferenceEffect)
note("baseline perseveration x MSIT interference effect: rho=%+.2f p=%.3f n=%d", s[1], s[2], s[3])
