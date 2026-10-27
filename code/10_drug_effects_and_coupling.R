# =============================================================================
# 10_drug_effects_and_coupling.R
#
# Manuscript: Results, "Drug effects across the battery", "Feedback-driven and cued shifting changed
#             in opposite directions", "Change toward stability accompanied drinking reduction";
#             Tables S1 and S2 (without Bayes factors; see 14_bayes_factors.R); Figures 2-4 and 6 data.
#
# For each of the 34 retest measures:
#   drug effect      ANCOVA post ~ pre + arm (anc() in config.R); BH-FDR across the 34 p values
#   coupling         Spearman rho between change in the measure and change in the proportion of
#                    drinking days (30 d after dosing - baseline), pooled across arms; arm-adjusted rank
#                    regression; BH-FDR across the 34 pooled p values. Same for heavy-drinking days.
#   PHQ coupling     Spearman rho between change in the measure and PHQ change, pooled across arms (BH-FDR across
#                    the 34); rank regression of PHQ change on change in the measure, adjusted for arm
#   baseline         Spearman rho between the baseline score and baseline drinking (Table S3; FDR in 12)
# Re-signing (17 measures with a pole): d_stab = d * pole (positive = drug moved the measure toward
# stability); cpl_stab = -rho * pole (positive = gain in stability accompanied drinking REDUCTION,
# because a reduction is a negative change in drinking).
#
# Outputs (OUT_DIR/results)
#   all_measures_results.csv        one row per measure (Tables S1, S2, S3)
#   stability_composite.csv         narrow stability composite per participant
#   stability_broad_composite.csv   broad (8-index) composite per participant
#   switch_repeat_effects.csv       learning-phase switch-trial and repeat-trial RT (drug effects)
#   10_drug_effects_and_coupling.txt  every number quoted in the text
# =============================================================================
source("config.R")
start_log("10_drug_effects_and_coupling.txt")

# ---- per-measure statistics --------------------------------------------------
R <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
  v <- meta$var[i]; w <- wide(v); a <- anc(w); wg <- w[!is.na(w$Group), ]
  c0 <- rankAssoc("pct_drink_chg1", "change", w)                 # pooled Spearman
  c1 <- rankAssoc("pct_drink_chg1", "change", wg, "Group")       # arm-adjusted
  h0 <- rankAssoc("pct_heavy_chg1", "change", w)
  m0 <- rankAssoc("PHQ_change", "change", wg, "Group"); d0 <- rankAssoc("PHQ_change", "change", w)
  b0 <- rankAssoc("pct_drink_v0", "pre", w); b1 <- rankAssoc("pct_heavy_v0", "pre", w); b2 <- rankAssoc("dpdd_v0", "pre", w)
  data.frame(meta[i, ], t(a),
             cpl_n = c0[["n"]], cpl_rho = c0[["rho"]], cpl_p = c0[["rho_p"]], cpl_adj_p = c1[["p"]],
             cplH_rho = h0[["rho"]], cplH_p = h0[["rho_p"]], phq_b = m0[["beta"]], phq_p = m0[["p"]],
             phq_n = d0[["n"]], phq_rho = d0[["rho"]], phq_rho_p = d0[["rho_p"]],
             base_n = b0[["n"]], base_drink_rho = b0[["rho"]], base_drink_p = b0[["rho_p"]],
             base_heavy_rho = b1[["rho"]], base_heavy_p = b1[["rho_p"]], base_dpdd_rho = b2[["rho"]], base_dpdd_p = b2[["rho_p"]])
}))
R$q_arm <- p.adjust(R$p, "BH")
R$q_cpl <- p.adjust(R$cpl_p, "BH")
R$q_phq <- p.adjust(R$phq_rho_p, "BH")
R$d_stab   <- R$d * R$pole
R$cpl_stab <- -R$cpl_rho * R$pole
R$phq_stab <- -R$phq_rho * R$pole      # positive = gain in stability accompanied a PHQ reduction
write.csv(R, file.path(RESULTS, "all_measures_results.csv"), row.names = FALSE)

# ---- drug effects -------------------------------------------------------------
note("== Drug effects (34 ANCOVAs)")
note("smallest q (BH-FDR, 34 tests) = %.2f", min(R$q_arm))
for (v in c("MSIT_InterferenceEffect", "MSIT_AccuracyDrop", "PRL_RewardAcc", "PRL_PunishAcc")) {
  r <- R[R$var == v, ]; note("%-28s d=%+.2f [%.2f, %.2f] p=%.3f", r$label, r$d, r$lo, r$hi, r$p) }
note("uncorrected p<.05:")
for (i in which(R$p < .05)) with(R[i, ], note("  %-36s d=%+.2f [%.2f, %.2f] p=%.3f  placebo %.3f->%.3f  psilocybin %.3f->%.3f  r(baseline,change)=%.2f  n=%d/%d",
                                              label, d, lo, hi, p, preP, postP, preS, postS, rbc, nP, nS))

# within-arm paired t statistics quoted in the text
paired <- function(v, g) { w <- wide(v); x <- w[w$Group %in% g & complete.cases(w[, c("pre", "post")]), ]
  t <- t.test(x$post, x$pre, paired = TRUE); note("  %-26s %-10s %.3f -> %.3f  t(%d)=%.2f p=%.3f n=%d",
                                                 v, g, mean(x$pre), mean(x$post), t$parameter, t$statistic, t$p.value, nrow(x)) }
note("within-arm change:")
for (v in c("PRL_Perseveration", "CF_SwitchCost_transfer")) for (g in c("Placebo", "Psilocybin")) paired(v, g)
for (v in c("CF_PostSwitchAcc_transfer", "CF_SwitchCost_transfer", "PRL_Perseveration", "ES_Learn_SwitchCost_RT", "ES_Learn_RT", "ES_Learn_Acc")) {
  r <- R[R$var == v, ]; note("%-36s d=%+.2f [%.2f, %.2f] p=%.3f n=%d/%d | placebo %.3f->%.3f, psilocybin %.3f->%.3f", r$label, r$d, r$lo, r$hi, r$p, r$nP, r$nS, r$preP, r$postP, r$preS, r$postS) }
# toward flexibility within the psilocybin arm = a pole measure whose within-arm change has sign -pole
tf <- R$pole != 0 & R$pS < .05 & sign(R$postS - R$preS) == -R$pole
note("pole measures moving significantly toward flexibility within the psilocybin arm: %d", sum(tf, na.rm = TRUE))

# ---- sensitivity: interval between drug administration and post-dose testing ----------------
# Post-dose testing took place 4-9 days after drug administration (01_allocation.py, session_dates.csv).
# The interval is added as a covariate to the ANCOVA for the drug effects quoted in the text.
SD <- read.csv(file.path(DERIVED, "session_dates.csv"))
note("\n== Sensitivity: days from drug administration to post-dose testing as a covariate")
for (g in c("Placebo", "Psilocybin")) { x <- SD$dose_to_post[SD$SubjectID %in% dose$SubjectID[dose$Group == g]]
  note("  interval %-10s median %.0f (range %d-%d), n=%d", g, median(x), min(x), max(x), length(x)) }
for (v in c("ES_Learn_RT", "ES_Learn_SwitchCost_RT", "ES_ChoiceProb_L1", "PRL_Perseveration")) {
  w <- merge(wide(v), SD[, c("SubjectID", "dose_to_post")], by = "SubjectID")
  d <- w[complete.cases(w[, c("pre", "post", "Group", "dose_to_post")]), ]
  p0 <- summary(lm(post ~ pre + Group, d))$coefficients["GroupPsilocybin", 4]
  co <- summary(lm(post ~ pre + Group + dose_to_post, d))$coefficients
  note("  %-24s arm p=%.3f -> %.3f with interval covariate (interval p=%.2f), n=%d", v, p0, co["GroupPsilocybin", 4], co["dose_to_post", 4], nrow(d)) }

# ---- re-signed drug effects and sign tests (descriptive) ---------------------------
S  <- R[R$pole != 0, ]                                   # 17 measures with a pole
fb <- S[S$task %in% c("PRL", "Rule switching"), ]          # 11 feedback-driven
cu <- S[S$task == "Demand selection", ]                    # 4 cued task-switch costs
note("\n== Position on the continuum (17 measures)")
note("feedback-driven measures moved toward flexibility: %d of %d, sign test p=%.3f", sum(fb$d_stab < 0), nrow(fb), binom.test(sum(fb$d_stab < 0), nrow(fb))$p.value)
note("cued switch costs moved toward stability: %d of %d; d from %.2f to %.2f", sum(cu$d_stab > 0), nrow(cu), min(cu$d), max(cu$d))

# ---- switch-trial and repeat-trial RT, demand selection learning phase ------------------
# effort_switch_repeat.csv already has chance-level baseline sessions masked (07_assemble_measures.py).
E  <- read.csv(file.path(DERIVED, "effort_switch_repeat.csv"))
LE <- add_measure(L, E[, c("SubjectID", "Session", "SwitchRT", "RepeatRT")])
SR <- do.call(rbind, lapply(c("SwitchRT", "RepeatRT"), function(v) data.frame(var = v, t(anc(wide(v, LE))))))
write.csv(SR, file.path(RESULTS, "switch_repeat_effects.csv"), row.names = FALSE)
for (i in 1:2) with(SR[i, ], note("%-9s d=%+.2f [%.2f, %.2f] p=%.3f n=%d/%d", var, d, lo, hi, p, nP, nS))

# ---- narrow stability composite ---------------------------------------------------------
# Mean of z-scored change in perseveration and in transfer-phase switch cost (one or both available).
W <- merge(wide("PRL_Perseveration")[, c("SubjectID", "change", "Group", "pct_drink_chg1", "pct_heavy_chg1", "PHQ_change")],
           wide("CF_SwitchCost_transfer")[, c("SubjectID", "change")], by = "SubjectID", all = TRUE, suffixes = c(".p", ".s"))
W$stab <- round(rowMeans(cbind(z(W$change.p), z(W$change.s)), na.rm = TRUE), 10)
W <- W[!is.nan(W$stab), ]; Wg <- W[!is.na(W$Group), ]
tt <- t.test(stab ~ Group, Wg)                       # Welch
cc <- Wg[!is.na(Wg$change.p) & !is.na(Wg$change.s), ]
note("\n== Narrow stability composite")
note("placebo mean z=%+.2f (n=%d); psilocybin %+.2f (n=%d); Welch t(%.1f)=%.2f (placebo - psilocybin), p=%.3f; complete cases p=%.3f (n=%d)",
     mean(Wg$stab[Wg$Group == "Placebo"]), sum(Wg$Group == "Placebo"), mean(Wg$stab[Wg$Group == "Psilocybin"]), sum(Wg$Group == "Psilocybin"),
     tt$parameter, tt$statistic, tt$p.value, t.test(stab ~ Group, cc)$p.value, nrow(cc))
# Cohen's d (psilocybin - placebo) with the usual large-sample CI (Table S7)
cohen <- function(x, g) { a <- x[g == "Placebo"]; b <- x[g == "Psilocybin"]; n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2)); d <- (mean(b) - mean(a)) / sp
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))); c(d = d, lo = d - 1.96 * se, hi = d + 1.96 * se) }
cd <- cohen(Wg$stab, Wg$Group); note("Cohen's d = %+.2f [%.2f, %.2f]", cd[1], cd[2], cd[3])
r0 <- rankAssoc("pct_drink_chg1", "stab", W); r1 <- rankAssoc("pct_drink_chg1", "stab", Wg, "Group")
ci <- bw_ci(r0[["rho"]], r0[["n"]])
note("x drinking-day change: rho=%+.2f [%.2f, %.2f], p=%.3f, n=%d; arm-adjusted p=%.3f", r0[4], ci[1], ci[2], r0[5], r0[1], r1[3])
dd <- Wg[complete.cases(Wg[, c("stab", "pct_drink_chg1")]), ]; dd$ry <- rank(dd$pct_drink_chg1); dd$rx <- rank(dd$stab)
note("arm x composite interaction p=%.2f", summary(lm(ry ~ rx * Group, dd))$coefficients["rx:GroupPsilocybin", 4])
for (g in c("Placebo", "Psilocybin")) { x <- Wg[Wg$Group == g, ]; s <- spear(x$stab, x$pct_drink_chg1)
  note("  within %-10s rho=%+.2f p=%.3f n=%d", g, s[1], s[2], s[3]) }
# participants with both component measures (the others contribute perseveration change only)
note("participants with both measures: %d of %d; perseveration only: %d; switch cost only: %d", sum(!is.na(W$change.p) & !is.na(W$change.s)), nrow(W),
     sum(!is.na(W$change.p) & is.na(W$change.s)), sum(is.na(W$change.p) & !is.na(W$change.s)))
Wc <- W[!is.na(W$change.p) & !is.na(W$change.s), ]
r0c <- rankAssoc("pct_drink_chg1", "stab", Wc); r1c <- rankAssoc("pct_drink_chg1", "stab", Wc[!is.na(Wc$Group), ], "Group"); cic <- bw_ci(r0c[["rho"]], r0c[["n"]])
note("complete cases x drinking-day change: rho=%+.2f [%.2f, %.2f], p=%.3f, n=%d; arm-adjusted p=%.3f", r0c[4], cic[1], cic[2], r0c[5], r0c[1], r1c[3])
s <- spear(W$stab, W$PHQ_change); r1 <- rankAssoc("PHQ_change", "stab", Wg, "Group")
note("x PHQ change: rho=%+.2f [%.2f, %.2f] p=%.3f n=%d; arm-adjusted p=%.3f", s[1], bw_ci(s[1], s[3])[1], bw_ci(s[1], s[3])[2], s[2], s[3], r1[3])
Wi <- merge(Wg, SD[, c("SubjectID", "dose_to_post")], by = "SubjectID")
note("x drinking-day change adjusted for arm and days from drug administration to retest: p=%.3f (n=%d)",
     rankAssoc("pct_drink_chg1", "stab", Wi, c("Group", "dose_to_post"))[3], sum(complete.cases(Wi[, c("stab", "pct_drink_chg1", "dose_to_post")])))
# Sensitivity: baseline-adjusted composite. The arms differed at baseline on transfer-phase switch cost, so change
# scores can carry regression to the mean. Each measure's post-dose score is residualized on its baseline score
# (pooled across arms) before z-scoring; the composite is the mean of the available residual z-scores.
resid_post <- function(w) { r <- rep(NA_real_, nrow(w)); ok <- complete.cases(w[, c("pre", "post")]); r[ok] <- residuals(lm(post ~ pre, w[ok, ])); r }
wp_ <- wide("PRL_Perseveration"); ws_ <- wide("CF_SwitchCost_transfer")
A <- merge(data.frame(SubjectID = wp_$SubjectID, adj.p = resid_post(wp_)), data.frame(SubjectID = ws_$SubjectID, adj.s = resid_post(ws_)), by = "SubjectID", all = TRUE)
W <- merge(W, A, by = "SubjectID", all.x = TRUE)
W$stab_adj <- round(rowMeans(cbind(z(W$adj.p), z(W$adj.s)), na.rm = TRUE), 10)
Wa <- W[!is.nan(W$stab_adj), ]; Wag <- Wa[!is.na(Wa$Group), ]
tt <- t.test(stab_adj ~ Group, Wag); s <- spear(Wa$stab_adj, Wa$pct_drink_chg1); r1 <- rankAssoc("pct_drink_chg1", "stab_adj", Wag, "Group")
note("\n== Baseline-adjusted stability composite (sensitivity)")
note("placebo %+.2f (n=%d) psilocybin %+.2f (n=%d); Welch t(%.1f)=%.2f, p=%.3f", mean(Wag$stab_adj[Wag$Group == "Placebo"]), sum(Wag$Group == "Placebo"),
     mean(Wag$stab_adj[Wag$Group == "Psilocybin"]), sum(Wag$Group == "Psilocybin"), tt$parameter, tt$statistic, tt$p.value)
note("x drinking-day change: rho=%+.2f [%.2f, %.2f], p=%.3f, n=%d; arm-adjusted p=%.3f", s[1], bw_ci(s[1], s[3])[1], bw_ci(s[1], s[3])[2], s[2], s[3], r1[3])
write.csv(W, file.path(RESULTS, "stability_composite.csv"), row.names = FALSE)

# ---- broad composite (sensitivity) -------------------------------------------------------
# One index per construct and task, each aligned so that positive change = toward stability;
# participants need at least 3 of the 8 indices.
bro <- list(PRL_Perseveration = 1, PRL_beta = 1, PRL_alpha = -1, CF_SwitchCost_transfer = 1, CF_SwitchCost_learning = 1,
            ES_Learn_SwitchCost_RT = 1, ES_Dec_SwitchCost_RT = 1, MSIT_InterferenceEffect = -1)
LB <- L; LB$PRL_beta <- (LB$PRL_betaGain + LB$PRL_betaLoss) / 2; LB$PRL_alpha <- (LB$PRL_alphaGain + LB$PRL_alphaLoss) / 2
Z <- Reduce(function(a, b) merge(a, b, by = "SubjectID", all = TRUE), lapply(names(bro), function(v) {
  w <- wide(v, LB); x <- data.frame(SubjectID = w$SubjectID, z = bro[[v]] * z(w$change)); names(x)[2] <- v; x }))
Z$broad <- rowMeans(Z[, -1], na.rm = TRUE); Z$nmeas <- rowSums(!is.na(Z[, names(bro)]))
Z <- merge(merge(Z, dose, by = "SubjectID", all.x = TRUE), dr, by = "SubjectID", all.x = TRUE)
Z <- merge(Z, phq, by = "SubjectID", all.x = TRUE); Z <- Z[Z$nmeas >= 3, ]; Zg <- Z[!is.na(Z$Group), ]
tt <- t.test(broad ~ Group, Zg); cd <- cohen(Zg$broad, Zg$Group)
note("\n== Broad composite (8 indices, >=3 per participant)")
note("placebo %+.2f (n=%d) psilocybin %+.2f (n=%d); Welch t(%.1f)=%.2f p=%.2f; Cohen's d=%+.2f [%.2f, %.2f]",
     mean(Zg$broad[Zg$Group == "Placebo"]), sum(Zg$Group == "Placebo"), mean(Zg$broad[Zg$Group == "Psilocybin"]), sum(Zg$Group == "Psilocybin"),
     tt$parameter, tt$statistic, tt$p.value, cd[1], cd[2], cd[3])
r0 <- rankAssoc("pct_drink_chg1", "broad", Z); ci <- bw_ci(r0[["rho"]], r0[["n"]])
note("x drinking-day change: rho=%+.2f [%.2f, %.2f] p=%.2f n=%d", r0[4], ci[1], ci[2], r0[5], r0[1])
write.csv(Z, file.path(RESULTS, "stability_broad_composite.csv"), row.names = FALSE)

# ---- coupling with drinking ------------------------------------------------------------------
note("\n== Coupling with drinking-day change (34 measures)")
note("smallest q (BH-FDR, 34 tests) = %.2f", min(R$q_cpl))
note("stability gain with drinking reduction: %d of %d measures, sign test p=%.3f", sum(S$cpl_stab > 0), nrow(S), binom.test(sum(S$cpl_stab > 0), nrow(S))$p.value)
note("  feedback-driven %d of %d; cued switch costs %d of %d; MSIT %d of %d", sum(fb$cpl_stab > 0), nrow(fb), sum(cu$cpl_stab > 0), nrow(cu),
     sum(S$cpl_stab[S$task == "MSIT"] > 0), sum(S$task == "MSIT"))
for (v in c("PRL_Perseveration", "PRL_RewardAcc")) { r <- R[R$var == v, ]; ci <- bw_ci(r$cpl_rho, r$cpl_n)
  note("%-26s rho=%+.2f [%.2f, %.2f] p=%.3f n=%d; arm-adjusted p=%.3f", r$label, r$cpl_rho, ci[1], ci[2], r$cpl_p, r$cpl_n, r$cpl_adj_p)
  w2 <- wide(v); for (g in c("Placebo", "Psilocybin")) { x <- w2[w2$Group %in% g, ]; s <- spear(x$change, x$pct_drink_chg1)
    note("    within %-10s rho=%+.2f n=%d", g, s[1], s[3]) } }
# retest perseveration and heavy drinking in the 30 days after dosing (concurrent)
wp <- wide("PRL_Perseveration"); r0 <- rankAssoc("pct_heavy_v1", "post", wp); r1 <- rankAssoc("pct_heavy_v1", "post", wp[!is.na(wp$Group), ], c("Group", "pre"))
note("retest perseveration x heavy-drinking days (30 d after dosing): rho=%+.2f p=%.3f n=%d; adjusted for arm and baseline perseveration p=%.3f", r0[4], r0[5], r0[1], r1[3])

note("\n== Couplings with uncorrected p<.05")
for (i in which(R$cpl_p < .05)) with(R[i, ], note("  %-30s rho=%+.2f p=%.3f n=%d adj p=%.3f q=%.2f", label, cpl_rho, cpl_p, cpl_n, cpl_adj_p, q_cpl))
note("PHQ couplings with p<.05: %s", paste(R$label[R$phq_p < .05], collapse = "; "))

# ---- coupling with depressive symptoms (PHQ change), parallel to the drinking analyses above -----------------
note("\n== Coupling with PHQ change (34 measures)")
note("pooled Spearman: smallest q (BH-FDR, 34 tests) = %.2f; uncorrected p<.05: %d; n=%d-%d", min(R$q_phq), sum(R$phq_rho_p < .05), min(R$phq_n), max(R$phq_n))
note("largest |rho|: %s", paste(with(R[order(-abs(R$phq_rho)), ][1:3, ], sprintf("%s rho=%+.2f p=%.2f n=%d", label, phq_rho, phq_rho_p, phq_n)), collapse = "; "))
note("stability gain with PHQ reduction: %d of %d measures, sign test p=%.2f; feedback-driven %d of %d; cued switch costs %d of %d",
     sum(S$phq_stab > 0), nrow(S), binom.test(sum(S$phq_stab > 0), nrow(S))$p.value, sum(fb$phq_stab > 0), nrow(fb), sum(cu$phq_stab > 0), nrow(cu))
qa <- p.adjust(R$phq_p, "BH")
note("arm-adjusted rank regression: smallest q (BH-FDR, 34 tests) = %.2f; p<.05: %s", min(qa),
     paste(sprintf("%s beta=%+.2f p=%.3f q=%.2f n=%d", R$label, R$phq_b, R$phq_p, qa, R$phq_n)[R$phq_p < .05], collapse = "; "))

# ---- change-score correlations between tasks (Discussion, "Placing the battery on the continuum") --------
note("\n== Correlations between change scores of different tasks")
chg <- function(v) { w <- wide(v); w[, c("SubjectID", "change")] }
for (pr in list(c("CF_SwitchCost_transfer", "ES_Learn_SwitchCost_RT"), c("CF_SwitchCost_transfer", "MSIT_InterferenceEffect"))) {
  x <- merge(chg(pr[1]), chg(pr[2]), by = "SubjectID"); s <- spear(x$change.x, x$change.y)
  note("change in %s x change in %s: rho=%+.2f p=%.3f n=%d", pr[1], pr[2], s[1], s[2], s[3]) }
