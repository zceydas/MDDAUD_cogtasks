# =============================================================================
# 16_test_statistics.R
#
# Test statistics (t, F, with degrees of freedom) for every p value quoted in the Results that the other
# scripts report as a p value only. Each test re-fits exactly the model of the script named in its
# section heading, and the p value is printed next to the statistic so the two can be checked against the
# manuscript and against that script's log. Correlations are reported in the text as Spearman r_s with n,
# which is their test statistic, and are not repeated here.
#
# Output: results/16_test_statistics.txt, results/drug_effect_statistics.csv
# =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
start_log("16_test_statistics.txt")

# rank-based (partial) regression, as rankAssoc() in config.R, returning the t statistic of the predictor
rk <- function(y, x, d, covs = character(0)) {
  dd <- d[complete.cases(d[, c(y, x, covs)]), ]
  dd$ry <- scale(rank(dd[[y]])); dd$rx <- scale(rank(dd[[x]]))
  m <- lm(as.formula(paste(c("ry ~ rx", covs), collapse = " + ")), dd); co <- summary(m)$coefficients
  c(beta = co["rx", 1], t = co["rx", 3], df = m$df.residual, p = co["rx", 4], n = nrow(dd))
}
fmt_t <- function(r) sprintf("beta=%+.2f, t(%d)=%.2f, p=%.4f, n=%d", r[["beta"]], r[["df"]], r[["t"]], r[["p"]], r[["n"]])
# ANCOVA post ~ pre + arm, as anc() in config.R
anc_t <- function(w, extra = character(0)) {
  d <- w[complete.cases(w[, c("pre", "post", "Group", extra)]), ]; d$Group <- factor(d$Group, levels = c("Placebo", "Psilocybin"))
  m <- lm(as.formula(paste(c("post ~ pre + Group", extra), collapse = " + ")), d); co <- summary(m)$coefficients["GroupPsilocybin", ]
  c(t = co[["t value"]], df = m$df.residual, p = co[["Pr(>|t|)"]], n = nrow(d))
}
fa <- function(a) sprintf("t(%d)=%.2f, p=%.4f, n=%d", a[["df"]], a[["t"]], a[["p"]], a[["n"]])
# interaction of arm with a rank predictor, ranks of all continuous variables (as in 10, 12, 13)
int_t <- function(y, x, d, cv = character(0)) {
  dd <- d[complete.cases(d[, c(y, x, "Group", cv)]), ]; dd$ry <- rank(dd[[y]]); dd$rx <- rank(dd[[x]])
  for (k in seq_along(cv)) dd[[paste0("rc", k)]] <- rank(dd[[cv[k]]])
  m <- lm(as.formula(paste(c("ry ~ rx * Group", if (length(cv)) paste0("rc", seq_along(cv))), collapse = " + ")), dd)
  co <- summary(m)$coefficients["rx:GroupPsilocybin", ]
  sprintf("t(%d)=%.2f, p=%.3f, n=%d", m$df.residual, co[["t value"]], co[["Pr(>|t|)"]], nrow(dd))
}

# ---- clinical course (11_sample_and_clinical.R) ------------------------------------------------
note("== Drinking and PHQ (11)")
long <- merge(tl[tl$SubjectID %in% COHORT & tl$visit %in% 0:2, ], dose, by = "SubjectID"); long$visitF <- factor(long$visit)
for (v in c("pct_drink", "pct_heavy")) {
  a <- anova(lmer(as.formula(paste(v, "~ Group * visitF + (1 | SubjectID)")), long))["Group:visitF", ]
  note("  arm x window, %-9s F(%d, %.1f)=%.2f, p=%.2f", v, a$NumDF, a$DenDF, a$`F value`, a$`Pr(>F)`) }
ph <- merge(phq, dose, by = "SubjectID"); ph$pre <- ph$PHQ_v6; ph$post <- ph$PHQ_v9
note("  PHQ ANCOVA (arm): %s", fa(anc_t(ph)))

# ---- drug effects (10_drug_effects_and_coupling.R) ------------------------------------------------
note("\n== Drug effects, ANCOVA t for every measure (10)")
DT <- do.call(rbind, lapply(meta$var, function(v) { a <- anc_t(wide(v)); data.frame(var = v, label = meta$label[meta$var == v], t = a[["t"]], df = a[["df"]], p = a[["p"]]) }))
write.csv(DT, file.path(RESULTS, "drug_effect_statistics.csv"), row.names = FALSE)
for (i in seq_len(nrow(DT))) with(DT[i, ], note("  %-40s t(%d)=%+.2f, p=%.3f", label, df, t, p))
E <- read.csv(file.path(DERIVED, "effort_switch_repeat.csv")); LE <- add_measure(L, E[, c("SubjectID", "Session", "SwitchRT", "RepeatRT")])
for (v in c("SwitchRT", "RepeatRT")) note("  learning-phase %-8s %s", v, fa(anc_t(wide(v, LE))))
note("within-arm paired t (post vs pre):")
for (v in c("PRL_Perseveration", "CF_SwitchCost_transfer")) for (g in c("Placebo", "Psilocybin")) {
  w <- wide(v); x <- w[w$Group %in% g & complete.cases(w[, c("pre", "post")]), ]; tt <- t.test(x$post, x$pre, paired = TRUE)
  note("  %-24s %-10s t(%d)=%.2f, p=%.3f", v, g, tt$parameter, tt$statistic, tt$p.value) }
SD <- read.csv(file.path(DERIVED, "session_dates.csv"))
note("drug effect adjusted for days from drug administration to retest:")
for (v in c("ES_Learn_RT", "ES_Learn_SwitchCost_RT", "ES_ChoiceProb_L1", "PRL_Perseveration"))
  note("  %-24s %s", v, fa(anc_t(merge(wide(v), SD[, c("SubjectID", "dose_to_post")], by = "SubjectID"), "dose_to_post")))

# ---- narrow and broad stability composites (10) ----------------------------------------------------
S2 <- read.csv(file.path(RESULTS, "stability_composite.csv")); S2g <- S2[!is.na(S2$Group), ]
cc <- S2g[!is.na(S2g$change.p) & !is.na(S2g$change.s), ]; tt <- t.test(stab ~ Group, cc)
note("\n== Stability composite (10)")
note("  arm difference, participants with both measures: Welch t(%.1f)=%.2f, p=%.3f, n=%d", tt$parameter, tt$statistic, tt$p.value, nrow(cc))
note("  x drinking-day change, arm-adjusted: %s", fmt_t(rk("pct_drink_chg1", "stab", S2g, "Group")))
Si <- merge(S2g, SD[, c("SubjectID", "dose_to_post")], by = "SubjectID")
note("  x drinking-day change, adjusted for arm and days to retest: %s", fmt_t(rk("pct_drink_chg1", "stab", Si, c("Group", "dose_to_post"))))
note("  arm x composite interaction: %s", int_t("pct_drink_chg1", "stab", S2g))
note("  x drinking-day change, both measures, arm-adjusted: %s", fmt_t(rk("pct_drink_chg1", "stab", cc, "Group")))
note("  x PHQ change, arm-adjusted: %s", fmt_t(rk("PHQ_change", "stab", S2g, "Group")))
S2a <- S2[!is.nan(S2$stab_adj) & !is.na(S2$stab_adj), ]
note("  baseline-adjusted composite x drinking-day change, arm-adjusted: %s", fmt_t(rk("pct_drink_chg1", "stab_adj", S2a[!is.na(S2a$Group), ], "Group")))
B2 <- read.csv(file.path(RESULTS, "stability_broad_composite.csv")); tt <- t.test(broad ~ Group, B2[!is.na(B2$Group), ])
note("  broad composite, arm difference: Welch t(%.1f)=%.2f, p=%.2f", tt$parameter, tt$statistic, tt$p.value)

# ---- couplings with drinking (10) -------------------------------------------------------------------
note("\n== Couplings, arm-adjusted rank regression (10)")
for (v in c("PRL_RewardAcc", "PRL_Perseveration", "CF_SwitchCost_transfer")) { w <- wide(v); wg <- w[!is.na(w$Group), ]
  note("  %-24s x drinking-day change: %s", v, fmt_t(rk("pct_drink_chg1", "change", wg, "Group"))) }
wp <- wide("PRL_Perseveration"); wpg <- wp[!is.na(wp$Group), ]
note("  retest perseveration x heavy-drinking days (30 d), adjusted for arm and baseline perseveration: %s",
     fmt_t(rk("pct_heavy_v1", "post", wpg, c("Group", "pre"))))

# ---- demand avoidance (13_demand_avoidance.R) --------------------------------------------------------
note("\n== Demand avoidance (13)")
LA <- add_measure(L, demand_avoidance()); w <- wide("DemandAvoid", LA); wg <- w[!is.na(w$Group), ]
note("  arm effect at retest (ANCOVA): %s", fa(anc_t(w)))
for (y in c("pct_drink_chg1", "pct_heavy_chg1", "PHQ_change")) {
  cv <- if (y == "PHQ_change") "PHQ_v6" else sub("chg1", "v0", y)
  dd <- wg[complete.cases(wg[, c(y, "pre", cv)]), ]
  note("  baseline avoidance -> %-14s %s; interaction with arm %s", y, fmt_t(rk(y, "pre", dd, c("Group", cv))), int_t(y, "pre", dd, cv)) }

# ---- prediction from baseline (12_prediction_and_wmc.R) ------------------------------------------------
note("\n== Prediction of drinking change from baseline (12)")
wo <- merge(merge(osp, dose, by = "SubjectID"), dr, by = "SubjectID"); wo$pre <- wo$OSpan
for (y in c("pct_heavy", "pct_drink")) {
  yy <- paste0(y, "_chg1"); b <- paste0(y, "_v0")
  note("  OSpan -> %-9s change: %s", y, fmt_t(rk(yy, "pre", wo, c("Group", b))))
  note("     90-day window: %s", fmt_t(rk(paste0(y, "_chg2"), "pre", wo, c("Group", b))))
  m <- lm(as.formula(paste(yy, "~ OSpan + Group +", b)), wo); co <- summary(m)$coefficients["OSpan", ]
  note("     Pearson (linear) model: t(%d)=%.2f, p=%.4f", m$df.residual, co[["t value"]], co[["Pr(>|t|)"]])
  lo <- sapply(seq_len(nrow(wo)), function(i) rk(yy, "pre", wo[-i, ], c("Group", b))[c("t", "p")])
  note("     leave-one-out: smallest t=%.2f (largest p=%.3f), df=%d", min(lo["t", ]), max(lo["p", ]), nrow(wo) - 1 - 4)
  note("     without 20040130: %s", fmt_t(rk(yy, "pre", wo[wo$SubjectID != 20040130, ], c("Group", b))))
  wc <- merge(wo, alloc[, c("SubjectID", "age", "sex", "QIDSR", "antidep")], by = "SubjectID")
  for (cv in c("age", "sex", "QIDSR", "antidep")) { m <- lm(as.formula(paste(yy, "~ OSpan + Group +", b, "+", cv)), wc); co <- summary(m)$coefficients["OSpan", ]
    note("     + %-8s (linear): t(%d)=%.2f, p=%.4f", cv, m$df.residual, co[["t value"]], co[["Pr(>|t|)"]]) }
  note("     processing accuracy >=85%%: %s", fmt_t(rk(yy, "pre", wo[wo$OSpan_math >= 85, ], c("Group", b))))
  note("     arm x span interaction: %s", int_t(yy, "pre", wo, b))
}
wq <- merge(merge(osp, dose, by = "SubjectID"), phq, by = "SubjectID"); wq$pre <- wq$OSpan
note("  OSpan -> PHQ change (arm + baseline PHQ adjusted): %s", fmt_t(rk("PHQ_change", "pre", wq, c("Group", "PHQ_v6"))))
w <- wide("ES_Dec_Acc"); wg <- w[!is.na(w$Group), ]
note("  decision-phase accuracy -> drinking-day change: %s", fmt_t(rk("pct_drink_chg1", "pre", wg, c("Group", "pct_drink_v0"))))

# ---- statistics for the supplementary tables (S1, S2) ---------------------------------------------------
# S1: within-arm paired t (post vs pre) and ANCOVA t; S2: arm-adjusted coupling t and PHQ-change rank-regression t.
paired_t <- function(w, g) { x <- w[w$Group %in% g & complete.cases(w[, c("pre", "post")]), ]
  if (nrow(x) < 3 || sd(x$post - x$pre) == 0) return(c(t = NA, df = NA))
  tt <- t.test(x$post, x$pre, paired = TRUE); c(t = unname(tt$statistic), df = unname(tt$parameter)) }
SS <- do.call(rbind, lapply(meta$var, function(v) { w <- wide(v); wg <- w[!is.na(w$Group), ]
  a <- anc_t(w); pP <- paired_t(w, "Placebo"); pS <- paired_t(w, "Psilocybin")
  c1 <- rk("pct_drink_chg1", "change", wg, "Group"); m0 <- rk("PHQ_change", "change", wg, "Group")
  data.frame(var = v, anc_t = a[["t"]], anc_df = a[["df"]], pairP_t = pP[["t"]], pairP_df = pP[["df"]], pairS_t = pS[["t"]], pairS_df = pS[["df"]],
             cpl_adj_t = c1[["t"]], cpl_adj_df = c1[["df"]], phq_t = m0[["t"]], phq_df = m0[["df"]]) }))
write.csv(SS, file.path(RESULTS, "supplement_statistics.csv"), row.names = FALSE)
note("\nwrote supplement_statistics.csv (S1 and S2 test statistics, %d measures)", nrow(SS))
