# =============================================================================
# 11_sample_and_clinical.R
#
# Manuscript: Methods "Participants and design" (sample sizes); Results "Sample, task validity and
#             clinical course" (Table 2, drinking and PHQ outcomes).
#
#   Table 2         baseline characteristics by arm for the 25 participants with a baseline session, with the
#                   arm comparison (Welch t-test or Fisher's exact test)
#   Drinking        linear mixed model, proportion of drinking (or heavy-drinking) days ~ arm x visit
#                   (baseline 90 d, 30 d after dosing, 90 d after dosing) + (1 | participant)
#   PHQ             ANCOVA PHQ at the post-dose visit ~ baseline PHQ + arm; d in pooled-SD units as in anc()
#   PHQ x drinking  Spearman correlation between PHQ change and drinking-day change
#
# Outputs: results/table2_baseline_characteristics.csv, results/11_sample_and_clinical.txt
# =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
start_log("11_sample_and_clinical.txt")

T1 <- merge(data.frame(SubjectID = COHORT), alloc, by = "SubjectID", all.x = TRUE)
T1 <- merge(merge(merge(T1, dr, by = "SubjectID", all.x = TRUE), phq, by = "SubjectID", all.x = TRUE), osp, by = "SubjectID", all.x = TRUE)
post_ids <- unique(L$SubjectID[L$Session == "post"])
note("cohort with a baseline session: n=%d; with post-dose task data: n=%d; baseline only: %s",
     length(COHORT), length(post_ids), paste(setdiff(COHORT, post_ids), dose$Group[match(setdiff(COHORT, post_ids), dose$SubjectID)], collapse = ", "))
msd <- function(x, dig = 1) sprintf(paste0("%.", dig, "f (%.", dig, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
tab <- do.call(cbind, lapply(c("Placebo", "Psilocybin"), function(g) { x <- T1[T1$alloc == g, ]
  setNames(data.frame(c(
    nrow(x), msd(x$age), sum(x$sex == "F"), sum(x$antidep == 1), msd(x$QIDSR),
    sprintf("%s [%d]", msd(x$PHQ_v6), sum(!is.na(x$PHQ_v6))),
    msd(x$pct_drink_v0, 2), msd(x$pct_heavy_v0, 2), msd(x$dpdd_v0),
    sprintf("%s [%d]", msd(x$OSpan), sum(!is.na(x$OSpan))))), g) }))
tab <- cbind(Characteristic = c("n", "Age, years, mean (SD)", "Female, n", "Taking an antidepressant, n", "QIDS-SR at screening, mean (SD)",
                                "PHQ at baseline visit, mean (SD) [n]", "Proportion of drinking days, baseline 90 days",
                                "Proportion of heavy-drinking days", "Drinks per drinking day", "Operation span, mean (SD) [n]"), tab)
# arm comparison for the Methods sentence "not statistically different": Welch t-test (continuous), Fisher's exact test (counts)
T1$Group <- factor(T1$alloc, levels = c("Placebo", "Psilocybin"))
welch <- function(v) t.test(as.formula(paste(v, "~ Group")), T1)$p.value
fish <- function(x) fisher.test(table(T1$Group, x))$p.value
tab$p <- c(NA, welch("age"), fish(T1$sex == "F"), fish(T1$antidep == 1), welch("QIDSR"), welch("PHQ_v6"),
           welch("pct_drink_v0"), welch("pct_heavy_v0"), welch("dpdd_v0"), welch("OSpan"))
tab$p <- ifelse(is.na(tab$p), "", sprintf("%.2f", tab$p))
write.csv(tab, file.path(RESULTS, "table2_baseline_characteristics.csv"), row.names = FALSE)
for (i in seq_len(nrow(tab))) note("%-48s %-18s %-18s %s", tab[i, 1], tab[i, 2], tab[i, 3], tab[i, 4])
note("smallest arm-comparison p (Welch t / Fisher): %s", min(tab$p[tab$p != ""]))

# ---- drinking -------------------------------------------------------------------
long <- merge(tl[tl$SubjectID %in% COHORT, ], dose, by = "SubjectID"); long$visitF <- factor(long$visit)
x <- merge(data.frame(SubjectID = COHORT), dr, by = "SubjectID")
note("\ndrinking days %.2f -> %.2f; heavy-drinking days %.2f -> %.2f (baseline -> 30 d after dosing; n=%d with both)",
     mean(x$pct_drink_v0), mean(x$pct_drink_v1, na.rm = TRUE), mean(x$pct_heavy_v0), mean(x$pct_heavy_v1, na.rm = TRUE), sum(!is.na(x$pct_drink_v1)))
# The TLFB file has five windows: 0 = baseline 90 d, 1 = 30 d after dosing, 2 = 90 d after dosing, and 3-4 = the
# corresponding windows after the second dosing session. The primary model uses the three windows described
# in the Methods (0-2). The model with all five windows is printed for reference only.
for (vs in list(0:2, 0:4)) for (v in c("pct_drink", "pct_heavy")) {
  lx <- long[long$visit %in% vs, ]; lx$visitF <- factor(lx$visit)
  m <- lmer(as.formula(paste(v, "~ Group * visitF + (1 | SubjectID)")), lx)
  a <- anova(m)["Group:visitF", ]
  note("%s, windows %d-%d%s: arm x visit F(%d, %.1f)=%.2f, p=%.2f", v, min(vs), max(vs), ifelse(max(vs) == 2, " (primary)", " (reference)"),
       a$NumDF, a$DenDF, a$`F value`, a$`Pr(>F)`) }
# arm effect on drinking at 30 d (ANCOVA d, Table S7)
for (v in c("pct_drink", "pct_heavy")) { w <- merge(x, dose, by = "SubjectID"); w$pre <- w[[paste0(v, "_v0")]]; w$post <- w[[paste0(v, "_v1")]]
  a <- anc(w); note("%s at 30 d: ANCOVA d=%+.2f [%.2f, %.2f] p=%.2f n=%d/%d", v, a["d"], a["lo"], a["hi"], a["p"], a["nP"], a["nS"]) }

# ---- PHQ ----------------------------------------------------------------------------
ph <- merge(phq, dose, by = "SubjectID")
pp <- ph[complete.cases(ph[, c("PHQ_v6", "PHQ_v9")]), ]
note("PHQ, all participants with both visits: %.1f -> %.1f (n=%d)", mean(pp$PHQ_v6), mean(pp$PHQ_v9), nrow(pp))
for (g in c("Placebo", "Psilocybin")) { y <- ph[ph$Group == g, ]; note("PHQ %-10s %.1f -> %.1f (n=%d)", g, mean(y$PHQ_v6), mean(y$PHQ_v9), nrow(y)) }
w <- ph; w$pre <- w$PHQ_v6; w$post <- w$PHQ_v9; a <- anc(w)
note("PHQ ANCOVA: d=%+.2f [%.2f, %.2f], p=%.3f, n=%d", a["d"], a["lo"], a["hi"], a["p"], a["nP"] + a["nS"])
y <- merge(phq, dr, by = "SubjectID"); s <- spear(y$PHQ_change, y$pct_drink_chg1); ci <- bw_ci(s[1], s[3])
note("PHQ change x drinking-day change: rho=%+.2f [%.2f, %.2f], p=%.2f, n=%d", s[1], ci[1], ci[2], s[2], s[3])
