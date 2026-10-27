# =============================================================================
# 17_supplementary_analyses.R
#
# Supplement, "Supplementary Results":
#   Table S8  post-dose operation span (Methods: "post-dose scores, available for a subset, are reported in the
#             Supplement"): pre and post means by arm, practice effect, drug effect (ANCOVA, as anc() in config.R)
#             and test-retest correlation.
#   Table S9  per-level demand selection models: learning-phase RT and accuracy and decision-phase choice
#             probability at each of the four demand levels, linear mixed model
#             value ~ arm x demand level x session + (1 + session | participant), demand level as the centred
#             task-switch probability (.13, .38, .63, .88); Satterthwaite F tests (lmerTest).
#             Baseline sessions whose learning phase was not above chance are excluded, as in the main analyses.
#
# Outputs: results/tableS8_ospan_postdose.csv, results/tableS9_demand_level_models.csv,
#          results/17_supplementary_analyses.txt
# =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
start_log("17_supplementary_analyses.txt")

# ---- Table S8: post-dose operation span --------------------------------------------------------------
O <- read.csv(file.path(DERIVED, "ospan_scores.csv"))
OL <- merge(L[, c("SubjectID", "Session")], O[, c("SubjectID", "Session", "OSpan", "OSpan_math")], by = c("SubjectID", "Session"), all = TRUE)
w <- wide("OSpan", OL); a <- anc(w); wc <- w[complete.cases(w[, c("pre", "post", "Group")]), ]
tt <- t.test(wc$post, wc$pre, paired = TRUE)
d <- wc; d$Group <- factor(d$Group, levels = c("Placebo", "Psilocybin")); m <- lm(post ~ pre + Group, d); co <- summary(m)$coefficients["GroupPsilocybin", ]
r <- spear(wc$pre, wc$post)
note("== Operation span at baseline and post-dose testing (participants with both sessions)")
note("n=%d (placebo %d, psilocybin %d)", nrow(wc), a[["nP"]], a[["nS"]])
note("placebo %.2f -> %.2f; psilocybin %.2f -> %.2f", a[["preP"]], a[["postP"]], a[["preS"]], a[["postS"]])
note("practice effect (both arms): %.2f -> %.2f, paired t(%d)=%.2f, p=%.3f", mean(wc$pre), mean(wc$post), tt$parameter, tt$statistic, tt$p.value)
note("drug effect (ANCOVA): d=%+.2f [%.2f, %.2f], t(%d)=%.2f, p=%.3f", a[["d"]], a[["lo"]], a[["hi"]], m$df.residual, co[["t value"]], co[["Pr(>|t|)"]])
note("baseline-post-dose correlation: rs=%+.2f, p=%.3f, n=%d", r[1], r[2], r[3])
S8 <- data.frame(
  Statistic = c("Participants with both sessions (placebo / psilocybin)", "Placebo, baseline → post-dose", "Psilocybin, baseline → post-dose",
                "Practice effect, both arms (paired t)", "Drug effect at post-dose testing (ANCOVA)", "Baseline–post-dose correlation (Spearman)"),
  Value = c(sprintf("%d (%d / %d)", nrow(wc), a[["nP"]], a[["nS"]]), sprintf("%.2f → %.2f", a[["preP"]], a[["postP"]]),
            sprintf("%.2f → %.2f", a[["preS"]], a[["postS"]]),
            sprintf("%.2f → %.2f; t(%d)=%.2f, p=%.3f", mean(wc$pre), mean(wc$post), tt$parameter, tt$statistic, tt$p.value),
            sprintf("d=%+.2f [%.2f, %.2f]; t(%d)=%.2f, p=%.3f", a[["d"]], a[["lo"]], a[["hi"]], m$df.residual, co[["t value"]], co[["Pr(>|t|)"]]),
            sprintf("rs=%+.2f, p=%.3f", r[1], r[2])))
write.csv(S8, file.path(RESULTS, "tableS8_ospan_postdose.csv"), row.names = FALSE)

# ---- Table S9: per-level demand selection models --------------------------------------------------------
E <- read.csv(file.path(DERIVED, "demand_selection_scores.csv"))
lvl <- c(.13, .38, .63, .88); lvc <- lvl - mean(lvl)
fam <- c(`Learning-phase RT (s)` = "ES_Learn_RT_L", `Learning-phase accuracy` = "ES_Learn_Acc_L", `P(choose deck when offered)` = "ES_ChoiceProb_L")
E[E$Session == "pre" & E$ES_chance_pre %in% c(TRUE, "True"), grep("^ES_", names(E))] <- NA   # chance-level baseline sessions
E <- merge(E, dose, by = "SubjectID")
terms <- c("Group", "Level_c", "Session", "Group:Level_c", "Group:Session", "Level_c:Session", "Group:Level_c:Session")
S9 <- do.call(rbind, lapply(names(fam), function(f) {
  long <- do.call(rbind, lapply(1:4, function(k) data.frame(SubjectID = E$SubjectID, Group = factor(E$Group, levels = c("Placebo", "Psilocybin")),
    Session = factor(E$Session, levels = c("pre", "post")), Level_c = lvc[k], Value = E[[paste0(fam[[f]], k)]])))
  long <- long[!is.na(long$Value), ]
  m <- suppressMessages(lmer(Value ~ Group * Level_c * Session + (1 + Session | SubjectID), long, control = lmerControl(optimizer = "bobyqa")))
  an <- anova(m)
  note("\n== %s: %d participants, %d observations%s", f, length(unique(long$SubjectID)), nrow(long), ifelse(isSingular(m), " (singular fit)", ""))
  for (t in terms) note("  %-24s F(%d, %.1f)=%.2f, p=%.3f", t, an[t, "NumDF"], an[t, "DenDF"], an[t, "F value"], an[t, "Pr(>F)"])
  data.frame(Measure = f, Term = terms, F = an[terms, "F value"], NumDF = an[terms, "NumDF"], DenDF = an[terms, "DenDF"], p = an[terms, "Pr(>F)"],
             n = length(unique(long$SubjectID)), singular = isSingular(m), row.names = NULL)
}))
write.csv(S9, file.path(RESULTS, "tableS9_demand_level_models.csv"), row.names = FALSE)
