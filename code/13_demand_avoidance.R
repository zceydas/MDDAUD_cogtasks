# =============================================================================
# 13_demand_avoidance.R
#
# Manuscript: Results "No evidence that psilocybin changed the willingness to exert control" (Figure 5,
#             Table S5); demand avoidance rows of Table 3.
#
# Demand avoidance = proportion of free choices of the less demanding deck in the decision phase
# (04_score_demand_selection.py). Baseline sessions whose learning phase was not above chance are
# excluded unless stated ("incl. chance-level sessions").
#
#   Baseline        mean vs .5 (indifference); relation to operation span, to every performance measure of the
#                   same task, to the other tasks' stability-flexibility measures, to PHQ and QIDS-SR.
#   Retest          drug effect (ANCOVA), test-retest correlation, retest avoidance x retest PHQ.
#   Clinical        change in avoidance x change in drinking / PHQ; baseline avoidance (and its interaction with
#                   arm) predicting drinking or PHQ change, rank regression adjusted for the baseline outcome.
#
# Outputs: results/tableS5_demand_avoidance.csv, results/13_demand_avoidance.txt
# =============================================================================
source("config.R")
start_log("13_demand_avoidance.txt")

LA  <- add_measure(L, demand_avoidance())              # chance-level baseline sessions excluded
LA2 <- add_measure(L, demand_avoidance(raw = TRUE))    # all sessions
w  <- merge(wide("DemandAvoid", LA), osp, by = "SubjectID", all.x = TRUE)
w  <- merge(w, alloc[, c("SubjectID", "QIDSR")], by = "SubjectID", all.x = TRUE)
w2 <- wide("DemandAvoid", LA2)
S5 <- list(); add <- function(row, n, txt) { S5[[length(S5) + 1]] <<- data.frame(row = row, n = n, result = txt); note("%-62s n=%-4s %s", row, n, txt) }
fmt_r <- function(s) sprintf("rho=%+.2f, p=%.3f", s[1], s[2])

# ---- baseline -----------------------------------------------------------------------------------
x <- w$pre[!is.na(w$pre)]; tt <- t.test(x, mu = .5)
add("Baseline mean vs .5", length(x), sprintf("mean=%.2f (SD %.2f); t(%d)=%.2f, p=%.2f; avoiders %d, seekers %d, indifferent %d",
    mean(x), sd(x), tt$parameter, tt$statistic, tt$p.value, sum(x > .5), sum(x < .5), sum(x == .5)))
s <- spear(w$pre, w$OSpan); ci <- bw_ci(s[1], s[3]); add("x operation span", s[3], sprintf("%s [%.2f, %.2f]", fmt_r(s), ci[1], ci[2]))
pre <- LA[LA$Session == "pre", ]
same <- c("ES_Learn_Acc", "ES_Learn_RT", "ES_Learn_SwitchCost_RT", "ES_Learn_SwitchCost_Acc", "ES_Dec_Acc", "ES_Dec_RT", "ES_Dec_SwitchCost_RT", "ES_Dec_SwitchCost_Acc")
ss <- sapply(same, function(v) spear(pre$DemandAvoid, pre[[v]]))
add("x 8 performance measures of the same task", paste(range(ss[3, ]), collapse = "-"),
    sprintf("rho from %+.2f to %+.2f; max |rho|=%.2f; min p=%.2f", min(ss[1, ]), max(ss[1, ]), max(abs(ss[1, ])), min(ss[2, ])))
other <- c("PRL_Perseveration", "PRL_RewardAcc", "CF_SwitchCost_transfer", "CF_SwitchCost_learning", "MSIT_InterferenceEffect")
for (v in other) { s <- spear(pre$DemandAvoid, pre[[v]]); add(paste("x", v), s[3], fmt_r(s)) }
s <- spear(w$pre, w$PHQ_v6); ci <- bw_ci(s[1], s[3]); add("x PHQ, same visit", s[3], sprintf("%s [%.2f, %.2f]", fmt_r(s), ci[1], ci[2]))
s <- spear(w2$pre, w2$PHQ_v6); add("x PHQ, same visit, incl. chance-level sessions", s[3], fmt_r(s))
s <- spear(w$pre, w$QIDSR); add("x QIDS-SR at screening", s[3], fmt_r(s))
for (y in c("pct_drink_v0", "pct_heavy_v0")) { s <- spear(w$pre, w[[y]]); add(paste("x baseline", y), s[3], fmt_r(s)) }

# ---- retest -------------------------------------------------------------------------------------
s <- spear(w$post, w$PHQ_v9); ci <- bw_ci(s[1], s[3]); add("Retest avoidance x retest PHQ", s[3], sprintf("%s [%.2f, %.2f]", fmt_r(s), ci[1], ci[2]))
for (lab in c("", ", incl. chance-level sessions")) { ww <- if (lab == "") w else w2
  # arm means from participants with both sessions, as in 21_construct_validation.py
  r <- ww[!is.na(ww$Group) & !is.na(ww$pre) & !is.na(ww$post), ]
  r$a <- r$pre - ave(r$pre, r$Group); r$b <- r$post - ave(r$post, r$Group)
  s <- spear(r$a, r$b); s0 <- spear(ww$pre, ww$post)
  add(paste0("Test-retest (within arm)", lab), s[3], sprintf("%s; without removing arm means rho=%+.2f, p=%.2f", fmt_r(s), s0[1], s0[2])) }
a <- anc(w)
add("Arm effect at retest (ANCOVA)", sprintf("%d/%d", a["nP"], a["nS"]),
    sprintf("d=%+.2f [%.2f, %.2f], p=%.2f; placebo %.2f->%.2f, psilocybin %.2f->%.2f; r(baseline, change)=%.2f",
            a["d"], a["lo"], a["hi"], a["p"], a["preP"], a["postP"], a["preS"], a["postS"], a["rbc"]))
a2 <- anc(w2); add("Arm effect, incl. chance-level sessions", sprintf("%d/%d", a2["nP"], a2["nS"]), sprintf("d=%+.2f, p=%.2f", a2["d"], a2["p"]))

# ---- clinical -----------------------------------------------------------------------------------
wg <- w[!is.na(w$Group), ]
for (y in c("pct_drink_chg1", "PHQ_change")) { r <- rankAssoc(y, "change", w); r1 <- rankAssoc(y, "change", wg, "Group")
  add(paste("Change x", y), r[1], sprintf("rho=%+.2f, p=%.2f; arm-adjusted p=%.2f", r[4], r[5], r1[3])) }
for (y in c("pct_drink_chg1", "pct_heavy_chg1", "PHQ_change")) {
  cv <- if (y == "PHQ_change") "PHQ_v6" else sub("chg1", "v0", y)
  dd <- wg[complete.cases(wg[, c(y, "pre", cv)]), ]
  # adjusted main effect: same scaled-rank model as rankAssoc(); interaction: ranks of all continuous variables
  r <- rankAssoc(y, "pre", dd, c("Group", cv))
  dd$ry <- rank(dd[[y]]); dd$rx <- rank(dd$pre); dd$rc <- rank(dd[[cv]])
  m <- summary(lm(ry ~ rx * Group + rc, dd))$coefficients
  s <- spear(dd$pre, dd[[y]])
  add(paste("Baseline avoidance ->", y), r[1], sprintf("adjusted beta=%+.2f, p=%.2f; interaction with arm p=%.2f; unadjusted rho=%+.2f, p=%.3f",
                                                       r[2], r[3], m["rx:GroupPsilocybin", 4], s[1], s[2]))
  for (g in c("Placebo", "Psilocybin")) { x <- dd[dd$Group == g, ]; s <- spear(x$pre, x[[y]]); add(paste("   within", g), s[3], fmt_r(s)) }
}
write.csv(do.call(rbind, S5), file.path(RESULTS, "tableS5_demand_avoidance.csv"), row.names = FALSE)
