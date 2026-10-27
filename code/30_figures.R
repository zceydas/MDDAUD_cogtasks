# =============================================================================
# 30_figures.R: Figures 2-7 of the manuscript, plus an export of Figure 1 (a schematic, assets/Figure1_framework_and_findings.svg;
#                Figure S1 is drawn by 21_construct_validation.py).
#
#   Figure 2  drug effects on all 34 retest measures, stability composites on top (results of 10)
#   Figure 3  stability-flexibility map: re-signed drug effects and correlations with drinking reduction,
#             stability composites on top                                     (results of 10)
#   Figure 4  baseline to post-dose change by arm for four measures           (ANCOVA as in 10)
#   Figure 5  demand avoidance: x PHQ, x operation span, by arm               (as in 13)
#   Figure 6  coupling of cognitive change with drinking change               (as in 10)
#   Figure 7  operation span predicts drinking change                         (as in 12)
# Run after 10_drug_effects_and_coupling.R. Every figure is drawn at its printed size (7.0-7.2 in wide, the
# full text width of a journal page) so that text is 7-9 pt in print, and written three ways: .pdf (vector,
# cairo), .tiff (600 dpi, LZW) for submission and .png (600 dpi) for the manuscript file. Bitmaps use ragg so
# that Greek letters render on every platform.
# Output: OUT_DIR/figures/Figure2_... to Figure7_... (.pdf, .tiff, .png), Figure3_data.csv
# =============================================================================
source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
save <- function(name, g, w, h) suppressWarnings({ f <- file.path(FIGURES, sub("\\.png$", "", name))
  ggsave(paste0(f, ".png"), g, width = w, height = h, dpi = 600, bg = "white", device = ragg::agg_png)
  ggsave(paste0(f, ".tiff"), g, width = w, height = h, dpi = 600, bg = "white", device = ragg::agg_tiff, compression = "lzw")
  ggsave(paste0(f, ".pdf"), g, width = w, height = h, bg = "white", device = cairo_pdf) })

COL  <- c(Placebo = "#4d4d4d", Psilocybin = "#9b27af", Unresolved = "grey70")
POLE <- c(Stability = "#0F6E56", Flexibility = "#D85A30", Neither = "#888780")
# Text with Greek letters and arrows is only drawn correctly under a UTF-8 locale; under the plain "C" locale
# (the default for Rscript on some servers) they are drawn as "..".
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
if (!l10n_info()[["UTF-8"]]) invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8")))
th <- theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 8.5),
                                           plot.title.position = "plot", plot.subtitle = element_text(size = 7),
                                           plot.caption = element_text(size = 6.5), legend.text = element_text(size = 7.5))
PT <- 1.7   # point size for the printed figures
R <- read.csv(file.path(RESULTS, "all_measures_results.csv"))
R$poleLab <- factor(c(`1` = "Stability", `-1` = "Flexibility", `0` = "Neither")[as.character(R$pole)], levels = names(POLE))
R$task <- factor(R$task, levels = c("PRL", "Rule switching", "Demand selection", "MSIT"))

# Spearman correlations are labelled r_s (italic r, subscript s) in every figure.
nz <- function(x, k) sub("^(-?)0\\.", "\\1.", sprintf(paste0("%.", k, "f"), x))          # .12 rather than 0.12
mn <- function(x) gsub("-", "\u2212", x)                                          # typographic minus
rs <- function(prefix, r, rest) bquote(.(prefix) * italic(r)[s] * .(paste0(" = ", mn(nz(r, 2)), rest)))
bottom <- plot_annotation(theme = theme(legend.position = "bottom"))   # legend position of collected guides

# ---- stability composites (10): arm difference (Cohen's d of change, psilocybin - placebo, Welch p) and
# correlation of composite change with drinking-day change (Spearman, approximate CI as in Figure 3B)
cohen <- function(x, g) { a <- x[g == "Placebo"]; b <- x[g == "Psilocybin"]; n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * var(a) + (n2 - 1) * var(b)) / (n1 + n2 - 2)); d <- (mean(b) - mean(a)) / sp
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2))); c(d = d, lo = d - 1.96 * se, hi = d + 1.96 * se, nP = n1, nS = n2) }
comp_row <- function(file, col, label) {
  X <- read.csv(file.path(RESULTS, file)); Xg <- X[!is.na(X$Group), ]
  cd <- cohen(Xg[[col]], Xg$Group); p <- t.test(as.formula(paste(col, "~ Group")), Xg)$p.value
  s <- spear(X[[col]], X$pct_drink_chg1)
  data.frame(var = col, label = label, task = "Composite", pole = 1, nP = cd[["nP"]], nS = cd[["nS"]], d = cd[["d"]], lo = cd[["lo"]],
             hi = cd[["hi"]], p = p, d_stab = cd[["d"]], cpl_n = s[["n"]], cpl_rho = s[["rho"]], cpl_p = s[["p"]], cpl_stab = -s[["rho"]])
}
COMP <- rbind(comp_row("stability_composite.csv", "stab", "Stability composite, change"),
              comp_row("stability_broad_composite.csv", "broad", "Broad composite (8 indices), change"))

# ---- Figure 2: drug effects on all measures ----------------------------------------------------------
F2 <- rbind(COMP[, c("label", "task", "pole", "nP", "nS", "d", "lo", "hi", "p")], R[, c("label", "task", "pole", "nP", "nS", "d", "lo", "hi", "p")])
F2$poleLab <- factor(c(`1` = "Stability", `-1` = "Flexibility", `0` = "Neither")[as.character(F2$pole)], levels = names(POLE))
F2$task <- factor(as.character(F2$task), levels = c("Composite", "PRL", "Rule switching", "Demand selection", "MSIT"))
levels(F2$task)[1] <- "Stability\ncomposites"
F2$lab2 <- factor(sprintf("%s (%d/%d)", F2$label, F2$nP, F2$nS), levels = rev(sprintf("%s (%d/%d)", F2$label, F2$nP, F2$nS)))
g2 <- ggplot(F2, aes(d, lab2)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = .25, orientation = "y", colour = "grey40") +
  geom_point(aes(fill = poleLab, shape = p < .05), size = PT) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), labels = c("p≥.05", "p<.05"), name = NULL) +
  scale_fill_manual(values = POLE, name = "Higher score:", labels = c("more stable", "more flexible", "neither")) +
  facet_grid(task ~ ., scales = "free_y", space = "free_y", switch = "y") + coord_cartesian(xlim = c(-2.6, 2.8)) +
  labs(x = "Arm difference at post-dose testing, SD units (psilocybin − placebo), 95% CI", y = NULL,
       caption = "Composites: Cohen's d of change (Welch test); measures: baseline-adjusted difference (ANCOVA).") + th +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        strip.background = element_blank(), panel.spacing = unit(5, "pt"), legend.box = "horizontal") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.2)))
save("Figure2_arm_effects_all_measures.png", g2, 7, 7.3)

# ---- Figure 4: change by arm ---------------------------------------------------------------------------
prepost <- function(w, ylab, ttl, hline = NULL, sub = "Arm \u00d7 session") {
  w <- w[complete.cases(w[, c("pre", "post", "Group")]), ]; a <- anc(w)
  long <- rbind(data.frame(SubjectID = w$SubjectID, Group = w$Group, Session = "Baseline", y = w$pre),
                data.frame(SubjectID = w$SubjectID, Group = w$Group, Session = "Post-dose", y = w$post))
  long$Session <- factor(long$Session, levels = c("Baseline", "Post-dose"))
  s <- aggregate(y ~ Group + Session, long, function(x) c(m = mean(x), se = sd(x) / sqrt(length(x))))
  s <- data.frame(s[, 1:2], m = s$y[, 1], se = s$y[, 2])
  g <- ggplot()
  if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linetype = 2, colour = "grey60")
  g + geom_line(data = long, aes(Session, y, group = SubjectID, colour = Group), alpha = .18) +
    geom_line(data = s, aes(Session, m, group = Group, colour = Group), linewidth = .8, position = position_dodge(.15)) +
    geom_errorbar(data = s, aes(Session, ymin = m - se, ymax = m + se, colour = Group), width = .1, position = position_dodge(.15)) +
    geom_point(data = s, aes(Session, m, colour = Group), size = PT, position = position_dodge(.15)) +
    scale_colour_manual(values = COL, name = NULL) +
    labs(x = NULL, y = ylab, title = ttl, subtitle = sprintf("%s: d=%+.2f, p=%.3f (n=%d/%d)", sub, a["d"], a["p"], a["nP"], a["nS"])) + th }
g3 <- (prepost(wide("PRL_Perseveration"), "Perseveration parameter", "A  Choice perseveration (PRL)") +
       prepost(wide("CF_SwitchCost_transfer"), "Accuracy drop after rule change", "B  Switch cost, feedback-driven") +
       prepost(wide("ES_Learn_SwitchCost_RT"), "Switch − repeat RT (s)", "C  Switch cost, cued (demand selection)") +
       prepost(wide("ES_Learn_RT"), "Learning-phase RT (s)", "D  Response time (demand selection)")) +
  plot_layout(ncol = 2, guides = "collect") + bottom
save("Figure4_change_by_arm.png", g3, 7, 6)

# ---- Figure 3: stability-flexibility map --------------------------------------------------------------
S <- rbind(COMP[, c("var", "label", "task", "pole", "d", "lo", "hi", "p", "d_stab", "cpl_n", "cpl_rho", "cpl_p", "cpl_stab")],
           R[R$pole != 0, c("var", "label", "task", "pole", "d", "lo", "hi", "p", "d_stab", "cpl_n", "cpl_rho", "cpl_p", "cpl_stab")])
S$family <- factor(ifelse(S$task == "Composite", "Stability\ncomposites", ifelse(S$task %in% c("PRL", "Rule switching"), "Feedback-driven\n(PRL, rule switching)",
                   ifelse(S$task == "Demand selection", "Cued task switching\n(demand selection)", "Interference\n(MSIT)"))),
                   levels = c("Stability\ncomposites", "Feedback-driven\n(PRL, rule switching)", "Cued task switching\n(demand selection)", "Interference\n(MSIT)"))
S$lab3 <- factor(sprintf("%s%s", S$label, ifelse(S$task == "Composite", "", ifelse(S$pole == 1, " [+S]", " [+F]"))),
                 levels = rev(sprintf("%s%s", S$label, ifelse(S$task == "Composite", "", ifelse(S$pole == 1, " [+S]", " [+F]")))))
S$dlo <- ifelse(S$pole == 1, S$lo, -S$hi); S$dhi <- ifelse(S$pole == 1, S$hi, -S$lo)
# approximate 95% CI for Spearman rho: Fisher z with the Fieller et al. (1957) variance 1.06/(n-3)
fz <- function(r, n, s) tanh(atanh(r) + s * 1.96 * sqrt(1.06 / (n - 3)))
S$c1 <- -S$pole * fz(S$cpl_rho, S$cpl_n, 1); S$c2 <- -S$pole * fz(S$cpl_rho, S$cpl_n, -1)
S$clo <- pmin(S$c1, S$c2); S$chi <- pmax(S$c1, S$c2)
bg <- list(annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#FAECE7", alpha = .6),
           annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#E1F5EE", alpha = .6), geom_vline(xintercept = 0, colour = "grey50"))
pa <- ggplot(S, aes(d_stab, lab3)) + bg + geom_errorbar(aes(xmin = dlo, xmax = dhi), width = .25, orientation = "y", colour = "grey40") +
  geom_point(aes(shape = p < .05), fill = "#9b27af", size = PT) + scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = expression(atop("Psilocybin vs placebo, " * italic(d) * " (95% CI)", "\u2190 flexibility      stability \u2192")), y = NULL, title = "A  Drug effect at post-dose testing") +
  coord_cartesian(xlim = c(-2.6, 2.6)) + th +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 7.5), strip.background = element_blank(), panel.spacing = unit(5, "pt"))
pb <- ggplot(S, aes(cpl_stab, lab3)) + bg + geom_errorbar(aes(xmin = clo, xmax = chi), width = .25, orientation = "y", colour = "grey40") +
  geom_point(aes(shape = cpl_p < .05), fill = "#0F6E56", size = PT) + scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression(atop(italic(r)[s] * " with drinking reduction (95% CI)", "\u2190 flexibility      stability \u2192")), y = NULL,
       title = "B  Correlation with drinking reduction") +
  coord_cartesian(xlim = c(-.95, .95)) + th +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(), strip.background = element_blank(),
        strip.placement = "outside", panel.spacing = unit(5, "pt"), plot.margin = margin(5.5, 14, 5.5, 5.5))
save("Figure3_stability_flexibility_map.png", pa + pb + plot_layout(widths = c(1, 1)), 7.2, 5.6)
write.csv(S[, c("var", "label", "task", "pole", "d", "lo", "hi", "p", "d_stab", "cpl_n", "cpl_rho", "cpl_p", "cpl_stab")], file.path(FIGURES, "Figure3_data.csv"), row.names = FALSE)

# ---- Figure 5: demand avoidance ------------------------------------------------------------------------
LA <- add_measure(L, demand_avoidance())
w5 <- merge(wide("DemandAvoid", LA), osp, by = "SubjectID", all.x = TRUE); w5$Group2 <- ifelse(is.na(w5$Group), "Unresolved", w5$Group)
sc5 <- function(x, y, xl, yl, ttl) { d <- w5[complete.cases(w5[, c(x, y)]), ]; s <- spear(d[[x]], d[[y]])
  ggplot(d, aes(.data[[x]], .data[[y]])) + geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey85", linewidth = .6) +
    geom_point(aes(colour = Group2), size = PT) + scale_colour_manual(values = COL, name = NULL) +
    labs(x = xl, y = yl, title = ttl, subtitle = rs("", s[1], sprintf(", p = %s, n = %d", nz(s[2], 3), s[3]))) + th }
g5 <- (sc5("PHQ_v6", "pre", "PHQ at baseline visit", "Demand avoidance (baseline)", "A  Avoidance and depressive symptoms") +
       sc5("OSpan", "pre", "Operation span (baseline)", "Demand avoidance (baseline)", "B  Avoidance and operation span") +
       (prepost(w5, "Demand avoidance", "C  Avoidance by arm", hline = .5, sub = "Baseline-adjusted") + guides(colour = "none"))) +
  plot_layout(guides = "collect") + bottom
save("Figure5_demand_avoidance.png", g5, 7.2, 3.1)

# ---- Figure 6: coupling with drinking -------------------------------------------------------------------
scat <- function(d, x, y, xl, yl, ttl, covs = "Group") {
  d <- d[complete.cases(d[, c(x, y)]), ]; r0 <- rankAssoc(y, x, d); r1 <- rankAssoc(y, x, d[!is.na(d$Group), ], covs)
  d$Group2 <- ifelse(is.na(d$Group), "Unresolved", d$Group)
  ggplot(d, aes(.data[[x]], .data[[y]])) + geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey85", linewidth = .6) +
    geom_smooth(data = d[!is.na(d$Group), ], aes(colour = Group), method = "lm", formula = y ~ x, se = FALSE, linetype = 2, linewidth = .45) +
    geom_point(aes(colour = Group2), size = PT * .9) + scale_colour_manual(values = COL, name = NULL) +
    labs(x = xl, y = yl, title = ttl, subtitle = rs("", r0[4], sprintf(", p = %s, n = %d; arm-adjusted p = %s", nz(r0[5], 3), r0[1], nz(r1[3], 3)))) + th }
S2 <- read.csv(file.path(RESULTS, "stability_composite.csv"))
g6 <- (scat(S2, "stab", "pct_drink_chg1", "Stability composite change (z)", "Change in drinking days", "A  Stability composite change") +
       scat(wide("PRL_RewardAcc"), "change", "pct_drink_chg1", "Change in reward-learning accuracy", "Change in drinking days", "B  Reward-learning accuracy change") +
       scat(wide("PRL_Perseveration"), "post", "pct_heavy_v1", "Perseveration after drug administration", "Heavy-drinking days, 30 d after dosing\n(proportion, not a change score)", "C  Post-dose perseveration", covs = c("Group", "pre"))) +
  plot_layout(guides = "collect") + bottom
save("Figure6_coupling_with_drinking.png", g6, 7.2, 3.3)

# ---- Figure 7: operation span and drinking change ---------------------------------------------------------
w7 <- merge(merge(osp, dose, by = "SubjectID"), dr, by = "SubjectID"); w7$pre <- w7$OSpan
w7$crit <- ifelse(w7$OSpan_math >= 85, "Processing accuracy ≥85%", "Below 85%")
wm <- function(y, b, yl, ttl) { r <- rankAssoc(y, "pre", w7, c("Group", b))
  ggplot(w7, aes(OSpan, .data[[y]])) + geom_hline(yintercept = 0, colour = "grey80") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey85", linewidth = .6) +
    geom_smooth(aes(colour = Group), method = "lm", formula = y ~ x, se = FALSE, linetype = 2, linewidth = .45) +
    geom_point(aes(colour = Group, shape = crit), size = PT) + scale_colour_manual(values = COL, name = NULL) +
    scale_shape_manual(values = c(`Processing accuracy ≥85%` = 16, `Below 85%` = 1), name = NULL) +
    labs(x = "Baseline operation span (partial-credit score)", y = yl, title = ttl,
         subtitle = sprintf("Rank-partial β=%+.2f, p=%.4f, n=%d (adjusted for arm and baseline)", r[2], r[3], r[1])) + th }
g7 <- (wm("pct_drink_chg1", "pct_drink_v0", "Change in drinking days", "A  Drinking days") +
       wm("pct_heavy_chg1", "pct_heavy_v0", "Change in heavy-drinking days", "B  Heavy-drinking days")) +
  plot_layout(guides = "collect") + bottom
save("Figure7_wmc_predicts_drinking_change.png", g7, 7, 3.6)
cat("figures written to", FIGURES, "\n")

# ---- Figure 1: the schematic is drawn by hand (assets/Figure1_framework_and_findings.svg); export it at the same
# printed width as the other figures when rsvg-convert is available
svg <- file.path("..", "assets", "Figure1_framework_and_findings.svg"); rc <- Sys.which("rsvg-convert")
if (nzchar(rc) && file.exists(svg)) {
  f1 <- file.path(FIGURES, "Figure1_framework_and_findings")
  system2(rc, c("-d", 600, "-p", 600, "-w", 7 * 600, "-b", "white", "-o", shQuote(paste0(f1, ".png")), shQuote(svg)))
  system2(rc, c("-f", "pdf", "-w", 7 * 72, "-o", shQuote(paste0(f1, ".pdf")), shQuote(svg)))
} else message("rsvg-convert not found: Figure 1 not exported (use the SVG in assets/)")
