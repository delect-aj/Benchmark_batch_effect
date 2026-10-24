# Figure 4 -- Q3: how far does batch-phenotype confounding carry a correction before it breaks?
#
#   a  performance along the simulated confounding gradient, both metric families
#   b  the overcorrection index: biological signal lost per unit of batch signal removed
#   c  with no biological signal at all, how much phenotype structure does a method invent
#   d  the same measures on real data, at the two ends of the axis
#   e  where each method crosses the point of keeping less than half the true signal
#
# Usage: Rscript figures/fig4.R <results dir> <out.pdf>
source("figures/_panels.R")

batch_metrics <- function() c("R2_batch_ait", "R2_batch_bc", "R2_batch_jac", "ASW_batch", "iLISI")
bio_metrics <- function() c("R2_pheno_ait", "R2_pheno_bc", "R2_pheno_jac", "ASW_pheno", "cLISI")

# Per dataset and method, the scaled score of each metric family, with the scenario parameters kept.
scored_by_scenario <- function(scored) {
  parameters <- scored[!duplicated(scored$dataset), c("dataset", "scenario", "family", "cramers_v",
                                                      "conf", "template", "da_prop", "affected")]
  families <- bind_rows(
    family_score(scored, batch_metrics(), by = c("dataset", "method")),
    family_score(scored, bio_metrics(), by = c("dataset", "method")))
  left_join(families, parameters, by = "dataset", suffix = c("", "_scenario"))
}

panel_gradient <- function(sweep_scores, methods_table = method_table()) {
  curves <- sweep_scores %>%
    group_by(method, family, template, cramers_v) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  ggplot(label_methods(curves, methods_table), aes(cramers_v, score, colour = label, group = label)) +
    geom_line(linewidth = 0.4) + geom_point(size = 0.7) +
    facet_grid(family ~ template, labeller = as_labeller(
      c(batch = "batch removal", bio = "signal retention", `16s.tsv` = "16S", `mgx.tsv` = "metagenome",
        `16s` = "16S", mgx = "metagenome"))) +
    scale_colour_viridis_d(option = "turbo", name = NULL) +
    labs(x = "batch-phenotype confounding (Cramer's V)", y = "scaled score",
         title = "Performance along the confounding gradient") +
    theme_bench() + theme(legend.position = "right")
}

# How much biological signal a method gives up for each unit of batch signal it removes. Taken as a
# median over scenarios: the ratio is unstable wherever a method barely moved the batch term.
panel_overcorrection <- function(simulated, methods_table = method_table()) {
  reference <- simulated[simulated$method == "raw", c("dataset", "R2_batch_ait", "R2_pheno_ait")]
  names(reference)[2:3] <- c("raw_batch", "raw_pheno")
  # oracle has no batch term to remove, so the ratio is not defined for it
  corrected <- inner_join(simulated[!simulated$method %in% c("raw", "oracle"), ], reference, by = "dataset")
  corrected$index <- (corrected$raw_pheno - corrected$R2_pheno_ait) /
                     (corrected$raw_batch - corrected$R2_batch_ait)
  summarised <- corrected %>%
    group_by(method, cramers_v) %>%
    summarise(index = median(index, na.rm = TRUE), .groups = "drop")
  ggplot(label_methods(summarised, methods_table), aes(cramers_v, index, colour = label)) +
    geom_hline(yintercept = 1, linetype = "22", linewidth = 0.3, colour = "grey50") +
    geom_line(linewidth = 0.4) + geom_point(size = 0.7) +
    scale_colour_viridis_d(option = "turbo", guide = "none") +
    scale_y_log10(labels = scales::label_number()) +
    coord_cartesian(ylim = c(0.01, 10)) +
    labs(x = "Cramer's V", y = "signal lost / batch removed",
         title = "Overcorrection index") +
    theme_bench()
}

panel_invented_signal <- function(null_scores, methods_table = method_table()) {
  ggplot(label_methods(null_scores, methods_table), aes(factor(round(cramers_v, 1)), R2_pheno_ait,
                                                        fill = track)) +
    geom_col(position = position_dodge(0.85), width = 0.8) +
    facet_wrap(~ label, nrow = 2) +
    scale_fill_manual(values = track_colours(), guide = "none") +
    labs(x = "Cramer's V", y = "R2 of a signal-free phenotype",
         title = "Invented structure in the null scenarios") +
    theme_bench()
}

panel_real_anchors <- function(real_scores, confounding, methods_table = method_table()) {
  anchored <- real_scores
  anchored$v <- confounding[sub("^real/", "", anchored$dataset)]
  anchored$end <- ifelse(anchored$v < 0.2, "balanced designs (V < 0.2)", "confounded designs (V >= 0.2)")
  summarised <- anchored %>%
    group_by(method, family, end) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  ggplot(label_methods(summarised, methods_table), aes(label, score, fill = family)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    facet_wrap(~ end, nrow = 1) +
    scale_fill_brewer(palette = "Set2", name = NULL,
                      labels = c(batch = "batch removal", bio = "signal retention")) +
    labs(x = NULL, y = "scaled score", title = "The same measures on real data") +
    theme_bench() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# The point on the gradient where a method keeps less than half of the signal the oracle still shows.
panel_breakpoint <- function(sweep_simulated, methods_table = method_table()) {
  oracle <- sweep_simulated[sweep_simulated$method == "oracle", c("dataset", "R2_pheno_ait")]
  names(oracle)[2] <- "oracle_pheno"
  retained <- inner_join(sweep_simulated[sweep_simulated$method != "oracle", ], oracle, by = "dataset")
  retained$retention <- retained$R2_pheno_ait / retained$oracle_pheno
  breaks <- retained %>%
    group_by(method, cramers_v) %>%
    summarise(retention = mean(retention, na.rm = TRUE), .groups = "drop") %>%
    group_by(method) %>%
    summarise(breakpoint = suppressWarnings(min(cramers_v[retention < 0.5])), .groups = "drop")
  breaks$breakpoint[!is.finite(breaks$breakpoint)] <- 1.1   # never fell below half across the gradient
  shown <- label_methods(breaks, methods_table)
  shown$label <- reorder(shown$label, -shown$breakpoint)
  ggplot(shown, aes(breakpoint, label, fill = track)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 1.05, linetype = "22", linewidth = 0.3, colour = "grey50") +
    scale_fill_manual(values = track_colours(), guide = "none") +
    scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.1), labels = c("0", "0.3", "0.6", "0.9", "never")) +
    labs(x = "Cramer's V where less than half the signal survives", y = NULL,
         title = "Breaking point") +
    theme_bench()
}

main <- function(args) {
  results_dir <- args[1]; out <- args[2]
  scored <- read_scores(file.path(results_dir, "scores.tsv"))
  simulated <- scored[scored$regime == "sim" & scored$status == "ok", ]
  by_scenario <- scored_by_scenario(simulated)
  real_by_dataset <- bind_rows(
    family_score(scored[scored$regime == "real" & scored$status == "ok", ], batch_metrics(),
                 by = c("dataset", "method")),
    family_score(scored[scored$regime == "real" & scored$status == "ok", ], bio_metrics(),
                 by = c("dataset", "method")))
  confounding <- dataset_confounding(results_dir, unique(scored$dataset[scored$regime == "real"]))
  names(confounding) <- sub("^real/", "", unique(scored$dataset[scored$regime == "real"]))

  figure <- (panel_gradient(by_scenario[by_scenario$family_scenario == "sweep", ]) /
             (panel_overcorrection(simulated[simulated$family == "sweep", ]) |
              panel_breakpoint(simulated[simulated$family == "sweep", ])) /
             panel_invented_signal(simulated[simulated$family == "null" & simulated$method != "oracle", ]) /
             panel_real_anchors(real_by_dataset, confounding)) +
    plot_layout(heights = c(1.4, 1, 1.2, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  save_figure(figure, out, 180, 215)
}

main(commandArgs(TRUE))
