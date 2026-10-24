# Figure 3 -- Q2: with the truth in hand, what did the methods actually do?
#
#   a  per-taxon bias against the composition that was pipetted, method by method
#   b  distance to that composition, for every designed mock dataset
#   c  false structure: one material, many protocols -- what is left between batches
#   d  fidelity to the batch-free simulated truth
#   e  systematic against non-systematic bias: does the multiplicative assumption decide the outcome
#   f  does the proxy metric move with the truth?
#
# Usage: Rscript figures/fig3.R <results dir> <out.pdf>
source("figures/_panels.R")

short_dataset <- function(paths) sub("jmbc_", "jmbc ", sub(".*/", "", paths))

# Each method's value minus the value `raw` got on the same dataset: the improvement correction bought.
change_from_raw <- function(scored, metric) {
  reference <- scored[scored$method == "raw", c("dataset", metric)]
  names(reference)[2] <- "raw_value"
  joined <- inner_join(scored[scored$method != "raw", c("dataset", "method", "scenario", metric)],
                       reference, by = "dataset")
  joined$change <- joined[[metric]] - joined$raw_value
  joined
}

panel_taxon_bias <- function(per_taxon, dataset, methods_table = method_table()) {
  shown <- per_taxon[per_taxon$dataset == dataset & !is.na(per_taxon$log2_ratio), ]
  shown$taxon <- factor(shown$taxon, unique(shown$taxon[order(shown$log2_ratio[shown$method == "raw"])]))
  ggplot(label_methods(shown, methods_table), aes(taxon, label, fill = log2_ratio)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-3, 3), oob = scales::squish, name = "log2\nobs/truth") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, title = paste("Per-taxon bias,", sub(".*/", "", dataset),
                                           "(taxa ordered by the bias before correction)")) +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 4.5), panel.grid = element_blank())
}

panel_distance_to_truth <- function(truth_summary, methods_table = method_table()) {
  shown <- label_methods(truth_summary, methods_table)
  shown$dataset <- short_dataset(shown$dataset)
  ggplot(shown, aes(aitchison_to_truth, label, fill = track)) +
    geom_col(width = 0.7) +
    facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = track_colours(), guide = "none") +
    scale_y_discrete(limits = rev) +
    labs(x = "Aitchison distance to truth", y = NULL, title = "Accuracy against the pipetted truth") +
    theme_bench()
}

panel_false_structure <- function(truth_summary, methods_table = method_table()) {
  shown <- label_methods(truth_summary, methods_table)
  shown$dataset <- short_dataset(shown$dataset)
  shown$false_da_fraction <- shown$false_da_n / shown$n_taxa
  ggplot(shown, aes(y = label)) +
    geom_col(aes(x = R2_batch, fill = track), width = 0.7) +
    geom_point(aes(x = false_da_fraction), shape = 21, fill = "white", size = 1.2) +
    facet_wrap(~ dataset, nrow = 1) +
    scale_fill_manual(values = track_colours(), guide = "none") +
    scale_y_discrete(limits = rev) +
    labs(x = "batch R2 among identical material", y = NULL,
         title = "False structure between protocols (points: taxa still differing)") +
    theme_bench()
}

panel_oracle_fidelity <- function(simulated, methods_table = method_table()) {
  long <- bind_rows(
    data.frame(method = simulated$method, value = simulated$oracle_mantel, measure = "distances (Mantel r)"),
    data.frame(method = simulated$method, value = simulated$oracle_taxon_rho, measure = "per taxon (rho)"))
  ggplot(label_methods(long[!is.na(long$value), ], methods_table), aes(value, label, fill = track)) +
    geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.7) +
    facet_wrap(~ measure, nrow = 1) +
    scale_fill_manual(values = track_colours(), guide = "none") +
    scale_y_discrete(limits = rev) +
    labs(x = "agreement with the batch-free truth", y = NULL,
         title = "Fidelity to the simulated truth") +
    theme_bench()
}

panel_systematic <- function(simulated, methods_table = method_table()) {
  shown <- simulated[!is.na(simulated$oracle_taxon_rho), ]
  shown$bias_kind <- cut(shown$affected, c(0, 0.4, 0.8, 1),
                         labels = c("10-40% of taxa", "40-80%", "80-100% (systematic)"))
  ggplot(label_methods(shown, methods_table), aes(oracle_taxon_rho, label, fill = bias_kind)) +
    stat_summary(fun = median, geom = "col", orientation = "y",
                 position = position_dodge(0.8), width = 0.75) +
    scale_fill_brewer(palette = "BuPu", name = "taxa hit by the bias") +
    scale_y_discrete(limits = rev) +
    labs(x = "median per-taxon rho", y = NULL,
         title = "Does the multiplicative assumption decide it?") +
    theme_bench()
}

panel_proxy_versus_truth <- function(simulated, methods_table = method_table()) {
  batch <- change_from_raw(simulated, "R2_batch_ait")
  truth <- change_from_raw(simulated, "oracle_mantel")
  joined <- inner_join(batch %>% select(dataset, method, batch_change = change),
                       truth %>% select(dataset, method, truth_change = change),
                       by = c("dataset", "method"))
  ggplot(label_methods(joined, methods_table), aes(-batch_change, truth_change, colour = track)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
    geom_point(size = 0.3, alpha = 0.35) +
    scale_colour_manual(values = track_colours(), name = NULL) +
    labs(x = "batch R2 removed (vs raw)", y = "fidelity gained (vs raw)",
         title = "Proxy against truth, one point per scenario") +
    theme_bench()
}

main <- function(args) {
  results_dir <- args[1]; out <- args[2]
  scored <- read_scores(file.path(results_dir, "scores.tsv"))
  simulated <- scored[scored$regime == "sim" & scored$status == "ok" & scored$method != "oracle", ]
  truth_summary <- read.delim(file.path(results_dir, "truth_summary.tsv"))
  truth_summary <- truth_summary[!is.na(truth_summary$aitchison_to_truth), ]
  per_taxon <- read.delim(file.path(results_dir, "truth_taxon.tsv"))

  figure <- (panel_taxon_bias(per_taxon, "real/jmbc_extraction") /
             (panel_distance_to_truth(truth_summary) | panel_false_structure(truth_summary)) /
             (panel_oracle_fidelity(simulated) | panel_systematic(simulated)) /
             panel_proxy_versus_truth(simulated)) +
    plot_layout(heights = c(1, 1.1, 1.1, 0.8)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  save_figure(figure, out, 180, 235)
}

main(commandArgs(TRUE))
