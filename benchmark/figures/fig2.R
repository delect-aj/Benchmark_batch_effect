# Figure 2 -- Q1: do the evaluation criteria agree with each other?
#
#   a  one dataset before and after correction, coloured by batch and by phenotype
#   b  every method against every metric, as a rank averaged over the evaluation datasets
#   c  the same methods ranked under three distances: does the order survive the choice?
#   d  how the metrics correlate with each other across all method x dataset results
#   e  batch removal against signal retention, with the undominated methods marked
#
# Usage: Rscript figures/fig2.R <results dir> <out.pdf>
source("figures/_panels.R")

batch_metrics <- function() c("R2_batch_ait", "R2_batch_bc", "R2_batch_jac", "ASW_batch", "iLISI")
bio_metrics <- function() c("R2_pheno_ait", "R2_pheno_bc", "R2_pheno_jac", "ASW_pheno", "cLISI")

# Principal coordinates of one corrected table, in the space its own metrics are computed in.
ordination <- function(results_dir, dataset, method, meta) {
  produced <- read_output(results_dir, dataset, method)
  coordinates <- cmdscale(dist(euclidean_space(produced$values, produced$kind)), k = 2)
  data.frame(method = method, dim1 = coordinates[, 1], dim2 = coordinates[, 2],
             batch = meta[rownames(coordinates), "batch"],
             phenotype = meta[rownames(coordinates), "phenotype"])
}

panel_ordination <- function(results_dir, dataset, methods, methods_table = method_table()) {
  meta <- read.delim(file.path(results_dir, "data", dataset, "meta.tsv"), row.names = 1)
  points <- label_methods(bind_rows(lapply(methods, function(m) ordination(results_dir, dataset, m, meta))),
                          methods_table)
  # The two rows colour by variables with different numbers of levels, so each is spread over the whole
  # palette by its own index rather than sharing one discrete scale that would squash the phenotypes.
  spread <- function(values) { f <- factor(values); (as.integer(f) - 1) / max(1, nlevels(f) - 1) }
  coloured <- bind_rows(transform(points, colour = spread(batch), facet = "by batch"),
                        transform(points, colour = spread(phenotype), facet = "by phenotype"))
  ggplot(coloured, aes(dim1, dim2, colour = colour)) +
    geom_point(size = 0.35, alpha = 0.7) +
    facet_grid(facet ~ label, scales = "free") +
    scale_colour_viridis_c(guide = "none", option = "turbo") +
    labs(x = "PCo 1", y = "PCo 2", title = paste("Ordination of", sub(".*/", "", dataset))) +
    theme_bench() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

panel_global <- function(scored, methods_table = method_table()) {
  metrics <- c(batch_metrics(), bio_metrics())
  scores <- bind_rows(lapply(metrics, function(m) scaled_score(scored, m))) %>%
    group_by(method, metric, family) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  scores$metric <- factor(scores$metric, metrics,
                          metric_catalog()$label[match(metrics, metric_catalog()$metric)])
  ggplot(label_methods(scores, methods_table), aes(metric, label)) +
    geom_tile(aes(fill = score), colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.na(score), "", sprintf("%.2f", score))), size = 1.6) +
    facet_grid(. ~ family, scales = "free_x", space = "free_x",
               labeller = as_labeller(c(batch = "batch removal", bio = "signal retention"))) +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1), name = "scaled\nscore") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, title = "Mean scaled score across the evaluation datasets") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
}

panel_flip <- function(scored, methods_table = method_table()) {
  distances <- c(R2_batch_bc = "Bray-Curtis", R2_batch_ait = "Aitchison", R2_batch_jac = "Jaccard")
  orders <- bind_rows(lapply(names(distances), function(m) ranked(scored, m))) %>%
    group_by(method, metric) %>%
    summarise(rank = mean(rank, na.rm = TRUE), .groups = "drop")
  orders$distance <- factor(distances[orders$metric], distances)
  ggplot(label_methods(orders, methods_table), aes(distance, rank, group = label, colour = track)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1) +
    ggrepel::geom_text_repel(data = ~ subset(.x, distance == tail(levels(distance), 1)),
                             aes(label = label), size = 1.8, direction = "y", hjust = 0,
                             nudge_x = 0.15, segment.size = 0.2, max.overlaps = 20) +
    scale_y_reverse(breaks = seq(2, 16, 2)) +
    scale_colour_manual(values = track_colours(), guide = "none") +
    expand_limits(x = 4.2) +
    labs(x = NULL, y = "mean rank", title = "Does the ranking survive the distance?") +
    theme_bench()
}

panel_metric_agreement <- function(scored) {
  metrics <- c(batch_metrics(), bio_metrics())
  values <- as.matrix(scored[scored$status == "ok", metrics])
  correlations <- cor(values, use = "pairwise.complete.obs", method = "spearman")
  long <- data.frame(expand.grid(x = metrics, y = metrics), rho = as.vector(correlations))
  labels <- metric_catalog()$label[match(metrics, metric_catalog()$metric)]
  long$x <- factor(long$x, metrics, labels); long$y <- factor(long$y, metrics, labels)
  ggplot(long, aes(x, y, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", limits = c(-1, 1),
                         name = "Spearman") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, title = "Agreement between metrics") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
}

# Methods that nothing else beats on both axes at once.
pareto_front <- function(points) {
  dominated <- vapply(seq_len(nrow(points)), function(i)
    any(points$batch >= points$batch[i] & points$bio >= points$bio[i] &
        (points$batch > points$batch[i] | points$bio > points$bio[i])), TRUE)
  points[!dominated, ][order(points$batch[!dominated]), ]
}

panel_tradeoff <- function(scored, methods_table = method_table()) {
  points <- full_join(family_score(scored, batch_metrics()) %>% select(method, batch = score),
                      family_score(scored, bio_metrics()) %>% select(method, bio = score),
                      by = "method")
  front <- pareto_front(points)
  ggplot(label_methods(points, methods_table), aes(batch, bio)) +
    geom_line(data = label_methods(front, methods_table), colour = "grey60", linewidth = 0.4) +
    geom_point(aes(colour = track), size = 1.3) +
    ggrepel::geom_text_repel(aes(label = label), size = 1.8, segment.size = 0.2, max.overlaps = 20) +
    scale_colour_manual(values = track_colours(), name = NULL) +
    labs(x = "batch removal", y = "signal retention", title = "The trade-off") +
    theme_bench()
}

main <- function(args) {
  results_dir <- args[1]; out <- args[2]
  scored <- read_scores(file.path(results_dir, "scores.tsv"))
  evaluated <- scored[scored$regime == "real" & scored$status == "ok", ]
  shown <- intersect(c("raw", "limma", "conqur", "mmuphin", "harmony"), unique(evaluated$method))

  figure <- (panel_ordination(results_dir, "real/hivrc_16s", shown) /
             (panel_global(evaluated) | (panel_flip(evaluated) / panel_metric_agreement(evaluated) /
                                         panel_tradeoff(evaluated)))) +
    plot_layout(heights = c(1, 3.2)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  save_figure(figure, out, 180, 200)
}

main(commandArgs(TRUE))
