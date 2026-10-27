# Figure 1 -- benchmark design and data landscape. No results, only what we put in.
#
#   a  the three data regimes and what each one can score
#   b  the datasets: size, number of batches, assay
#   c  the confounding axis, with the real datasets placed on the simulated gradient
#   d  the simulation parameter space actually sampled
#   e  the method inventory and what each method needs
#
# Usage: Rscript figures/fig1.R <results dir> <out.pdf>
source("figures/_panels.R")
suppressMessages(library(ggrepel))

dataset_confounding_table <- function(datasets, results_dir)
  dataset_confounding(results_dir, sub("^results/data/", "", datasets$path))

panel_design <- function() {
  boxes <- data.frame(
    x = c(1, 2, 3, 2, 1, 2, 3, 4), y = c(3, 3, 3, 2, 1, 1, 1, 1),
    fill = c(rep("regime", 3), "interface", rep("metric", 4)),
    text = c("curated real data\n10 sets, 7,115 samples\nproxy metrics only",
             "designed protocol runs\n4 sets, known composition\nor repeated specimens",
             "simulation\n213 scenarios\nfull ground truth",
             "one wrapper interface\ncounts + batch + phenotype  ->  corrected table or embedding",
             "batch removal\nR2, ASW, LISI", "signal retention\nR2, cLISI, oracle",
             "downstream\nDA, cross-study", "practicality\ntime, memory, failures"))
  boxes$width <- ifelse(boxes$fill == "interface", 2.94, 0.92)
  ggplot(boxes, aes(x, y)) +
    geom_tile(aes(fill = fill, width = width), height = 0.82, colour = "grey30", linewidth = 0.3) +
    geom_text(aes(label = text), size = 1.9, lineheight = 1.1) +
    annotate("segment", x = c(1, 2, 3), xend = c(1.9, 2, 2.1), y = 2.58, yend = 2.42,
             arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.3) +
    annotate("segment", x = 2, xend = c(1, 2, 3, 4), y = 1.58, yend = 1.42,
             arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.3) +
    scale_fill_manual(values = c(regime = "#deebf7", interface = "grey90", metric = "#fee6ce"),
                      guide = "none") +
    scale_x_continuous(expand = expansion(0.02)) +
    coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0.5, 3.5)) +
    labs(title = "Three data regimes, one interface, four metric families") +
    theme_void(base_size = 7) + theme(plot.title = element_text(size = 7, face = "bold"))
}

panel_landscape <- function(datasets) {
  ggplot(datasets, aes(samples, n_batch)) +
    geom_point(aes(size = features, shape = assay, colour = evaluated), alpha = 0.85) +
    geom_text_repel(aes(label = dataset), size = 1.9, min.segment.length = 0, max.overlaps = Inf,
                    force = 12, box.padding = 0.45, point.padding = 0.25, seed = 1,
                    segment.size = 0.2, segment.colour = "grey50") +
    scale_x_log10() + scale_size_area(max_size = 5, breaks = c(10, 100, 700)) +
    scale_colour_manual(values = c(yes = "#2c7fb8", no = "#d95f0e"),
                        labels = c(yes = "evaluation", no = "calibration")) +
    labs(x = "samples (log scale)", y = "batches", size = "features",
         shape = "assay", colour = "role", title = "Dataset landscape") +
    guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2),
           size = guide_legend(order = 3)) +
    theme_bench() +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 6),
          legend.spacing.x = unit(1, "mm"))
}

panel_confounding <- function(datasets, confounding, simulated_levels) {
  placed <- data.frame(dataset = datasets$dataset, v = confounding, assay = datasets$assay,
                       evaluated = datasets$evaluated)
  placed$dataset <- factor(placed$dataset, placed$dataset[order(placed$v)])
  ggplot(placed, aes(v, dataset)) +
    geom_vline(xintercept = simulated_levels, colour = "grey70", linetype = "22", linewidth = 0.3) +
    annotate("text", x = simulated_levels, y = nlevels(placed$dataset) + 1.1,
             label = format(simulated_levels), size = 1.9, colour = "grey40") +
    annotate("text", x = 0, y = nlevels(placed$dataset) + 2.3, label = "simulated gradient",
             size = 2, colour = "grey40", hjust = 0) +
    geom_segment(aes(x = 0, xend = v, yend = dataset), colour = "grey60", linewidth = 0.3) +
    geom_point(aes(colour = assay), size = 1.6) +
    scale_colour_brewer(palette = "Dark2") +
    # declared so the numeric y of the two annotations does not type the axis as continuous first
    scale_y_discrete(expand = expansion(add = c(0.6, 3.2))) +
    scale_x_continuous(limits = c(0, 1.02), expand = expansion(0.01)) +
    labs(x = "batch-phenotype confounding (Cramer's V)", y = NULL, colour = NULL,
         title = "Where the real datasets sit on the confounding axis") +
    theme_bench() + theme(panel.grid.major.y = element_blank())
}

panel_simulation_space <- function(scenarios) {
  axes <- c(conf = "confounding (V)", n_batch = "batches", n_per = "samples per batch",
            imbalance = "batch size ratio", bias_sd = "multiplicative bias (sd)",
            affected = "taxa hit by bias", scale_sd = "variance shift (sd)",
            dropout = "extra dropout", depth_fold = "depth ratio",
            da_prop = "DA taxa", da_log2fc = "DA effect (log2)")
  sampled <- bind_rows(lapply(names(axes), function(k)
    data.frame(parameter = axes[[k]], value = scenarios[[k]])))
  sampled$parameter <- factor(sampled$parameter, axes)
  ggplot(sampled, aes(value, 1)) +
    geom_point(size = 0.5, alpha = 0.25, position = position_jitter(height = 0.3)) +
    facet_wrap(~ parameter, ncol = 4, scales = "free_x") +
    scale_x_continuous(n.breaks = 4) +
    labs(x = NULL, y = NULL, title = "Simulation parameter space sampled (213 scenarios)") +
    theme_bench() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())
}

panel_inventory <- function(methods) {
  properties <- c(needs_phenotype = "needs phenotype", needs_controls = "needs per-batch controls",
                  binary_only = "binary phenotype only", transductive = "fits on all samples")
  long <- bind_rows(lapply(names(properties), function(k)
    data.frame(method = methods$method, label = methods$label, track = methods$track,
               property = properties[[k]], value = methods[[k]])))
  long$label <- factor(long$label, rev(methods$label[match(method_order(methods), methods$method)]))
  long$property <- factor(long$property, properties)
  output <- data.frame(label = long$label[!duplicated(long$label)], property = "output",
                       value = methods$output[match(levels(long$label)[as.integer(long$label[!duplicated(long$label)])],
                                                    methods$label)])
  ggplot(long, aes(property, label)) +
    geom_tile(aes(fill = value), colour = "white", linewidth = 0.4) +
    geom_text(data = output, aes(label = value), size = 1.8) +
    scale_fill_manual(values = c(yes = "#2c7fb8", no = "grey92", na = "grey75"), name = NULL,
                      breaks = c("yes", "no", "na"),
                      labels = c(yes = "yes", no = "no", na = "not applicable")) +
    labs(x = NULL, y = NULL, title = "Method inventory") +
    theme_bench() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
          legend.position = "bottom")
}

main <- function(args) {
  results_dir <- args[1]; out <- args[2]
  datasets <- read.delim("data/DATASETS.tsv")
  datasets <- datasets[datasets$track %in% c("real", "mock"), ]
  scenarios <- read.delim("config/scenarios.tsv")
  methods <- method_table()
  methods <- methods[methods$role != "sensitivity", ]

  figure <- (panel_design() /
             (panel_landscape(datasets) | panel_confounding(datasets, dataset_confounding_table(datasets, results_dir),
                                                            c(0, 0.3, 0.6, 0.9, 1))) /
             panel_simulation_space(scenarios) /
             panel_inventory(methods)) +
    plot_layout(heights = c(1, 1.1, 1.1, 1.3)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
  save_figure(figure, out, 180, 230)
}

main(commandArgs(TRUE))
