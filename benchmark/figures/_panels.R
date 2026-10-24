# Shared pieces for the five main figures: the method inventory, the metric catalogue, the two
# normalisations used for ranking, and one place where figure size and type are set.
#
# Every figure script sources this file, defines one function per panel and assembles them in main();
# nothing here depends on state defined by the caller.
suppressMessages({ library(ggplot2); library(patchwork); library(dplyr) })

# --- Inventories ---------------------------------------------------------------------------------

method_table <- function(root = ".") read.delim(file.path(root, "config/methods.tsv"))

# Axis order used in every panel: baselines first, then the methods that return a table, then the
# methods that return an embedding, keeping the inventory's own order inside each group.
method_order <- function(methods = method_table())
  methods$method[order(factor(methods$track, c("baseline", "A", "B")))]

# Attach label/track to any table with a `method` column and order it for plotting.
label_methods <- function(scored, methods = method_table()) {
  joined <- left_join(scored, methods[, c("method", "label", "track", "output", "role")], by = "method")
  joined$label <- factor(joined$label, methods$label[match(method_order(methods), methods$method)])
  joined$track <- factor(joined$track, c("baseline", "A", "B"))
  joined
}

# Which columns of results/scores.tsv are metrics, what they mean and which way is better.
# `transform = "abs"` marks a metric whose distance from zero is what matters (batch silhouette).
metric_catalog <- function() read.delim(text = "metric\tfamily\tdirection\ttransform\tlabel
R2_batch_ait\tbatch\t-1\tidentity\tR2 batch (Aitchison)
R2_batch_bc\tbatch\t-1\tidentity\tR2 batch (Bray-Curtis)
R2_batch_jac\tbatch\t-1\tidentity\tR2 batch (Jaccard)
ASW_batch\tbatch\t-1\tabs\tbatch ASW
iLISI\tbatch\t1\tidentity\tiLISI
R2_pheno_ait\tbio\t1\tidentity\tR2 phenotype (Aitchison)
R2_pheno_bc\tbio\t1\tidentity\tR2 phenotype (Bray-Curtis)
R2_pheno_jac\tbio\t1\tidentity\tR2 phenotype (Jaccard)
ASW_pheno\tbio\t1\tidentity\tphenotype ASW
cLISI\tbio\t1\tidentity\tcLISI
oracle_mantel\ttruth\t1\tidentity\tdistance fidelity (Mantel r)
oracle_taxon_rho\ttruth\t1\tidentity\tper-taxon fidelity (rho)
DA_FDR\tda\t-1\tidentity\tDA false discovery rate
DA_power\tda\t1\tidentity\tDA power
DA_AP\tda\t1\tidentity\tDA average precision", sep = "\t")

# --- Reading results -----------------------------------------------------------------------------

# results/scores.tsv, with the dataset path split into the parts the figures group by.
read_scores <- function(path) {
  scored <- read.delim(path, check.names = FALSE)
  parts <- strsplit(scored$dataset, "/")
  scored$regime <- vapply(parts, `[`, "", 1)
  scored$scenario <- vapply(parts, function(p) p[2], "")
  scored$rep <- vapply(parts, function(p) if (length(p) > 2) p[3] else NA_character_, "")
  scored$family <- sub("_(16s|mgx)_.*$", "", scored$scenario)   # sweep / null / lhs
  scored
}

# --- Normalisation -------------------------------------------------------------------------------

# One metric, oriented so that larger is better, with the transform the catalogue asks for.
oriented <- function(values, metric, catalog = metric_catalog()) {
  spec <- catalog[match(metric, catalog$metric), ]
  if (spec$transform == "abs") values <- abs(values)
  spec$direction * values
}

# scIB-style score: min-max inside each dataset so metrics on different scales can be averaged.
# Datasets where every method ties (or the metric does not apply) contribute nothing.
scaled_score <- function(scored, metric, group = "dataset", catalog = metric_catalog()) {
  scored$value <- oriented(scored[[metric]], metric, catalog)
  scored %>%
    group_by(across(all_of(group))) %>%
    mutate(score = (value - min(value, na.rm = TRUE)) /
                   (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(metric = metric, family = catalog$family[match(metric, catalog$metric)])
}

# Plain competition rank, 1 = best, used where the comparison is about order rather than distance.
ranked <- function(scored, metric, group = "dataset", catalog = metric_catalog()) {
  scored$value <- oriented(scored[[metric]], metric, catalog)
  scored %>%
    group_by(across(all_of(group))) %>%
    mutate(rank = rank(-value, na.last = "keep", ties.method = "min")) %>%
    ungroup() %>%
    mutate(metric = metric)
}

# Mean scaled score over several metrics, per method -- the summary used by the heatmaps and the
# trade-off and recommendation panels.
family_score <- function(scored, metrics, group = "dataset", catalog = metric_catalog()) {
  bind_rows(lapply(metrics, function(m) scaled_score(scored, m, group, catalog))) %>%
    group_by(method, family) %>%
    summarise(score = mean(score, na.rm = TRUE), n = sum(!is.na(score)), .groups = "drop")
}

# --- Output --------------------------------------------------------------------------------------

theme_bench <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(size = base_size, margin = margin(1.5, 1.5, 1.5, 1.5, "mm")),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
          plot.tag = element_text(size = base_size + 2, face = "bold"),
          plot.title = element_text(size = base_size, face = "bold"),
          axis.text = element_text(colour = "black"))
}

track_colours <- function() c(baseline = "grey40", A = "#2c7fb8", B = "#d95f0e")

# Nature-style double column is 180 mm; single column 88 mm.
save_figure <- function(figure, path, width_mm, height_mm) {
  ggsave(path, figure, width = width_mm, height = height_mm, units = "mm", device = "pdf")
  ggsave(sub("\\.pdf$", ".png", path), figure, width = width_mm, height = height_mm,
         units = "mm", dpi = 300)
  cat("wrote", path, sprintf("(%g x %g mm)\n", width_mm, height_mm))
}
