# =============================================================================
# Gene tree topology & support analysis
# Species: N. graecorum (outgroup), D. muscipula, D. regia, D. capensis
#
# Pipeline:
#   1. Parse all tree copy-number modalities
#   2. Analyse strictly single-copy (1:1:1:1) trees
#   3. Expand: draw one tip/species from multi-copy trees ("synthetic")
#   4. Filter both sets for strong support (aLRT >= 70, UFBoot >= 70)
#   5. Plot support distributions and topology frequency diagrams
#   6. Add inference:
#        - bootstrap 95% CIs on topology proportions
#        - formal goodness-of-fit test against equal-topology null
#   7. Assess robustness of synthetic set:
#        - repeated draw-one replicates per orthogroup
#        - per-orthogroup stability summaries
#        - replicate-level topology proportion summaries
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

# ---- Constants ---------------------------------------------------------------

TREE_DIR      <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/nucleotide_trees/og_iqtree"
SUPPORT_THRES <- 70

# Number of repeated synthetic draw-one replicates per multi-copy orthogroup
N_SYNTHETIC_REPS <- 100

# Seeds for reproducibility
SEED_SINGLE_SYNTHETIC <- 1
SEED_BOOTSTRAP        <- 1
SEED_REPLICATES       <- 1

species_order <- c("N_gra_dom", "Dio_muscipula", "D_regia", "D_capensis")

ingroup_pairs <- list(
  Dio_Reg = c("Dio_muscipula", "D_regia"),
  Dio_Cap = c("Dio_muscipula", "D_capensis"),
  Reg_Cap = c("D_regia",       "D_capensis")
)

topology_levels <- names(ingroup_pairs)

sp_cols <- c(
  N_gra_dom     = "#1b9e77",
  Dio_muscipula = "#d95f02",
  D_regia       = "#7570b3",
  D_capensis    = "#e7298a",
  mixed         = "grey60"
)

topology_plot_cols <- c(
  Dio_Reg = "red",
  Dio_Cap = "green3",
  Reg_Cap = "blue"
)

# Canonical 4-taxon tree strings for each topology
topology_strings <- list(
  Dio_Reg = "(N_gra_dom,(D_capensis,(Dio_muscipula,D_regia)));",
  Dio_Cap = "(N_gra_dom,(D_regia,(Dio_muscipula,D_capensis)));",
  Reg_Cap = "(N_gra_dom,(Dio_muscipula,(D_regia,D_capensis)));"
)

# ---- Helper functions --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Map tip labels to species names
tip_to_species <- function(x) {
  case_when(
    str_detect(x, "^Dcp_")      ~ "D_capensis",
    str_detect(x, "^DRG__")     ~ "D_regia",
    str_detect(x, "^_scaffold") ~ "Dio_muscipula",
    str_detect(x, "^Nepgr")     ~ "N_gra_dom",
    TRUE                        ~ NA_character_
  )
}

# Root tree on N_gra_dom (single tip or MRCA of multiple copies)
root_on_nepgr <- function(tr) {
  nep_tips <- tr$tip.label[tip_to_species(tr$tip.label) == "N_gra_dom"]
  if (length(nep_tips) == 0) return(tr)
  outgroup <- if (length(nep_tips) == 1) nep_tips else getMRCA(tr, nep_tips)
  root(tr, outgroup = outgroup, resolve.root = TRUE)
}

# Parse IQ-TREE node labels "aLRT/UFBoot" and attach as tree attributes
attach_support_vectors <- function(tr) {
  raw <- tr$node.label %||% rep(NA_character_, tr$Nnode)
  raw <- trimws(raw)
  raw[raw == "" | tolower(raw) == "root"] <- NA_character_

  parts <- strsplit(ifelse(is.na(raw), NA_character_, raw), "/", fixed = TRUE)

  parse_idx <- function(z, i) {
    if (length(z) < i || all(is.na(z))) {
      NA_real_
    } else {
      suppressWarnings(as.numeric(z[i]))
    }
  }

  attr(tr, "alrt")   <- vapply(parts, parse_idx, numeric(1), i = 1)
  attr(tr, "ufboot") <- vapply(parts, parse_idx, numeric(1), i = 2)
  tr
}

# Get aLRT/UFBoot for a specific node (NA-safe)
node_support <- function(tr, node) {
  if (is.na(node)) return(c(alrt = NA_real_, ufboot = NA_real_))
  idx <- node - length(tr$tip.label)
  if (idx < 1 || idx > tr$Nnode) return(c(alrt = NA_real_, ufboot = NA_real_))
  c(alrt = attr(tr, "alrt")[idx], ufboot = attr(tr, "ufboot")[idx])
}

# Get the single tip for a species (returns NA if not exactly 1)
tip_for_species <- function(tr, sp_name) {
  hits <- tr$tip.label[tip_to_species(tr$tip.label) == sp_name]
  if (length(hits) != 1) NA_character_ else hits
}

# Infer which sister pair is monophyletic within the ingroup
infer_topology <- function(tr, dio, reg, cap) {
  flags <- c(
    Dio_Reg = is.monophyletic(tr, c(dio, reg)),
    Dio_Cap = is.monophyletic(tr, c(dio, cap)),
    Reg_Cap = is.monophyletic(tr, c(reg, cap))
  )
  if (sum(flags, na.rm = TRUE) != 1) return(NA_character_)
  names(which(flags))
}

# Compact copy-number pattern string for a tree file
tree_pattern_key <- function(tr) {
  sp  <- tip_to_species(tr$tip.label)
  tab <- table(factor(sp[!is.na(sp)], levels = species_order))
  paste0(names(tab), "=", as.integer(tab), collapse = ";")
}

# Build canonical tree for a named topology mode
make_topology_tree <- function(mode) {
  read.tree(text = topology_strings[[mode]])
}

# Create one random synthetic single-copy tree from a multi-copy tree
make_random_single_copy <- function(tr) {
  sp <- tip_to_species(tr$tip.label)
  if (!all(species_order %in% na.omit(sp))) return(NULL)

  keep <- sapply(species_order, function(s) sample(tr$tip.label[sp == s], 1))
  attach_support_vectors(keep.tip(tr, keep))
}

# ---- Load all trees ----------------------------------------------------------

tree_files <- list.files(TREE_DIR, pattern = "\\.treefile$", full.names = TRUE)
stopifnot(length(tree_files) > 0)

trees_all <- setNames(
  lapply(tree_files, function(f) root_on_nepgr(read.tree(f))),
  tree_files
)

keys <- vapply(trees_all, tree_pattern_key, character(1))

cat("Total trees:", length(trees_all), "\n")

# ---- 1. Copy-number modality summary -----------------------------------------

modality_tab <- as.data.frame(table(keys), stringsAsFactors = FALSE) %>%
  rename(pattern = keys, n_trees = Freq) %>%
  mutate(pct = round(100 * n_trees / sum(n_trees), 2)) %>%
  arrange(desc(n_trees))

print(head(modality_tab, 30))
write.csv(modality_tab, "copy_number_modalities.csv", row.names = FALSE)
cat("Wrote: copy_number_modalities.csv\n")

# ---- 2. Single-copy (1:1:1:1) tree analysis ----------------------------------

SINGLE_KEY <- paste0(species_order, "=1", collapse = ";")

single_copy_trees <- lapply(trees_all[keys == SINGLE_KEY], attach_support_vectors)
multi_copy_trees  <- trees_all[keys != SINGLE_KEY]

cat("Single-copy trees:", length(single_copy_trees), "\n")
cat("Multi-copy trees: ", length(multi_copy_trees), "\n")

# ---- 3. Synthetic single-copy trees (one random draw from multi-copy) --------

set.seed(SEED_SINGLE_SYNTHETIC)
synthetic_trees <- discard(
  lapply(multi_copy_trees, make_random_single_copy),
  is.null
)

cat("Synthetic single-copy trees:", length(synthetic_trees), "\n")

# ---- 4. Build results table from a list of trees -----------------------------

extract_topology_support <- function(tree_list) {
  imap_dfr(tree_list, function(tr, nm) {
    nep <- tip_for_species(tr, "N_gra_dom")
    dio <- tip_for_species(tr, "Dio_muscipula")
    reg <- tip_for_species(tr, "D_regia")
    cap <- tip_for_species(tr, "D_capensis")

    if (any(is.na(c(nep, dio, reg, cap)))) return(NULL)

    topo <- infer_topology(tr, dio, reg, cap)

    # support on the sister-pair clade that defines the topology
    pair_tips <- switch(
      topo,
      Dio_Reg = c(dio, reg),
      Dio_Cap = c(dio, cap),
      Reg_Cap = c(reg, cap),
      NULL
    )

    pair_node <- if (!is.null(pair_tips) && is.monophyletic(tr, pair_tips)) {
      getMRCA(tr, pair_tips)
    } else {
      NA_integer_
    }

    sup <- node_support(tr, pair_node)

    tibble(
      og       = sub("\\.treefile$", "", basename(nm)),
      topology = topo,
      alrt     = unname(sup["alrt"]),
      ufboot   = unname(sup["ufboot"])
    )
  })
}

# Build tables for both sets
res_single   <- extract_topology_support(single_copy_trees)
res_expanded <- extract_topology_support(c(single_copy_trees, synthetic_trees))

write.csv(res_single,   "single_copy_supports.csv", row.names = FALSE)
write.csv(res_expanded, "expanded_supports.csv",    row.names = FALSE)
cat("Wrote: single_copy_supports.csv, expanded_supports.csv\n")

# Strong-support subsets
res_single_strong <- filter(res_single, alrt >= SUPPORT_THRES, ufboot >= SUPPORT_THRES)
res_expanded_strong <- filter(res_expanded, alrt >= SUPPORT_THRES, ufboot >= SUPPORT_THRES)

write.csv(res_single_strong,   "single_copy_supports_strong.csv", row.names = FALSE)
write.csv(res_expanded_strong, "expanded_supports_strong.csv",    row.names = FALSE)
cat("Wrote: single_copy_supports_strong.csv, expanded_supports_strong.csv\n")

# ---- 5. Plotting helpers -----------------------------------------------------

# Faceted UFBoot / aLRT histograms by topology
plot_support_histograms <- function(res, title, outfile_stem) {
  plot_df <- res %>%
    select(og, topology, alrt, ufboot) %>%
    pivot_longer(c(alrt, ufboot), names_to = "metric", values_to = "support") %>%
    filter(!is.na(support), !is.na(topology)) %>%
    mutate(
      metric = recode(
        metric,
        alrt   = "SH-aLRT (1000 reps)",
        ufboot = "UFBoot (1000 reps)"
      ),
      topology = factor(topology, levels = topology_levels)
    )

  if (nrow(plot_df) == 0) {
    cat("Skipping histogram plot for", outfile_stem, "(no data)\n")
    return(invisible(NULL))
  }

  p <- ggplot(plot_df, aes(x = support, fill = topology)) +
    geom_histogram(binwidth = 100 / 30, boundary = 0, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = topology_plot_cols, drop = FALSE) +
    facet_grid(metric ~ topology, scales = "fixed") +
    labs(x = "Support value", y = "Number of gene trees", title = title) +
    theme_bw() +
    theme(legend.position = "none")

  ggsave(paste0(outfile_stem, ".pdf"), p, width = 11, height = 6)
  ggsave(paste0(outfile_stem, ".png"), p, width = 11, height = 6, dpi = 200)
  cat("Wrote:", outfile_stem, ".pdf/.png\n")
}

# Stacked cladogram plot (edge width proportional to topology frequency)
plot_topology_modes <- function(res, outfile, title_suffix = "", n_boot = 10000, seed = 1) {
  ci_df <- bootstrap_topology_props(res, n_boot = n_boot, seed = seed)

  freq <- ci_df$n
  names(freq) <- ci_df$topology

  if (sum(freq, na.rm = TRUE) == 0) {
    cat("Skipping topology mode plot for", outfile, "(no resolved topologies)\n")
    return(invisible(NULL))
  }

  maxf <- max(freq, na.rm = TRUE)

  pdf(outfile, width = 10, height = 8.8)
  par(mfrow = c(3, 1), mar = c(1, 1, 3, 1), oma = c(0, 0, 3.5, 0))

  for (m in topology_levels) {
    tr <- make_topology_tree(m)

    row_m <- ci_df %>% filter(topology == m)

    label_text <- sprintf(
      "%s  (%d / %d = %.1f%%, 95%% CI %.1f–%.1f%%)%s",
      m,
      row_m$n,
      row_m$total,
      100 * row_m$prop,
      100 * row_m$conf.low,
      100 * row_m$conf.high,
      title_suffix
    )

    plot.phylo(
      tr,
      type = "cladogram",
      use.edge.length = FALSE,
      edge.width = 1 + 6 * (row_m$n / maxf),
      cex = 1.3,
      tip.color = sp_cols[tr$tip.label],
      main = label_text
    )

    tiplabels(pch = 19, col = sp_cols[tr$tip.label], cex = 1.3)
    legend(
      "topleft",
      legend = names(sp_cols)[1:4],
      col = sp_cols[1:4],
      pch = 19,
      bty = "n",
      cex = 0.9
    )
  }

mtext(
  sprintf(
    paste(
      "Topology frequencies across resolved gene trees.",
      "95%% confidence intervals were estimated by nonparametric bootstrap",
      "resampling of orthogroups with replacement (%s bootstrap replicates).",
      sep = "\n"
    ),
    format(n_boot, big.mark = ",")
  ),
  side = 3,
  outer = TRUE,
  line = 1.2,
  cex = 0.95
)

  dev.off()
  cat("Wrote:", outfile, "\n")

  invisible(ci_df)
}

plot_synthetic_stability <- function(per_og_df, outfile_stem) {
  if (nrow(per_og_df) == 0) {
    cat("Skipping synthetic stability plot for", outfile_stem, "(no data)\n")
    return(invisible(NULL))
  }

  p <- ggplot(per_og_df, aes(x = modal_prop)) +
    geom_histogram(binwidth = 0.05, boundary = 0, colour = "white") +
    labs(
      x = "Fraction of draw-one replicates supporting the modal topology",
      y = "Number of orthogroups",
      title = "Stability of topology assignment across repeated synthetic draws"
    ) +
    theme_bw()

  ggsave(paste0(outfile_stem, ".pdf"), p, width = 8, height = 5)
  ggsave(paste0(outfile_stem, ".png"), p, width = 8, height = 5, dpi = 200)
  cat("Wrote:", outfile_stem, ".pdf/.png\n")
}

plot_synthetic_stability_bars <- function(rep_res, outfile_stem) {
  plot_df <- rep_res %>%
    filter(!is.na(topology)) %>%
    count(og, topology, name = "n") %>%
    tidyr::complete(og, topology = topology_levels, fill = list(n = 0)) %>%
    group_by(og) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  if (nrow(plot_df) == 0) {
    cat("Skipping synthetic stability bar plot for", outfile_stem, "(no resolved topologies)\n")
    return(invisible(NULL))
  }

  og_meta <- plot_df %>%
    group_by(og) %>%
    summarise(
      modal_topology = topology[which.max(prop)][1],
      modal_prop = max(prop),
      .groups = "drop"
    ) %>%
    mutate(
      modal_topology = factor(modal_topology, levels = topology_levels)
    ) %>%
    arrange(modal_topology, desc(modal_prop), og)

  plot_df <- plot_df %>%
    left_join(og_meta, by = "og") %>%
    mutate(
      og = factor(og, levels = og_meta$og),
      topology = factor(topology, levels = topology_levels),
      modal_topology = factor(modal_topology, levels = topology_levels)
    )

  p <- ggplot(plot_df, aes(x = og, y = prop, fill = topology)) +
    geom_col(width = 0.95) +
    facet_grid(. ~ modal_topology, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = topology_plot_cols, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      x = "Orthogroup",
      y = "Fraction of synthetic replicates",
      fill = "Topology",
      title = "Topology assignment across repeated synthetic draws"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold")
    )

  ggsave(paste0(outfile_stem, ".pdf"), p, width = 12, height = 5.5)
  ggsave(paste0(outfile_stem, ".png"), p, width = 12, height = 5.5, dpi = 200)
  cat("Wrote:", outfile_stem, ".pdf/.png\n")

  invisible(plot_df)
}
plot_synthetic_stability_heatmap <- function(rep_res, outfile_stem) {
  plot_df <- rep_res %>%
    filter(!is.na(topology)) %>%
    count(og, topology, name = "n") %>%
    tidyr::complete(og, topology = topology_levels, fill = list(n = 0)) %>%
    group_by(og) %>%
    mutate(
      total = sum(n),
      prop = ifelse(total > 0, n / total, NA_real_)
    ) %>%
    ungroup() %>%
    filter(!is.na(prop))

  if (nrow(plot_df) == 0) {
    cat("Skipping synthetic stability heatmap for", outfile_stem, "(no resolved topologies after filtering)\n")
    return(invisible(NULL))
  }

  og_meta <- plot_df %>%
    group_by(og) %>%
    summarise(
      modal_topology = topology[which.max(prop)][1],
      modal_prop = max(prop),
      .groups = "drop"
    ) %>%
    mutate(
      modal_topology = factor(modal_topology, levels = topology_levels)
    ) %>%
    arrange(modal_topology, desc(modal_prop), og)

  plot_df <- plot_df %>%
    mutate(
      og = factor(og, levels = rev(og_meta$og)),
      topology = factor(topology, levels = topology_levels)
    )

  p <- ggplot(plot_df, aes(x = topology, y = og, fill = prop)) +
    geom_tile() +
    scale_fill_gradient(
      limits = c(0, 1),
      low = "white",
      high = "steelblue",
      na.value = "grey90"
    ) +
    labs(
      x = "Topology",
      y = "Orthogroup",
      fill = "Replicate\nfraction",
      title = "Topology assignment across repeated synthetic draws"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )

  ggsave(paste0(outfile_stem, ".pdf"), p, width = 7, height = 10)
  ggsave(paste0(outfile_stem, ".png"), p, width = 7, height = 10, dpi = 200)
  cat("Wrote:", outfile_stem, ".pdf/.png\n")

  invisible(plot_df)
}
# ---- 6. Inferential helpers --------------------------------------------------

bootstrap_topology_props <- function(res, n_boot = 10000, seed = 1) {
  dat <- res %>%
    filter(!is.na(topology)) %>%
    mutate(topology = factor(topology, levels = topology_levels))

  if (nrow(dat) == 0) {
    return(tibble(
      topology  = topology_levels,
      n         = 0L,
      total     = 0L,
      prop      = NA_real_,
      conf.low  = NA_real_,
      conf.high = NA_real_
    ))
  }

  obs_counts <- table(dat$topology)
  obs_props  <- as.numeric(obs_counts) / sum(obs_counts)

  set.seed(seed)
  boot_mat <- replicate(n_boot, {
    idx <- sample.int(nrow(dat), size = nrow(dat), replace = TRUE)
    sampled <- factor(dat$topology[idx], levels = topology_levels)
    as.numeric(table(sampled)) / sum(table(sampled))
  })

  tibble(
    topology  = topology_levels,
    n         = as.integer(obs_counts),
    total     = sum(obs_counts),
    prop      = obs_props,
    conf.low  = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    conf.high = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}


# ---- 7. Repeated synthetic draw-one replicates -------------------------------

generate_synthetic_replicates <- function(tree_list, n_reps = 100, seed = 1) {
  nm_vec <- names(tree_list)

  map2_dfr(tree_list, seq_along(tree_list), function(tr, i) {
    nm <- nm_vec[i]

    map_dfr(seq_len(n_reps), function(rep_id) {
      # deterministic seed per orthogroup x replicate
      set.seed(seed + i * 100000L + rep_id)

      syn_tr <- make_random_single_copy(tr)
      if (is.null(syn_tr)) return(NULL)

      out <- extract_topology_support(setNames(list(syn_tr), nm))
      if (nrow(out) == 0) return(NULL)

      out %>%
        mutate(replicate = rep_id)
    })
  })
}

summarize_synthetic_robustness <- function(rep_res, support_thres = SUPPORT_THRES) {
  if (nrow(rep_res) == 0) {
    empty_per_og <- tibble(
      og = character(),
      modal_topology = character(),
      modal_prop = numeric(),
      n_topologies_seen = integer(),
      n_reps = integer(),
      strong_prop = numeric()
    )
    empty_overall <- tibble(
      n_ogs = integer(),
      mean_modal_prop = numeric(),
      median_modal_prop = numeric(),
      frac_perfectly_stable = numeric(),
      frac_ge_95pct_stable = numeric(),
      frac_switching = numeric(),
      mean_strong_prop = numeric()
    )
    return(list(per_og = empty_per_og, overall = empty_overall))
  }

  strong_df <- rep_res %>%
    mutate(
      strong = !is.na(alrt) & !is.na(ufboot) &
        alrt >= support_thres & ufboot >= support_thres
    )

  per_og <- strong_df %>%
    filter(!is.na(topology)) %>%
    count(og, topology, name = "n") %>%
    group_by(og) %>%
    mutate(prop = n / sum(n)) %>%
    summarise(
      modal_topology = topology[which.max(prop)][1],
      modal_prop = max(prop),
      n_topologies_seen = n(),
      n_reps = sum(n),
      .groups = "drop"
    ) %>%
    left_join(
      strong_df %>%
        group_by(og) %>%
        summarise(
          strong_prop = mean(strong, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "og"
    )

  overall <- per_og %>%
    summarise(
      n_ogs = n(),
      mean_modal_prop = mean(modal_prop, na.rm = TRUE),
      median_modal_prop = median(modal_prop, na.rm = TRUE),
      frac_perfectly_stable = mean(modal_prop == 1, na.rm = TRUE),
      frac_ge_95pct_stable = mean(modal_prop >= 0.95, na.rm = TRUE),
      frac_switching = mean(n_topologies_seen > 1, na.rm = TRUE),
      mean_strong_prop = mean(strong_prop, na.rm = TRUE)
    )

  list(per_og = per_og, overall = overall)
}

summarize_replicate_topology_props <- function(rep_res) {
  if (nrow(rep_res) == 0) {
    return(tibble(
      topology = character(),
      mean_prop = numeric(),
      conf.low = numeric(),
      conf.high = numeric()
    ))
  }

  rep_res %>%
    filter(!is.na(topology)) %>%
    count(replicate, topology, name = "n") %>%
    tidyr::complete(replicate, topology = topology_levels, fill = list(n = 0)) %>%
    group_by(replicate) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    group_by(topology) %>%
    summarise(
      mean_prop = mean(prop),
      conf.low = quantile(prop, 0.025),
      conf.high = quantile(prop, 0.975),
      .groups = "drop"
    )
}

# ---- 8. Generate core plots --------------------------------------------------

# All single-copy trees
plot_support_histograms(
  res_single,
  "Support distributions – single-copy trees",
  "support_distributions_single_copy"
)
plot_topology_modes(res_single, "topology_modes_all.pdf")

# Strong-support single-copy
plot_topology_modes(
  res_single_strong,
  "topology_modes_strong_support.pdf",
  " – strong support only"
)

# Expanded (single + synthetic)
plot_support_histograms(
  res_expanded,
  "Support distributions – single-copy + synthetic",
  "support_distributions_expanded"
)
plot_topology_modes(res_expanded, "topology_modes_expanded_all.pdf")

# Strong-support expanded
plot_topology_modes(
  res_expanded_strong,
  "topology_modes_expanded_strong_support.pdf",
  " – strong support only"
)

# ---- 9. Inferential summaries ------------------------------------------------

# ---- 9. Bootstrap CI summaries -----------------------------------------------

ci_single <- bootstrap_topology_props(
  res_single,
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

ci_single_strong <- bootstrap_topology_props(
  res_single_strong,
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

ci_expanded <- bootstrap_topology_props(
  res_expanded,
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

ci_expanded_strong <- bootstrap_topology_props(
  res_expanded_strong,
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

write.csv(ci_single,          "topology_bootstrap_ci_single.csv", row.names = FALSE)
write.csv(ci_single_strong,   "topology_bootstrap_ci_single_strong.csv", row.names = FALSE)
write.csv(ci_expanded,        "topology_bootstrap_ci_expanded.csv", row.names = FALSE)
write.csv(ci_expanded_strong, "topology_bootstrap_ci_expanded_strong.csv", row.names = FALSE)

cat("Wrote bootstrap CI tables.\n")


# ---- 10. Repeated synthetic replicate analysis -------------------------------

synthetic_reps <- generate_synthetic_replicates(
  multi_copy_trees,
  n_reps = N_SYNTHETIC_REPS,
  seed = SEED_REPLICATES
)

write.csv(synthetic_reps, "synthetic_replicate_topologies.csv", row.names = FALSE)

syn_robust <- summarize_synthetic_robustness(
  synthetic_reps,
  support_thres = SUPPORT_THRES
)

syn_rep_props <- summarize_replicate_topology_props(synthetic_reps)

write.csv(syn_robust$per_og,   "synthetic_per_og_robustness.csv", row.names = FALSE)
write.csv(syn_robust$overall,  "synthetic_overall_robustness.csv", row.names = FALSE)
write.csv(syn_rep_props,       "synthetic_replicate_topology_props.csv", row.names = FALSE)

cat("Wrote repeated synthetic replicate summaries.\n")

plot_synthetic_stability(
  syn_robust$per_og,
  "synthetic_topology_stability"
)

plot_synthetic_stability_bars(
  synthetic_reps,
  "synthetic_topology_stability_bars"
)
plot_synthetic_stability_heatmap(
  synthetic_reps,
  "synthetic_topology_stability_heatmap"
)
# All single-copy trees
plot_support_histograms(
  res_single,
  "Support distributions – single-copy trees",
  "support_distributions_single_copy"
)
plot_topology_modes(
  res_single,
  "topology_modes_all.pdf",
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

# Strong-support single-copy
plot_topology_modes(
  res_single_strong,
  "topology_modes_strong_support.pdf",
  title_suffix = " – strong support only",
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

# Expanded (single + synthetic)
plot_support_histograms(
  res_expanded,
  "Support distributions – single-copy + synthetic",
  "support_distributions_expanded"
)
plot_topology_modes(
  res_expanded,
  "topology_modes_expanded_all.pdf",
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)

# Strong-support expanded
plot_topology_modes(
  res_expanded_strong,
  "topology_modes_expanded_strong_support.pdf",
  title_suffix = " – strong support only",
  n_boot = 10000,
  seed = SEED_BOOTSTRAP
)
# ---- 11. Console summaries ---------------------------------------------------

cat("\n=== Topology counts: single-copy ===\n")
print(table(factor(res_single$topology, levels = topology_levels), useNA = "ifany"))

cat("\n=== Topology counts: single-copy, strong support ===\n")
print(table(factor(res_single_strong$topology, levels = topology_levels), useNA = "ifany"))

cat("\n=== Topology counts: expanded ===\n")
print(table(factor(res_expanded$topology, levels = topology_levels), useNA = "ifany"))

cat("\n=== Topology counts: expanded, strong support ===\n")
print(table(factor(res_expanded_strong$topology, levels = topology_levels), useNA = "ifany"))

cat("\n=== Bootstrap CIs: single-copy ===\n")
print(ci_single)

cat("\n=== Bootstrap CIs: single-copy, strong support ===\n")
print(ci_single_strong)

cat("\n=== Bootstrap CIs: expanded ===\n")
print(ci_expanded)

cat("\n=== Bootstrap CIs: expanded, strong support ===\n")
print(ci_expanded_strong)


cat("\n=== Synthetic robustness: overall ===\n")
print(syn_robust$overall)

cat("\n=== Synthetic topology stability (sanity check for bar plot) ===\n")

stab_summary <- syn_robust$per_og %>%
  mutate(
    modal_topology = factor(modal_topology, levels = topology_levels)
  )

cat("\nOrthogroups by modal topology:\n")
print(table(stab_summary$modal_topology))

cat("\nModal topology stability (fraction of replicates supporting modal topology):\n")
print(summary(stab_summary$modal_prop))

cat("\nOrthogroups switching topology across replicates:\n")
print(
  stab_summary %>%
    summarise(
      total_ogs = n(),
      switching_ogs = sum(n_topologies_seen > 1),
      frac_switching = mean(n_topologies_seen > 1)
    )
)

cat("\nStability by modal topology:\n")
print(
  stab_summary %>%
    group_by(modal_topology) %>%
    summarise(
      n_ogs = n(),
      median_modal_prop = median(modal_prop),
      mean_modal_prop = mean(modal_prop),
      frac_perfect = mean(modal_prop == 1),
      frac_ge_0.95 = mean(modal_prop >= 0.95),
      .groups = "drop"
    )
)
cat("\nDone.\n")
