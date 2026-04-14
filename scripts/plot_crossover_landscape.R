#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript plot_crossover_landscape.R <project_root>")
}

project_root <- normalizePath(args[1], mustWork = TRUE)

co_intervals_file <- file.path(project_root, "results/07_crossovers/co_intervals.bed")
co_summary_file   <- file.path(project_root, "results/07_crossovers/co_summary.tsv")

# set this explicitly to avoid path guessing problems
fai_file <- "/biodata/dep_mercier/grp_marques/marques/Cuscuta-genomes/final_assemblies/C_epithymum_hap1_chr.fasta.fai"

plot_dir <- file.path(project_root, "results/08_plots/crossover_landscape")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, filename, width = 10, height = 7) {
  ggsave(file.path(plot_dir, filename), p, width = width, height = height, dpi = 300)
}

co_intervals <- read_tsv(
  co_intervals_file,
  show_col_types = FALSE,
  col_names = c("chrom", "start", "end", "sample"),
  col_types = cols(
    chrom = col_character(),
    start = col_double(),
    end = col_double(),
    sample = col_character()
  )
) %>%
  mutate(
    width_bp = end - start + 1,
    midpoint_bp = (start + end) / 2,
    midpoint_mbp = midpoint_bp / 1e6
  )

co_summary <- read_tsv(
  co_summary_file,
  show_col_types = FALSE,
  col_types = cols(
    prefix = col_character(),
    informative_markers = col_double(),
    nonambiguous_markers = col_double(),
    valid_segments = col_double(),
    n_crossovers = col_double(),
    min_markers_cell = col_double(),
    min_segment_bp = col_double(),
    min_markers_segment = col_double(),
    status = col_character()
  )
)

n_cells <- nrow(co_summary)

chrom_sizes <- read_tsv(
  fai_file,
  show_col_types = FALSE,
  col_names = c("chrom", "length_bp", "offset", "line_bases", "line_width"),
  col_types = cols(
    chrom = col_character(),
    length_bp = col_double(),
    offset = col_double(),
    line_bases = col_double(),
    line_width = col_double()
  )
) %>%
  select(chrom, length_bp) %>%
  mutate(length_mbp = length_bp / 1e6)

# keep only chromosomes actually present in the CO table
chrom_sizes <- chrom_sizes %>%
  filter(chrom %in% unique(co_intervals$chrom))

# ------------------------------------------------------------
# Plot 1: recombination landscape
# uncertainty-aware sliding window version
# ------------------------------------------------------------

winsize_bp <- 1e6
stepsize_bp <- 1e5
n_random <- 500
set.seed(123)

get_sliding_windows <- function(chr_len, winsize, stepsize) {
  if (chr_len <= winsize) {
    tibble(
      window_start = 1,
      window_end = chr_len
    )
  } else {
    starts <- seq(1, chr_len - winsize + 1, by = stepsize)
    tibble(
      window_start = starts,
      window_end = pmin(starts + winsize - 1, chr_len)
    ) %>%
      distinct()
  }
}

sample_positions_in_interval <- function(start, end, n = 500) {
  if (is.na(start) || is.na(end) || end < start) {
    rep(NA_real_, n)
  } else if (start == end) {
    rep(start, n)
  } else {
    sample(seq.int(start, end), size = n, replace = TRUE)
  }
}

compute_landscape_one_chr <- function(chr_name, chr_len, intervals_chr, n_cells,
                                      winsize_bp, stepsize_bp, n_random) {
  windows <- get_sliding_windows(chr_len, winsize_bp, stepsize_bp)

  if (nrow(intervals_chr) == 0) {
    return(
      windows %>%
        mutate(
          chrom = chr_name,
          window_mid_bp = (window_start + window_end) / 2,
          window_mid_mbp = window_mid_bp / 1e6,
          rate_mean = 0,
          rate_low = 0,
          rate_high = 0
        )
    )
  }

sampled_pos_list <- purrr::map2(
  intervals_chr$start,
  intervals_chr$end,
  ~ sample_positions_in_interval(.x, .y, n = n_random)
)

sampled_pos_mat <- do.call(cbind, sampled_pos_list)

# ensure matrix even if only one interval
sampled_pos_mat <- as.matrix(sampled_pos_mat)

  window_results <- purrr::pmap_dfr(
    windows,
    function(window_start, window_end) {
      counts_per_draw <- apply(
        sampled_pos_mat,
        1,
        function(pos_vec) sum(pos_vec >= window_start & pos_vec <= window_end, na.rm = TRUE)
      )

      rate_vec <- 100 * (counts_per_draw / n_cells) / ((window_end - window_start + 1) / 1e6)

      tibble(
        chrom = chr_name,
        window_start = window_start,
        window_end = window_end,
        window_mid_bp = (window_start + window_end) / 2,
        window_mid_mbp = ((window_start + window_end) / 2) / 1e6,
        rate_mean = mean(rate_vec, na.rm = TRUE),
      rate_low = as.numeric(quantile(rate_vec, probs = 0.025, na.rm = TRUE)),
rate_high = as.numeric(quantile(rate_vec, probs = 0.975, na.rm = TRUE))
      )
    }
  )

  window_results
}

landscape_df <- purrr::pmap_dfr(
  chrom_sizes %>% select(chrom, length_bp),
  function(chrom, length_bp) {
    chr_name <- chrom
    chr_len <- length_bp

    intervals_chr <- co_intervals %>%
      filter(.data$chrom == .env$chr_name) %>%
      select(start, end)

    compute_landscape_one_chr(
      chr_name = chr_name,
      chr_len = chr_len,
      intervals_chr = intervals_chr,
      n_cells = n_cells,
      winsize_bp = winsize_bp,
      stepsize_bp = stepsize_bp,
      n_random = n_random
    )
  }
)
write_tsv(landscape_df, file.path(plot_dir, "recombination_landscape_table.tsv"))


p1 <- ggplot(landscape_df, aes(x = window_mid_mbp, y = rate_mean)) +
  geom_ribbon(aes(ymin = rate_low, ymax = rate_high), fill = "grey70", alpha = 0.5) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ chrom, scales = "fixed", ncol = 1) +
  theme_bw() +
  labs(
    title = "Recombination landscape",
    subtitle = paste0(
      "Sliding windows: ", winsize_bp / 1e6, " Mbp, step ", stepsize_bp / 1e6,
      " Mbp | Monte Carlo n = ", n_random, " | cells = ", n_cells
    ),
    x = "Physical position (Mbp)",
    y = "Crossover rate (cM/Mbp)"
  )

save_plot(p1, "01_recombination_landscape.png", width = 12, height = 16)

# ------------------------------------------------------------
# Plot 2: Marey map
# ------------------------------------------------------------

marey <- co_intervals %>%
  arrange(chrom, midpoint_bp) %>%
  group_by(chrom) %>%
  mutate(
    co_index = row_number(),
    cumulative_cM = co_index * (100 / n_cells),
    midpoint_mbp = midpoint_bp / 1e6
  ) %>%
  ungroup()

p2 <- ggplot(marey, aes(midpoint_mbp, cumulative_cM, color = chrom)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  theme_bw() +
  labs(
    title = "Marey map",
    subtitle = paste0("All chromosomes | cells = ", n_cells),
    x = "Physical position (Mbp)",
    y = "Genetic position (cM)",
    color = "Chromosome"
  )

save_plot(p2, "02_marey_map.png", width = 10, height = 7)

# ------------------------------------------------------------
# Plot 3: per-chromosome CO class proportions
# ------------------------------------------------------------

all_samples <- tibble(sample = co_summary$prefix)
all_sample_chrom <- tidyr::crossing(sample = all_samples$sample, chrom = chrom_sizes$chrom)

co_per_sample_chr <- co_intervals %>%
  count(sample, chrom, name = "n_co")

co_per_sample_chr_full <- all_sample_chrom %>%
  left_join(co_per_sample_chr, by = c("sample", "chrom")) %>%
  mutate(n_co = replace_na(n_co, 0L)) %>%
  mutate(class = case_when(
    n_co == 0 ~ "0",
    n_co == 1 ~ "1",
    n_co >= 2 ~ "2+"
  ))

prop_df <- co_per_sample_chr_full %>%
  count(chrom, class, name = "n") %>%
  group_by(chrom) %>%
  mutate(prop = 100 * n / sum(n)) %>%
  ungroup()

prop_df <- prop_df %>%
  mutate(class = factor(class, levels = c("0", "1", "2+")))

p3 <- ggplot(prop_df, aes(chrom, prop, fill = class)) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.1f%%", prop)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  theme_bw() +
  labs(
    title = "Per-chromosome CO class proportions",
    subtitle = "Classes: 0, 1, or 2+ crossovers per chromosome per cell",
    x = "Chromosome",
    y = "Percentage of cells",
    fill = "COs"
  )

save_plot(p3, "03_chromosome_co_proportions.png", width = 11, height = 7)

# ------------------------------------------------------------
# Plot 4: coefficient of coincidence curve
# ------------------------------------------------------------

n_bins_coc <- 15

interval_tbl <- purrr::pmap_dfr(
  chrom_sizes %>% select(chrom, length_bp),
  function(chrom, length_bp) {
    tibble(
      chrom = chrom,
      bin = seq_len(n_bins_coc),
      bin_start = floor((seq_len(n_bins_coc) - 1) * length_bp / n_bins_coc) + 1,
      bin_end = floor(seq_len(n_bins_coc) * length_bp / n_bins_coc)
    )
  }
)

co_binned <- interval_tbl %>%
  left_join(co_intervals %>% select(chrom, sample, midpoint_bp), by = "chrom") %>%
  filter(midpoint_bp >= bin_start, midpoint_bp <= bin_end | is.na(midpoint_bp)) %>%
  filter(!is.na(sample)) %>%
  distinct(chrom, sample, bin)

co_presence <- tidyr::crossing(
  sample = all_samples$sample,
  chrom = chrom_sizes$chrom,
  bin = seq_len(n_bins_coc)
) %>%
  left_join(interval_tbl, by = c("chrom", "bin")) %>%
  left_join(co_binned %>% mutate(has_co = 1L), by = c("sample", "chrom", "bin")) %>%
  mutate(has_co = replace_na(has_co, 0L))

pair_results <- list()
idx <- 1L

for (chr in unique(co_presence$chrom)) {
  chr_dat <- co_presence %>% filter(chrom == chr)

  for (i in 1:(n_bins_coc - 1)) {
    for (j in (i + 1):n_bins_coc) {
      d1 <- chr_dat %>% filter(bin == i) %>% arrange(sample)
      d2 <- chr_dat %>% filter(bin == j) %>% arrange(sample)

      merged <- d1 %>%
        select(sample, has_co_i = has_co) %>%
        inner_join(d2 %>% select(sample, has_co_j = has_co), by = "sample")

      p_i <- mean(merged$has_co_i)
      p_j <- mean(merged$has_co_j)
      obs_double <- mean((merged$has_co_i == 1) & (merged$has_co_j == 1))
      exp_double <- p_i * p_j
      coc <- ifelse(exp_double > 0, obs_double / exp_double, NA_real_)

      pair_results[[idx]] <- tibble(
        chrom = chr,
        bin_i = i,
        bin_j = j,
        separation = j - i,
        p_i = p_i,
        p_j = p_j,
        obs_double = obs_double,
        exp_double = exp_double,
        coc = coc
      )
      idx <- idx + 1L
    }
  }
}

coc_tbl <- bind_rows(pair_results)

coc_curve <- coc_tbl %>%
  group_by(separation) %>%
  summarise(
    mean_coc = mean(coc, na.rm = TRUE),
    n_pairs = sum(!is.na(coc)),
    .groups = "drop"
  )

write_tsv(coc_tbl, file.path(plot_dir, "coc_pairwise_table.tsv"))
write_tsv(coc_curve, file.path(plot_dir, "coc_curve_summary.tsv"))

p4 <- ggplot(coc_curve, aes(separation, mean_coc)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title = "Coefficient of coincidence curve",
    subtitle = paste0("Chromosomes divided into ", n_bins_coc, " intervals | cells = ", n_cells),
    x = "Interval separation",
    y = "Mean CoC"
  )

save_plot(p4, "04_coefficient_of_coincidence_curve.png", width = 9, height = 6)

message("Crossover landscape plots written to: ", plot_dir)
