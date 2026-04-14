#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: plot_co_qc_distributions.R <project_root>")
}

project_root <- normalizePath(args[1], mustWork = TRUE)

marker_dir    <- file.path(project_root, "results", "06_markers")
switches_file <- file.path(marker_dir, "switches.stats.tsv")
selected_file <- file.path(marker_dir, "selected_for_co.txt")
plot_dir      <- file.path(project_root, "results", "08_plots", "co_qc")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# current thresholds
THRESH_C <- 500
THRESH_S <- 3000000
THRESH_N <- 10
THRESH_SWITCH_RATE <- 0.10

message("Reading switch summary...")
switches <- read_tsv(
  switches_file,
  show_col_types = FALSE,
  col_types = cols(
    barcode = col_character(),
    marker_num = col_double(),
    switch_num = col_double(),
    switch_rate = col_double()
  )
)

selected <- read_lines(selected_file)

save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

#------------------------------------------------------------
# 1. Per-cell distributions already available
#------------------------------------------------------------

p1 <- ggplot(switches, aes(x = marker_num)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_C, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Marker count per cell",
    subtitle = paste0(
      "Current threshold -c = ", THRESH_C,
      " | passing = ", sum(switches$marker_num >= THRESH_C),
      " / ", nrow(switches)
    ),
    x = "marker_num (log10 scale)",
    y = "Number of cells"
  ) +
  theme_bw()

save_plot(p1, "01_marker_num_hist.png")

p2 <- ggplot(switches, aes(x = switch_num)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  labs(
    title = "Switch count per cell",
    x = "switch_num (log10 scale)",
    y = "Number of cells"
  ) +
  theme_bw()

save_plot(p2, "02_switch_num_hist.png")

p3 <- ggplot(switches, aes(x = switch_rate)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_SWITCH_RATE, linetype = "dashed") +
  labs(
    title = "Switch rate per cell",
    subtitle = paste0(
      "Current threshold = ", THRESH_SWITCH_RATE,
      " | passing = ", sum(switches$switch_rate <= THRESH_SWITCH_RATE),
      " / ", nrow(switches)
    ),
    x = "switch_rate",
    y = "Number of cells"
  ) +
  theme_bw()

save_plot(p3, "03_switch_rate_hist.png")

p4 <- ggplot(switches, aes(x = marker_num, y = switch_rate)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_vline(xintercept = THRESH_C, linetype = "dashed") +
  geom_hline(yintercept = THRESH_SWITCH_RATE, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Marker count vs switch rate",
    x = "marker_num (log10 scale)",
    y = "switch_rate"
  ) +
  theme_bw()

save_plot(p4, "04_marker_num_vs_switch_rate.png")

#------------------------------------------------------------
# 2. Build genotype block summaries from input_corrected.txt
#------------------------------------------------------------

message("Building genotype block summaries from selected_for_co cells...")

state_from_counts <- function(ref_count, alt_count) {
  denom <- ref_count + alt_count
  frac <- ifelse(denom > 0, ref_count / denom, NA_real_)
  case_when(
    is.na(frac) ~ NA_real_,
    frac < 0.2 ~ 0,
    frac > 0.8 ~ 1,
    TRUE ~ 0.5
  )
}

block_summaries <- vector("list", length(selected))

for (i in seq_along(selected)) {
  bc <- selected[[i]]
  f <- file.path(marker_dir, bc, "input_corrected.txt")

  if (!file.exists(f) || file.info(f)$size == 0) {
    next
  }

  dat <- read_table(
    f,
    col_names = c("chrom", "pos", "ref", "ref_count", "alt", "alt_count"),
    show_col_types = FALSE
  )

  dat <- dat %>%
    filter(ref_count > 0 | alt_count > 0) %>%
    mutate(
      state = state_from_counts(ref_count, alt_count)
    )

  if (nrow(dat) == 0) next

  # keep chromosome-aware block splitting
  dat <- dat %>%
    arrange(chrom, pos) %>%
    mutate(
      new_block = if_else(
        row_number() == 1 |
          chrom != lag(chrom) |
          state != lag(state),
        TRUE, FALSE
      ),
      block_id = cumsum(new_block)
    )

  blocks <- dat %>%
    group_by(barcode = bc, chrom, block_id, state) %>%
    summarise(
      n_markers = n(),
      start_pos = min(pos),
      end_pos = max(pos),
      span_bp = end_pos - start_pos + 1,
      .groups = "drop"
    )

  block_summaries[[i]] <- blocks
}

blocks_all <- bind_rows(block_summaries)

if (nrow(blocks_all) == 0) {
  stop("No genotype blocks were built; cannot plot CO QC distributions.")
}

write_tsv(blocks_all, file.path(plot_dir, "genotype_block_summary.tsv"))

#------------------------------------------------------------
# 3. Distributions relevant to -n and -s
#------------------------------------------------------------

p5 <- ggplot(blocks_all, aes(x = n_markers)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_N, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Genotype block size in markers",
    subtitle = paste0(
      "Current threshold -n = ", THRESH_N,
      " | passing blocks = ", sum(blocks_all$n_markers >= THRESH_N),
      " / ", nrow(blocks_all)
    ),
    x = "Markers per genotype block (log10 scale)",
    y = "Number of blocks"
  ) +
  theme_bw()

save_plot(p5, "05_block_marker_count_hist.png")

p6 <- ggplot(blocks_all, aes(x = span_bp)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_S, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Genotype block span in bp",
    subtitle = paste0(
      "Current threshold -s = ", THRESH_S,
      " bp | passing blocks = ", sum(blocks_all$span_bp >= THRESH_S),
      " / ", nrow(blocks_all)
    ),
    x = "Block span (bp, log10 scale)",
    y = "Number of blocks"
  ) +
  theme_bw()

save_plot(p6, "06_block_span_bp_hist.png")

p7 <- ggplot(blocks_all, aes(x = n_markers, y = span_bp)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_vline(xintercept = THRESH_N, linetype = "dashed") +
  geom_hline(yintercept = THRESH_S, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Genotype block markers vs genomic span",
    x = "Markers per block (log10 scale)",
    y = "Block span in bp (log10 scale)"
  ) +
  theme_bw()

save_plot(p7, "07_block_markers_vs_span.png")

#------------------------------------------------------------
# 4. Threshold summaries
#------------------------------------------------------------

threshold_summary <- tibble(
  metric = c(
    "marker_num_per_cell",
    "switch_rate_per_cell",
    "markers_per_block",
    "block_span_bp"
  ),
  threshold = c(
    THRESH_C,
    THRESH_SWITCH_RATE,
    THRESH_N,
    THRESH_S
  ),
  n_total = c(
    nrow(switches),
    nrow(switches),
    nrow(blocks_all),
    nrow(blocks_all)
  ),
  n_pass = c(
    sum(switches$marker_num >= THRESH_C, na.rm = TRUE),
    sum(switches$switch_rate <= THRESH_SWITCH_RATE, na.rm = TRUE),
    sum(blocks_all$n_markers >= THRESH_N, na.rm = TRUE),
    sum(blocks_all$span_bp >= THRESH_S, na.rm = TRUE)
  )
) %>%
  mutate(frac_pass = n_pass / n_total)

write_tsv(threshold_summary, file.path(plot_dir, "co_qc_threshold_summary.tsv"))

message("CO QC plots written to: ", plot_dir)
