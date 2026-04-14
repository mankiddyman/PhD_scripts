#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: plot_qc_distributions.R <project_root>")
}

project_root <- normalizePath(args[1], mustWork = TRUE)

barcode_qc_file   <- file.path(project_root, "results/03_barcode_qc/barcode_read_qc.tsv")
demux_summary_file <- file.path(project_root, "results/03_barcode_qc/demux_summary.tsv")
switches_file     <- file.path(project_root, "results/06_markers/switches.stats.tsv")
plots_dir         <- file.path(project_root, "results/08_plots/qc")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# thresholds currently in use
THRESH_TOTAL_READS <- 10000
THRESH_MAPQ_UNIQ   <- 3
THRESH_UNIQ_READS  <- 5000
THRESH_UNIQ_RATIO  <- 0.2
THRESH_MARKERS     <- 500
THRESH_SWITCH_RATE <- 0.10

message("Reading input tables...")

barcode_qc <- read_tsv(
  barcode_qc_file,
  show_col_types = FALSE,
  col_types = cols(
    barcode = col_character(),
    total_reads = col_double()
  )
)

demux_summary <- read_tsv(
  demux_summary_file,
  show_col_types = FALSE,
  col_types = cols(
    barcode = col_character(),
    total_reads = col_double(),
    uniq_map_reads = col_double(),
    uniq_ratio = col_double()
  )
)

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

# Join stage-03 and stage-04 summaries
qc_all <- demux_summary %>%
  left_join(switches, by = "barcode")

# helper to save plots
save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(plots_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

# helper for threshold summary
threshold_summary <- tibble(
  metric = c(
    "total_reads",
    "uniq_map_reads",
    "uniq_ratio",
    "marker_num",
    "switch_rate"
  ),
  threshold = c(
    THRESH_TOTAL_READS,
    THRESH_UNIQ_READS,
    THRESH_UNIQ_RATIO,
    THRESH_MARKERS,
    THRESH_SWITCH_RATE
  ),
  n_total = c(
    nrow(barcode_qc),
    nrow(demux_summary),
    nrow(demux_summary),
    nrow(switches),
    nrow(switches)
  ),
  n_pass = c(
    sum(barcode_qc$total_reads >= THRESH_TOTAL_READS, na.rm = TRUE),
    sum(demux_summary$uniq_map_reads >= THRESH_UNIQ_READS, na.rm = TRUE),
    sum(demux_summary$uniq_ratio >= THRESH_UNIQ_RATIO, na.rm = TRUE),
    sum(switches$marker_num >= THRESH_MARKERS, na.rm = TRUE),
    sum(switches$switch_rate <= THRESH_SWITCH_RATE, na.rm = TRUE)
  )
) %>%
  mutate(frac_pass = n_pass / n_total)

write_tsv(threshold_summary, file.path(plots_dir, "qc_threshold_summary.tsv"))

# 1) total reads / barcode
p1 <- ggplot(barcode_qc, aes(x = total_reads)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_TOTAL_READS, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Total reads per barcode",
    subtitle = paste0(
      "Threshold = ", THRESH_TOTAL_READS,
      " | kept = ", sum(barcode_qc$total_reads >= THRESH_TOTAL_READS),
      " / ", nrow(barcode_qc)
    ),
    x = "Total reads per barcode (log10 scale)",
    y = "Number of barcodes"
  ) +
  theme_bw()

save_plot(p1, "01_total_reads_per_barcode_hist.png")

# zoomed version
p1z <- ggplot(barcode_qc %>% filter(total_reads <= 100000), aes(x = total_reads)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_TOTAL_READS, linetype = "dashed") +
  labs(
    title = "Total reads per barcode (zoomed)",
    subtitle = paste0(
      "Threshold = ", THRESH_TOTAL_READS,
      " | kept = ", sum(barcode_qc$total_reads >= THRESH_TOTAL_READS),
      " / ", nrow(barcode_qc)
    ),
    x = "Total reads per barcode",
    y = "Number of barcodes"
  ) +
  theme_bw()

save_plot(p1z, "02_total_reads_per_barcode_hist_zoom.png")

# barcode rank plot
barcode_rank <- barcode_qc %>%
  arrange(desc(total_reads)) %>%
  mutate(rank = row_number())

p_rank <- ggplot(barcode_rank, aes(x = rank, y = total_reads)) +
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = THRESH_TOTAL_READS, linetype = "dashed") +
  labs(
    title = "Barcode rank plot",
    x = "Barcode rank (log10 scale)",
    y = "Total reads (log10 scale)"
  ) +
  theme_bw()

save_plot(p_rank, "03_barcode_rank_plot.png")

# 2) uniq_map_reads distribution
p2 <- ggplot(demux_summary, aes(x = uniq_map_reads)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_UNIQ_READS, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Uniquely mapped reads per selected barcode",
    subtitle = paste0(
      "MAPQ unique definition currently: MAPQ > ", THRESH_MAPQ_UNIQ,
      " | threshold = ", THRESH_UNIQ_READS,
      " | kept = ", sum(demux_summary$uniq_map_reads >= THRESH_UNIQ_READS),
      " / ", nrow(demux_summary)
    ),
    x = "Uniquely mapped reads (log10 scale)",
    y = "Number of barcodes"
  ) +
  theme_bw()

save_plot(p2, "04_uniq_map_reads_hist.png")

# 3) uniq ratio distribution
p3 <- ggplot(demux_summary, aes(x = uniq_ratio)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_UNIQ_RATIO, linetype = "dashed") +
  labs(
    title = "Unique mapping ratio",
    subtitle = paste0(
      "Threshold = ", THRESH_UNIQ_RATIO,
      " | kept = ", sum(demux_summary$uniq_ratio >= THRESH_UNIQ_RATIO),
      " / ", nrow(demux_summary)
    ),
    x = "uniq_map_reads / total_reads",
    y = "Number of barcodes"
  ) +
  theme_bw()

save_plot(p3, "05_uniq_ratio_hist.png")

# 4) marker_num distribution
p4 <- ggplot(switches, aes(x = marker_num)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_MARKERS, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Marker count per barcode",
    subtitle = paste0(
      "Threshold = ", THRESH_MARKERS,
      " | kept = ", sum(switches$marker_num >= THRESH_MARKERS),
      " / ", nrow(switches)
    ),
    x = "Marker count (log10 scale)",
    y = "Number of barcodes"
  ) +
  theme_bw()

save_plot(p4, "06_marker_num_hist.png")

# 5) switch_rate distribution
p5 <- ggplot(switches, aes(x = switch_rate)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_SWITCH_RATE, linetype = "dashed") +
  labs(
    title = "Switch rate per barcode",
    subtitle = paste0(
      "Threshold = ", THRESH_SWITCH_RATE,
      " | kept = ", sum(switches$switch_rate <= THRESH_SWITCH_RATE),
      " / ", nrow(switches)
    ),
    x = "switch_num / marker_num",
    y = "Number of barcodes"
  ) +
  theme_bw()

save_plot(p5, "07_switch_rate_hist.png")

# 6) marker_num vs switch_rate
p6 <- ggplot(switches, aes(x = marker_num, y = switch_rate)) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_vline(xintercept = THRESH_MARKERS, linetype = "dashed") +
  geom_hline(yintercept = THRESH_SWITCH_RATE, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Marker count vs switch rate",
    x = "Marker count (log10 scale)",
    y = "Switch rate"
  ) +
  theme_bw()

save_plot(p6, "08_marker_num_vs_switch_rate.png")

# QC funnel
funnel <- tibble(
  step = c(
    "All barcodes",
    "Reads >= 1e4",
    "uniq_map_reads >= 5e3",
    "uniq_ratio >= 0.2",
    "marker_num >= 500",
    "switch_rate <= 0.10"
  ),
  n = c(
    nrow(barcode_qc),
    sum(barcode_qc$total_reads >= THRESH_TOTAL_READS, na.rm = TRUE),
    sum(demux_summary$uniq_map_reads >= THRESH_UNIQ_READS, na.rm = TRUE),
    sum(demux_summary$uniq_ratio >= THRESH_UNIQ_RATIO, na.rm = TRUE),
    sum(switches$marker_num >= THRESH_MARKERS, na.rm = TRUE),
    sum(switches$switch_rate <= THRESH_SWITCH_RATE, na.rm = TRUE)
  )
)

write_tsv(funnel, file.path(plots_dir, "qc_funnel.tsv"))

p_funnel <- ggplot(funnel, aes(x = reorder(step, -n), y = n)) +
  geom_col() +
  labs(
    title = "QC retention summary",
    x = "",
    y = "Number retained"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p_funnel, "09_qc_funnel.png", width = 9, height = 6)

message("QC plots written to: ", plots_dir)
