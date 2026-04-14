#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript plot_final_qc.R <project_root>")
}

project_root <- normalizePath(args[1], mustWork = TRUE)

barcode_qc_file    <- file.path(project_root, "results/03_barcode_qc/barcode_read_qc.tsv")
demux_summary_file <- file.path(project_root, "results/03_barcode_qc/demux_summary.tsv")
switches_file      <- file.path(project_root, "results/06_markers/switches.stats.tsv")
co_summary_file    <- file.path(project_root, "results/07_crossovers/co_summary.tsv")
co_intervals_file  <- file.path(project_root, "results/07_crossovers/co_intervals.bed")

plot_dir <- file.path(project_root, "results/08_plots/qc_final")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

THRESH_TOTAL_READS <- 10000
THRESH_UNIQ_READS  <- 5000
THRESH_UNIQ_RATIO  <- 0.2
THRESH_MARKERS     <- 500
THRESH_SWITCH_RATE <- 0.10

save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(file.path(plot_dir, filename), p, width = width, height = height, dpi = 300)
}

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
  mutate(width_bp = end - start + 1)

co_summary_fixed <- co_summary %>%
  mutate(barcode = paste0(prefix, "-1"))

qc_join <- demux_summary %>%
  left_join(switches, by = "barcode") %>%
  left_join(co_summary_fixed, by = "barcode")

# 1 total reads
p1 <- ggplot(barcode_qc, aes(total_reads)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_TOTAL_READS, linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Total reads per barcode",
    subtitle = paste0("Threshold = ", THRESH_TOTAL_READS,
                      " | kept = ", sum(barcode_qc$total_reads >= THRESH_TOTAL_READS), "/", nrow(barcode_qc)),
    x = "Total reads per barcode (log10)",
    y = "Count"
  )
save_plot(p1, "01_total_reads_per_barcode.png")

# 2 barcode rank plot
barcode_rank <- barcode_qc %>%
  arrange(desc(total_reads)) %>%
  mutate(rank = row_number())

p2 <- ggplot(barcode_rank, aes(rank, total_reads)) +
  geom_point(size = 0.4) +
  geom_hline(yintercept = THRESH_TOTAL_READS, linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    title = "Barcode rank plot",
    x = "Barcode rank (log10)",
    y = "Total reads (log10)"
  )
save_plot(p2, "02_barcode_rank_plot.png")

# 3 uniq mapped reads
p3 <- ggplot(demux_summary, aes(uniq_map_reads)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_UNIQ_READS, linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Uniquely mapped reads per barcode",
    subtitle = paste0("Threshold = ", THRESH_UNIQ_READS,
                      " | kept = ", sum(demux_summary$uniq_map_reads >= THRESH_UNIQ_READS), "/", nrow(demux_summary)),
    x = "Uniquely mapped reads (log10)",
    y = "Count"
  )
save_plot(p3, "03_uniq_mapped_reads.png")

# 4 uniq ratio
p4 <- ggplot(demux_summary, aes(uniq_ratio)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_UNIQ_RATIO, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Unique mapping ratio",
    subtitle = paste0("Threshold = ", THRESH_UNIQ_RATIO,
                      " | kept = ", sum(demux_summary$uniq_ratio >= THRESH_UNIQ_RATIO), "/", nrow(demux_summary)),
    x = "uniq_map_reads / total_reads",
    y = "Count"
  )
save_plot(p4, "04_unique_mapping_ratio.png")

# 5 marker num
p5 <- ggplot(switches, aes(marker_num)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_MARKERS, linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Marker count per barcode",
    subtitle = paste0("Threshold = ", THRESH_MARKERS,
                      " | kept = ", sum(switches$marker_num >= THRESH_MARKERS), "/", nrow(switches)),
    x = "Marker count (log10)",
    y = "Count"
  )
save_plot(p5, "05_marker_count.png")

# 6 switch num
p6 <- ggplot(switches, aes(switch_num)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Switch count per barcode",
    x = "Switch count (log10)",
    y = "Count"
  )
save_plot(p6, "06_switch_count.png")

# 7 switch rate
p7 <- ggplot(switches, aes(switch_rate)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = THRESH_SWITCH_RATE, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Switch rate per barcode",
    subtitle = paste0("Threshold = ", THRESH_SWITCH_RATE,
                      " | kept = ", sum(switches$switch_rate <= THRESH_SWITCH_RATE), "/", nrow(switches)),
    x = "Switch rate",
    y = "Count"
  )
save_plot(p7, "07_switch_rate.png")

# 8 marker vs switch rate
p8 <- ggplot(switches, aes(marker_num, switch_rate)) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_vline(xintercept = THRESH_MARKERS, linetype = "dashed") +
  geom_hline(yintercept = THRESH_SWITCH_RATE, linetype = "dashed") +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Marker count vs switch rate",
    x = "Marker count (log10)",
    y = "Switch rate"
  )
save_plot(p8, "08_marker_vs_switch_rate.png")

# 9 COs per cell
p9 <- ggplot(co_summary, aes(n_crossovers)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  theme_bw() +
  labs(
    title = "Crossovers per selected cell",
    subtitle = paste0("Mean = ", round(mean(co_summary$n_crossovers), 2),
                      " | n = ", nrow(co_summary)),
    x = "Number of crossovers",
    y = "Count"
  )
save_plot(p9, "09_crossovers_per_cell.png")

# 10 CO interval widths
p10 <- ggplot(co_intervals, aes(width_bp / 1e6)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(
    title = "CO interval widths",
    x = "Interval width (Mbp)",
    y = "Count"
  )
save_plot(p10, "10_co_interval_widths.png")

# 11 marker count vs COs
p11 <- ggplot(qc_join %>% filter(!is.na(n_crossovers)), aes(marker_num, n_crossovers)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Marker count vs crossovers",
    x = "Marker count (log10)",
    y = "Number of crossovers"
  )
save_plot(p11, "11_marker_count_vs_crossovers.png")

# 12 switch rate vs COs
p12 <- ggplot(qc_join %>% filter(!is.na(n_crossovers)), aes(switch_rate, n_crossovers)) +
  geom_point(alpha = 0.4, size = 0.8) +
  theme_bw() +
  labs(
    title = "Switch rate vs crossovers",
    x = "Switch rate",
    y = "Number of crossovers"
  )
save_plot(p12, "12_switch_rate_vs_crossovers.png")

# Summary TSV
qc_summary <- tibble(
  metric = c("all_barcodes", "reads_ge_1e4", "uniq_reads_ge_5e3", "uniq_ratio_ge_0.2",
             "marker_num_ge_500", "switch_rate_le_0.10", "selected_for_co"),
  n = c(
    nrow(barcode_qc),
    sum(barcode_qc$total_reads >= THRESH_TOTAL_READS),
    sum(demux_summary$uniq_map_reads >= THRESH_UNIQ_READS),
    sum(demux_summary$uniq_ratio >= THRESH_UNIQ_RATIO),
    sum(switches$marker_num >= THRESH_MARKERS),
    sum(switches$switch_rate <= THRESH_SWITCH_RATE),
    nrow(co_summary)
  )
)

write_tsv(qc_summary, file.path(plot_dir, "qc_summary.tsv"))

message("QC plots written to: ", plot_dir)
