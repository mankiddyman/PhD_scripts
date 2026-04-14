#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------
# argument parsing
# ------------------------------------------------------------

parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "-")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    if (i == length(x)) {
      stop(sprintf("Missing value for argument: %s", key))
    }
    val <- x[[i + 1]]
    out[[sub("^-+", "", key)]] <- val
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("i", "p", "g", "o", "c", "s", "n")
missing_args <- setdiff(required, names(args))
if (length(missing_args) > 0) {
  stop(
    paste0(
      "Missing required arguments: ",
      paste(missing_args, collapse = ", "),
      "\nUsage: Rscript hapCO_identification.R ",
      "-i <input_corrected.txt> ",
      "-p <prefix> ",
      "-g <genome.fai> ",
      "-o <output_dir> ",
      "-c <min_markers_per_cell> ",
      "-s <min_segment_bp> ",
      "-n <min_markers_per_segment>"
    )
  )
}

input_file <- normalizePath(args$i, mustWork = TRUE)
prefix <- args$p
genome_fai <- normalizePath(args$g, mustWork = TRUE)
output_dir <- args$o
min_markers_cell <- as.integer(args$c)
min_segment_bp <- as.integer(args$s)
min_markers_segment <- as.integer(args$n)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# output paths
# ------------------------------------------------------------

marker_out <- file.path(output_dir, paste0(prefix, "_markers_annotated.tsv"))
segment_out <- file.path(output_dir, paste0(prefix, "_segments.tsv"))
co_out <- file.path(output_dir, paste0(prefix, "_allele_cnts_at_markers_sorted_co_pred.txt"))
summary_out <- file.path(output_dir, paste0(prefix, "_co_summary.tsv"))

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

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

empty_outputs <- function(reason, informative_markers = 0L, valid_segments = 0L, n_cos = 0L) {
  readr::write_tsv(
    tibble(
      prefix = prefix,
      informative_markers = informative_markers,
      valid_segments = valid_segments,
      n_crossovers = n_cos,
      min_markers_cell = min_markers_cell,
      min_segment_bp = min_segment_bp,
      min_markers_segment = min_markers_segment,
      status = reason
    ),
    summary_out
  )

  readr::write_tsv(
    tibble(
      chrom = character(),
      pos = integer(),
      ref = character(),
      ref_count = integer(),
      alt = character(),
      alt_count = integer(),
      total_count = integer(),
      ref_fraction = double(),
      state = double(),
      informative = logical()
    ),
    marker_out
  )

  readr::write_tsv(
    tibble(
      chrom = character(),
      segment_id = integer(),
      state = double(),
      start_pos = integer(),
      end_pos = integer(),
      n_markers = integer(),
      span_bp = integer(),
      keep_segment = logical()
    ),
    segment_out
  )

  readr::write_tsv(
    tibble(
      chrom = character(),
      co_start = integer(),
      co_end = integer(),
      co_midpoint = double(),
      left_segment_id = integer(),
      right_segment_id = integer(),
      left_state = double(),
      right_state = double(),
      left_n_markers = integer(),
      right_n_markers = integer(),
      left_span_bp = integer(),
      right_span_bp = integer(),
      inter_segment_gap_bp = integer()
    ),
    co_out
  )
}

# ------------------------------------------------------------
# read inputs
# ------------------------------------------------------------

genome_tbl <- read_tsv(
  genome_fai,
  col_names = c("chrom", "length", "offset", "line_bases", "line_width"),
  show_col_types = FALSE
) %>%
  select(chrom, length)

dat <- read_table(
  input_file,
  col_names = c("chrom", "pos", "ref", "ref_count", "alt", "alt_count"),
  show_col_types = FALSE
) %>%
  mutate(
    pos = as.integer(pos),
    ref_count = as.integer(ref_count),
    alt_count = as.integer(alt_count)
  ) %>%
  inner_join(genome_tbl, by = "chrom") %>%
  arrange(chrom, pos) %>%
  mutate(
    total_count = ref_count + alt_count,
    ref_fraction = if_else(total_count > 0, ref_count / total_count, NA_real_),
    informative = total_count > 0,
    state = state_from_counts(ref_count, alt_count)
  ) %>%
  select(-length)

# always write annotated marker table
write_tsv(dat, marker_out)

informative_dat <- dat %>%
  filter(informative)

if (nrow(informative_dat) < min_markers_cell) {
  empty_outputs(
    reason = "too_few_informative_markers",
    informative_markers = nrow(informative_dat),
    valid_segments = 0L,
    n_cos = 0L
  )
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------
# build segments
# ------------------------------------------------------------

# for segment construction, use only confident haplotype states (0 or 1)
seg_input <- informative_dat %>%
  filter(!is.na(state), state %in% c(0, 1)) %>%
  arrange(chrom, pos)

if (nrow(seg_input) == 0) {
  empty_outputs(
    reason = "no_nonambiguous_markers",
    informative_markers = nrow(informative_dat),
    valid_segments = 0L,
    n_cos = 0L
  )
  quit(save = "no", status = 0)
}

seg_input <- seg_input %>%
  mutate(
    is_new_segment = if_else(
      row_number() == 1L |
        chrom != lag(chrom) |
        state != lag(state),
      TRUE, FALSE
    ),
    segment_id = cumsum(is_new_segment)
  )

segments <- seg_input %>%
  group_by(chrom, segment_id, state) %>%
  summarise(
    start_pos = min(pos),
    end_pos = max(pos),
    n_markers = n(),
    span_bp = end_pos - start_pos + 1L,
    .groups = "drop"
  ) %>%
  mutate(
    keep_segment = n_markers >= min_markers_segment & span_bp >= min_segment_bp
  )

write_tsv(segments, segment_out)

valid_segments <- segments %>%
  filter(keep_segment) %>%
  arrange(chrom, start_pos)

if (nrow(valid_segments) < 2) {
  empty_outputs(
    reason = "fewer_than_two_valid_segments",
    informative_markers = nrow(informative_dat),
    valid_segments = nrow(valid_segments),
    n_cos = 0L
  )
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------
# call COs from boundaries between adjacent valid segments
# ------------------------------------------------------------

left_seg <- valid_segments %>%
  mutate(
    left_segment_id = segment_id,
    left_state = state,
    left_start_pos = start_pos,
    left_end_pos = end_pos,
    left_n_markers = n_markers,
    left_span_bp = span_bp
  ) %>%
  select(
    chrom,
    left_segment_id,
    left_state,
    left_start_pos,
    left_end_pos,
    left_n_markers,
    left_span_bp
  )

right_seg <- valid_segments %>%
  mutate(
    right_segment_id = segment_id,
    right_state = state,
    right_start_pos = start_pos,
    right_end_pos = end_pos,
    right_n_markers = n_markers,
    right_span_bp = span_bp
  ) %>%
  select(
    chrom,
    right_segment_id,
    right_state,
    right_start_pos,
    right_end_pos,
    right_n_markers,
    right_span_bp
  )

co_tbl <- tibble(
  chrom = valid_segments$chrom[-nrow(valid_segments)],
  left_segment_id = valid_segments$segment_id[-nrow(valid_segments)],
  right_segment_id = valid_segments$segment_id[-1],
  left_state = valid_segments$state[-nrow(valid_segments)],
  right_state = valid_segments$state[-1],
  left_end_pos = valid_segments$end_pos[-nrow(valid_segments)],
  right_start_pos = valid_segments$start_pos[-1],
  left_n_markers = valid_segments$n_markers[-nrow(valid_segments)],
  right_n_markers = valid_segments$n_markers[-1],
  left_span_bp = valid_segments$span_bp[-nrow(valid_segments)],
  right_span_bp = valid_segments$span_bp[-1]
) %>%
  filter(
    chrom == lead(chrom, default = NA_character_) | TRUE
  )

# keep only true adjacent same-chromosome state flips
co_tbl <- valid_segments %>%
  arrange(chrom, start_pos) %>%
  mutate(next_chrom = lead(chrom),
         next_segment_id = lead(segment_id),
         next_state = lead(state),
         next_start_pos = lead(start_pos),
         next_n_markers = lead(n_markers),
         next_span_bp = lead(span_bp)) %>%
  filter(
    !is.na(next_segment_id),
    chrom == next_chrom,
    state != next_state
  ) %>%
  transmute(
    chrom = chrom,
    co_start = end_pos,
    co_end = next_start_pos,
    co_midpoint = (end_pos + next_start_pos) / 2,
    left_segment_id = segment_id,
    right_segment_id = next_segment_id,
    left_state = state,
    right_state = next_state,
    left_n_markers = n_markers,
    right_n_markers = next_n_markers,
    left_span_bp = span_bp,
    right_span_bp = next_span_bp,
    inter_segment_gap_bp = next_start_pos - end_pos
  )

write_tsv(co_tbl, co_out)

# ------------------------------------------------------------
# summary
# ------------------------------------------------------------

write_tsv(
  tibble(
    prefix = prefix,
    informative_markers = nrow(informative_dat),
    nonambiguous_markers = nrow(seg_input),
    valid_segments = nrow(valid_segments),
    n_crossovers = nrow(co_tbl),
    min_markers_cell = min_markers_cell,
    min_segment_bp = min_segment_bp,
    min_markers_segment = min_markers_segment,
    status = "ok"
  ),
  summary_out
)
