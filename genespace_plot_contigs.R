# Redoing haplotype phased assembly for C_epithymum
# Module note: load genespace before running R

library(GENESPACE)
library(Biostrings)

# =========================
# User settings
# =========================

CHR_N <- 12   # number of chromosomes / scaffolds to plot
TELO_KMERS <- c("TTTAGGG", "CCCTAAA")

PARAMS <- list(
  minContigGapSize = 100,
  maxDistBtwTelo   = 20,
  minTeloSize      = 200,
  minTeloDens      = 0.75,
  minChrSize       = 0,
  maxDist2end      = 10000,
  verbose          = TRUE
)

PLOT_PALETTE <- colorRampPalette(c("blue", "green"))
PLOT_NCOLORS <- 4

# =========================
# Helper functions
# =========================

run_telo_analysis <- function(asm, chr_n = CHR_N,
                              telo_kmers = TELO_KMERS,
                              params = PARAMS,
                              palette = PLOT_PALETTE,
                              nColors = PLOT_NCOLORS) {
  
  message("Reading assembly: ", asm)
  dnass <- readDNAStringSet(asm)
  
  result <- find_contigsGapsTelos(
    dnass = dnass,
    teloKmers = telo_kmers,
    minContigGapSize = params$minContigGapSize,
    maxDistBtwTelo   = params$maxDistBtwTelo,
    minTeloSize      = params$minTeloSize,
    minTeloDens      = params$minTeloDens,
    minChrSize       = params$minChrSize,
    maxDist2end      = params$maxDist2end,
    verbose          = params$verbose
  )
  
  valid_scaffolds <- paste0("scaffold_", seq_len(chr_n))
  
  filtered_result <- lapply(result, function(gr) {
    gr[seqnames(gr) %in% valid_scaffolds]
  })
  
  plot_contigs(
    cgt = filtered_result,
    nColors = nColors,
    palette = palette
  )
  
  return(result)
}

get_large_telomeres <- function(result, min_width = 1000) {
  telo_df <- as.data.frame(result$telo)
  telo_df[telo_df$width > min_width, ]
  #sort by highest width, then sort by chromosome
  telo_df <- telo_df[order(-telo_df$width, telo_df$seqnames), ]
}

get_scaffold_telomeres <- function(result, scaffold_name) {
  telo_df <- as.data.frame(result$telo)
  telo_df[telo_df$seqnames == scaffold_name, ]
}

collapse_nearby_telomeres <- function(df, max_gap = 1000) {
  # Expect columns: seqnames, start, end, width, strand, position
  
  if (nrow(df) == 0) return(df)
  
  # sort first
  df <- df[order(df$seqnames, df$start, df$end), ]
  rownames(df) <- NULL
  
  # gap from previous row within same scaffold
  same_seq <- c(FALSE, df$seqnames[-1] == df$seqnames[-nrow(df)])
  gap_from_prev <- c(Inf, df$start[-1] - df$end[-nrow(df)] - 1)
  gap_from_prev[!same_seq] <- Inf
  
  # start a new group when gap is too large or scaffold changes
  group_id <- cumsum(c(1, gap_from_prev[-1] > max_gap))
  
  # collapse each group
  out <- do.call(rbind, lapply(split(df, group_id), function(x) {
    data.frame(
      seqnames = x$seqnames[1],
      start = min(x$start),
      end = max(x$end),
      width = sum(x$width),              # sum of individual widths
      span = max(x$end) - min(x$start) + 1,  # genomic span of merged block
      strand = if ("strand" %in% names(x)) x$strand[1] else NA,
      position = if ("position" %in% names(x)) {
        if (length(unique(x$position)) == 1) unique(x$position) else "mixed"
      } else NA,
      n_merged = nrow(x)
    )
  }))
  
  rownames(out) <- NULL
  out
}

# =========================
# Single runs
# =========================
asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/assemblies/C_epithymum_haphic_mar2026/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_v0.FINAL.fa"
result <- run_telo_analysis(asm,chr_n = 14)

#check large telomeres
df <- get_large_telomeres(result, min_width = 2e3)
get_scaffold_telomeres(result, "scaffold_1")
collapse_nearby_telomeres(get_scaffold_telomeres(result, "scaffold_2"), max_gap = 5e2)
