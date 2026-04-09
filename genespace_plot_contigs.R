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
}

get_scaffold_telomeres <- function(result, scaffold_name) {
  telo_df <- as.data.frame(result$telo)
  telo_df[telo_df$seqnames == scaffold_name, ]
}

# =========================
# Single runs
# =========================

# Example 1: paradoxa older assembly
asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/assemblies/C_epithymum_haphic_mar2026/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_v0.FINAL.fa"
result <- run_telo_analysis(asm, chr_n = 14)


# Telomeres over 1 kb
large_telos <- get_large_telomeres(result, min_width = 1000)
print(large_telos)

# =========================
# Multiple paradoxa versions
# =========================

paradoxa_versions <- c(
  "V1.2" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.2.FINAL.fa",
  "V1.3" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.3.FINAL.fa",
  "V1.4" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.4.FINAL.fa",
  "V1.5" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.5.FINAL.fa",
  "V1.6" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.6.FINAL.fa",
  "V1.7" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.7.FINAL.fa",
  "V1.8" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V1.8.FINAL.fa",
  "V2.0" = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/results/assemblies/Drosera_paradoxa_assembly_hic_6181Ndeeplight_6986Adeeplight/HiC_scaffolding/haphic/05.post_juicebox/out_JBAT_V2.FINAL.fa"
)

results_list <- list()

for (ver in names(paradoxa_versions)) {
  cat("\n====================\n")
  cat("Running:", ver, "\n")
  cat("====================\n")
  
  results_list[[ver]] <- run_telo_analysis(
    asm = paradoxa_versions[[ver]],
    chr_n = CHR_N
  )
}

# =========================
# Optional scaffold checks
# =========================

get_scaffold_telomeres(results_list[["V1.5"]], "scaffold_7")
get_scaffold_telomeres(results_list[["V1.6"]], "scaffold_9")
get_scaffold_telomeres(results_list[["V1.7"]], "scaffold_11")
get_scaffold_telomeres(results_list[["V1.8"]], "scaffold_12")

# =========================
# Scorpioides assembly
# =========================

asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/assemblies/Drosera_scorpioides_HiC_2025_04_29_saurab_no_hom_cov/newres/HiC_scaffolding/ragtag_putg_hap1hap2/ragtag_output/juicebox_ragtag/HapHiC/haphic_with_params/05.post_juicebox/out_JBAT.FINAL.fa"
result_scorpioides <- run_telo_analysis(asm, chr_n = 16)
