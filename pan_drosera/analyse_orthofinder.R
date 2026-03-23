# analysing the orthofinder results from the pandrosera work, objective is to  the copy number and use it to predict holocentricity


plots <- list()

cap <- function(expr) {
  val <- eval.parent(substitute(expr))

  # If it is a ggplot-like object, explicitly draw it
  if (inherits(val, c("gg", "ggplot", "grob", "gtable"))) {
    print(val)
  }

  # Record whatever is now on the current graphics device
  plots[[length(plots) + 1]] <<- recordPlot()

  invisible(val)
}
# importing
wd <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/pan_drosera/peptide/OrthoFinder/modelling"
# make and cd
if(!dir.exists(wd)) {
  dir.create(wd)
}
setwd(wd)


#ok now importing, because the fuckin orthofinder run failed I must be a bit creative, i have an orthogroups.txt which has the name of the orthogroup and just the name of the gene, so i need to manually first create vectors of gene names for each species
# then mimick the creation of the copy number table, then i can begin modelling

library(data.table)

## =========================================================
## 0) Paths
## =========================================================

peptide_dir <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/pan_drosera/peptide"
orthogroups_file <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/pan_drosera/peptide/OrthoFinder/Results_Mar17/Orthogroups/Orthogroups.txt"

library(data.table)

## paths
orthofinder_base <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/pan_drosera/peptide/OrthoFinder/Results_Mar17"
wd_dir <- file.path(orthofinder_base, "WorkingDirectory")

species_names <- c(
  "D_aliciae", "D_binata", "D_capensis",
  "D_paradoxa", "D_regia", "D_scorpioides"
)

## species id -> species
species_id_dt <- fread(
  file.path(wd_dir, "SpeciesIDs.txt"),
  sep = ":",
  header = FALSE,
  col.names = c("species_id", "species_file"),
  strip.white = TRUE
)[
  , .(
    species_id = as.integer(species_id),
    species = sub("\\.from_gff3\\.pep\\.fa$", "", species_file)
  )
]

## seq id -> gene + species
seq_map_dt <- fread(
  file.path(wd_dir, "SequenceIDs.txt"),
  sep = ":",
  header = FALSE,
  col.names = c("seq_id", "gene"),
  strip.white = TRUE
)[
  , species_id := as.integer(sub("_.*$", "", seq_id))
][
  species_id_dt, on = "species_id"
][
  , species_gene := paste(species, gene, sep = "::")
]

## parse MCL cluster file
cluster_lines <- readLines(
  file.path(wd_dir, "clusters_OrthoFinder_I1.2.txt_id_pairs.txt"),
  warn = FALSE
)

cluster_lines <- trimws(cluster_lines)
cluster_lines <- cluster_lines[(match("begin", cluster_lines) + 1):length(cluster_lines)]
cluster_lines <- cluster_lines[nzchar(cluster_lines) & cluster_lines != ")"]

tokens <- unlist(strsplit(paste(cluster_lines, collapse = " "), "[[:space:]]+"))
end_idx <- which(tokens == "$")
start_idx <- c(1, head(end_idx + 1, -1))

cluster_long_dt <- rbindlist(
  Map(function(s, e) {
    x <- tokens[s:(e - 1)]
    data.table(cluster_idx = x[1], seq_id = x[-1])
  }, start_idx, end_idx)
)

cluster_ids <- unique(cluster_long_dt$cluster_idx)
cluster_long_dt[, Orthogroup := sprintf("OG%07d", match(cluster_idx, cluster_ids) - 1L)]

## final rescued gene-level orthogroup membership
gene_long <- cluster_long_dt[
  seq_map_dt,
  on = "seq_id",
  nomatch = 0
][
  , .(Orthogroup, seq_id, gene, species, species_gene)
]

## wide orthogroup copy-number table
wide_dt <- dcast(
  gene_long[, .N, by = .(Orthogroup, species)],
  Orthogroup ~ species,
  value.var = "N",
  fill = 0
)

for (sp in species_names) {
  if (!sp %in% names(wide_dt)) wide_dt[, (sp) := 0L]
}

wide_dt[, Total := rowSums(.SD), .SDcols = species_names]
setcolorder(wide_dt, c("Orthogroup", species_names, "Total"))
setkey(wide_dt, Orthogroup)

## quick assertions
stopifnot(nrow(gene_long) == nrow(seq_map_dt))
stopifnot(sum(is.na(gene_long$gene)) == 0)
stopifnot(sum(is.na(gene_long$species)) == 0)
stopifnot(nrow(wide_dt) == uniqueN(gene_long$Orthogroup))

head(wide_dt)

# now make a df for holocentricity

holo_df <- data.frame(
  species = c("D_aliciae","D_binata","D_capensis","D_paradoxa","D_regia","D_scorpioides"),
  holocentric = c(0,0,0,1,1,1)
)

ploidy_df <- data.frame(
  species = c("D_aliciae","D_binata","D_capensis","D_paradoxa","D_regia","D_scorpioides"),
  ploidy = c(1,2,2,2,1,4)
)

# ============================================================
# Holocentricity screen from orthogroup copy number
# Ploidy-adjusted, exact permutation test, fast/vectorized
# ============================================================

library(data.table)
library(matrixStats)
library(ggplot2)
library(pheatmap)

# ------------------------------------------------------------
# 0) INPUTS
# wide_dt: data.table with columns:
#   Orthogroup, D_aliciae, D_binata, D_capensis, D_paradoxa, D_regia, D_scorpioides, Total
#
# holo_df:
# species, holocentric
#
# ploidy_df:
# species, ploidy
# ------------------------------------------------------------

# Example metadata (use yours as given)
holo_df <- data.frame(
  species = c("D_aliciae","D_binata","D_capensis","D_paradoxa","D_regia","D_scorpioides"),
  holocentric = c(0,0,0,1,1,1)
)

ploidy_df <- data.frame(
  species = c("D_aliciae","D_binata","D_capensis","D_paradoxa","D_regia","D_scorpioides"),
  ploidy = c(1,2,2,2,1,4)
)

# ------------------------------------------------------------
# 1) ALIGN SPECIES + BASIC CHECKS
# ------------------------------------------------------------

meta_df <- merge(holo_df, ploidy_df, by = "species", all = FALSE)

species_cols <- intersect(meta_df$species, colnames(wide_dt))
meta_df <- meta_df[match(species_cols, meta_df$species), ]

cat("========== BASIC CHECKS ==========\n")
cat("Species found in wide_dt:", paste(species_cols, collapse = ", "), "\n")
cat("Number of species:", length(species_cols), "\n")
cat("Holocentric counts:\n")
print(table(meta_df$holocentric))
cat("Ploidy values:\n")
print(meta_df[, c("species", "ploidy")])
cat("\n")

stopifnot(length(species_cols) == 6)
stopifnot(all(meta_df$species == species_cols))

# ------------------------------------------------------------
# 2) BUILD COPY NUMBER MATRIX
# ------------------------------------------------------------

# Keep Orthogroup IDs and species matrix only
og_ids <- wide_dt$Orthogroup
X_raw <- as.matrix(wide_dt[, ..species_cols])
storage.mode(X_raw) <- "numeric"

cat("========== MATRIX SUMMARY ==========\n")
cat("Orthogroups:", nrow(X_raw), "\n")
cat("Species columns:", ncol(X_raw), "\n")
cat("Raw copy number range:", paste(range(X_raw, na.rm = TRUE), collapse = " to "), "\n")
cat("\n")

# ------------------------------------------------------------
# 3) PLOIDY ADJUSTMENT + LOG TRANSFORM
# ------------------------------------------------------------

ploidy_vec <- meta_df$ploidy
names(ploidy_vec) <- meta_df$species

# Divide each species column by its ploidy
X_ploidy_adj <- sweep(X_raw, 2, ploidy_vec, "/")

# Log transform to stabilize huge values
X <- log1p(X_ploidy_adj)

#histogram of copy numbers to visualise skew faceted with unlogged
a <-ggplot(data.table(value = as.vector(X_raw), type = "raw"), aes(x = value)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(title = "Raw copy number distribution", x = "Copy number", y = "Count")
b<-ggplot(data.table(value = as.vector(X), type = "adjusted"), aes(x = value)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(title = "Ploidy-adjusted log copy number distribution", x = "log1p(copy/ploidy)", y = "Count")

library(gridExtra)
grid.arrange(a, b, nrow = 1)
cap(grid.arrange(a, b, nrow = 1))
cat("========== AFTER PLOIDY ADJUSTMENT ==========\n")
cat("Adjusted+log range:", paste(round(range(X, na.rm = TRUE), 4), collapse = " to "), "\n")
cat("\n")

# ------------------------------------------------------------
# 4) FILTER LOW-INFORMATION ORTHOGROUPS
# ------------------------------------------------------------

# useful filters:
# - total adjusted copy number > 0
# - variance > 0
# - present in at least 2 species
row_sum <- rowSums(X_ploidy_adj, na.rm = TRUE)
row_var <- rowVars(X, na.rm = TRUE)
row_present <- rowSums(X_raw > 0, na.rm = TRUE)

keep <- (row_sum > 0) & (row_var > 0) & (row_present >= 2)

Xf <- X[keep, , drop = FALSE]
Xf_raw <- X_raw[keep, , drop = FALSE]
og_ids_f <- og_ids[keep]

cat("========== FILTERING ==========\n")
cat("Filtering criteria:\n")
cat("- Total ploidy-adjusted copy number > 0\n")
cat("- Variance of log-adjusted copy number > 0\n")
cat("- Present in at least 2 species (raw count > 0)\n")

cat("Orthogroups before filtering:", nrow(X_raw), "\n")
cat("Orthogroups after filtering :", nrow(Xf), "\n")
cat("Removed:", nrow(X_raw) - nrow(Xf), "\n")
cat("\n")

# ------------------------------------------------------------
# 5) OBSERVED EFFECT SIZE
#    Delta = mean(holo) - mean(nonholo)
# ------------------------------------------------------------

holo_idx <- which(meta_df$holocentric == 1)
non_idx  <- which(meta_df$holocentric == 0)

obs_effect <- rowMeans(Xf[, holo_idx, drop = FALSE]) -
              rowMeans(Xf[, non_idx,  drop = FALSE])

# Also compute simple standardized effect for intuition
# (difference in means divided by SD across all species)
row_sd <- rowSds(Xf)
std_effect <- obs_effect / ifelse(row_sd == 0, NA, row_sd)

cat("========== EFFECT SIZE SUMMARY ==========\n")
print(summary(obs_effect))
cat("\nTop positive raw effects:\n")
print(head(sort(obs_effect, decreasing = TRUE), 10))
cat("\nTop negative raw effects:\n")
print(head(sort(obs_effect, decreasing = FALSE), 10))
cat("\n")

# ------------------------------------------------------------
# 6) EXACT PERMUTATION TEST
#    There are choose(6,3)=20 possible 3-vs-3 labelings
# ------------------------------------------------------------

perm_combos <- combn(seq_len(ncol(Xf)), length(holo_idx), simplify = FALSE)
n_perm <- length(perm_combos)

cat("========== EXACT PERMUTATION TEST ==========\n")
cat("Number of exact labelings considered:", n_perm, "\n")
cat("This is the full null space for 3-vs-3 grouping.\n\n")

# For speed, compute all permutation effects into a matrix
perm_effect_mat <- matrix(NA_real_, nrow = nrow(Xf), ncol = n_perm)

for (j in seq_along(perm_combos)) {
  h_idx <- perm_combos[[j]]
  nh_idx <- setdiff(seq_len(ncol(Xf)), h_idx)
  perm_effect_mat[, j] <- rowMeans(Xf[, h_idx, drop = FALSE]) -
                          rowMeans(Xf[, nh_idx, drop = FALSE])
}

# Two-sided exact p-value
# proportion of null effects at least as extreme as observed
abs_obs <- abs(obs_effect)
abs_null <- abs(perm_effect_mat)

p_exact <- rowMeans(abs_null >= abs_obs)

# BH FDR
q_bh <- p.adjust(p_exact, method = "BH")

# ------------------------------------------------------------
# 7) ASSEMBLE RESULTS TABLE
# ------------------------------------------------------------

results_dt <- data.table(
  Orthogroup = og_ids_f,
  effect = obs_effect,
  abs_effect = abs(obs_effect),
  std_effect = std_effect,
  p_exact = p_exact,
  q_bh = q_bh,
  mean_holo = rowMeans(Xf[, holo_idx, drop = FALSE]),
  mean_nonholo = rowMeans(Xf[, non_idx, drop = FALSE]),
  present_species = row_present[keep],
  raw_total = rowSums(Xf_raw)
)

setorder(results_dt, p_exact, -abs_effect)

cat("========== RESULTS SUMMARY ==========\n")
cat("Min exact p-value possible here is limited by only 20 labelings.\n")
cat("Number of orthogroups with p_exact <= 0.05:", sum(results_dt$p_exact <= 0.05), "\n")
cat("Number of orthogroups with q_bh <= 0.10 :", sum(results_dt$q_bh <= 0.10), "\n")
cat("Number of orthogroups with q_bh <= 0.20 :", sum(results_dt$q_bh <= 0.20), "\n")
cat("\nTop 20 hits:\n")
print(results_dt[1:20])
cat("\n")

# Save results
fwrite(results_dt, "orthogroup_holocentricity_screen.tsv", sep = "\t")
cat("Wrote results to: orthogroup_holocentricity_screen.tsv\n\n")

# ------------------------------------------------------------
# 8) VOLCANO-LIKE PLOT
# ------------------------------------------------------------

volcano_df <- as.data.frame(results_dt)
volcano_df$neglog10_p <- -log10(volcano_df$p_exact + 1e-12)

p1 <- ggplot(volcano_df, aes(x = effect, y = neglog10_p)) +
  geom_point(alpha = 0.35, size = 1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Orthogroup association with holocentricity",
    subtitle = "Effect = mean(log1p(copy/ploidy)) in holo - nonholo",
    x = "Effect size",
    y = "-log10(exact p-value)"
  )

print(p1)
cap(p1)
# ------------------------------------------------------------
# 9) HISTOGRAM OF EFFECT SIZES
# ------------------------------------------------------------

p2 <- ggplot(volcano_df, aes(x = effect)) +
  geom_histogram(bins = 100) +
  theme_bw(base_size = 12) +
  labs(
    title = "Distribution of orthogroup effect sizes",
    x = "Effect size",
    y = "Count"
  )

print(p2)
cap(p2)
# ------------------------------------------------------------
# 10) HEATMAP OF TOP HITS
# ------------------------------------------------------------

top_n <- 30
top_hits <- results_dt[order(p_exact, -abs_effect)][1:top_n]$Orthogroup
top_idx <- match(top_hits, og_ids_f)

heat_mat <- Xf[top_idx, , drop = FALSE]
rownames(heat_mat) <- top_hits
colnames(heat_mat) <- species_cols

annotation_col <- data.frame(
  holocentric = factor(meta_df$holocentric, levels = c(0,1)),
  ploidy = factor(meta_df$ploidy)
)
rownames(annotation_col) <- species_cols

cat("Drawing heatmap for top", top_n, "hits...\n")
pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = paste("Top", top_n, "candidate orthogroups \n Colors = log1p(copy number / ploidy)")
)
cap(pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = paste("Top", top_n, "candidate orthogroups \n Colors = log1p(copy number / ploidy)")
))
# ------------------------------------------------------------
# 11) BOXPLOTS FOR TOP 6 HITS
# ------------------------------------------------------------

top_show <- 6
top6 <- results_dt[order(p_exact, -abs_effect)][1:top_show]$Orthogroup
top6_idx <- match(top6, og_ids_f)

plot_long <- data.table(
  Orthogroup = rep(top6, each = ncol(Xf)),
  species = rep(species_cols, times = top_show),
  value = as.vector(t(Xf[top6_idx, , drop = FALSE]))
)

plot_long <- merge(plot_long, meta_df, by = "species", all.x = TRUE)

p3 <- ggplot(plot_long, aes(x = factor(holocentric), y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.08), size = 2, alpha = 0.9) +
  facet_wrap(~ Orthogroup, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Top orthogroups by holocentricity association",
    x = "Holocentricity (0 = nonholo, 1 = holo)",
    y = "log1p(copy number / ploidy)"
  )

print(p3)
cap(p3)
# ------------------------------------------------------------
# 12) VERY SIMPLE HUMAN-READABLE HIT SUMMARY
# ------------------------------------------------------------

cat("========== HUMAN-READABLE TOP HIT SUMMARY ==========\n")
top_summary <- results_dt[1:15, .(
  Orthogroup,
  effect = round(effect, 4),
  mean_holo = round(mean_holo, 4),
  mean_nonholo = round(mean_nonholo, 4),
  p_exact = round(p_exact, 4),
  q_bh = round(q_bh, 4),
  present_species,
  raw_total
)]
print(top_summary)
cat("\n")

# ------------------------------------------------------------
# 13) OPTIONAL: POSITIVE VS NEGATIVE HIT TABLES
# ------------------------------------------------------------

top_positive <- results_dt[order(-effect, p_exact)][1:20]
top_negative <- results_dt[order(effect, p_exact)][1:20]

cat("Top 20 orthogroups enriched in holocentric species:\n")
print(top_positive)
cat("\n")

cat("Top 20 orthogroups depleted in holocentric species:\n")
print(top_negative)
cat("\n")

# ------------------------------------------------------------
# 14) INTERPRETATION NOTES
# ------------------------------------------------------------

cat("========== INTERPRETATION NOTES ==========\n")
cat("* Positive effect: orthogroup tends to have higher ploidy-adjusted copy number in holocentric species.\n")
cat("* Negative effect: orthogroup tends to have lower ploidy-adjusted copy number in holocentric species.\n")
cat("* Because n=6, treat these as candidate orthogroups, not definitive biomarkers.\n")
cat("* FDR may be very weak because the sample size is tiny and there are many orthogroups.\n")
cat("* Strong candidates are those with: low p_exact, large |effect|, and visually clean separation in the boxplots/heatmap.\n")
cat("\n")



library(data.table)
library(matrixStats)

# ------------------------------------------------------------
# JACKKNIFE STABILITY ANALYSIS
# ------------------------------------------------------------

cat("========== JACKKNIFE STABILITY ANALYSIS ==========\n")

species <- meta_df$species
n_species <- length(species)

loo_effect_mat <- matrix(
  NA_real_,
  nrow = nrow(Xf),
  ncol = n_species,
  dimnames = list(og_ids_f, species)
)

for (leave_out in seq_len(n_species)) {

  keep_species <- setdiff(seq_len(n_species), leave_out)

  holo_idx_sub <- which(meta_df$holocentric[keep_species] == 1)
  non_idx_sub  <- which(meta_df$holocentric[keep_species] == 0)

  X_sub <- Xf[, keep_species, drop = FALSE]

  effect_sub <- rowMeans(X_sub[, holo_idx_sub, drop = FALSE]) -
                rowMeans(X_sub[, non_idx_sub, drop = FALSE])

  loo_effect_mat[, leave_out] <- effect_sub
}

loo_abs_dev <- abs(loo_effect_mat - obs_effect)

jackknife_max_influence <- rowMaxs(loo_abs_dev, na.rm = TRUE)
jackknife_median_influence <- rowMedians(loo_abs_dev, na.rm = TRUE)
jackknife_mean_influence <- rowMeans(loo_abs_dev, na.rm = TRUE)

full_sign <- sign(obs_effect)

same_sign_mat <- sweep(
  sign(loo_effect_mat),
  1,
  full_sign,
  FUN = "=="
)

jackknife_sign_stability <- rowMeans(same_sign_mat, na.rm = TRUE)
jackknife_sign_stability[full_sign == 0] <- NA_real_

jackknife_stability <- 1 / (1 + jackknife_median_influence)

cat("Jackknife effects computed.\n\n")

# ------------------------------------------------------------
# CONSISTENCY SCORE
# ------------------------------------------------------------

cat("========== CONSISTENCY ANALYSIS ==========\n")

pairwise_consistency <- numeric(nrow(Xf))

for (i in seq_len(nrow(Xf))) {

  vals <- Xf[i, ]
  holo_vals <- vals[holo_idx]
  non_vals  <- vals[non_idx]

  pair_diff <- outer(holo_vals, non_vals, "-")

  if (obs_effect[i] > 0) {
    pair_scores <- ifelse(pair_diff > 0, 1,
                          ifelse(pair_diff == 0, 0.5, 0))
  } else if (obs_effect[i] < 0) {
    pair_scores <- ifelse(pair_diff < 0, 1,
                          ifelse(pair_diff == 0, 0.5, 0))
  } else {
    pair_scores <- rep(0.5, length(pair_diff))
  }

  pairwise_consistency[i] <- mean(pair_scores)
}

cat("Pairwise consistency scoring done.\n\n")

# ------------------------------------------------------------
# SEPARATION METRICS
# ------------------------------------------------------------

cat("========== SEPARATION ANALYSIS ==========\n")

separation_gap <- numeric(nrow(Xf))

for (i in seq_len(nrow(Xf))) {
  vals <- Xf[i, ]
  holo_vals <- vals[holo_idx]
  non_vals  <- vals[non_idx]

  if (obs_effect[i] > 0) {
    separation_gap[i] <- min(holo_vals) - max(non_vals)
  } else if (obs_effect[i] < 0) {
    separation_gap[i] <- min(non_vals) - max(holo_vals)
  } else {
    separation_gap[i] <- 0
  }
}

# positive gap only; overlap gets 0
gap_score <- pmax(separation_gap, 0)

if (max(gap_score, na.rm = TRUE) > 0) {
  gap_score <- gap_score / max(gap_score, na.rm = TRUE)
} else {
  gap_score <- rep(0, length(gap_score))
}

present_mat <- (Xf_raw > 0) * 1

presence_holo <- rowMeans(present_mat[, holo_idx, drop = FALSE], na.rm = TRUE)
presence_non  <- rowMeans(present_mat[, non_idx, drop = FALSE], na.rm = TRUE)

presence_separation <- ifelse(
  obs_effect > 0,
  presence_holo - presence_non,
  ifelse(obs_effect < 0,
         presence_non - presence_holo,
         0)
)

presence_separation <- pmax(presence_separation, 0)

separation_score <- 0.5 * gap_score + 0.5 * presence_separation

cat("Separation metrics computed.\n\n")

# ------------------------------------------------------------
# BUILD METRIC TABLE IN Xf / og_ids_f ORDER
# ------------------------------------------------------------

metric_dt <- data.table(
  Orthogroup = og_ids_f,
  stability = jackknife_stability,
  consistency = pairwise_consistency,
  jackknife_max_influence = jackknife_max_influence,
  jackknife_median_influence = jackknife_median_influence,
  jackknife_mean_influence = jackknife_mean_influence,
  jackknife_sign_stability = jackknife_sign_stability,
  separation_gap = separation_gap,
  gap_score = gap_score,
  presence_separation = presence_separation,
  separation_score = separation_score
)

# ------------------------------------------------------------
# MERGE METRICS BACK TO RESULTS BY ORTHOGROUP
# ------------------------------------------------------------

cat("========== COMBINING SCORES ==========\n")

final_dt <- merge(
  copy(results_dt),
  metric_dt,
  by = "Orthogroup",
  all.x = TRUE,
  sort = FALSE
)

final_dt[, effect_scaled := abs_effect / max(abs_effect, na.rm = TRUE)]
final_dt[is.na(jackknife_sign_stability), jackknife_sign_stability := 0.5]

# recommended additive score:
# prioritize clean separation first, then consistency, then robustness
final_dt[, final_score :=
           (0.35 * separation_score) +
           (0.30 * consistency) +
           (0.20 * stability) +
           (0.15 * effect_scaled)
]

# combine separation + consistency into unified signal
final_dt[, signal_score := 0.5 * consistency + 0.5 * separation_score]

final_dt[, final_score :=
           (0.45  * signal_score) +
           (0.25 * stability) +
           (0.30 * effect_scaled)
]

# separation is a BONUS, not required
final_dt[, signal_score := consistency + 0.3 * separation_score]

# normalize to 0–1
final_dt[, signal_score := signal_score / max(signal_score, na.rm = TRUE)]

final_dt[, final_score :=
           (0.50 * signal_score) +
           (0.25 * stability) +
           (0.25 * effect_scaled)
]
setorder(final_dt, -final_score)

cat("Final scoring complete.\n\n")

cat("========== FINAL SCORE SUMMARY ==========\n")
print(summary(final_dt$final_score))
cat("\nTop 20 candidates:\n")
print(
  final_dt[1:20, .(
    Orthogroup,
    final_score = round(final_score, 4),
    separation_score = round(separation_score, 4),
    gap_score = round(gap_score, 4),
    presence_separation = round(presence_separation, 4),
    consistency = round(consistency, 4),
    stability = round(stability, 4),
    effect = round(effect, 4),
    abs_effect = round(abs_effect, 4),
    p_exact = round(p_exact, 4),
    q_bh = round(q_bh, 4)
  )]
)
cat("\n")
#------------------------------------------------------------
# OUTPUT TOP CANDIDATES
# ------------------------------------------------------------

cat("========== TOP HIGH-CONFIDENCE CANDIDATES ==========\n")

top_candidates <- final_dt[1:30, .(
  Orthogroup,
  final_score = round(final_score, 3),
  stability = round(stability, 3),
  consistency,
  effect = round(effect, 3),
  mean_holo = round(mean_holo, 3),
  mean_nonholo = round(mean_nonholo, 3),
  present_species,
  raw_total
)]

print(top_candidates)

fwrite(final_dt, "orthogroup_high_confidence.tsv", sep = "\t")

cat("\nSaved full ranked table: orthogroup_high_confidence.tsv\n\n")

# hist of final scores
p_final <- ggplot(final_dt, aes(x = final_score)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(
    title = "Distribution of final composite scores",
    x = "Final score",
    y = "Count"
  )
# redraw on values above 0.5
p_final_zoom <- ggplot(final_dt[final_score > 0.7], aes(x = final_score)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(
    title = "Zoomed distribution of final scores > 0.7",
    x = "Final score",
    y = "Count"
  )
grid.arrange(p_final, p_final_zoom, nrow = 1)
cap(grid.arrange(p_final, p_final_zoom, nrow = 1))
# ------------------------------------------------------------

# VISUAL 1: Stability vs Effect
# ------------------------------------------------------------

library(ggplot2)

p1 <- ggplot(final_dt, aes(x = stability, y = abs_effect)) +
  geom_point(alpha = 0.4, size = 1) +
  theme_bw() +
  labs(
    title = "Stability vs Effect Size",
    x = "Stability (fraction of leave-one-out runs)",
    y = "Absolute effect size"
  )

print(p1)
cap(p1)
# visual 2 3d graph of stability, consistency, effect
library(plotly)
p2 <- plot_ly(
  final_dt,
  x = ~stability,
  y = ~consistency,
  z = ~effect_scaled,
  text = ~Orthogroup,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = ~final_score, colorscale = "Viridis", showscale = TRUE)
) %>%
  layout(
    title = "3D view of stability, consistency, and effect",
    scene = list(
      xaxis = list(title = "Stability"),
      yaxis = list(title = "Consistency"),
      zaxis = list(title = "Scaled Effect")
    )
  )
htmlwidgets::saveWidget(p2, "p2.html", selfcontained = TRUE)
# colored by fiunal score
cat("Number of orthogroups with final_score > 0.8:", sum(final_dt$final_score > 0.8), "\n")

# ------------------------------------------------------------
# VISUAL 3: Heatmap of BEST candidates only
# ------------------------------------------------------------

library(pheatmap)

best_n <- nrow(final_dt[final_score > 0.7])


best_ids <- final_dt$Orthogroup[1:best_n]
best_idx <- match(best_ids, og_ids_f)

heat_mat <- Xf[best_idx, , drop = FALSE]
rownames(heat_mat) <- best_ids

annotation_col <- data.frame(
  holocentric = factor(meta_df$holocentric),
  ploidy = factor(meta_df$ploidy)
)
rownames(annotation_col) <- species

pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Best high-confidence orthogroups"
)

cap(pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Best high-confidence orthogroups"
))


# ok fine now we can bring in the  functional annotation of the orthogroups and see if any of the top candidates have interesting annotations, but that is for another day












library(data.table)

cat("=====================================================\n")
cat("BUILDING MASTER GENE TABLE\n")
cat("=====================================================\n\n")

## =========================================================
## 0) Absolute paths
## =========================================================

fantasia_base <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/pan_drosera/fantasia_output"

species_names <- c(
  "D_aliciae",
  "D_binata",
  "D_capensis",
  "D_paradoxa",
  "D_regia",
  "D_scorpioides"
)

## =========================================================
## 1) Gene -> Orthogroup table
## =========================================================

cat("========== STEP 1: GENE -> ORTHOGROUP ==========\n")

# gene_long already has Orthogroup and gene from your earlier code
gene_og_dt <- unique(gene_long[, .(gene, Orthogroup)])

cat("Rows in gene_og_dt:", nrow(gene_og_dt), "\n")
cat("Unique genes with orthogroup assignment:", uniqueN(gene_og_dt$gene), "\n")
cat("Unique orthogroups:", uniqueN(gene_og_dt$Orthogroup), "\n\n")

## =========================================================
## 2) Load all annotation files with wildcard discovery
## =========================================================

cat("========== STEP 2: LOADING FUNCTIONAL ANNOTATIONS ==========\n")

anno_list <- vector("list", length(species_names))
names(anno_list) <- species_names

for (sp in species_names) {
  sp_dir <- file.path(fantasia_base, sp)

  if (!dir.exists(sp_dir)) {
    stop(sprintf("Species directory not found: %s", sp_dir))
  }

  # detect run dirs dynamically
  run_dirs <- list.dirs(sp_dir, recursive = FALSE, full.names = TRUE)

  if (length(run_dirs) == 0) {
    warning(sprintf("No run directories found for %s", sp))
    next
  }

  # look for results*.csv in all run dirs
  candidate_files <- unlist(
    lapply(run_dirs, function(x) {
      list.files(x, pattern = "^results.*\\.csv$", full.names = TRUE)
    }),
    use.names = FALSE
  )

  if (length(candidate_files) == 0) {
    warning(sprintf("No results CSV found for %s", sp))
    next
  }

  # choose the newest by modification time
  finfo <- file.info(candidate_files)
  file_path <- rownames(finfo)[which.max(finfo$mtime)]

  cat("Reading", sp, "from:\n  ", file_path, "\n")

  dt <- fread(file_path)

  dt[, species_from_file := sp]
  dt[, source_file := file_path]

  anno_list[[sp]] <- dt
}

anno_dt <- rbindlist(anno_list, fill = TRUE, use.names = TRUE)

cat("\nAnnotation rows loaded:", nrow(anno_dt), "\n")
cat("Unique annotated genes:", uniqueN(anno_dt$query_accession), "\n")
cat("Unique GO IDs:", uniqueN(anno_dt$go_id), "\n\n")

if (nrow(anno_dt) == 0) {
  stop("No annotation data loaded.")
}

## =========================================================
## 3) Clean annotation table
## =========================================================

cat("========== STEP 3: CLEANING ANNOTATIONS ==========\n")

# Standardize gene column name
setnames(anno_dt, old = "query_accession", new = "gene", skip_absent = TRUE)

# Keep only columns likely to be useful
keep_cols <- c(
  "gene",
  "species_from_file",
  "reference_id",
  "model_key",
  "distance",
  "reliability_index",
  "distance_metric",
  "uniprot_accession",
  "go_id",
  "category",
  "go_description",
  "evidence_codes",
  "source_file"
)

keep_cols <- intersect(keep_cols, names(anno_dt))
anno_dt <- anno_dt[, ..keep_cols]

# Remove exact duplicate rows if any
anno_dt <- unique(anno_dt)

cat("Rows after deduplication:", nrow(anno_dt), "\n")
cat("Unique genes after deduplication:", uniqueN(anno_dt$gene), "\n\n")

cat("Reliability summary:\n")
print(summary(anno_dt$reliability_index))
cat("\n")

## =========================================================
## 4) Finalize annotation table
## =========================================================

anno_dt[, species := species_from_file]
anno_dt[, species_from_file := NULL]
anno_dt[, species_gene := paste(species, gene, sep = "::")]

setcolorder(
  anno_dt,
  c("species", "gene", "species_gene",
    setdiff(names(anno_dt), c("species", "gene", "species_gene")))
)

stopifnot(!anyNA(anno_dt$species))
stopifnot(!anyNA(anno_dt$gene))

## =========================================================
## 5) Basic checks: did import work?
## =========================================================

cat("========== STEP 4: ANNOTATION IMPORT CHECKS ==========\n")

cat("Rows in anno_dt:", nrow(anno_dt), "\n")
cat("Unique species:", uniqueN(anno_dt$species), "\n")
cat("Unique genes:", uniqueN(anno_dt$species_gene), "\n")
cat("Unique GO IDs:", uniqueN(anno_dt$go_id), "\n\n")

cat("Rows per species:\n")
print(anno_dt[, .N, by = species][order(species)])
cat("\n")

cat("Unique annotated genes per species:\n")
print(anno_dt[, .(annotated_genes = uniqueN(gene)), by = species][order(species)])
cat("\n")

cat("GO category counts:\n")
print(anno_dt[, .N, by = category][order(-N)])
cat("\n")

cat("Reliability summary overall:\n")
print(summary(anno_dt$reliability_index))
cat("\n")

cat("Reliability summary by species:\n")
print(
  anno_dt[, .(
    n_rows = .N,
    n_genes = uniqueN(gene),
    median_rel = median(reliability_index, na.rm = TRUE),
    mean_rel = mean(reliability_index, na.rm = TRUE)
  ), by = species][order(species)]
)
cat("\n")

cat("Missingness:\n")
print(anno_dt[, .(
  missing_go_id = sum(is.na(go_id) | go_id == ""),
  missing_go_desc = sum(is.na(go_description) | go_description == ""),
  missing_reliability = sum(is.na(reliability_index)),
  missing_uniprot = sum(is.na(uniprot_accession) | uniprot_accession == "")
)])
cat("\n")

## optional: high-confidence subset
anno_hc_dt <- anno_dt[!is.na(reliability_index) & reliability_index >= 0.90]

cat("High-confidence rows (reliability >= 0.90):", nrow(anno_hc_dt), "\n")
cat("High-confidence unique genes:", uniqueN(anno_hc_dt$species_gene), "\n\n")

## =========================================================
## 6) Collapse to one row per gene
## =========================================================

collapse_unique <- function(x) {
  x <- sort(unique(na.omit(x)))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(x, collapse = "; ")
}

gene_func_dt <- anno_dt[
  ,
  .(
    n_rows = .N,
    n_go_terms = uniqueN(go_id[!is.na(go_id) & go_id != ""]),
    max_reliability = suppressWarnings(max(reliability_index, na.rm = TRUE)),
    mean_reliability = suppressWarnings(mean(reliability_index, na.rm = TRUE)),
    uniprot_accessions = collapse_unique(uniprot_accession),
    go_ids = collapse_unique(go_id),
    go_descriptions = collapse_unique(go_description),
    go_categories = collapse_unique(category),
    evidence_codes = collapse_unique(evidence_codes),
    model_keys = collapse_unique(model_key)
  ),
  by = .(species, gene, species_gene)
]

gene_func_dt[!is.finite(max_reliability), max_reliability := NA_real_]
gene_func_dt[!is.finite(mean_reliability), mean_reliability := NA_real_]

cat("========== STEP 5: COLLAPSED GENE-LEVEL TABLE ==========\n")
cat("Rows in gene_func_dt:", nrow(gene_func_dt), "\n")
cat("Unique genes in gene_func_dt:", uniqueN(gene_func_dt$species_gene), "\n\n")

cat("Per-species collapsed gene counts:\n")
print(gene_func_dt[, .N, by = species][order(species)])
cat("\n")

cat("Preview anno_dt:\n")
print(anno_dt[1:10])
cat("\n")

cat("Preview gene_func_dt:\n")
print(gene_func_dt[1:10])
cat("\n")
fwrite(anno_dt, file.path(wd, "annotation_long.tsv"), sep = "\t")
fwrite(gene_func_dt, file.path(wd, "gene_function_collapsed.tsv"), sep = "\t")


gene_og_dt <- unique(gene_long[, .(species, gene, species_gene, Orthogroup)])
gene_master_dt <- merge(gene_og_dt, gene_func_dt, by = c("species", "gene", "species_gene"), all = TRUE)

## merge gene → OG with annotations
og_anno_dt <- merge(
  gene_long[, .(Orthogroup, species, species_gene)],
  anno_dt[, .(species, species_gene, go_id)],
  by = c("species", "species_gene"),
  allow.cartesian = TRUE
)

## remove missing GO
og_anno_dt <- og_anno_dt[!is.na(go_id) & go_id != ""]
## number of annotated genes per OG
og_gene_counts <- unique(og_anno_dt[, .(Orthogroup, species_gene)])[
  , .(n_genes_annotated = .N), by = Orthogroup
]

## GO counts per OG
og_go_counts <- og_anno_dt[
  , .(
    n_genes_with_go = uniqueN(species_gene),
    n_species_with_go = uniqueN(species)
  ),
  by = .(Orthogroup, go_id)
]
## merge + compute support
og_go_dt <- merge(og_go_counts, og_gene_counts, by = "Orthogroup")

og_go_dt[, gene_support := n_genes_with_go / n_genes_annotated]

## species support (normalize by species present in OG)
og_species_counts <- unique(gene_long[, .(Orthogroup, species)])[
  , .(n_species = .N), by = Orthogroup
]

og_go_dt <- merge(og_go_dt, og_species_counts, by = "Orthogroup")

og_go_dt[, species_support := n_species_with_go / n_species]

library(ggplot2)

a <- ggplot(og_go_dt, aes(x = gene_support)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(limits = c(0,1)) +
  theme_minimal() +
  ggtitle("Distribution of gene-level GO support")
# save it
ggsave(filename = file.path(wd, "gene_support_histogram.png"), plot = a, width = 6, height = 4, dpi = 300)
cap(a)
a <- ggplot(og_go_dt, aes(x = n_genes_with_go)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_minimal() +
  ggtitle("Number of genes supporting each GO term (log scale)")
ggsave(filename = file.path(wd, "n_genes_with_go_histogram.png"), plot = a, width = 6, height = 4, dpi = 300)
cap(a)
a<- ggplot(og_go_dt, aes(x = species_support)) +
  geom_histogram(bins = 20) +
  theme_minimal() +
  ggtitle("Distribution of species support")
ggsave(filename = file.path(wd, "species_support_histogram.png"), plot = a, width = 6, height = 4, dpi = 300)
cap(a)
a <- ggplot(og_go_dt, aes(x = gene_support, y = species_support)) +
  geom_bin2d() +
  theme_minimal() +
  ggtitle("Gene vs Species support")
ggsave(filename = file.path(wd, "gene_vs_species_support.png"), plot = a, width = 6, height = 4, dpi = 300)
cap(a)
## =========================================================
## Neutral orthogroup function assignment
## Rule:
##   n_genes_with_go >= 2
##   gene_support    >= 0.4
##   species_support >= 0.5
## =========================================================

cat("========== BUILDING NEUTRAL ORTHOGROUP FUNCTION ASSIGNMENTS ==========\n")

## keep only GO-annotated rows we need
anno_go_dt <- unique(
  anno_dt[!is.na(go_id) & go_id != "",
          .(species, species_gene, go_id, go_description, category)]
)

## merge orthogroup membership with GO annotations
og_anno_dt <- merge(
  gene_long[, .(Orthogroup, species, species_gene)],
  anno_go_dt,
  by = c("species", "species_gene"),
  allow.cartesian = TRUE
)

cat("Rows in og_anno_dt:", nrow(og_anno_dt), "\n")
cat("Unique orthogroups with any GO annotation:", uniqueN(og_anno_dt$Orthogroup), "\n")
cat("Unique GO terms:", uniqueN(og_anno_dt$go_id), "\n\n")

## denominator 1: total annotated genes per orthogroup
og_annot_gene_dt <- unique(
  og_anno_dt[, .(Orthogroup, species_gene)]
)[
  , .(n_annotated_genes = .N), by = Orthogroup
]

## denominator 2: total species represented in the orthogroup
og_species_dt <- unique(
  gene_long[, .(Orthogroup, species)]
)[
  , .(n_species_in_og = .N), by = Orthogroup
]

## numerator tables: support for each GO term within each orthogroup
og_go_dt <- og_anno_dt[
  ,
  .(
    n_genes_with_go = uniqueN(species_gene),
    n_species_with_go = uniqueN(species),
    go_description = go_description[1],
    category = category[1]
  ),
  by = .(Orthogroup, go_id)
]

## add support denominators
og_go_dt <- merge(og_go_dt, og_annot_gene_dt, by = "Orthogroup")
og_go_dt <- merge(og_go_dt, og_species_dt, by = "Orthogroup")

## support metrics
og_go_dt[, gene_support := n_genes_with_go / n_annotated_genes]
og_go_dt[, species_support := n_species_with_go / n_species_in_og]

## strict neutral assignment
og_go_neutral_dt <- og_go_dt[
  n_genes_with_go >= 2 &
  gene_support >= 0.4 &
  species_support >= 0.5
][
  order(Orthogroup, -species_support, -gene_support, -n_genes_with_go, go_id)
]

cat("Orthogroup x GO rows total:", nrow(og_go_dt), "\n")
cat("Rows passing neutral thresholds:", nrow(og_go_neutral_dt), "\n")
cat("Orthogroups with >=1 assigned GO:", uniqueN(og_go_neutral_dt$Orthogroup), "\n\n")

cat("Threshold summary:\n")
cat("  n_genes_with_go >= 2\n")
cat("  gene_support    >= 0.4\n")
cat("  species_support >= 0.5\n\n")

## top GO term per orthogroup
og_top_func_dt <- og_go_neutral_dt[, .SD[1], by = Orthogroup]

cat("Orthogroups with a top assigned GO:", nrow(og_top_func_dt), "\n\n")

## collapsed summary: top few GO terms per orthogroup
collapse_unique <- function(x) {
  x <- sort(unique(na.omit(x)))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(x, collapse = "; ")
}

og_func_summary_dt <- og_go_neutral_dt[
  ,
  .(
    n_assigned_go = .N,
    top_go_id = go_id[1],
    top_go_description = go_description[1],
    top_category = category[1],
    assigned_go_ids = collapse_unique(go_id),
    assigned_go_descriptions = collapse_unique(go_description),
    assigned_categories = collapse_unique(category),
    max_gene_support = max(gene_support),
    max_species_support = max(species_support)
  ),
  by = Orthogroup
]

cat("Rows in og_func_summary_dt:", nrow(og_func_summary_dt), "\n\n")

## =========================================================
## Diagnostics
## =========================================================

cat("========== DIAGNOSTICS ==========\n")

cat("Assigned GO count per orthogroup:\n")
print(summary(og_func_summary_dt$n_assigned_go))
cat("\n")

cat("Gene support among retained OG-GO assignments:\n")
print(summary(og_go_neutral_dt$gene_support))
cat("\n")

cat("Species support among retained OG-GO assignments:\n")
print(summary(og_go_neutral_dt$species_support))
cat("\n")

cat("Top GO categories among retained assignments:\n")
print(og_go_neutral_dt[, .N, by = category][order(-N)])
cat("\n")

cat("Preview of top assigned functions:\n")
print(og_top_func_dt[1:10])
cat("\n")

## =========================================================
## Save
## =========================================================

fwrite(og_go_dt, file.path(wd, "og_go_support_table.tsv"), sep = "\t")
fwrite(og_go_neutral_dt, file.path(wd, "og_go_neutral_assignments.tsv"), sep = "\t")
fwrite(og_top_func_dt, file.path(wd, "og_top_function.tsv"), sep = "\t")
fwrite(og_func_summary_dt, file.path(wd, "og_function_summary.tsv"), sep = "\t")

cat("Saved:\n")
cat("  ", file.path(wd, "og_go_support_table.tsv"), "\n")
cat("  ", file.path(wd, "og_go_neutral_assignments.tsv"), "\n")
cat("  ", file.path(wd, "og_top_function.tsv"), "\n")
cat("  ", file.path(wd, "og_function_summary.tsv"), "\n\n")



final_annotated_dt <- merge(
  final_dt,
  og_func_summary_dt,
  by = "Orthogroup",
  all.x = TRUE
)



library(data.table)
library(ggplot2)

run_go_enrichment <- function(final_dt,
                              og_go_neutral_dt,
                              top_n = 200,
                              min_top_hits = 3,
                              max_terms_plot = 20) {
  stopifnot(all(c("Orthogroup", "final_score") %in% names(final_dt)))
  stopifnot(all(c("Orthogroup", "go_id", "go_description", "category") %in% names(og_go_neutral_dt)))

  # rank orthogroups by final score
  ranked_dt <- copy(final_dt)[order(-final_score)]

  # use only orthogroups that actually have at least one assigned GO term
  background_ogs <- unique(og_go_neutral_dt$Orthogroup)
  ranked_dt <- ranked_dt[Orthogroup %in% background_ogs]

  if (nrow(ranked_dt) < top_n) {
    warning(sprintf("Only %d annotated orthogroups available; using all of them.", nrow(ranked_dt)))
    top_n <- nrow(ranked_dt)
  }

  top_ogs <- ranked_dt[1:top_n, unique(Orthogroup)]
  bg_ogs  <- unique(ranked_dt$Orthogroup)

  # binary OG x GO presence table
  og_go_bin <- unique(
    og_go_neutral_dt[, .(Orthogroup, go_id, go_description, category)]
  )

  # enrichment per GO term
  enrich_dt <- og_go_bin[
    ,
    {
      ogs_with_go <- unique(Orthogroup)

      a <- sum(ogs_with_go %in% top_ogs)                  # top OGs with GO
      b <- sum(ogs_with_go %in% bg_ogs) - a              # other OGs with GO
      c <- length(top_ogs) - a                           # top OGs without GO
      d <- length(bg_ogs) - length(top_ogs) - b          # other OGs without GO

      # avoid invalid tests
      if ((a + b) == 0 || (c + d) == 0) {
        p <- NA_real_
        or <- NA_real_
      } else {
        ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
        p <- ft$p.value
        or <- unname(ft$estimate)
      }

      .(
        n_top_with_go = a,
        n_bg_with_go = a + b,
        p_value = p,
        odds_ratio = or
      )
    },
    by = .(go_id, go_description, category)
  ]

  enrich_dt[, p_adj := p.adjust(p_value, method = "BH")]
  enrich_dt[, neglog10_padj := -log10(p_adj)]
  enrich_dt[, fold_in_top := n_top_with_go / length(top_ogs)]
  enrich_dt[, fold_in_bg := n_bg_with_go / length(bg_ogs)]

  setorder(enrich_dt, p_adj, -odds_ratio, -n_top_with_go)

  # plotting subset
  plot_dt <- enrich_dt[
    is.finite(neglog10_padj) &
      !is.na(odds_ratio) &
      n_top_with_go >= min_top_hits
  ][1:min(.N, max_terms_plot)]

  p <- ggplot(plot_dt,
              aes(x = neglog10_padj,
                  y = reorder(go_description, neglog10_padj))) +
    geom_point(aes(size = n_top_with_go, color = category)) +
    theme_bw(base_size = 12) +
    labs(
      title = sprintf("GO enrichment in top %d orthogroups", top_n),
      subtitle = "Ranked by final_score; background = all GO-annotated orthogroups",
      x = "-log10(BH-adjusted p-value)",
      y = NULL,
      size = "Top OGs with GO",
      color = "GO category"
    )

  list(
    top_ogs = top_ogs,
    background_ogs = bg_ogs,
    enrichment_table = enrich_dt,
    plot_data = plot_dt,
    plot = p
  )
}



res500 <- run_go_enrichment(
  final_dt = final_dt,
  og_go_neutral_dt = og_go_neutral_dt,
  top_n = 500,
  min_top_hits = 3,
  max_terms_plot = 20
)

print(res500$plot)
cap(res500$plot)
head(res500$enrichment_table)


library(data.table)
library(ggplot2)

run_topN_go_enrichment <- function(final_dt,
                                   og_go_dt,
                                   top_n = 200,
                                   min_genes_with_go = 2,
                                   min_bg_with_go = 10,
                                   max_terms_plot = 30) {
  stopifnot(all(c("Orthogroup", "final_score") %in% names(final_dt)))
  stopifnot(all(c("Orthogroup", "go_id", "go_description", "category", "n_genes_with_go") %in% names(og_go_dt)))

  ## rank orthogroups
  rank_dt <- unique(final_dt[, .(Orthogroup, final_score)])
  rank_dt <- rank_dt[!is.na(final_score)]
  setorder(rank_dt, -final_score)

  ## GO presence table for enrichment
  og_go_bin <- unique(
    og_go_dt[
      n_genes_with_go >= min_genes_with_go &
        !is.na(go_id) & go_id != "",
      .(Orthogroup, go_id, go_description, category)
    ]
  )

  ## background = all orthogroups with at least one GO in this table
  background_ogs <- unique(og_go_bin$Orthogroup)
  rank_dt <- rank_dt[Orthogroup %in% background_ogs]

  if (nrow(rank_dt) == 0) stop("No ranked orthogroups overlap GO-annotated orthogroups.")

  top_n <- min(top_n, nrow(rank_dt))
  top_ogs <- rank_dt[1:top_n, Orthogroup]
  bg_ogs <- rank_dt$Orthogroup

  ## Fisher test per GO term
  enrich_dt <- og_go_bin[
    ,
    {
      ogs_with_go <- unique(Orthogroup)

      a <- sum(ogs_with_go %in% top_ogs)                # top OGs with GO
      b <- sum(ogs_with_go %in% bg_ogs) - a            # non-top OGs with GO
      c <- length(top_ogs) - a                         # top OGs without GO
      d <- length(bg_ogs) - length(top_ogs) - b        # non-top OGs without GO

      if (a + b == 0 || c + d == 0) {
        p <- NA_real_
        or <- NA_real_
      } else {
        ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
        p <- ft$p.value
        or <- unname(ft$estimate)
      }

      .(
        n_top_with_go = a,
        n_bg_with_go = a + b,
        odds_ratio = or,
        p_value = p
      )
    },
    by = .(go_id, go_description, category)
  ]

  enrich_dt <- enrich_dt[n_bg_with_go >= min_bg_with_go]
  enrich_dt[, p_adj := p.adjust(p_value, method = "BH")]
  enrich_dt[, neglog10_padj := -log10(p_adj)]
  enrich_dt[, top_fraction := n_top_with_go / length(top_ogs)]
  enrich_dt[, bg_fraction := n_bg_with_go / length(bg_ogs)]

  setorder(enrich_dt, p_adj, -odds_ratio, -n_top_with_go)

  plot_dt <- enrich_dt[
    is.finite(neglog10_padj) &
      !is.na(odds_ratio)
  ][1:min(.N, max_terms_plot)]

  p <- ggplot(
    plot_dt,
    aes(x = neglog10_padj,
        y = reorder(go_description, neglog10_padj),
        size = n_top_with_go,
        color = category)
  ) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(
      title = sprintf("GO enrichment in top %d orthogroups", top_n),
      subtitle = "Ranked by final_score",
      x = "-log10(BH-adjusted p-value)",
      y = NULL,
      size = "Top OGs with GO",
      color = "GO category"
    )

  list(
    top_ogs = top_ogs,
    background_ogs = bg_ogs,
    enrichment_table = enrich_dt,
    plot_data = plot_dt,
    plot = p
  )
}

res200 <- run_topN_go_enrichment(
  final_dt = final_dt,
  og_go_dt = og_go_dt,
  top_n = 200,
  min_genes_with_go = 2,
  min_bg_with_go = 10,
  max_terms_plot = 30
)

print(res200$plot)
cap(res200$plot)
head(res200$enrichment_table, 20)

res500 <- run_topN_go_enrichment(final_dt, og_go_dt, top_n = 500)
print(res500$plot)
cap(res500$plot)
res1000 <- run_topN_go_enrichment(final_dt, og_go_dt, top_n = 1000)
print(res1000$plot)
cap(res1000$plot)

number <- nrow(final_dt[final_score > 0.7])
res_top <- run_topN_go_enrichment(final_dt, og_go_dt, top_n = number)
print(res_top$plot)
cap(res_top$plot)

## install once if needed
## BiocManager::install(c("rrvgo", "GOSemSim", "org.At.tair.db"))

library(data.table)
library(rrvgo)
library(ggplot2)

reduce_go_semantic <- function(enrich_dt,
                               ontology = "P",
                               padj_cutoff = 0.25,
                               max_terms = 100,
                               sim_threshold = 0.7,
                               orgdb = "org.At.tair.db") {
  dt <- copy(enrich_dt)[
    category == ontology &
      !is.na(go_id) & go_id != "" &
      !is.na(p_adj) &
      p_adj <= padj_cutoff
  ]

  if (nrow(dt) == 0) {
    stop(sprintf("No GO terms passed filters for ontology %s", ontology))
  }

  dt <- dt[order(p_adj, -odds_ratio)][1:min(.N, max_terms)]

  ## scores for representative-term selection
  scores <- dt$neglog10_padj
  names(scores) <- dt$go_id

  ## semantic similarity matrix
  simMatrix <- rrvgo::calculateSimMatrix(
    x = dt$go_id,
    orgdb = orgdb,
    ont = switch(ontology, P = "BP", F = "MF", C = "CC"),
    method = "Wang"
  )

  ## reduce terms
  reduced <- rrvgo::reduceSimMatrix(
    simMatrix,
    scores = scores,
    threshold = sim_threshold,
    orgdb = orgdb
  )

  reduced_dt <- as.data.table(reduced)
  setnames(reduced_dt, c("go", "parent", "score", "size", "term", "parentTerm"),
                      c("go_id", "cluster_rep_go_id", "semantic_score", "cluster_size",
                        "go_description", "cluster_rep_description"))

  reduced_dt <- merge(
    reduced_dt,
    dt[, .(go_id, p_adj, p_value, neglog10_padj, odds_ratio, n_top_with_go, n_bg_with_go)],
    by = "go_id",
    all.x = TRUE
  )

  ## representative terms only
  reps_dt <- reduced_dt[go_id == cluster_rep_go_id][order(p_adj, -odds_ratio)]

  p <- ggplot(
    reps_dt,
    aes(x = neglog10_padj,
        y = reorder(cluster_rep_description, neglog10_padj),
        size = n_top_with_go)
  ) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(
      title = sprintf("Semantic-reduced GO enrichment (%s)", ontology),
      subtitle = sprintf("Similarity threshold = %.2f", sim_threshold),
      x = "-log10(BH-adjusted p-value)",
      y = NULL,
      size = "Top OGs with GO"
    )

  list(
    input_terms = dt,
    reduced_table = reduced_dt,
    representative_terms = reps_dt,
    plot = p
  )
}

number <- nrow(final_dt[final_score > 0.7])
res_top <- run_topN_go_enrichment(final_dt, og_go_dt, top_n = number,max_terms_plot = 100)

res_good <- run_topN_go_enrichment(
  final_dt = final_dt,
  og_go_dt = og_go_dt,
  top_n = number,
  min_genes_with_go = 2,
  min_bg_with_go = 10,
  max_terms_plot = 50)
cap(res_top$plot)
cap(res_good$plot)

bp_red <- reduce_go_semantic(
  enrich_dt = res_top$enrichment_table,
  ontology = "P",
  padj_cutoff = 1,
  max_terms = 100,
  sim_threshold = 0.7,
  orgdb = "org.At.tair.db"
)
print(bp_red$plot)
cap(bp_red$plot)
bp_red$representative_terms


mf_red <- reduce_go_semantic(
  enrich_dt = res_top$enrichment_table,
  ontology = "F",
  padj_cutoff = 1,
  max_terms = 100,
  sim_threshold = 0.7,
  orgdb = "org.At.tair.db"
)
print(mf_red$plot)
cap(mf_red$plot)

cl_red <- reduce_go_semantic(
  enrich_dt = res_top$enrichment_table,
  ontology = "C",
  padj_cutoff = 1,
  max_terms = 100,
  sim_threshold = 0.7,
  orgdb = "org.At.tair.db"
)
print(cl_red$plot)
cap(cl_red$plot)
library(data.table)

target_term <- "meiotic cell cycle"

## 1) representative semantic cluster
target_rep_dt <- bp_red$representative_terms[
  cluster_rep_description == target_term | go_description == target_term
]

if (nrow(target_rep_dt) == 0) {
  stop(sprintf("Could not find term: %s", target_term))
}

## 2) all GO IDs in that semantic cluster
target_cluster_ids <- unique(target_rep_dt$cluster_rep_go_id)

target_go_dt <- bp_red$reduced_table[
  cluster_rep_go_id %in% target_cluster_ids
]

cat("GO terms in selected semantic cluster:\n")
print(unique(target_go_dt[, .(go_id, go_description, cluster_rep_go_id, cluster_rep_description)]))

target_go_ids <- unique(target_go_dt$go_id)

## 3) all OGs carrying those GO IDs
target_og_dt <- unique(
  og_go_dt[
    go_id %in% target_go_ids & n_genes_with_go >= 2,
    .(
      Orthogroup, go_id, go_description, category,
      n_genes_with_go, n_annotated_genes, gene_support,
      n_species_with_go, n_species_in_og, species_support
    )
  ]
)

cat("Number of orthogroups linked to target term:", uniqueN(target_og_dt$Orthogroup), "\n")


# find intersection of target OGs with top-N set used for enrichment
target_og_top_dt <- target_og_dt[Orthogroup %in% res_top$top_ogs]
wide_dt[Orthogroup %in% target_og_top_dt$Orthogroup]
gene_long[Orthogroup %in% target_og_top_dt$Orthogroup]


# save this
fwrite(target_og_dt, file.path(wd, "meiotic_cell_cycle_og_go_rows.tsv"), sep = "\t")





## 4) subset final annotated and wide tables
target_final_dt <- merge(
  final_annotated_dt,
  unique(target_og_dt[, .(Orthogroup)]),
  by = "Orthogroup"
)

target_wide_dt <- merge(
  wide_dt,
  unique(target_og_dt[, .(Orthogroup)]),
  by = "Orthogroup"
)

## 5) attach GO rows
target_final_go_dt <- merge(
  target_final_dt,
  target_og_dt,
  by = "Orthogroup",
  allow.cartesian = TRUE
)

## 6) sort
if ("abs_effect" %in% names(target_final_go_dt)) {
  setorder(target_final_go_dt, -final_score, -abs_effect)
} else {
  target_final_go_dt[, abs_effect_tmp := abs(effect)]
  setorder(target_final_go_dt, -final_score, -abs_effect_tmp)
}
setorder(target_wide_dt, Orthogroup)

cat("\nPreview: final_annotated rows\n")
print(target_final_go_dt[1:20])

cat("\nPreview: wide rows\n")
print(target_wide_dt[1:20])

## 7) save
fwrite(target_final_go_dt, file.path(wd, "meiotic_cell_cycle_final_annotated.tsv"), sep = "\t")
fwrite(target_wide_dt, file.path(wd, "meiotic_cell_cycle_wide.tsv"), sep = "\t")
fwrite(target_og_dt, file.path(wd, "meiotic_cell_cycle_og_go_rows.tsv"), sep = "\t")



library(data.table)

target_term <- "meiotic cell cycle"

## 1) identify semantic cluster / GO IDs
target_rep_dt <- bp_red$representative_terms[
  cluster_rep_description == target_term | go_description == target_term
]

if (nrow(target_rep_dt) == 0) {
  stop(sprintf("Could not find term: %s", target_term))
}

target_cluster_ids <- unique(target_rep_dt$cluster_rep_go_id)

target_go_dt <- bp_red$reduced_table[
  cluster_rep_go_id %in% target_cluster_ids
]

target_go_ids <- unique(target_go_dt$go_id)

cat("GO terms in selected semantic cluster:\n")
print(unique(target_go_dt[, .(go_id, go_description, cluster_rep_go_id, cluster_rep_description)]))

## 2) OGs carrying those GO IDs (same rule as enrichment input)
target_og_all_dt <- unique(
  og_go_dt[
    go_id %in% target_go_ids & n_genes_with_go >= 2,
    .(
      Orthogroup, go_id, go_description, category,
      n_genes_with_go, n_annotated_genes, gene_support,
      n_species_with_go, n_species_in_og, species_support
    )
  ]
)

cat("All orthogroups linked to target term:", uniqueN(target_og_all_dt$Orthogroup), "\n")

## 3) restrict to the actual top-N enrichment set
## IMPORTANT: res200$top_ogs must come from the same enrichment run that gave the meiotic term
target_og_top_dt <- target_og_all_dt[Orthogroup %in% res200$top_ogs]

cat("Top-N orthogroups linked to target term:", uniqueN(target_og_top_dt$Orthogroup), "\n")

## sanity check: should match n_top_with_go from enrichment table
meiotic_enrich_row <- res200$enrichment_table[
  go_id %in% target_go_ids
][order(p_adj)]

cat("Relevant enrichment rows:\n")
print(meiotic_enrich_row[, .(go_id, go_description, n_top_with_go, n_bg_with_go, odds_ratio, p_value, p_adj)])

## 4) define one shared OG set
shared_ogs <- sort(unique(target_og_top_dt$Orthogroup))

cat("Shared OG set size used for ALL output tables:", length(shared_ogs), "\n")
print(shared_ogs)

## 5) subset all tables using exactly the same OG set
target_final_dt <- final_annotated_dt[Orthogroup %in% shared_ogs]
target_wide_dt  <- wide_dt[Orthogroup %in% shared_ogs]
target_go_rows_dt <- target_og_top_dt[Orthogroup %in% shared_ogs]

## 6) attach GO rows to final table
target_final_go_dt <- merge(
  target_final_dt,
  target_go_rows_dt,
  by = "Orthogroup",
  allow.cartesian = TRUE
)

## 7) force same ordering everywhere
og_order_dt <- data.table(Orthogroup = shared_ogs, og_order = seq_along(shared_ogs))

target_final_go_dt <- merge(target_final_go_dt, og_order_dt, by = "Orthogroup")
target_wide_dt     <- merge(target_wide_dt, og_order_dt, by = "Orthogroup")
target_go_rows_dt  <- merge(target_go_rows_dt, og_order_dt, by = "Orthogroup")

setorder(target_final_go_dt, og_order, -final_score)
setorder(target_wide_dt, og_order)
setorder(target_go_rows_dt, og_order)

## 8) print matching Orthogroup vectors to verify
cat("\nOrthogroups in final table:\n")
print(unique(target_final_go_dt$Orthogroup))

cat("\nOrthogroups in wide table:\n")
print(unique(target_wide_dt$Orthogroup))

cat("\nOrthogroups in GO rows table:\n")
print(unique(target_go_rows_dt$Orthogroup))

## 9) preview
cat("\nPreview: final rows\n")
print(target_final_go_dt[, .(
  Orthogroup, final_score, effect, abs_effect, stability, consistency,
  go_id, go_description, n_genes_with_go, gene_support, species_support
)][1:20])

cat("\nPreview: wide rows\n")
print(target_wide_dt[1:20])

cat("\nPreview: GO support rows\n")
print(target_go_rows_dt[1:20])

## 10) save
fwrite(target_final_go_dt, file.path(wd, "meiotic_cell_cycle_TOPN_final_annotated.tsv"), sep = "\t")
fwrite(target_wide_dt, file.path(wd, "meiotic_cell_cycle_TOPN_wide.tsv"), sep = "\t")
fwrite(target_go_rows_dt, file.path(wd, "meiotic_cell_cycle_TOPN_og_go_rows.tsv"), sep = "\t")

cat("\nSaved:\n")
cat(file.path(wd, "meiotic_cell_cycle_TOPN_final_annotated.tsv"), "\n")
cat(file.path(wd, "meiotic_cell_cycle_TOPN_wide.tsv"), "\n")
cat(file.path(wd, "meiotic_cell_cycle_TOPN_og_go_rows.tsv"), "\n")



library(data.table)
library(ggplot2)

plot_go_vs_effect <- function(final_dt,
                             og_go_dt,
                             enrich_dt,
                             top_n_terms = 10,
                             min_genes_with_go = 2) {

  ## 1) take top GO terms from enrichment
  top_terms_dt <- enrich_dt[
    !is.na(p_adj)
  ][order(p_adj)][1:top_n_terms]

  top_go_ids <- unique(top_terms_dt$go_id)

  ## 2) map OGs to those GO terms
  og_go_sub <- unique(
    og_go_dt[
      go_id %in% top_go_ids &
        n_genes_with_go >= min_genes_with_go,
      .(Orthogroup, go_id, go_description, category)
    ]
  )

  ## 3) attach effect scores
  dt <- merge(
    og_go_sub,
    final_dt[, .(Orthogroup, effect, final_score)],
    by = "Orthogroup",
    all.x = TRUE
  )

  ## 4) summarise per GO term
  summary_dt <- dt[
    ,
    .(
      mean_effect = mean(effect, na.rm = TRUE),
      median_effect = median(effect, na.rm = TRUE),
      weighted_effect = weighted.mean(effect, w = final_score, na.rm = TRUE),
      n_ogs = uniqueN(Orthogroup)
    ),
    by = .(go_id, go_description, category)
  ]

  ## 5) order terms by weighted effect
  setorder(summary_dt, -weighted_effect)

  ## 6) plot
  p <- ggplot(
    summary_dt,
    aes(
      x = weighted_effect,
      y = reorder(go_description, weighted_effect),
      fill = category
    )
  ) +
    geom_col() +
    theme_bw(base_size = 12) +
    labs(
      title = sprintf("Top %d GO terms vs effect size", top_n_terms),
      subtitle = "Orthogroups weighted by final_score",
      x = "Weighted mean effect",
      y = NULL,
      fill = "GO category"
    )

  list(
    summary = summary_dt,
    plot = p
  )
}
res_plot <- plot_go_vs_effect(
  final_dt = final_dt,
  og_go_dt = og_go_dt,
  enrich_dt = res_top$enrichment_table,
  top_n_terms = 30
)

print(res_plot$plot)
cap(res_plot$plot)
res_plot$summary



top25_ogs <- final_dt[
  order(-final_score)
][1:28, Orthogroup]


og_label_dt <- og_go_dt[
  n_genes_with_go >= 2
][
  , .SD[which.max(gene_support * species_support)],
  by = Orthogroup
][
  , .(
    Orthogroup,
    label = go_description,
    category
  )
]

og_label_dt <- merge(
  data.table(Orthogroup = top25_ogs),
  og_label_dt,
  by = "Orthogroup",
  all.x = TRUE
)

og_label_dt <- merge(
  og_label_dt,
  final_annotated_dt[, .(
    Orthogroup,
    top_go_description,
    top_category
  )],
  by = "Orthogroup",
  all.x = TRUE
)

## fallback logic
og_label_dt[, label := fifelse(
  !is.na(label), label,
  top_go_description
)]

og_label_dt[, category := fifelse(
  !is.na(category), category,
  top_category
)]

## final fallback
og_label_dt[is.na(label), label := "unannotated"]
og_label_dt[is.na(category), category := "unknown"]  



og_label_dt[, label := sub("^(.{40}).*$", "\\1...", label)]



annotation_row <- og_label_dt[
  match(rownames(heat_mat), Orthogroup)
][
  , .(Function = label, Category = category)
]

rownames(annotation_row) <- rownames(heat_mat)



## build labels aligned to best_ids
og_label_vec <- og_go_dt[
  n_genes_with_go >= 2
][
  order(Orthogroup, -species_support, -gene_support, -n_genes_with_go)
][
  , .SD[1], by = Orthogroup
][
  , .(Orthogroup, Function = go_description)
]

## fallback to neutral assignment if available
if ("og_func_summary_dt" %in% ls()) {
  og_label_vec <- merge(
    og_label_vec,
    og_func_summary_dt[, .(Orthogroup, top_go_description)],
    by = "Orthogroup",
    all = TRUE
  )

  og_label_vec[is.na(Function), Function := top_go_description]
}

## final fallback
og_label_vec[is.na(Function), Function := "unannotated"]

## shorten labels (important for readability)
og_label_vec[, Function := fifelse(
  nchar(Function) > 40,
  paste0(substr(Function, 1, 37), "..."),
  Function
)]

## align to heatmap rows
label_vec <- og_label_vec[
  match(best_ids, Orthogroup),
  Function
]

## combine with OG IDs
nice_rownames <- paste0(best_ids, " | ", label_vec)

## apply
rownames(heat_mat) <- nice_rownames



pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Best high-confidence orthogroups",
  fontsize_row = 6   # reduce for readability
)

cap(pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Best high-confidence orthogroups",
  fontsize_row = 6
))

## ensure sorted by final_score (descending)
setorder(final_dt, -final_score)

best_n <- 28
best_ids <- final_dt$Orthogroup[1:best_n]

## match to matrix
best_idx <- match(best_ids, og_ids_f)

heat_mat <- Xf[best_idx, , drop = FALSE]
rownames(heat_mat) <- best_ids

og_label_vec <- og_go_dt[
  n_genes_with_go >= 2
][
  order(Orthogroup, -species_support, -gene_support, -n_genes_with_go)
][
  , .SD[1], by = Orthogroup
][
  , .(Orthogroup, Function = go_description)
]

if ("og_func_summary_dt" %in% ls()) {
  og_label_vec <- merge(
    og_label_vec,
    og_func_summary_dt[, .(Orthogroup, top_go_description)],
    by = "Orthogroup",
    all = TRUE
  )
  og_label_vec[is.na(Function), Function := top_go_description]
}

og_label_vec[is.na(Function), Function := "unannotated"]

og_label_vec[, Function := fifelse(
  nchar(Function) > 40,
  paste0(substr(Function, 1, 37), "..."),
  Function
)]

label_vec <- og_label_vec[
  match(best_ids, Orthogroup),
  Function
]

rownames(heat_mat) <- paste0(best_ids, " | ", label_vec)

annotation_col <- data.frame(
  holocentric = factor(meta_df$holocentric),
  ploidy = factor(meta_df$ploidy)
)
rownames(annotation_col) <- species

library(pheatmap)

pheatmap(
  heat_mat,

  ## biological normalization
  scale = "row",

  ## keep clustering (structure discovery)
  cluster_rows = TRUE,
  cluster_cols = TRUE,

  ## BUT preserve initial ordering influence
  clustering_method = "ward.D2",

  ## stronger separation of modules
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",

  ## annotation
  annotation_col = annotation_col,

  ## readability
  show_rownames = TRUE,
  fontsize_row = 6,

  main = "Top orthogroups (sorted by score + clustered)"
)

cap(pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 6,
  main = "Top orthogroups (sorted by score + clustered)"
))


pdf("all_plots.pdf", width = 10, height = 8, onefile = TRUE)
for (p in plots) replayPlot(p)
dev.off()
cat("Saved all plots to all_plots.pdf\n")
