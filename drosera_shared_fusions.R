suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ---------------------------
# 1) Read and standardize
# ---------------------------
blocks <- fread("/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/results/syntenicBlock_coordinates.csv")

# Helper: make sure your names match the file (edit these)
ANCESTOR <- "N_gra_dom"
REGIA    <- "D_regia"
CAPENSIS <- "D_capensis"

# Genespace file has genome1/genome2 with chr1/chr2 referring to those genomes.
# We'll build a "mapping table" that always has:
#   ext_genome, ext_chr, anc_chr, ext_pos (ordering coordinate along ext chr)
# We use startOrd/endOrd if available; otherwise startBp/endBp.

make_ext_anc_map <- function(dt, ancestor_name, ext_name) {
  dt_sub <- dt[
    (genome1 == ancestor_name & genome2 == ext_name) |
      (genome2 == ancestor_name & genome1 == ext_name)
  ]
  if (nrow(dt_sub) == 0) stop("No rows found for ancestor/ext pair: ", ancestor_name, " vs ", ext_name)

  # Determine which side is ext vs ancestor row-wise
  dt_sub[, ext_genome := ext_name]

  dt_sub[, ext_chr := ifelse(genome1 == ext_name, chr1, chr2)]
  dt_sub[, anc_chr := ifelse(genome1 == ancestor_name, chr1, chr2)]

  # Pick an ordering coordinate along ext chromosome
  # Prefer Ord; fall back to Bp. Use midpoint for robustness.
  has_ord <- all(c("startOrd1","endOrd1","startOrd2","endOrd2") %in% names(dt_sub))
  has_bp  <- all(c("startBp1","endBp1","startBp2","endBp2") %in% names(dt_sub))
 # behaviour if both are present is
  if (has_ord) {
    dt_sub[, ext_mid := ifelse(genome1 == ext_name,
                               (startOrd1 + endOrd1)/2,
                               (startOrd2 + endOrd2)/2)]
  } else if (has_bp) {
    dt_sub[, ext_mid := ifelse(genome1 == ext_name,
                               (startBp1 + endBp1)/2,
                               (startBp2 + endBp2)/2)]
  } else {
    stop("Cannot find ordering columns (Ord or Bp) to order blocks along ext chromosome.")
  }

  dt_sub[, .(ext_genome, ext_chr, anc_chr, ext_mid)]
}

# ---------------------------
# 2) Convert blocks -> segments
#    (collapse adjacent blocks from same ancestor chr)
# ---------------------------
blocks_to_segments <- function(map_dt) {
  setorder(map_dt, ext_chr, ext_mid)

  # For each ext chr, collapse consecutive identical anc_chr into one segment "run"
  map_dt[, anc_chr := as.character(anc_chr)]
  map_dt[, run_id := rleid(anc_chr), by = ext_chr]

  segs <- map_dt[, .(
    anc_chr = anc_chr[1],
    seg_mid = mean(ext_mid),
    n_blocks = .N
  ), by = .(ext_chr, run_id)]

  # segment index = order along chromosome
  setorder(segs, ext_chr, seg_mid)
  segs[, seg_index := seq_len(.N), by = ext_chr]

  # keep only what we need for permutations
  segs[, .(ext_chr, seg_index, anc_chr)]
}

# ---------------------------
# 3) Extract adjacency pairs (fusion edges)
# ---------------------------
segments_to_adjacencies <- function(segs, unique_pairs = TRUE) {
  setorder(segs, ext_chr, seg_index)
  segs[, next_anc := shift(anc_chr, type = "lead"), by = ext_chr]

  adj <- segs[!is.na(next_anc) & anc_chr != next_anc,
              .(a1 = anc_chr, a2 = next_anc)]

  # Make adjacency undirected: AK1|AK5 same as AK5|AK1
  adj[, pair := ifelse(a1 < a2, paste0(a1, "|", a2), paste0(a2, "|", a1))]

  if (unique_pairs) {
    unique(adj[, .(pair)])
  } else {
    # multiset with counts
    adj[, .N, by = pair]
  }
}

# ---------------------------
# 4) Permutation: shuffle anc labels across segments
#    preserving:
#      - number of segments per ext chr (fixed by seg table)
#      - global frequency of each anc_chr (counts across all segments)
# ---------------------------
permute_segments <- function(segs) {
  perm <- copy(segs)
  perm[, anc_chr := sample(anc_chr, size = .N, replace = FALSE)]
  perm
}

# ---------------------------
# 5) Overlap metric
# ---------------------------
overlap_count <- function(adj1, adj2) {
  # adj tables expected: one column 'pair' for unique case
  length(intersect(adj1$pair, adj2$pair))
}

# ---------------------------
# 6) Run end-to-end
# ---------------------------
map_reg <- make_ext_anc_map(blocks, ANCESTOR, REGIA)
map_cap <- make_ext_anc_map(blocks, ANCESTOR, CAPENSIS)

segs_reg <- blocks_to_segments(map_reg)
segs_cap <- blocks_to_segments(map_cap)

adj_reg_obs <- segments_to_adjacencies(segs_reg, unique_pairs = TRUE)
adj_cap_obs <- segments_to_adjacencies(segs_cap, unique_pairs = TRUE)

obs_shared <- overlap_count(adj_reg_obs, adj_cap_obs)

cat("Observed shared adjacency pairs (unique, undirected):", obs_shared, "\n")
cat("Regia unique adjacencies:", nrow(adj_reg_obs), "\n")
cat("Capensis unique adjacencies:", nrow(adj_cap_obs), "\n")



# =========================================================
# EXTRA NULL MODEL: chromosome-specific permutation
# =========================================================

# Shuffle anc_chr labels only within each ext_chr
permute_segments_chr_specific <- function(segs) {
  perm <- copy(segs)
  perm[, anc_chr := sample(anc_chr, size = .N, replace = FALSE), by = ext_chr]
  perm
}

# =========================================================
# SANITY CHECK HELPERS
# =========================================================

# Count anc_chr frequencies globally
count_global_anc <- function(segs) {
  out <- segs[, .N, by = anc_chr][order(anc_chr)]
  out[]
}

# Count anc_chr frequencies within chromosome
count_chr_anc <- function(segs) {
  out <- segs[, .N, by = .(ext_chr, anc_chr)][order(ext_chr, anc_chr)]
  out[]
}

# Check whether global counts are preserved
check_global_counts_preserved <- function(original, permuted) {
  orig <- count_global_anc(original)
  perm <- count_global_anc(permuted)
  setnames(orig, "N", "N_orig")
  setnames(perm, "N", "N_perm")
  chk <- merge(orig, perm, by = "anc_chr", all = TRUE)
  chk[is.na(N_orig), N_orig := 0L]
  chk[is.na(N_perm), N_perm := 0L]
  chk[, identical_count := (N_orig == N_perm)]
  chk[]
}

# Check whether chromosome-specific counts are preserved
check_chr_counts_preserved <- function(original, permuted) {
  orig <- count_chr_anc(original)
  perm <- count_chr_anc(permuted)
  setnames(orig, "N", "N_orig")
  setnames(perm, "N", "N_perm")
  chk <- merge(orig, perm, by = c("ext_chr", "anc_chr"), all = TRUE)
  chk[is.na(N_orig), N_orig := 0L]
  chk[is.na(N_perm), N_perm := 0L]
  chk[, identical_count := (N_orig == N_perm)]
  chk[]
}


# Quick sanity report for a permutation function
run_permutation_sanity_checks <- function(segs, perm_fun, label = "perm") {
  cat("\n============================\n")
  cat("Sanity checks for", label, "\n")
  cat("============================\n")

  perm <- perm_fun(segs)

  cat("\nOriginal first 10 rows:\n")
  print(segs[1:min(10, .N)])

  cat("\nPermuted first 10 rows:\n")
  print(perm[1:min(10, .N)])

  cat("\nGlobal anc_chr counts preserved?\n")
  gchk <- check_global_counts_preserved(segs, perm)
  print(gchk)
  cat("All global counts preserved:", all(gchk$identical_count), "\n")

  cat("\nPer-chromosome anc_chr counts preserved?\n")
  cchk <- check_chr_counts_preserved(segs, perm)
  print(cchk)
  cat("All chromosome-specific counts preserved:", all(cchk$identical_count), "\n")

  invisible(list(global = gchk, chr = cchk, perm = perm))
}

# =========================================================
# RUN SANITY CHECKS ON BOTH NULL MODELS
# =========================================================

# For GLOBAL null:
# - global counts SHOULD be preserved
# - per-chromosome counts usually should NOT all be preserved
sanity_reg_global <- run_permutation_sanity_checks(
  segs_reg, permute_segments, "GLOBAL permutation on D. regia"
)

sanity_cap_global <- run_permutation_sanity_checks(
  segs_cap, permute_segments, "GLOBAL permutation on D. capensis"
)

# For CHROMOSOME-SPECIFIC null:
# - global counts SHOULD be preserved
# - per-chromosome counts SHOULD also be preserved
sanity_reg_chr <- run_permutation_sanity_checks(
  segs_reg, permute_segments_chr_specific, "CHR-SPECIFIC permutation on D. regia"
)

sanity_cap_chr <- run_permutation_sanity_checks(
  segs_cap, permute_segments_chr_specific, "CHR-SPECIFIC permutation on D. capensis"
)


# =========================================================
# PERMUTATION TESTS: global and chromosome-specific
# =========================================================

set.seed(1)
n_perm <- 5000

shared_null_global <- integer(n_perm)
shared_null_chr    <- integer(n_perm)

for (i in seq_len(n_perm)) {

  # -------------------------
  # GLOBAL null
  # -------------------------
  segs_reg_p_g <- permute_segments(segs_reg)
  segs_cap_p_g <- permute_segments(segs_cap)

  adj_reg_p_g <- segments_to_adjacencies(segs_reg_p_g, unique_pairs = TRUE)
  adj_cap_p_g <- segments_to_adjacencies(segs_cap_p_g, unique_pairs = TRUE)

  shared_null_global[i] <- overlap_count(adj_reg_p_g, adj_cap_p_g)

  # -------------------------
  # CHROMOSOME-SPECIFIC null
  # -------------------------
  segs_reg_p_c <- permute_segments_chr_specific(segs_reg)
  segs_cap_p_c <- permute_segments_chr_specific(segs_cap)

  adj_reg_p_c <- segments_to_adjacencies(segs_reg_p_c, unique_pairs = TRUE)
  adj_cap_p_c <- segments_to_adjacencies(segs_cap_p_c, unique_pairs = TRUE)

  shared_null_chr[i] <- overlap_count(adj_reg_p_c, adj_cap_p_c)

  if (i %% 500 == 0) cat("Permutation", i, "of", n_perm, "\n")
}

# =========================================================
# SUMMARY STATS
# =========================================================

summarize_null <- function(null_vec, obs, label) {
  p_emp <- (sum(null_vec >= obs) + 1) / (length(null_vec) + 1)
  z <- (obs - mean(null_vec)) / sd(null_vec)

  cat("\n----------------------------------\n")
  cat(label, "\n")
  cat("----------------------------------\n")
  cat("Observed shared adjacency pairs:", obs, "\n")
  cat("Null mean:", mean(null_vec), "\n")
  cat("Null sd:", sd(null_vec), "\n")
  cat("Empirical p-value (>= observed):", p_emp, "\n")
  cat("Z-score:", z, "\n")

  invisible(list(
    mean = mean(null_vec),
    sd = sd(null_vec),
    p = p_emp,
    z = z
  ))
}

res_global <- summarize_null(shared_null_global, obs_shared, "GLOBAL null")
res_chr    <- summarize_null(shared_null_chr,    obs_shared, "CHROMOSOME-SPECIFIC null")


# =========================================================
# EXTRA SANITY CHECKS ON WHAT IS BEING PLOTTED
# =========================================================

cat("\n============================\n")
cat("Histogram sanity checks\n")
cat("============================\n")
cat("Observed shared value used in BOTH histograms:", obs_shared, "\n")
cat("Range of global null:", min(shared_null_global), "to", max(shared_null_global), "\n")
cat("Range of chr-specific null:", min(shared_null_chr), "to", max(shared_null_chr), "\n")
cat("Observed exceeds global null mean by:", obs_shared - mean(shared_null_global), "\n")
cat("Observed exceeds chr-specific null mean by:", obs_shared - mean(shared_null_chr), "\n")

# optional tabulation of most frequent null outcomes
cat("\nMost frequent values in GLOBAL null:\n")
print(sort(table(shared_null_global), decreasing = TRUE)[1:min(10, length(table(shared_null_global)))])

cat("\nMost frequent values in CHR-SPECIFIC null:\n")
print(sort(table(shared_null_chr), decreasing = TRUE)[1:min(10, length(table(shared_null_chr)))])

# =========================================================
# PLOTTING FUNCTIONS
# =========================================================

hist_plot_null <- function(null_vec, obs, title_txt, subtitle_txt = NULL) {
  hist(null_vec,
       breaks = "FD",
       main = title_txt,
       xlab = "Number of shared adjacency pairs",
       border = "grey30")

  axis(1, at = seq(floor(min(null_vec)), ceiling(max(null_vec)), by = 1))
  abline(v = obs, lwd = 2)

  p_emp <- (sum(null_vec >= obs) + 1) / (length(null_vec) + 1)
  z <- (obs - mean(null_vec)) / sd(null_vec)

  ann <- sprintf(
    "Observed = %d\nE(X) = %.2f\nSD = %.2f\np(X >= Xobs) = %.3g\nz = %.2f\nn = %d",
    obs, mean(null_vec), sd(null_vec), p_emp, z, length(null_vec)
  )

  usr <- par("usr")
  x_pos <- usr[2] - 0.02 * (usr[2] - usr[1])
  y_pos <- usr[4] - 0.05 * (usr[4] - usr[3])
  text(x_pos, y_pos, labels = ann, adj = c(1, 1))

  if (!is.null(subtitle_txt)) {
    mtext(subtitle_txt, side = 3, line = 0.2, cex = 0.8)
  }
}

# record global histogram
hist_plot_global <- function() {
  hist_plot_null(
    null_vec = shared_null_global,
    obs = obs_shared,
    title_txt = "Global permutation null",
    subtitle_txt = "Ancestral labels can move between chromosomes"
  )
}
hist_plot_global()
hist_plot_global_rec <- recordPlot()

# record chromosome-specific histogram
hist_plot_chr <- function() {
  hist_plot_null(
    null_vec = shared_null_chr,
    obs = obs_shared,
    title_txt = "Chromosome-specific permutation null",
    subtitle_txt = "Ancestral labels shuffled only within each extant chromosome"
  )
}
hist_plot_chr()
hist_plot_chr_rec <- recordPlot()


# making a visualisation
library(igraph)


# Your fixed known colors (in your stated order)
fixed_colors <- c(
  scaffold3  = "#D61F27",
  scaffold5  = "#F9AC60",
  scaffold6  = "#FCF7BF",
  scaffold7  = "#75c4fd",
  scaffold9  = "#7a63f1",
  scaffold10 = "#f3bbff"  # this would be the last interpolation color if 6 were generated
)

extra_colors <- c(
  scaffold2 = "#f16393",  # deep purple
  scaffold4 = "#1B9E77"   # teal green
)
node_colors <- c(fixed_colors, extra_colors)
node_colors <- node_colors[order(names(node_colors))]


# adj_reg_obs and adj_cap_obs have column 'pair' like "AK1|AK5"
split_pair <- function(pairs) {
  tstrsplit(pairs, "\\|")
}

reg_pairs <- adj_reg_obs$pair
cap_pairs <- adj_cap_obs$pair

all_pairs <- sort(unique(c(reg_pairs, cap_pairs)))
A <- split_pair(all_pairs)
edges <- data.frame(from = A[[1]], to = A[[2]], stringsAsFactors = FALSE)

edges$in_reg <- all_pairs %in% reg_pairs
edges$in_cap <- all_pairs %in% cap_pairs
edges$class <- ifelse(edges$in_reg & edges$in_cap, "shared",
                      ifelse(edges$in_reg, "regia_only", "capensis_only"))

g <- graph_from_data_frame(edges, directed = FALSE)
layout_circle <- layout_in_circle(g)

# assign colors to nodes
V(g)$color <- node_colors[V(g)$name]

# safety fallback if any unexpected names appear
V(g)$color[is.na(V(g)$color)] <- "grey80"


plot(g,
     layout = layout_circle,
     vertex.size = 35,
     vertex.label.cex = 1,
     vertex.label.family = "sans",
     vertex.frame.color = "black",
     edge.width = ifelse(E(g)$class == "shared", 2, 1),
     edge.lty = ifelse(E(g)$class == "shared", 1, 2)
)

legend("topleft",
       legend = names(node_colors),
       col = node_colors,
       pch = 19,
       pt.cex = 1.5,
       bty = "n")
#add a legend for edge types
legend("topright",
  # make legend have italics
       legend = c("Shared", expression(italic("D. regia")~"only"), expression(italic("D. capensis")~"only")),
       col = "black",
       lty = c(1, 2, 2),
       lwd = c(2, 1, 1),
       bty = "n")
# add main
title_expr <- bquote(
  "Shared vs unique adjacency pairs of " *
  italic(Nepenthes) *
  " ancestral chromosomes in extant " *
  italic(Drosera) *
  " genomes"
)

title(main = title_expr, cex.main = 1.2)


# --- Network plot (recorded) ---
net_plot <- function() {

  # draw network bigger (less internal whitespace)
  plot(
    g,
    layout = layout_circle,
    vertex.size = 50,
    vertex.label.cex = 0.8,
    vertex.label.family = "sans",
    vertex.frame.color = "black",
    edge.width = ifelse(E(g)$class == "shared", 2, 1),
    edge.lty   = ifelse(E(g)$class == "shared", 1, 2),
    margin = 0.02   # igraph plot margin (big win)
  )

  usr <- par("usr")

  # Allow drawing outside plot region
  op <- par(xpd = NA)
  on.exit(par(op), add = TRUE)

  # Left legend: push slightly into the left margin
  legend(
    x = usr[1]-1, y = usr[4],
    legend = names(node_colors),
    col = node_colors,
    pch = 19,
    pt.cex = 1.5,
    bty = "n",
    xjust = 0, yjust = 1,
    inset = c(-0.20, 0)   # negative = into margin
  )

  # Right legend: push slightly into the right margin
  legend(
    x = usr[2]+1, y = usr[4],
    legend = c("Shared",
               expression(italic("D. regia")~"only"),
               expression(italic("D. capensis")~"only")),
    col = "black",
    lty = c(1, 2, 2),
    lwd = c(2, 1, 1),
    bty = "n",
    xjust = 1, yjust = 1,
    inset = c(-0.20, 0)
  )

  title_expr <- bquote(
    "Shared vs unique adjacency pairs of " *
      italic(Nepenthes) *
      " ancestral chromosomes in extant " *
      italic(Drosera) *
      " genomes"
  )
  title(main = title_expr, cex.main = 1.2)
}# draw once and record it
net_plot()
net_plot_rec <- recordPlot()

# 3-panel figure: network + both nulls
pdf("adjacency_network_global_and_chrSpecific_nulls.pdf", width = 10, height = 5)
par(mfrow = c(1, 3),
    mar = c(3.2, 3.2, 3.0, 1.2),
    oma = c(0, 0, 0, 0))
replayPlot(net_plot_rec)
replayPlot(hist_plot_global_rec)
replayPlot(hist_plot_chr_rec)
dev.off()
