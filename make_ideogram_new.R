library(data.table)

synteny_file <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/results/syntenicBlock_coordinates.csv"

REF <- "N_gra_dom"
TARGETS <- c("D_capensis", "D_regia")
GENOMES <- TARGETS

dt <- fread(synteny_file)

# ---------- Filter rows (exact python logic) ----------
filtered <- dt[
  (genome1 == REF & genome2 %in% TARGETS) |
    (genome2 == REF & genome1 %in% TARGETS)
]

cat(sprintf("Total filtered rows %d\n", nrow(filtered)))

# quick sanity breakdown (tuple(sorted([genome1, genome2])))
filtered[, pair := ifelse(genome1 < genome2,
                         paste(genome1, genome2, sep="|"),
                         paste(genome2, genome1, sep="|"))]
counts <- filtered[, .N, by=pair][order(pair)]

cat("\nCounts by genome pair:\n")
for(i in seq_len(nrow(counts))) {
  p <- strsplit(counts$pair[i], "\\|")[[1]]
  cat(sprintf("(%s, %s) %d\n", p[1], p[2], counts$N[i]))
}

# ---------- Step 1: collect REF-side (chr_ref, startOrd_ref, endOrd_ref) ----------
# Vectorized (FIX for your error)
ref_intervals <- filtered[, .(
  chr_ref = fifelse(genome1 == REF, chr1, chr2),
  s_ref   = as.integer(fifelse(genome1 == REF, startOrd1, startOrd2)),
  e_ref   = as.integer(fifelse(genome1 == REF, endOrd1,   endOrd2))
)][s_ref < e_ref]

ref_chr_list <- split(ref_intervals[, .(s_ref, e_ref)], ref_intervals$chr_ref)
cat(sprintf("Chromosomes with reference blocks: %d\n", length(ref_chr_list)))

# ---------- Step 2: build atomic segments from unique endpoints per REF chromosome ----------
atoms_by_chr <- lapply(ref_chr_list, function(df) {
  endpoints <- sort(unique(c(df$s_ref, df$e_ref)))
  if(length(endpoints) < 2) return(data.table(s=integer(), e=integer()))
  data.table(
    s = endpoints[-length(endpoints)],
    e = endpoints[-1]
  )[s < e]
})

total_atoms <- sum(vapply(atoms_by_chr, nrow, integer(1)))
cat(sprintf("Total atomic reference segments: %d\n", total_atoms))

cat("\nAtoms per chromosome (first 10 chromosomes):\n")
for(chr_ref in head(names(atoms_by_chr), 10)) {
  cat(chr_ref, nrow(atoms_by_chr[[chr_ref]]), "\n")
}

sizes <- unlist(lapply(atoms_by_chr, function(segs) segs$e - segs$s), use.names=FALSE)
sizes_sorted <- sort(sizes)
median_size <- sizes_sorted[ceiling(length(sizes_sorted)/2)]
cat("Atom size summary (ord units):\n")
cat(sprintf("  n = %d\n", length(sizes)))
cat(sprintf("  min = %d median = %d max = %d\n", min(sizes), median_size, max(sizes)))

# ---------- stable atom id mapping: (chr, start, end) -> integer ----------
atom_by_chr <- list()
atom_key_to_id <- new.env(parent=emptyenv())
aid <- 0L

for(chr_ref in names(atoms_by_chr)) {
  segs <- atoms_by_chr[[chr_ref]][order(s, e)]
  if(nrow(segs) == 0) next
  atom_by_chr[[chr_ref]] <- data.table(
    s = segs$s,
    e = segs$e,
    atomID = aid + seq_len(nrow(segs)) - 1L
  )
  for(i in seq_len(nrow(segs))) {
    key <- paste(chr_ref, segs$s[i], segs$e[i], sep="|")
    atom_key_to_id[[key]] <- atom_by_chr[[chr_ref]]$atomID[i]
  }
  aid <- aid + nrow(segs)
}

cat(sprintf("Total atoms: %d\n", aid))

# ---------- Reference adjacencies A_ref ----------
edge_key <- function(a,b) if(a < b) paste0(a,"|",b) else paste0(b,"|",a)

A_ref <- new.env(parent=emptyenv())
M <- 0L

for(chr_ref in names(atom_by_chr)) {
  segs <- atom_by_chr[[chr_ref]][order(s)]
  if(nrow(segs) < 2) next
  for(i in 1:(nrow(segs)-1)) {
    k <- edge_key(segs$atomID[i], segs$atomID[i+1])
    if(is.null(A_ref[[k]])) {
      A_ref[[k]] <- TRUE
      M <- M + 1L
    }
  }
}
cat(sprintf("M (reference adjacencies): %d\n", M))

# ---------- Precompute atom midpoints per chr ----------
atom_mids <- lapply(atom_by_chr, function(df) {
  data.table(mid = (df$s + df$e)/2, atomID=df$atomID)
})

# ---------- map_ref_to_other_pos (python) ----------
map_ref_to_other_pos <- function(mid_ref, s_ref, e_ref, s_o, e_o, orient) {
  if(e_ref == s_ref) return(NA_real_)
  t <- (mid_ref - s_ref) / (e_ref - s_ref)
  if(orient == "+") s_o + t*(e_o - s_o) else e_o - t*(e_o - s_o)
}

# ---------- Project atoms into other genomes ----------
proj <- list(D_regia = list(), D_capensis = list())

for(r in 1:nrow(filtered)) {
  row <- filtered[r]

  # reference side interval (chr_ref, s_ref, e_ref)
  if(row$genome1 == REF) {
    chr_ref <- row$chr1
    s_ref <- as.integer(row$startOrd1); e_ref <- as.integer(row$endOrd1)
    other <- row$genome2
    chr_other <- row$chr2
    s_o <- as.integer(row$startOrd2); e_o <- as.integer(row$endOrd2)
  } else {
    chr_ref <- row$chr2
    s_ref <- as.integer(row$startOrd2); e_ref <- as.integer(row$endOrd2)
    other <- row$genome1
    chr_other <- row$chr1
    s_o <- as.integer(row$startOrd1); e_o <- as.integer(row$endOrd1)
  }

  if(!(other %in% GENOMES)) next

  orient <- trimws(as.character(row$orient))
  if(!(orient %in% c("+","-"))) orient <- "+"

  # normalize intervals (intended swap; if you want the python bug, tell me)
  if(s_ref > e_ref) { tmp <- s_ref; s_ref <- e_ref; e_ref <- tmp }
  if(s_o > e_o)     { tmp <- s_o;   s_o   <- e_o;   e_o   <- tmp }

  mids <- atom_mids[[chr_ref]]
  if(is.null(mids) || nrow(mids) == 0) next

  inside <- mids[mid >= s_ref & mid <= e_ref]
  if(nrow(inside) == 0) next

  pos <- vapply(inside$mid, function(m) map_ref_to_other_pos(m, s_ref, e_ref, s_o, e_o, orient), numeric(1))
  keep <- is.finite(pos)
  if(!any(keep)) next

  rec <- data.table(pos=pos[keep], atomID=inside$atomID[keep])

  if(is.null(proj[[other]][[chr_other]])) proj[[other]][[chr_other]] <- rec
  else proj[[other]][[chr_other]] <- rbind(proj[[other]][[chr_other]], rec)
}

# ---------- Build adjacency sets A_obs ----------
A_obs <- list(D_regia=new.env(parent=emptyenv()), D_capensis=new.env(parent=emptyenv()))
A_obs_size <- c(D_regia=0L, D_capensis=0L)

for(g in GENOMES) {
  for(chr_other in names(proj[[g]])) {
    lst <- proj[[g]][[chr_other]]
    if(is.null(lst) || nrow(lst) < 2) next
    lst <- lst[order(pos)]
    for(i in 1:(nrow(lst)-1)) {
      a <- lst$atomID[i]; b <- lst$atomID[i+1]
      if(a == b) next
      k <- edge_key(a,b)
      if(is.null(A_obs[[g]][[k]])) {
        A_obs[[g]][[k]] <- TRUE
        A_obs_size[[g]] <- A_obs_size[[g]] + 1L
      }
    }
  }
}

cat("Adjacency sizes:\n")
for(g in GENOMES) cat(g, A_obs_size[[g]], "\n")

# ---------- Set ops + summaries ----------
env_keys <- function(env) ls(env, all.names=TRUE)

Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia <- setdiff(Aref, Areg)
B_cap   <- setdiff(Aref, Acap)
X_shared <- length(intersect(B_regia, B_cap))

D_regia <- setdiff(Areg, Aref)
D_cap   <- setdiff(Acap, Aref)

cat(sprintf("Derived regia: %d\n", length(D_regia)))
cat(sprintf("Derived cap  : %d\n", length(D_cap)))

summarize_against_reference <- function(Aref, Ag, name) {
  preserved <- length(intersect(Aref, Ag))
  broken <- length(setdiff(Aref, Ag))
  derived <- length(setdiff(Ag, Aref))
  list(
    genome=name,
    A_ref=length(Aref),
    A_g=length(Ag),
    preserved=preserved,
    broken=broken,
    derived=derived,
    prop_broken=if(length(Aref)>0) broken/length(Aref) else NaN
  )
}

cat("\nSummary vs reference:\n")
print(summarize_against_reference(Aref, Areg, "D_regia"))
print(summarize_against_reference(Aref, Acap, "D_capensis"))

cat("\nShared breaks:\n")
cat("K_regia =", length(B_regia), "\n")
cat("K_cap   =", length(B_cap), "\n")
cat("X (shared broken ancestral adjacencies) =", X_shared, "\n")

cat("\nDerived adjacency counts:\n")
cat("|D_regia| =", length(D_regia), "\n")
cat("|D_cap|   =", length(D_cap), "\n")
cat("|D_shared|=", length(intersect(D_regia, D_cap)), "\n")

E <- (length(B_regia) * length(B_cap)) / M
cat("Expected shared breaks E[X] =", E, "\n")

# ---------- (Extra) breaks per reference chromosome ----------
# Map each reference adjacency edge to its reference chromosome
edge_to_chr <- new.env(parent=emptyenv())
for(chr_ref in names(atom_by_chr)) {
  segs <- atom_by_chr[[chr_ref]][order(s)]
  if(nrow(segs) < 2) next
  for(i in 1:(nrow(segs)-1)) {
    k <- edge_key(segs$atomID[i], segs$atomID[i+1])
    edge_to_chr[[k]] <- chr_ref
  }
}

breaks_by_chr <- function(B, genome_name) {
  if(length(B) == 0) return(data.table(genome=genome_name, chr=character(), n_breaks=integer()))
  chr <- vapply(B, function(k) {
    v <- edge_to_chr[[k]]
    if(is.null(v)) NA_character_ else v
  }, character(1))
  data.table(genome=genome_name, chr=chr)[!is.na(chr), .N, by=.(genome, chr)][order(-N)]
}

cat("\nBreaks per chromosome (top 10) - D_regia:\n")
print(head(breaks_by_chr(B_regia, "D_regia"), 10))
cat("\nBreaks per chromosome (top 10) - D_capensis:\n")
print(head(breaks_by_chr(B_cap, "D_capensis"), 10))
cat("\nShared breaks per chromosome (top 10):\n")
print(head(breaks_by_chr(intersect(B_regia, B_cap), "shared"), 10))



library(data.table)
library(karyoploteR)
library(GenomicRanges)

# ---- helper from earlier script ----
edge_key <- function(a,b) if(a < b) paste0(a,"|",b) else paste0(b,"|",a)
env_keys <- function(env) ls(env, all.names=TRUE)

# A_ref is an env; D_regia/D_cap are character vectors in the earlier script
Aref_set <- env_keys(A_ref)
D_regia_set <- D_regia
D_cap_set   <- D_cap
D_shared_set <- intersect(D_regia_set, D_cap_set)

# ---- Build breakpoint table for a genome from proj ----
# Breakpoint = boundary between consecutive atoms on the SAME chromosome that forms a DERIVED adjacency (edge not in A_ref)
# Position = midpoint between consecutive projected positions
breakpoints_from_proj <- function(proj_for_genome, genome_name, Aref_set, D_shared_set,
                                 window = 1) {
  out <- list()

  for(chr in names(proj_for_genome)) {
    lst <- proj_for_genome[[chr]]
    if(is.null(lst) || nrow(lst) < 2) next

    lst <- as.data.table(lst)
    setorder(lst, pos)

    # boundaries between consecutive entries
    a <- lst$atomID[-nrow(lst)]
    b <- lst$atomID[-1]
    p1 <- lst$pos[-nrow(lst)]
    p2 <- lst$pos[-1]

    # skip trivial self-adjacency
    keep <- a != b
    if(!any(keep)) next

    a <- a[keep]; b <- b[keep]; p1 <- p1[keep]; p2 <- p2[keep]
    edges <- mapply(edge_key, a, b, USE.NAMES = FALSE)

    # derived if edge not in reference
    is_derived <- !(edges %in% Aref_set)
    if(!any(is_derived)) next

    edges <- edges[is_derived]
    midpos <- (p1[is_derived] + p2[is_derived]) / 2

    # classify shared vs unique using shared-derived edge set
    is_shared <- edges %in% D_shared_set

    # represent breakpoints as tiny intervals for plotting
    # (karyoploteR needs ranges; "window" is in your coordinate units)
    dt_chr <- data.table(
      genome = genome_name,
      chr = chr,
      bp = midpos,
      start = pmax(1, floor(midpos - window)),
      end   = ceiling(midpos + window),
      shared = is_shared,
      edge = edges
    )

    out[[chr]] <- dt_chr
  }

  if(length(out) == 0) {
    return(data.table(genome=genome_name, chr=character(), bp=numeric(),
                      start=integer(), end=integer(), shared=logical(), edge=character()))
  }
  rbindlist(out, use.names=TRUE)
}

# Build breakpoint tables for both genomes
bp_regia <- breakpoints_from_proj(proj[["D_regia"]], "D_regia", Aref_set, D_shared_set, window=1)
bp_cap   <- breakpoints_from_proj(proj[["D_capensis"]], "D_capensis", Aref_set, D_shared_set, window=1)

# ---- choose "best" chromosome per genome ----
# Criteria: maximize shared breakpoints first, then total breakpoints
pick_best_chr <- function(bp_dt) {
  if(nrow(bp_dt) == 0) return(NA_character_)
  score <- bp_dt[, .(
    n_shared = sum(shared, na.rm=TRUE),
    n_total  = .N
  ), by=chr][order(-n_shared, -n_total)]
  score$chr[1]
}

best_reg_chr <- pick_best_chr(bp_regia)
best_cap_chr <- pick_best_chr(bp_cap)

cat("\nChosen chromosomes:\n")
cat("  D_regia    :", best_reg_chr, "\n")
cat("  D_capensis :", best_cap_chr, "\n\n")

# Print per-chromosome summary so you can sanity-check selection
cat("Top chromosomes D_regia (shared, total):\n")
print(head(bp_regia[, .(shared=sum(shared), total=.N), by=chr][order(-shared, -total)], 15))

cat("\nTop chromosomes D_capensis (shared, total):\n")
print(head(bp_cap[, .(shared=sum(shared), total=.N), by=chr][order(-shared, -total)], 15))

# Subset to chosen chromosomes
bp_reg_plot <- bp_regia[chr == best_reg_chr]
bp_cap_plot <- bp_cap[chr == best_cap_chr]

# ---- Build custom genomes for karyoploteR (coordinate scale = projected order-units) ----
# Use max observed breakpoint end as chromosome length, with padding
make_genome_one_chr <- function(chr, length_val) {
  gr <- GRanges(seqnames = chr, ranges = IRanges(start=1, end=as.integer(ceiling(length_val))))
  seqlengths(gr) <- as.integer(ceiling(length_val))
  gr
}

len_reg <- if(nrow(bp_reg_plot) > 0) max(bp_reg_plot$end) else 1
len_cap <- if(nrow(bp_cap_plot) > 0) max(bp_cap_plot$end) else 1
pad <- 10

gen_reg <- make_genome_one_chr(best_reg_chr, len_reg + pad)
gen_cap <- make_genome_one_chr(best_cap_chr, len_cap + pad)

# Convert breakpoint tables to GRanges
dt_to_gr <- function(dt) {
  if(nrow(dt) == 0) return(GRanges())
  GRanges(seqnames=dt$chr, ranges=IRanges(start=dt$start, end=dt$end))
}

reg_shared_gr   <- dt_to_gr(bp_reg_plot[shared==TRUE])
reg_unique_gr   <- dt_to_gr(bp_reg_plot[shared==FALSE])
cap_shared_gr   <- dt_to_gr(bp_cap_plot[shared==TRUE])
cap_unique_gr   <- dt_to_gr(bp_cap_plot[shared==FALSE])

# ---- Plot stacked: capensis on top, regia on bottom ----
library(data.table)
library(karyoploteR)
library(GenomicRanges)

# ---- assumes you already created these from the previous code ----
# best_cap_chr, best_reg_chr
# bp_cap_plot, bp_reg_plot   (data.tables with chr/start/end/shared)
# dt_to_gr(), make_genome_one_chr()

# Recreate GRanges for plotting
dt_to_gr <- function(dt) {
  if(nrow(dt) == 0) return(GRanges())
  GRanges(seqnames=dt$chr, ranges=IRanges(start=dt$start, end=dt$end))
}

cap_shared_gr <- dt_to_gr(bp_cap_plot[shared==TRUE])
cap_unique_gr <- dt_to_gr(bp_cap_plot[shared==FALSE])

reg_shared_gr <- dt_to_gr(bp_reg_plot[shared==TRUE])
reg_unique_gr <- dt_to_gr(bp_reg_plot[shared==FALSE])

# Build a "combined genome" with TWO chromosomes so both can appear in one plot
# (each chromosome is one track row in plot.type=2? actually plot.type=2 gives 2 data panels)
# We will use a single 'genome' containing both chromosomes, and plot both chromosomes.
make_genome <- function(chrs, lens) {
  gr <- GRanges(seqnames=chrs, ranges=IRanges(start=1, end=as.integer(lens)))
  seqlengths(gr) <- as.integer(lens)
  gr
}

# chromosome lengths from breakpoint coordinates (projected ord units)
cap_len <- if(nrow(bp_cap_plot) > 0) max(bp_cap_plot$end) else 1
reg_len <- if(nrow(bp_reg_plot) > 0) max(bp_reg_plot$end) else 1
pad <- 10

gen_both <- make_genome(
  chrs = c(best_cap_chr, best_reg_chr),
  lens = c(cap_len + pad, reg_len + pad)
)

# IMPORTANT: cap GRanges must be on cap chr, reg GRanges on reg chr
# (they already are, because bp_cap_plot$chr == best_cap_chr etc.)

# One plot, two data panels (tracks)
kp <- plotKaryotype(genome=gen_both,
                    chromosomes=c(best_cap_chr, best_reg_chr),
                    plot.type=2)

# Add "genome labels" as main title + panel labels
kpAddMainTitle(kp, "Breakpoints (grey = unique, purple = shared)")

# Label each panel explicitly (top = panel 1, bottom = panel 2)
kpAddLabels(kp, labels=paste0("D_capensis — ", best_cap_chr),
            srt=0, pos=3, cex=0.9, r0=0.92, r1=0.98, data.panel=1)

kpAddLabels(kp, labels=paste0("D_regia — ", best_reg_chr),
            srt=0, pos=3, cex=0.9, r0=0.92, r1=0.98, data.panel=2)

# Base numbers for both panels (optional)
kpAddBaseNumbers(kp, data.panel=1, add.units=FALSE)
kpAddBaseNumbers(kp, data.panel=2, add.units=FALSE)

# Plot Capensis in panel 1 (top track)
if(length(cap_unique_gr) > 0) kpPlotRegions(kp, cap_unique_gr, col="grey70", border=NA,
                                            r0=0.15, r1=0.85, data.panel=1)
if(length(cap_shared_gr) > 0) kpPlotRegions(kp, cap_shared_gr, col="#7B2CBF", border=NA,
                                            r0=0.15, r1=0.85, data.panel=1)

# Plot Regia in panel 2 (bottom track)
if(length(reg_unique_gr) > 0) kpPlotRegions(kp, reg_unique_gr, col="grey70", border=NA,
                                            r0=0.15, r1=0.85, data.panel=2)
if(length(reg_shared_gr) > 0) kpPlotRegions(kp, reg_shared_gr, col="#7B2CBF", border=NA,
                                            r0=0.15, r1=0.85, data.panel=2)

# Legend (base graphics overlay)
legend("topright",
       legend=c("Unique breakpoint", "Shared breakpoint"),
       fill=c("grey70", "#7B2CBF"),
       border=NA,
       bty="n",
       cex=0.9)



# Ground truth from edge sets
K_regia <- length(B_regia)
K_cap   <- length(B_cap)
X_edge  <- length(intersect(B_regia, B_cap))

cat("\n=== EDGE-SET (GROUND TRUTH) ===\n")
cat("K_regia =", K_regia, "\n")
cat("K_cap   =", K_cap, "\n")
cat("X_shared =", X_edge, "\n")


# Helper: ensure unique edges (avoid duplicates due to repeated projections)
unique_edges <- function(bp_dt) {
  if(is.null(bp_dt) || nrow(bp_dt) == 0) return(character())
  unique(bp_dt$edge)
}

edges_regia_all <- unique_edges(bp_regia)
edges_cap_all   <- unique_edges(bp_cap)

# K recovered from breakpoint tables (all derived edges represented as breakpoints)
K_regia_df <- length(edges_regia_all)
K_cap_df   <- length(edges_cap_all)
X_df       <- length(intersect(edges_regia_all, edges_cap_all))

cat("\n=== FROM BREAKPOINT DATAFRAMES (ALL CHROMS) ===\n")
cat("K_regia_df =", K_regia_df, "\n")
cat("K_cap_df   =", K_cap_df, "\n")
cat("X_df       =", X_df, "\n")



library(data.table)

# --- helpers ---
env_keys <- function(env) ls(env, all.names=TRUE)
edge_key <- function(a,b) if(a < b) paste0(a,"|",b) else paste0(b,"|",a)

Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia <- setdiff(Aref, Areg)               # K_regia broken ancestral adjacencies
B_cap   <- setdiff(Aref, Acap)               # K_cap
X_true  <- length(intersect(B_regia, B_cap)) # X_shared (true, from sets)

K_regia <- length(B_regia)
K_cap   <- length(B_cap)

# --- extract breakpoints (derived edges) from proj for all chromosomes ---
breakpoints_from_proj <- function(proj_for_genome, Aref) {
  out <- list()
  for(chr in names(proj_for_genome)) {
    lst <- as.data.table(proj_for_genome[[chr]])
    if(nrow(lst) < 2) next
    setorder(lst, pos)

    a <- lst$atomID[-.N]
    b <- lst$atomID[-1L]
    p1 <- lst$pos[-.N]
    p2 <- lst$pos[-1L]

    keep <- a != b
    if(!any(keep)) next

    a <- a[keep]; b <- b[keep]; p1 <- p1[keep]; p2 <- p2[keep]
    e <- mapply(edge_key, a, b, USE.NAMES = FALSE)

    # derived adjacency = not in reference adjacency set
    is_derived <- !(e %in% Aref)
    if(!any(is_derived)) next

    out[[chr]] <- data.table(
      chr = chr,
      edge = e[is_derived],
      bp   = (p1[is_derived] + p2[is_derived]) / 2
    )
  }
  if(length(out) == 0) return(data.table(chr=character(), edge=character(), bp=numeric()))
  rbindlist(out)
}

bp_regia <- breakpoints_from_proj(proj[["D_regia"]], Aref)
bp_cap   <- breakpoints_from_proj(proj[["D_capensis"]], Aref)

# --- ASSERT: extracted edges equal broken edges (K) and intersection equals X ---
edges_regia <- unique(bp_regia$edge)
edges_cap   <- unique(bp_cap$edge)

stopifnot(
  K_regia == length(edges_regia),
  K_cap   == length(edges_cap),
  X_true  == length(intersect(edges_regia, edges_cap)),
  setequal(edges_regia, B_regia),
  setequal(edges_cap,   B_cap)
)

cat(sprintf("\n✅ Breakpoint extraction matches edge-set truth: K_regia=%d K_cap=%d X_shared=%d\n",
            K_regia, K_cap, X_true))

# --- shared flag (for plotting) ---
shared_edges <- intersect(edges_regia, edges_cap)
bp_regia[, shared := edge %in% shared_edges]
bp_cap[,   shared := edge %in% shared_edges]

# --- chromosome choice table + GLOBAL totals check printed alongside ---
score_chr <- function(bp_dt, genome_label) {
  bp_dt[, .(
    genome = genome_label,
    n_total  = uniqueN(edge),
    n_shared = uniqueN(edge[shared == TRUE])
  ), by=chr][order(-n_shared, -n_total)]
}

score_reg <- score_chr(bp_regia, "D_regia")
score_cap <- score_chr(bp_cap,   "D_capensis")

best_reg_chr <- score_reg$chr[1]
best_cap_chr <- score_cap$chr[1]

cat("\n=== Chromosome scoring (top 10) ===\n")
print(rbindlist(list(head(score_cap,10), head(score_reg,10))))

# global totals check (THIS is what tells you if selection/extraction is losing edges)
cat("\n=== GLOBAL totals check (should equal K and X) ===\n")
cat(sprintf("Regia: unique edges in bp_regia = %d (K_regia=%d)\n", uniqueN(bp_regia$edge), K_regia))
cat(sprintf("Cap  : unique edges in bp_cap   = %d (K_cap=%d)\n",   uniqueN(bp_cap$edge),   K_cap))
cat(sprintf("Shared edges (intersection)     = %d (X_shared=%d)\n",
            length(intersect(unique(bp_regia$edge), unique(bp_cap$edge))), X_true))

cat("\nChosen chromosomes:\n")
cat("  D_capensis —", best_cap_chr, "\n")
cat("  D_regia    —", best_reg_chr, "\n")

# OPTIONAL: subset for plotting later
bp_cap_plot <- bp_cap[chr == best_cap_chr]
bp_reg_plot <- bp_regia[chr == best_reg_chr]



library(data.table)

# helpers
env_keys <- function(env) ls(env, all.names=TRUE)
edge_key <- function(a,b) if(a < b) paste0(a,"|",b) else paste0(b,"|",a)

# edge sets from your python-mirroring script
Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia <- setdiff(Aref, Areg)               # broken ancestral adjacencies
B_cap   <- setdiff(Aref, Acap)
K_regia <- length(B_regia)
K_cap   <- length(B_cap)
X_true  <- length(intersect(B_regia, B_cap))

cat(sprintf("\nEDGE-SET TRUTH: K_regia=%d K_cap=%d X_shared=%d\n", K_regia, K_cap, X_true))

# --- find a broken edge (a|b) as an ADJACENT pair in projected order, per chromosome ---
# Returns data.table with one row per occurrence; we will keep the first occurrence.
locate_edge_in_chr <- function(lst_dt, edge) {
  # lst_dt columns: pos, atomID (already sorted by pos)
  ab <- strsplit(edge, "\\|", fixed=FALSE)[[1]]
  a <- as.integer(ab[1]); b <- as.integer(ab[2])

  if(nrow(lst_dt) < 2) return(NULL)

  idx <- seq_len(nrow(lst_dt) - 1L)
  left  <- lst_dt$atomID[idx]
  right <- lst_dt$atomID[idx + 1L]

  hit <- (left == a & right == b) | (left == b & right == a)
  if(!any(hit)) return(NULL)

  j <- idx[which(hit)[1]]  # take first hit
  data.table(
    edge = edge,
    bp = (lst_dt$pos[j] + lst_dt$pos[j + 1L]) / 2
  )
}

# --- build breakpoint table for a genome from its broken-edge set B_g ---
breakpoints_from_broken <- function(proj_for_genome, B_set, genome_label) {
  out <- list()
  missing <- character()

  # Pre-sort each chromosome once
  proj_sorted <- lapply(proj_for_genome, function(x) {
    dt <- as.data.table(x)
    if(nrow(dt) > 1) setorder(dt, pos)
    dt
  })

  for(edge in B_set) {
    found <- FALSE
    for(chr in names(proj_sorted)) {
      lst <- proj_sorted[[chr]]
      rec <- locate_edge_in_chr(lst, edge)
      if(!is.null(rec)) {
        out[[length(out) + 1L]] <- cbind(data.table(genome=genome_label, chr=chr), rec)
        found <- TRUE
        break
      }
    }
    if(!found) missing <- c(missing, edge)
  }

  bp_dt <- if(length(out) == 0) data.table(genome=genome_label, chr=character(), edge=character(), bp=numeric())
           else rbindlist(out, use.names=TRUE)

  # Hard assert: every broken edge should be found somewhere in projections
  if(length(missing) > 0) {
    cat(sprintf("\nERROR: %d/%d broken edges not found as adjacent pairs in projections for %s.\n",
                length(missing), length(B_set), genome_label))
    cat("First 20 missing edges:\n")
    print(head(missing, 20))
    stop("Cannot proceed until missing broken edges are explained (projection duplicates / mapping gaps).")
  }

  # De-duplicate by edge just in case an edge appears multiple times; keep first
  setorder(bp_dt, chr, bp)
  bp_dt <- bp_dt[!duplicated(edge)]

  bp_dt
}

bp_regia <- breakpoints_from_broken(proj[["D_regia"]], B_regia, "D_regia")
bp_cap   <- breakpoints_from_broken(proj[["D_capensis"]], B_cap, "D_capensis")

# --- shared classification is by broken-edge intersection ---
shared_edges <- intersect(B_regia, B_cap)
bp_regia[, shared := edge %in% shared_edges]
bp_cap[,   shared := edge %in% shared_edges]

# --- ASSERT totals match K and X ---
stopifnot(
  uniqueN(bp_regia$edge) == K_regia,
  uniqueN(bp_cap$edge)   == K_cap,
  uniqueN(bp_regia[shared==TRUE]$edge) == X_true,
  uniqueN(bp_cap[shared==TRUE]$edge)   == X_true
)

cat("✅ Breakpoint tables match K and X exactly (broken-edge definition).\n")

# --- chromosome scoring + include global totals check right there ---
score_chr <- function(bp_dt) {
  bp_dt[, .(
    n_total  = uniqueN(edge),
    n_shared = uniqueN(edge[shared==TRUE])
  ), by=chr][order(-n_shared, -n_total)]
}

score_reg <- score_chr(bp_regia)
score_cap <- score_chr(bp_cap)

best_reg_chr <- score_reg$chr[1]
best_cap_chr <- score_cap$chr[1]

cat("\n=== Chromosome scoring (top 10) ===\n")
print(rbindlist(list(
  cbind(genome="D_capensis", head(score_cap,10)),
  cbind(genome="D_regia",    head(score_reg,10))
), use.names=TRUE, fill=TRUE))

cat("\n=== GLOBAL totals check (must equal K and X) ===\n")
cat(sprintf("Regia: total=%d shared=%d  (K=%d X=%d)\n",
            uniqueN(bp_regia$edge), uniqueN(bp_regia[shared==TRUE]$edge), K_regia, X_true))
cat(sprintf("Cap  : total=%d shared=%d  (K=%d X=%d)\n",
            uniqueN(bp_cap$edge), uniqueN(bp_cap[shared==TRUE]$edge), K_cap, X_true))

cat("\nChosen chromosomes:\n")
cat("  D_capensis —", best_cap_chr, "\n")
cat("  D_regia    —", best_reg_chr, "\n")

# For plotting later:
bp_cap_plot <- bp_cap[chr == best_cap_chr]
bp_reg_plot <- bp_regia[chr == best_reg_chr]





library(data.table)

env_keys <- function(env) ls(env, all.names=TRUE)

# --- edge sets truth from your python-mirroring objects ---
Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia <- setdiff(Aref, Areg)
B_cap   <- setdiff(Aref, Acap)

K_regia <- length(B_regia)
K_cap   <- length(B_cap)
X_true  <- length(intersect(B_regia, B_cap))

cat(sprintf("\nEDGE-SET TRUTH: K_regia=%d K_cap=%d X_shared=%d\n", K_regia, K_cap, X_true))

# --- build fast lookup: atomID -> positions by chromosome in a target genome ---
# proj[[genome]][[chr]] has (pos, atomID)
build_atom_pos_index <- function(proj_for_genome) {
  # returns list: chr -> data.table(atomID, pos) sorted; also an atom->rows index via split
  chr_tables <- lapply(proj_for_genome, function(x) {
    dt <- as.data.table(x)
    if(nrow(dt) == 0) return(dt)
    setorder(dt, atomID, pos)
    dt
  })
  chr_tables
}

# For one chromosome table, get all positions for atom a and atom b, choose closest pair
closest_pair_midpoint <- function(dt_chr, a, b) {
  pa <- dt_chr[atomID == a, pos]
  pb <- dt_chr[atomID == b, pos]
  if(length(pa) == 0 || length(pb) == 0) return(NULL)

  # find closest pair by absolute difference (O(n*m), but usually small; safe here)
  best_d <- Inf
  best_mid <- NA_real_
  for(x in pa) {
    d <- abs(pb - x)
    j <- which.min(d)
    if(d[j] < best_d) {
      best_d <- d[j]
      best_mid <- (x + pb[j]) / 2
    }
  }
  list(mid=best_mid, dist=best_d)
}

# Map a set of broken edges B_set into breakpoint rows on the target genome
# One row per edge (if mappable). Picks best chromosome = smallest dist.
map_broken_edges_to_breakpoints <- function(B_set, proj_for_genome, genome_label) {
  chr_tables <- build_atom_pos_index(proj_for_genome)

  out <- vector("list", length(B_set))
  missing <- character()

  for(i in seq_along(B_set)) {
    edge <- B_set[i]
    ab <- strsplit(edge, "\\|")[[1]]
    a <- as.integer(ab[1]); b <- as.integer(ab[2])

    best <- NULL
    best_chr <- NA_character_

    for(chr in names(chr_tables)) {
      dt_chr <- chr_tables[[chr]]
      if(is.null(dt_chr) || nrow(dt_chr) == 0) next
      cand <- closest_pair_midpoint(dt_chr, a, b)
      if(is.null(cand)) next
      if(is.null(best) || cand$dist < best$dist) {
        best <- cand
        best_chr <- chr
      }
    }

    if(is.null(best)) {
      missing <- c(missing, edge)
    } else {
      out[[i]] <- data.table(genome=genome_label, chr=best_chr, edge=edge, bp=best$mid, dist=best$dist)
    }
  }

  bp <- rbindlist(out, use.names=TRUE, fill=TRUE)
  if(nrow(bp) == 0) bp <- data.table(genome=genome_label, chr=character(), edge=character(), bp=numeric(), dist=numeric())

  # Report (don’t silently drop!)
  if(length(missing) > 0) {
    cat(sprintf("\nWARNING: %d/%d broken edges could not be mapped to %s (atoms never co-occur in projections).\n",
                length(missing), length(B_set), genome_label))
    cat("First 20 missing:\n")
    print(head(missing, 20))
  }

  # ensure 1 row per edge
  bp <- bp[!duplicated(edge)]
  setorder(bp, chr, bp)
  bp
}

bp_regia <- map_broken_edges_to_breakpoints(B_regia, proj[["D_regia"]], "D_regia")
bp_cap   <- map_broken_edges_to_breakpoints(B_cap,   proj[["D_capensis"]], "D_capensis")

# shared = shared broken edges (this MUST give X_true in both)
shared_edges <- intersect(B_regia, B_cap)
bp_regia[, shared := edge %in% shared_edges]
bp_cap[,   shared := edge %in% shared_edges]

# --- ASSERT the mapping preserved edge counts (allowing unmapped edges if warning happened) ---
if(uniqueN(bp_regia$edge) != K_regia || uniqueN(bp_cap$edge) != K_cap) {
  cat("\nNOTE: bp tables do not contain all edges (see WARNING above). Selection/plot will be based on mapped edges only.\n")
}

# These assertions are strict *only if* no missing happened:
if(uniqueN(bp_regia$edge) == K_regia && uniqueN(bp_cap$edge) == K_cap) {
  stopifnot(
    uniqueN(bp_regia[shared==TRUE]$edge) == X_true,
    uniqueN(bp_cap[shared==TRUE]$edge)   == X_true
  )
  cat("✅ bp tables include all broken edges; shared counts match X.\n")
}

# --- chromosome scoring + global check ---
score_chr <- function(bp_dt, genome_label) {
  bp_dt[, .(
    genome = genome_label,
    n_total  = uniqueN(edge),
    n_shared = uniqueN(edge[shared==TRUE])
  ), by=chr][order(-n_shared, -n_total)]
}

score_reg <- score_chr(bp_regia, "D_regia")
score_cap <- score_chr(bp_cap,   "D_capensis")

best_reg_chr <- score_reg$chr[1]
best_cap_chr <- score_cap$chr[1]

cat("\n=== Chromosome scoring (top 10) ===\n")
print(rbindlist(list(head(score_cap,10), head(score_reg,10)), use.names=TRUE, fill=TRUE))

cat("\n=== GLOBAL totals check (mapped edges) ===\n")
cat(sprintf("Regia mapped: total=%d shared=%d  | target K=%d X=%d\n",
            uniqueN(bp_regia$edge), uniqueN(bp_regia[shared==TRUE]$edge), K_regia, X_true))
cat(sprintf("Cap   mapped: total=%d shared=%d  | target K=%d X=%d\n",
            uniqueN(bp_cap$edge), uniqueN(bp_cap[shared==TRUE]$edge), K_cap, X_true))

cat("\nChosen chromosomes:\n")
cat("  D_capensis —", best_cap_chr, "\n")
cat("  D_regia    —", best_reg_chr, "\n")

# For plotting later
bp_cap_plot <- bp_cap[chr == best_cap_chr]
bp_reg_plot <- bp_regia[chr == best_reg_chr]



library(data.table)

env_keys <- function(env) ls(env, all.names=TRUE)

Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia <- setdiff(Aref, Areg)
B_cap   <- setdiff(Aref, Acap)

K_regia <- length(B_regia)
K_cap   <- length(B_cap)
X_true  <- length(intersect(B_regia, B_cap))

cat(sprintf("\nPYTHON-STYLE (set-based) truth:\nK_regia=%d  K_cap=%d  X_shared=%d\n",
            K_regia, K_cap, X_true))

# atoms that actually have positions in proj (so we can plot anything)
atoms_present <- function(proj_for_genome) {
  unique(unlist(lapply(proj_for_genome, function(x) as.data.table(x)$atomID), use.names=FALSE))
}

reg_atoms <- atoms_present(proj[["D_regia"]])
cap_atoms <- atoms_present(proj[["D_capensis"]])

edge_atoms <- function(edge_vec) {
  m <- t(vapply(strsplit(edge_vec, "\\|"), function(z) as.integer(z[1:2]), integer(2)))
  list(a=m[,1], b=m[,2])
}

# keep only broken edges where BOTH atoms exist in that genome's projection
filter_callable <- function(B_set, present_atoms) {
  ab <- edge_atoms(B_set)
  keep <- (ab$a %in% present_atoms) & (ab$b %in% present_atoms)
  B_set[keep]
}

B_regia_call <- filter_callable(B_regia, reg_atoms)
B_cap_call   <- filter_callable(B_cap,   cap_atoms)

K_regia_call <- length(B_regia_call)
K_cap_call   <- length(B_cap_call)
X_call       <- length(intersect(B_regia_call, B_cap_call))

cat(sprintf("\nCALLABLE (mappable to target coords):\nK_regia_call=%d  K_cap_call=%d  X_call=%d\n",
            K_regia_call, K_cap_call, X_call))


library(data.table)
library(karyoploteR)
library(GenomicRanges)

env_keys <- function(env) ls(env, all.names=TRUE)

# --- edge sets (same as your python logic) ---
Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia <- setdiff(Aref, Areg)  # broken ancestral adjacencies in regia
B_cap   <- setdiff(Aref, Acap)  # broken ancestral adjacencies in capensis

B_shared <- intersect(B_regia, B_cap)
cat(sprintf("K_regia=%d  K_cap=%d  X_shared=%d\n", length(B_regia), length(B_cap), length(B_shared)))

# --- Build a lookup: reference edge -> (chr, breakpoint position in ORD units) ---
# Since A_ref edges are only between consecutive atoms on each chr, we can locate their boundary.
edge_key <- function(a,b) if(a < b) paste0(a,"|",b) else paste0(b,"|",a)

edge_to_refpos <- new.env(parent=emptyenv())

for(chr in names(atom_by_chr)) {
  segs <- as.data.table(atom_by_chr[[chr]])
  setorder(segs, s, e)
  if(nrow(segs) < 2) next

  for(i in 1:(nrow(segs)-1)) {
    a <- segs$atomID[i]
    b <- segs$atomID[i+1]
    k <- edge_key(a,b)

    # breakpoint is the boundary between consecutive atoms in reference order
    # For atoms built from endpoints, this boundary is typically segs$e[i] == segs$s[i+1]
    bp <- (segs$e[i] + segs$s[i+1]) / 2

    edge_to_refpos[[k]] <- list(chr=chr, bp=bp)
  }
}

# --- Convert a set of broken edges into GRanges on the reference ---
edges_to_gr <- function(edge_set, window=0.4) {
  # window is in ORD units; keep small so it looks like a tick/point
  rows <- vector("list", length(edge_set))
  j <- 0L
  for(k in edge_set) {
    rec <- edge_to_refpos[[k]]
    if(is.null(rec)) next  # should not happen if edge_set ⊆ A_ref
    j <- j + 1L
    rows[[j]] <- data.table(chr=rec$chr, bp=rec$bp)
  }
  if(j == 0) return(GRanges())
  dt <- rbindlist(rows[1:j])
  dt[, start := pmax(1, floor(bp - window))]
  dt[, end   := ceiling(bp + window)]
  GRanges(seqnames=dt$chr, ranges=IRanges(start=dt$start, end=dt$end))
}

# Shared + unshared
B_any      <- union(B_regia, B_cap)
B_unshared <- setdiff(B_any, B_shared)

gr_unshared <- edges_to_gr(B_unshared)
gr_shared   <- edges_to_gr(B_shared)

# --- Build a custom karyoploteR genome for N_gra_dom in ORD units ---
# Chromosome length = max endpoint of last atom on that chr
ref_chr <- names(atom_by_chr)
ref_len <- vapply(atom_by_chr, function(segs) max(as.data.table(segs)$e), numeric(1))

ref_genome <- GRanges(seqnames=ref_chr, ranges=IRanges(start=1, end=as.integer(ref_len)))
seqlengths(ref_genome) <- as.integer(ref_len)

# --- Plot ---
kp <- plotKaryotype(genome=ref_genome, plot.type=1)
kpAddMainTitle(kp, "N_gra_dom reference: broken ancestral adjacencies (grey=unshared, purple=shared)")
kpAddBaseNumbers(kp, add.units=FALSE)

# unshared (grey) then shared (purple) on top
if(length(gr_unshared) > 0) kpPlotRegions(kp, gr_unshared, col="grey70", border=NA, r0=0.05, r1=0.95)
if(length(gr_shared)   > 0) kpPlotRegions(kp, gr_shared,   col="#7B2CBF", border=NA, r0=0.05, r1=0.95)

legend("topright",
       legend=c("Broken in only one (unshared)", "Broken in both (shared)"),
       fill=c("grey70", "#7B2CBF"),
       border=NA, bty="n", cex=0.9)

# Optional: print counts actually plotted
cat(sprintf("Plotted: unshared=%d shared=%d\n", length(B_unshared), length(B_shared)))


library(data.table)
library(karyoploteR)
library(GenomicRanges)

# helper
edge_key <- function(a,b) if(a < b) paste0(a,"|",b) else paste0(b,"|",a)

# Build reference bp-boundary lookup for each ancestral adjacency
edge_to_refbp <- new.env(parent=emptyenv())

for(chr in names(atom_by_chr)) {

  segs <- as.data.table(atom_by_chr[[chr]])
  setorder(segs, s, e)

  if(nrow(segs) < 2) next

  for(i in 1:(nrow(segs)-1)) {

    a <- segs$atomID[i]
    b <- segs$atomID[i+1]
    k <- edge_key(a,b)

    # Recover BP boundary:
    # find a synteny block row that spans this atom boundary on the reference
    hits <- filtered[
      (genome1 == REF & chr1 == chr &
       startOrd1 <= segs$e[i] & endOrd1 >= segs$s[i+1]) |
      (genome2 == REF & chr2 == chr &
       startOrd2 <= segs$e[i] & endOrd2 >= segs$s[i+1])
    ]

    if(nrow(hits) == 0) next

    # Use midpoint of reference BP interval as boundary
    if(hits$genome1[1] == REF) {
      bp <- (hits$endBp1[1] + hits$startBp1[1]) / 2
    } else {
      bp <- (hits$endBp2[1] + hits$startBp2[1]) / 2
    }

    edge_to_refbp[[k]] <- list(chr=chr, bp=bp)
  }
}



# Edge sets
B_shared <- intersect(B_regia, B_cap)
B_any    <- union(B_regia, B_cap)
B_unshared <- setdiff(B_any, B_shared)

edges_to_gr <- function(edge_set, window_bp=5e4) {
  rows <- vector("list", length(edge_set))
  j <- 0L

  for(k in edge_set) {
    rec <- edge_to_refbp[[k]]
    if(is.null(rec)) next
    j <- j + 1L
    rows[[j]] <- data.table(
      chr=rec$chr,
      start=as.integer(rec$bp - window_bp),
      end  =as.integer(rec$bp + window_bp)
    )
  }

  if(j == 0) return(GRanges())
  dt <- rbindlist(rows[1:j])
  dt$start <- pmax(dt$start, 1)
  GRanges(seqnames=dt$chr, ranges=IRanges(start=dt$start, end=dt$end))
}

gr_shared   <- edges_to_gr(B_shared)
gr_unshared <- edges_to_gr(B_unshared)


# Reference chromosome lengths from synteny table
ref_chr_len <- filtered[
  genome1 == REF,
  .(len = max(endBp1)),
  by=chr1
]

ref_genome <- GRanges(
  seqnames = ref_chr_len$chr1,
  ranges   = IRanges(start=1, end=ref_chr_len$len)
)
seqlengths(ref_genome) <- ref_chr_len$len

pdf("ideogram_breakpoints_fig2.pdf", width=8, height=4)
kp <- plotKaryotype(genome=ref_genome, plot.type=1)
kpAddMainTitle(kp,
  "N_gra_dom reference (bp): broken ancestral synteny\npurple = shared (Regia + Capensis), grey = unshared"
)

kpAddBaseNumbers(kp, add.units=TRUE)
colunshared="#000000"
colshared="#7B2CBF"
# plot unshared first
if(length(gr_unshared) > 0)
  kpPlotRegions(kp, gr_unshared, col=colunshared, border=NA, r0=0.05, r1=0.95)

# plot shared on top
if(length(gr_shared) > 0)
  kpPlotRegions(kp, gr_shared, col=colshared, border=NA, r0=0.05, r1=0.95)

legend("topright",
       legend=c("Broken in one genome", "Broken in both genomes"),
       fill=c(colunshared, colshared),
       border=NA, bty="n", cex=0.9)


# save that to pdf with a particular size 20cm wide 10cm tall
dev.off()
# not quite right


library(karyoploteR)
library(GenomicRanges)

# 20cm x 10cm in inches
pdf("ideogram_breakpoints_fig2.pdf", width = 20/2.54, height = 10/2.54, onefile = TRUE)

# Give karyoploteR extra room
pp <- getDefaultPlotParams(plot.type = 1)
pp$leftmargin   <- 0.10
pp$rightmargin  <- 0.32
pp$topmargin    <- 180   # try 180–300 (default is 120)
pp$bottommargin <- 100   # keep default-ish unless you need more room

pp$data1height  <- 160   # optional: shrink the main panel a bit (default is 200)

# IMPORTANT: avoid the canonical-chromosome filtering warning
# Use chromosomes="all" to keep everything in your custom genome
kp <- plotKaryotype(genome = ref_genome,
                    plot.type = 1,
                    plot.params = pp,
                    chromosomes = "all")

# Title: do NOT use y= (not supported in your version)
kpAddMainTitle(kp, "N_gra_dom reference (bp): broken ancestral synteny", cex = 1.0,
line=2
)

kpAddBaseNumbers(kp, add.units = TRUE)

# Plot breaks
if(length(gr_unshared) > 0)
  kpPlotRegions(kp, gr_unshared, col = "grey70", border = NA, r0 = 0.05, r1 = 0.95)

if(length(gr_shared) > 0)
  kpPlotRegions(kp, gr_shared, col = "#7B2CBF", border = NA, r0 = 0.05, r1 = 0.95)

# Legend: allow drawing outside plot region
par(xpd = NA)
legend("topright",
       inset = c(0, 0),   # push into right margin
       legend = c("Broken in one genome", "Broken in both (shared)"),
       fill = c("grey70", "#7B2CBF"),
       border = NA,
       bty = "n",
       cex = 0.9)

dev.off()


# now figure 1
suppressPackageStartupMessages({
  library(karyoploteR)
  library(GenomicRanges)
})

syn_path <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/results/syntenicBlock_coordinates.csv"

# ---- settings ----
nep_scaffold <- "scaffold2"        # name as it appears in the synteny file
nep_chr_name <- "N.gracilis_scaffold2"  # name you want to appear in karyoploteR
cap_col <- "#e4542d"
reg_col <- "#50ac72"
atom_col <- "#4D4D4D"
# ------------------

syn <- read.csv(syn_path, stringsAsFactors = FALSE, check.names = FALSE)

# blocks on Nep coords where genome2 is capensis or regia
sub <- syn[
  trimws(syn$genome1) == "N_gra_dom" &
    trimws(syn$chr1) == nep_scaffold &
    trimws(syn$genome2) %in% c("D_capensis", "D_regia"),
]
stopifnot(nrow(sub) > 0)

blk_df <- data.frame(
  chr     = trimws(sub$chr1),          # Nep scaffold2
  start   = as.numeric(sub$startBp1),  # Nep coords
  end     = as.numeric(sub$endBp1),
  genome2 = trimws(sub$genome2),       # D_capensis / D_regia
  blkID   = sub$blkID,
  g2_chr   = trimws(sub$chr2),         # genome2 scaffold (regia OR capensis)
  g2_start = as.numeric(sub$startBp2),
  g2_end   = as.numeric(sub$endBp2),
  stringsAsFactors = FALSE
)

# clean Nep coords
blk_df <- blk_df[!is.na(blk_df$start) & !is.na(blk_df$end), ]
swap <- blk_df$start > blk_df$end
if (any(swap)) {
  tmp <- blk_df$start[swap]
  blk_df$start[swap] <- blk_df$end[swap]
  blk_df$end[swap] <- tmp
}
blk_df <- blk_df[blk_df$end > blk_df$start, ]

blocks_nep_gr <- makeGRangesFromDataFrame(
  blk_df,
  seqnames.field = "chr",
  start.field    = "start",
  end.field      = "end",
  keep.extra.columns = TRUE
)

# color by genome2
blocks_nep_gr$col <- ifelse(blocks_nep_gr$genome2 == "D_capensis", cap_col, reg_col)

# Define Nep genome (only scaffold2)
nep_len <- max(end(blocks_nep_gr), na.rm = TRUE)
genome_nep_gr <- toGRanges(data.frame(chr = nep_scaffold, start = 0, end = nep_len))

# ---------------- RENAME scaffold2 -> N_gracilis_chr2 (so karyoploteR prints that) ----------------
old_chr <- nep_scaffold
new_chr <- nep_chr_name

seqlevels(genome_nep_gr)[seqlevels(genome_nep_gr) == old_chr] <- new_chr
seqlevels(blocks_nep_gr)[seqlevels(blocks_nep_gr) == old_chr] <- new_chr

# update variable used in plotting calls
nep_scaffold <- new_chr
# --------------------------------------------------------------------------------------------------

# split tracks for Figure 1
cap_gr <- blocks_nep_gr[blocks_nep_gr$genome2 == "D_capensis"]
reg_gr <- blocks_nep_gr[blocks_nep_gr$genome2 == "D_regia"]

# ---------------- Figure 2: build atoms (on Nep coords) ----------------
breaks <- sort(unique(c(0, nep_len, start(blocks_nep_gr), end(blocks_nep_gr))))
atom_starts <- breaks[-length(breaks)]
atom_ends   <- breaks[-1]

keep_nonempty <- atom_ends > atom_starts
atoms_gr <- GRanges(
  nep_scaffold,
  IRanges(atom_starts[keep_nonempty], atom_ends[keep_nonempty])
)

# keep only atoms covered by at least one synteny block
atoms_gr <- atoms_gr[countOverlaps(atoms_gr, blocks_nep_gr) > 0]
stopifnot(length(atoms_gr) > 0)

# ensure atoms seqname matches renamed chromosome (usually already does, but safe)
seqlevels(atoms_gr)[seqlevels(atoms_gr) == old_chr] <- new_chr
# ----------------------------------------------------------------------

# ---------------- PLOT (Figure 1 + Figure 2 in one plot) ----------------

# match the "second one" PDF size (cm -> inches)
pdf("ideogram_breakpoints_fig1.pdf",
    width  = 20/2.54,
    height = 10/2.54,
    onefile = TRUE)

# match the "second one" plot params (margins/heights)
pp <- getDefaultPlotParams(plot.type = 1)
pp$cex <- 0.2  # default is usually 1
pp$leftmargin   <- 0.18
pp$rightmargin  <- 0.02

pp$topmargin    <- 50
pp$bottommargin <- 20

pp$ideogramheight <- 12
pp$data1height    <- 240

kp <- plotKaryotype(genome = genome_nep_gr,
                    chromosomes = "all",
                    plot.type = 1,
                    plot.params = pp)

# ideogram grey overlay (like second one)
kpRect(kp,
       chr = nep_scaffold,
       x0 = 1,
       x1 = nep_len,
       y0 = 0,
       y1 = 1,
       col = "grey80",
       border = "grey50",
       data.panel = "ideogram")

# Figure 1: capensis (upper)
kpPlotRegions(kp, data = cap_gr, col = cap_col, border = NA, r0 = 0.55, r1 = 0.95)

# Figure 1: regia (lower)
kpPlotRegions(kp, data = reg_gr, col = reg_col, border = NA, r0 = 0.15, r1 = 0.50)

# Figure 2: atoms as a top band (adjusted to match second one aesthetics)
kpPlotRegions(kp, data = atoms_gr, col = atom_col, border = NA, r0 = 0.92, r1 = 0.95)

# separators + '*' markers at atom boundaries
sep_x <- sort(unique(start(atoms_gr)))
sep_x <- sep_x[sep_x > 0 & sep_x < nep_len]

if (length(sep_x) > 0) {
  kpSegments(
    kp,
    chr = nep_scaffold,
    x0  = sep_x, x1 = sep_x,
    y0  = rep(0.05, length(sep_x)),
    y1  = rep(0.95, length(sep_x)),
    col = "grey85",
    lwd = 1
  )

}

# legend placement like second one (single legend, outside plot region)
par(xpd = NA)

legend(
  x = par("usr")[2],
  y = par("usr")[4] + 0.01,
  legend = c("Atomic segments (breakpoint-derived)",
             "Synteny to D_capensis",
             "Synteny to D_regia"),
  fill   = c(atom_col, cap_col, reg_col),
  border = NA,
  bty = "n",
  xjust = 1,
  yjust = 0,
  cex = 0.9
)

dev.off()



library(data.table)
library(karyoploteR)
library(GenomicRanges)
library(grid)
library(gridBase)
library(VennDiagram)

# -----------------------------
# User settings (EDIT HERE)
# -----------------------------
OUT_PDF   <- "ideogram_breakpoints_fig2.pdf"
PDF_W_CM  <- 20
PDF_H_CM  <- 10

# All plotted breakpoints will have the same "tick thickness" (window) and same height (r0/r1)
WINDOW_BP <- 5e4     # half-window in bp (tick thickness)
R0 <- 0.05
R1 <- 0.80

# Colors (placeholders for now, as requested)
COL_SHARED   <- "#7B2CBF"  # shared (Regia + Capensis)
COL_REGIA    <- "#50ac72"  # regia-specific (user defines)
COL_CAPENSIS <- "#e4542d"  # capensis-specific (user defines)

# -----------------------------
# Helpers
# -----------------------------
env_keys <- function(env) ls(env, all.names = TRUE)
edge_key <- function(a, b) if(a < b) paste0(a, "|", b) else paste0(b, "|", a)

# Build reference adjacency -> reference breakpoint (bp) lookup
# Requires: atom_by_chr, filtered, REF
build_edge_to_refbp <- function(atom_by_chr, filtered, REF) {
  edge_to_refbp <- new.env(parent = emptyenv())

  for(chr in names(atom_by_chr)) {
    segs <- as.data.table(atom_by_chr[[chr]])
    setorder(segs, s, e)
    if(nrow(segs) < 2) next

    for(i in 1:(nrow(segs) - 1)) {
      a <- segs$atomID[i]
      b <- segs$atomID[i + 1]
      k <- edge_key(a, b)

      # find a synteny block row that spans this atom boundary on the reference
      hits <- filtered[
        (genome1 == REF & chr1 == chr &
           startOrd1 <= segs$e[i] & endOrd1 >= segs$s[i + 1]) |
        (genome2 == REF & chr2 == chr &
           startOrd2 <= segs$e[i] & endOrd2 >= segs$s[i + 1])
      ]

      if(nrow(hits) == 0) next

      # midpoint of reference BP interval
      if(hits$genome1[1] == REF) {
        bp <- (hits$endBp1[1] + hits$startBp1[1]) / 2
      } else {
        bp <- (hits$endBp2[1] + hits$startBp2[1]) / 2
      }

      edge_to_refbp[[k]] <- list(chr = chr, bp = bp)
    }
  }

  edge_to_refbp
}

edges_to_gr <- function(edge_set, edge_to_refbp, window_bp) {
  rows <- vector("list", length(edge_set))
  j <- 0L

  for(k in edge_set) {
    rec <- edge_to_refbp[[k]]
    if(is.null(rec)) next
    j <- j + 1L
    rows[[j]] <- data.table(
      chr = rec$chr,
      start = as.integer(rec$bp - window_bp),
      end   = as.integer(rec$bp + window_bp)
    )
  }

  if(j == 0L) return(GRanges())
  dt <- rbindlist(rows[1:j])
  dt[, start := pmax(start, 1L)]
  GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$start, end = dt$end))
}

build_ref_genome_from_filtered <- function(filtered, REF) {
  ref_chr_len <- filtered[
    genome1 == REF,
    .(len = max(endBp1)),
    by = chr1
  ]

  ref_genome <- GRanges(
    seqnames = ref_chr_len$chr1,
    ranges   = IRanges(start = 1, end = ref_chr_len$len)
  )
  seqlengths(ref_genome) <- ref_chr_len$len
  ref_genome
}

# -----------------------------
# 1) Edge sets: broken adjacencies
# -----------------------------
Aref <- env_keys(A_ref)
Areg <- env_keys(A_obs[["D_regia"]])
Acap <- env_keys(A_obs[["D_capensis"]])

B_regia    <- setdiff(Aref, Areg)  # broken ancestral adjacencies in regia
B_capensis <- setdiff(Aref, Acap)  # broken ancestral adjacencies in capensis

B_shared     <- intersect(B_regia, B_capensis)
B_regia_only <- setdiff(B_regia, B_shared)
B_cap_only   <- setdiff(B_capensis, B_shared)

cat(sprintf("K_regia=%d  K_cap=%d  shared=%d  regia_only=%d  cap_only=%d\n",
            length(B_regia), length(B_capensis), length(B_shared),
            length(B_regia_only), length(B_cap_only)))

# -----------------------------
# 2) Map adjacency edges -> reference bp positions
# -----------------------------
edge_to_refbp <- build_edge_to_refbp(atom_by_chr = atom_by_chr, filtered = filtered, REF = REF)

gr_shared     <- edges_to_gr(B_shared,     edge_to_refbp, window_bp = WINDOW_BP)
gr_regia_only <- edges_to_gr(B_regia_only, edge_to_refbp, window_bp = WINDOW_BP)
gr_cap_only   <- edges_to_gr(B_cap_only,   edge_to_refbp, window_bp = WINDOW_BP)

# -----------------------------
# 3) Custom genome + plotting
# -----------------------------
ref_genome <- build_ref_genome_from_filtered(filtered, REF)

# PDF size in cm -> inches
pdf(OUT_PDF, width = PDF_W_CM / 2.54, height = PDF_H_CM / 2.54, onefile = TRUE)

pp <- getDefaultPlotParams(plot.type = 1)
pp$leftmargin   <- 0.10
pp$rightmargin  <- 0.32
pp$topmargin    <- 180
pp$bottommargin <- 100
pp$data1height <- 320

kp <- plotKaryotype(
  genome      = ref_genome,
  plot.type   = 1,
  plot.params = pp,
  chromosomes = "all"
)

kpAddMainTitle(
  kp,
  "N_gra_dom reference (bp): broken ancestral synteny",
  cex = 1.0,
  line = 2
)

kpAddBaseNumbers(kp, add.units = TRUE)

# Plot in a stable order so overlaps look consistent:
# draw specific-only first, shared last (shared on top)
#DP <- 1  # choose 1 or 2, but use the same for all
#
#kpPlotRegions(kp, gr_regia_only, col=COL_REGIA,    border=NA, r0=R0, r1=R1, data.panel=DP)
#kpPlotRegions(kp, gr_cap_only,   col=COL_CAPENSIS, border=NA, r0=R0, r1=R1, data.panel=DP)
#kpPlotRegions(kp, gr_shared,     col=COL_SHARED,   border=NA, r0=R0, r1=R1, data.panel=DP)
plot_gr_as_rects <- function(kp, gr, col, r0, r1, data.panel=1) {
  if(length(gr)==0) return(invisible())
  kpRect(kp,
         chr=as.character(seqnames(gr)),
         x0=start(gr), x1=end(gr),
         y0=r0, y1=r1,
         col=col, border=NA,
         data.panel=data.panel)
}

DP <- 1
plot_gr_as_rects(kp, gr_regia_only, COL_REGIA,    R0, R1, DP)
plot_gr_as_rects(kp, gr_cap_only,   COL_CAPENSIS, R0, R1, DP)
plot_gr_as_rects(kp, gr_shared,     COL_SHARED,   R0, R1, DP)


par(xpd = NA)

 Draw the legend
legend(
  "topright",
  inset  = c(0, 0),
  legend = c("Broken in Regia only", "Broken in Capensis only", "Broken in both (shared)"),
  fill   = c(COL_REGIA, COL_CAPENSIS, COL_SHARED),
  border = NA,
  bty    = "n",
  cex    = 0.9
)

# Draw the base graphics legend FIRST (before entering grid viewports)
# Create euler diagram and convert to grob
library(eulerr)
library(grid)

# Create the euler fit
fit <- euler(c(
  "Regia" = n_regia_only,
  "Capensis" = n_cap_only,
  "Regia&Capensis" = n_shared
))

# Save to temporary device to capture as grob
tmp_file <- tempfile(fileext = ".png")
png(tmp_file, width = 400, height = 400, bg = "transparent")
plot(fit,
     fills = list(fill = c(COL_REGIA, COL_CAPENSIS, COL_SHARED),
                  alpha = 0.4),
     edges = list(col = c(COL_REGIA, COL_CAPENSIS), lwd = 2.5),
     labels = FALSE,
     quantities = list(cex = 4, fontface = "bold"))
dev.off()

# Read it back and draw as raster
library(png)
venn_img <- readPNG(tmp_file)

# Position it in a viewport
vp <- viewport(x = 0.85, y = 0.12, width = 0.38, height = 0.38, 
               just = c("center", "bottom"))
pushViewport(vp)
grid.raster(venn_img)
popViewport()

# Clean up
unlink(tmp_file)

dev.off()

cat(sprintf("Plotted: regia_only=%d cap_only=%d shared=%d\n",
            length(gr_regia_only), length(gr_cap_only), length(gr_shared)))
