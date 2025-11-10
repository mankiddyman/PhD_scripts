#!/bin/bash
# mash_cluster.sh
# Usage: ./mash_cluster.sh <genome.fa> <output_prefix> <num_chromosomes>

GENOME="$1"
OUT_PREFIX="$2"
NUM_CHRS="$3"

# Folders
CHRS_DIR="${OUT_PREFIX}_chroms"
SKETCH_DIR="${OUT_PREFIX}_sketches"
mkdir -p "$CHRS_DIR" "$SKETCH_DIR"

echo "Step 1: Splitting genome into chromosomes (first $NUM_CHRS scaffolds)..."
# Check if chromosomes already exist
if [ -z "$(ls -A $CHRS_DIR 2>/dev/null)" ]; then
    i=0
    awk -v max="$NUM_CHRS" '
      /^>/ {i++; if(i > max) exit; f=sprintf("'"$CHRS_DIR"'/%02d.fa", i)}
      {if(i>0 && i<=max) print > f}
    ' "$GENOME"
    echo "Chromosomes split into $CHRS_DIR/"
else
    echo "Chromosomes already exist in $CHRS_DIR/, skipping split."
fi

# Get all chromosome files
CHRS=("$CHRS_DIR"/*.fa)

echo "Step 2: Computing mash sketches (skipping existing sketches)..."
for chr in "${CHRS[@]}"; do
    SK="$SKETCH_DIR/$(basename $chr).msh"
    if [ ! -f "$SK" ]; then
        echo "Sketching $chr..."
        mash sketch -o "$SKETCH_DIR/$(basename $chr)" "$chr"
    else
        echo "Skipping $chr, sketch exists."
    fi
done

echo "Step 3: Computing pairwise mash distances..."
# Create a list file with all sketch paths
SKETCH_LIST="${OUT_PREFIX}.sketch_list.txt"
ls "$SKETCH_DIR"/*.msh > "$SKETCH_LIST"

# Use mash triangle for all-vs-all comparison, or paste the list
mash paste "${OUT_PREFIX}.combined" "$SKETCH_DIR"/*.msh
mash dist "${OUT_PREFIX}.combined.msh" "${OUT_PREFIX}.combined.msh" > "${OUT_PREFIX}.dist.tsv"

echo "Step 4: Generating heatmap plot..."
micromamba run -n r_lang Rscript -e '
library(tidyverse)
library(pheatmap)

# Read Mash distance table
dist_df <- read_tsv("'"${OUT_PREFIX}.dist.tsv"'", col_names=c("chr1","chr2","dist","pval","shared"))

# Extract just the chromosome number from the full path
dist_df <- dist_df %>%
  mutate(chr1 = basename(chr1),
         chr2 = basename(chr2)) %>%
  mutate(chr1 = gsub("\\\\.fa\\\\.msh$", "", chr1),
         chr2 = gsub("\\\\.fa\\\\.msh$", "", chr2))

# Build symmetric distance matrix
chromosomes <- sort(unique(c(dist_df$chr1, dist_df$chr2)))
dist_matrix <- matrix(0, nrow=length(chromosomes), ncol=length(chromosomes),
                      dimnames=list(chromosomes, chromosomes))

for(i in 1:nrow(dist_df)) {
  dist_matrix[dist_df$chr1[i], dist_df$chr2[i]] <- dist_df$dist[i]
  dist_matrix[dist_df$chr2[i], dist_df$chr1[i]] <- dist_df$dist[i]
}

# Generate heatmap PDF
pheatmap(dist_matrix, cluster_rows=TRUE, cluster_cols=TRUE,
         main="Chromosome similarity (Mash distance)",
         filename=paste0("'"${OUT_PREFIX}.heatmap.pdf"'" ))
'

echo "Done! Outputs:"
echo " - Chromosomes: $CHRS_DIR/"
echo " - Mash sketches: $SKETCH_DIR/"
echo " - Distance table: ${OUT_PREFIX}.dist.tsv"
echo " - Heatmap PDF: ${OUT_PREFIX}.heatmap.pdf"
