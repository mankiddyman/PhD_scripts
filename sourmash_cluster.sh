#!/bin/bash
# sourmash_cluster.sh
# Usage: ./sourmash_cluster.sh <genome.fa> <output_prefix> <num_chromosomes>

GENOME="$1"
OUT_PREFIX="$2"
NUM_CHRS="$3"

# Folders
CHRS_DIR="${OUT_PREFIX}_chroms"
SIG_DIR="${OUT_PREFIX}_sigs"
mkdir -p "$CHRS_DIR" "$SIG_DIR"

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

echo "Step 2: Computing sourmash signatures (sensitive)..."
for chr in "${CHRS[@]}"; do
    SIG="$SIG_DIR/$(basename $chr).sig"
    if [ ! -f "$SIG" ]; then
        echo "Computing signature for $chr..."
        sourmash compute -k 15 --scaled 1000 "$chr" -o "$SIG"
    else
        echo "Skipping $chr, signature exists."
    fi
done

echo "Step 3: Computing pairwise distances..."
DIST_MATRIX="${OUT_PREFIX}.dist.tsv"

# Sourmash compare produces a symmetric distance matrix
sourmash compare "$SIG_DIR"/*.sig -o "${OUT_PREFIX}.compare.npy"
micromamba activate r_lang
# Convert .npy matrix to TSV for R heatmap
python3 - <<EOF
import numpy as np
import pandas as pd
import os

compare_file = "${OUT_PREFIX}.compare.npy"
sigs = sorted([os.path.basename(f).replace(".fa.sig","") for f in os.listdir("${SIG_DIR}") if f.endswith(".sig")])

# Load distance matrix
matrix = np.load(compare_file)
df = pd.DataFrame(matrix, index=sigs, columns=sigs)
df.to_csv("${DIST_MATRIX}", sep="\t")
EOF

echo "Step 4: Generating heatmap..."
micromamba run -n r_lang Rscript -e '
library(tidyverse)
library(pheatmap)

# Read the distance TSV
dist_df <- read_tsv("'"${DIST_MATRIX}"'")

# Convert to matrix for heatmap
dist_matrix <- as.matrix(dist_df[,-1])
rownames(dist_matrix) <- dist_df[[1]]

# Generate heatmap PDF
pheatmap(dist_matrix, cluster_rows=TRUE, cluster_cols=TRUE,
         main="Chromosome similarity (Sourmash distance)",
         filename=paste0("'"${OUT_PREFIX}.heatmap.pdf"'"))
'

echo "Done! Outputs:"
echo " - Chromosomes: $CHRS_DIR/"
echo " - Sourmash signatures: $SIG_DIR/"
echo " - Distance table: ${DIST_MATRIX}"
echo " - Heatmap PDF: ${OUT_PREFIX}.heatmap.pdf"

