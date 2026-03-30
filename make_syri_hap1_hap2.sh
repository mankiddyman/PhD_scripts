#!/bin/bash
set -euo pipefail

# Activate environment
micromamba activate syri_env

# -----------------------------
# Paths
# -----------------------------
BASE="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa"
INPUT_DIR="$BASE/Files/genomes/derived/haplotypes"
OUT_DIR="$BASE/results/syri_hap1_hap2"

mkdir -p "$OUT_DIR"

# Input FASTA files
ref="$INPUT_DIR/D_paradoxa_hap1_chr.fasta"
qry="$INPUT_DIR/D_paradoxa_hap2_chr.fasta"

# Threads
threads=40
sort_threads=16
syri_threads=20

# Work in results directory
cd "$OUT_DIR"

# -----------------------------
# 1) Align hap2 to hap1
# -----------------------------
minimap2 -ax asm5 --eqx -t "$threads" "$ref" "$qry" | \
    samtools sort -@ "$sort_threads" -o hap2_vs_hap1.sorted.bam

samtools index hap2_vs_hap1.sorted.bam

# -----------------------------
# 2) Run SyRI
# -----------------------------
syri \
    -c hap2_vs_hap1.sorted.bam \
    -r "$ref" \
    -q "$qry" \
    -F B \
    -k \
    --nc "$syri_threads"

# -----------------------------
# 3) Prepare genomes.txt for plotsr
# -----------------------------
cat > genomes.txt <<EOF
$ref	hap1	lw:1.5
$qry	hap2	lw:1.5
EOF

# -----------------------------
# 4) Plot
# -----------------------------
plotsr --sr syri.out --genomes genomes.txt -o syri.pdf -S 0.6

echo "✅ DONE: Results in $OUT_DIR"
