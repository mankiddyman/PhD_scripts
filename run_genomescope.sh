#!/usr/bin/env bash
set -euo pipefail

# Reusable GenomeScope pipeline for FASTA/FASTQ reads
# Author: Aaryan Bhatia

usage() {
    cat <<EOF
Usage:
  $(basename "$0") -i READS.fa -o OUTDIR [options]

Required:
  -i FILE          Input reads file (FASTA or FASTQ)
  -o OUTDIR        Output directory

Optional:
  -k INT           K-mer size (default: 21)
  -t INT           Threads for jellyfish (default: 16)
  -s SIZE          Jellyfish hash size (default: 10G)
  -m INT           Max histogram coverage for jellyfish histo (default: 1000000)
  -p INT           GenomeScope ploidy (default: 2)
  -n STR           Sample name override
  --jf FILE        Use existing jellyfish .jf database instead of recounting
  -h               Show help

Example:
  bash $(basename "$0") \
    -i /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/Files/HiFi_reads/hifi_reads.fasta \
    -o /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/genomescope \
    -k 21 -t 32 -s 20G -p 2
EOF
}

READS=""
OUTDIR=""
K=21
THREADS=16
HASH_SIZE="10G"
HISTO_MAX=1000000
PLOIDY=2
SAMPLE=""
EXISTING_JF=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) READS="$2"; shift 2 ;;
        -o) OUTDIR="$2"; shift 2 ;;
        -k) K="$2"; shift 2 ;;
        -t) THREADS="$2"; shift 2 ;;
        -s) HASH_SIZE="$2"; shift 2 ;;
        -m) HISTO_MAX="$2"; shift 2 ;;
        -p) PLOIDY="$2"; shift 2 ;;
        -n) SAMPLE="$2"; shift 2 ;;
        --jf) EXISTING_JF="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    echo "Error: -o OUTDIR is required" >&2
    usage
    exit 1
fi

if [[ -z "$READS" && -z "$EXISTING_JF" ]]; then
    echo "Error: provide either -i READS or --jf existing.jf" >&2
    usage
    exit 1
fi

if [[ -n "$READS" && ! -f "$READS" ]]; then
    echo "Error: input reads file not found: $READS" >&2
    exit 1
fi

if [[ -n "$EXISTING_JF" && ! -f "$EXISTING_JF" ]]; then
    echo "Error: jellyfish db not found: $EXISTING_JF" >&2
    exit 1
fi

module load jellyfish || true
module load R || true

command -v jellyfish >/dev/null 2>&1 || { echo "Error: jellyfish not found in PATH"; exit 1; }
command -v genomescope.R >/dev/null 2>&1 || command -v Rscript >/dev/null 2>&1 || {
    echo "Error: neither genomescope.R nor Rscript found in PATH"
    exit 1
}

mkdir -p "$OUTDIR"

if [[ -z "$SAMPLE" ]]; then
    if [[ -n "$READS" ]]; then
        b=$(basename "$READS")
        SAMPLE="${b%%.*}"
    else
        b=$(basename "$EXISTING_JF")
        SAMPLE="${b%%.*}"
    fi
fi

WORKDIR="${OUTDIR}/${SAMPLE}_k${K}"
mkdir -p "$WORKDIR"

JF="${WORKDIR}/${SAMPLE}.jf"
HISTO="${WORKDIR}/${SAMPLE}.histo"
GS_OUT="${WORKDIR}/genomescope"

echo "Sample:      $SAMPLE"
echo "Workdir:     $WORKDIR"
echo "k-mer size:  $K"
echo "Threads:     $THREADS"
echo "Hash size:   $HASH_SIZE"
echo "Ploidy:      $PLOIDY"

if [[ -n "$EXISTING_JF" ]]; then
    echo "Using existing jellyfish db: $EXISTING_JF"
    JF="$EXISTING_JF"
else
    echo "Counting k-mers with jellyfish..."
    jellyfish count \
        -m "$K" \
        -s "$HASH_SIZE" \
        -t "$THREADS" \
        -C \
        -o "$JF" \
        "$READS"
fi

echo "Making k-mer histogram..."
jellyfish histo \
    -t "$THREADS" \
    -h "$HISTO_MAX" \
    "$JF" > "$HISTO"

echo "Running GenomeScope..."
if command -v genomescope.R >/dev/null 2>&1; then
    genomescope.R \
        -i "$HISTO" \
        -k "$K" \
        -p "$PLOIDY" \
        -o "$GS_OUT"
else
    # Fallback if genomescope is installed as an R script you run manually
    echo "genomescope.R not in PATH."
    echo "Histogram is ready here: $HISTO"
    echo "Run GenomeScope manually with something like:"
    echo "Rscript /path/to/genomescope.R -i $HISTO -k $K -p $PLOIDY -o $GS_OUT"
fi

echo
echo "Done."
echo "Histogram:   $HISTO"
echo "GenomeScope: $GS_OUT"
