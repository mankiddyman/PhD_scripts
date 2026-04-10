#!/usr/bin/env bash
set -euo pipefail

# Author: Aaryan Bhatia
# Reusable Helixer wrapper for one or more genome FASTA files

usage() {
    cat <<EOF
Usage:
  $(basename "$0") -n PRIMARY_COUNT -o OUTDIR -s HELIXER_SIF [options] genome1.fa [genome2.fa ...]

Required:
  -n PRIMARY_COUNT   Number of primary scaffolds/chromosomes to split individually
  -o OUTDIR          Parent output directory
  -s HELIXER_SIF     Path to Helixer singularity image

Optional:
  -l LINEAGE         Helixer lineage (default: land_plant)
  -t THREADS         Reserved for future use
  --skip-existing    Skip Helixer runs if output GFF3 already exists
  -h                 Show this help

Example:
  $(basename "$0") -n 7 \
    -o /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/Helixer \
    -s /netscratch/dep_mercier/grp_marques/Aaryan/methods/helixer_docker/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif \
    --skip-existing \
    /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/Files/aaryan_asm_2026/hap1.fasta \
    /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/Files/aaryan_asm_2026/hap2.fasta
EOF
}

PRIMARY_COUNT=""
OUTDIR=""
HELIXER_SIF=""
LINEAGE="land_plant"
SKIP_EXISTING=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n)
            PRIMARY_COUNT="$2"
            shift 2
            ;;
        -o)
            OUTDIR="$2"
            shift 2
            ;;
        -s)
            HELIXER_SIF="$2"
            shift 2
            ;;
        -l)
            LINEAGE="$2"
            shift 2
            ;;
        --skip-existing)
            SKIP_EXISTING=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [[ -z "$PRIMARY_COUNT" || -z "$OUTDIR" || -z "$HELIXER_SIF" || $# -lt 1 ]]; then
    usage
    exit 1
fi

module load singularity

[[ -f "$HELIXER_SIF" ]] || { echo "Helixer SIF not found: $HELIXER_SIF" >&2; exit 1; }

run_helixer_one_genome() {
    local genome="$1"
    [[ -f "$genome" ]] || { echo "Genome file not found: $genome" >&2; return 1; }

    local sample base workdir split_dir results_dir merged_gff
    base="$(basename "$genome")"
    sample="${base%.*}"

    workdir="${OUTDIR}/${sample}"
    split_dir="${workdir}/split_genome"
    results_dir="${workdir}/annotations"
    merged_gff="${results_dir}/${sample}_helixer_merged.gff3"

    mkdir -p "$split_dir" "$results_dir"

    echo
    echo "========================================"
    echo "Genome:   $genome"
    echo "Sample:   $sample"
    echo "Workdir:  $workdir"
    echo "========================================"

    echo "Splitting genome..."
    awk -v split_dir="$split_dir" -v n="$PRIMARY_COUNT" '
        /^>/ {
            count++
            if (count <= n) {
                file = sprintf("%s/chr_%d.fa", split_dir, count)
            } else {
                file = sprintf("%s/debris.fa", split_dir)
            }
        }
        { print >> file }
    ' "$genome"

    echo "Running Helixer on primary scaffolds..."
    shopt -s nullglob
    for chr_fasta in "$split_dir"/chr_*.fa; do
        local chr_name out_gff
        chr_name="$(basename "$chr_fasta" .fa)"
        out_gff="${results_dir}/${chr_name}_helixer.gff3"

        if [[ $SKIP_EXISTING -eq 1 && -s "$out_gff" ]]; then
            echo "Skipping existing: $out_gff"
            continue
        fi

        echo "Processing: $chr_name"
        singularity run --nv \
            -B "$workdir:$workdir" \
            "$HELIXER_SIF" Helixer.py \
            --fasta-path "$chr_fasta" \
            --lineage "$LINEAGE" \
            --gff-output-path "$out_gff"
    done
    shopt -u nullglob

    if [[ -s "$split_dir/debris.fa" ]]; then
        local debris_gff="${results_dir}/debris_helixer.gff3"
        if [[ $SKIP_EXISTING -eq 1 && -s "$debris_gff" ]]; then
            echo "Skipping existing: $debris_gff"
        else
            echo "Processing debris..."
            singularity run --nv \
                -B "$workdir:$workdir" \
                "$HELIXER_SIF" Helixer.py \
                --fasta-path "$split_dir/debris.fa" \
                --lineage "$LINEAGE" \
                --gff-output-path "$debris_gff"
        fi
    else
        echo "No debris sequences found."
    fi

    echo "Merging GFF3 files..."
    {
        echo "##gff-version 3"
        shopt -s nullglob
        for gff in "$results_dir"/chr_*_helixer.gff3 "$results_dir"/debris_helixer.gff3; do
            awk 'BEGIN{skip=1} /^##gff-version 3/ {next} {print}' "$gff"
        done
        shopt -u nullglob
    } > "$merged_gff"

    echo "Done: $merged_gff"
}

for genome in "$@"; do
    run_helixer_one_genome "$genome"
done

echo
echo "All genomes processed."
