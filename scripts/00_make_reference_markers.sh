#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# setup
###############################################################################

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

stage_start "00_make_reference_markers"

safe_mkdir "$MARKER_REF_DIR"

###############################################################################
# define paths
###############################################################################

ALIGN_BAM="${MARKER_REF_DIR}/hap2_vs_hap1.sorted.bam"
RAW_VCF="${MARKER_REF_DIR}/hap1_vs_hap2.raw.vcf.gz"
FILTERED_VCF="${MARKER_REF_DIR}/hap1_vs_hap2.filtered.vcf.gz"

###############################################################################
# step 1: align hap2 to hap1
###############################################################################

if ! file_exists_and_nonempty "$ALIGN_BAM"; then
    log "Aligning hap2 to hap1"

    TMP_BAM=$(make_temp_path "$ALIGN_BAM")

    run_in_env minimap2 -ax asm5 -t "$MINIMAP2_THREADS" "$HAP1_FASTA" "$HAP2_FASTA" \
    | run_in_env samtools sort -@ "$SAMTOOLS_THREADS" -o "$TMP_BAM"

    move_if_success "$TMP_BAM" "$ALIGN_BAM"

    run_in_env samtools index "$ALIGN_BAM"
else
    log "Alignment BAM already exists, skipping"
fi

###############################################################################
# step 2: variant calling
###############################################################################

if ! file_exists_and_nonempty "$RAW_VCF"; then
    log "Calling variants"

    TMP_VCF=$(make_temp_path "$RAW_VCF")

    run_in_env bcftools mpileup \
        -f "$HAP1_FASTA" \
        -Ou \
        "$ALIGN_BAM" \
    | run_in_env bcftools call \
        -mv \
        -Oz \
        -o "$TMP_VCF"

    move_if_success "$TMP_VCF" "$RAW_VCF"

    run_in_env bcftools index "$RAW_VCF"
else
    log "Raw VCF exists, skipping"
fi

###############################################################################
# step 3: filtering
###############################################################################

if ! file_exists_and_nonempty "$FILTERED_VCF"; then
    log "Filtering variants (SNPs, biallelic)"

    TMP_VCF=$(make_temp_path "$FILTERED_VCF")

    run_in_env bcftools view \
        -v snps \
        -m2 -M2 \
        "$RAW_VCF" \
    | run_in_env bcftools filter \
        -e 'QUAL<20' \
        -Oz \
        -o "$TMP_VCF"

    move_if_success "$TMP_VCF" "$FILTERED_VCF"

    run_in_env bcftools index "$FILTERED_VCF"
else
    log "Filtered VCF exists, skipping"
fi

###############################################################################
# step 4: convert to marker table
###############################################################################

if ! file_exists_and_nonempty "$MARKER_FILE"; then
    log "Converting VCF to marker table"

    TMP_TSV=$(make_temp_path "$MARKER_FILE")

    run_in_env bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        "$FILTERED_VCF" \
        > "$TMP_TSV"

    move_if_success "$TMP_TSV" "$MARKER_FILE"

    assert_nonempty_file "$MARKER_FILE"
else
    log "Marker TSV exists, skipping"
fi

###############################################################################
# step 5: create BED file
###############################################################################

if ! file_exists_and_nonempty "$MARKER_BED"; then
    log "Creating BED file"

    TMP_BED=$(make_temp_path "$MARKER_BED")

    awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' "$MARKER_FILE" > "$TMP_BED"

    move_if_success "$TMP_BED" "$MARKER_BED"

    assert_nonempty_file "$MARKER_BED"
else
    log "Marker BED exists, skipping"
fi

###############################################################################
# final checks
###############################################################################

assert_nonempty_file "$FILTERED_VCF"
assert_nonempty_file "$MARKER_FILE"
assert_nonempty_file "$MARKER_BED"

log "Marker generation complete"

stage_end "00_make_reference_markers"
