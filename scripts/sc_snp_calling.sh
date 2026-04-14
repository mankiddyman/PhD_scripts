#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <barcode_id> <reference_fasta> <input_bam> <output_dir>" >&2
    exit 1
fi

BARCODE_ID="$1"
REFERENCE_FASTA="$2"
INPUT_BAM="$3"
OUTPUT_DIR="$4"

[[ -n "$REFERENCE_FASTA" ]] || { echo "ERROR: reference_fasta argument is empty" >&2; exit 1; }
[[ -n "$INPUT_BAM" ]] || { echo "ERROR: input_bam argument is empty" >&2; exit 1; }
[[ -n "$OUTPUT_DIR" ]] || { echo "ERROR: output_dir argument is empty" >&2; exit 1; }


PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "sc_snp_calling.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

require_nonempty_file "$REFERENCE_FASTA"
require_nonempty_file "$INPUT_BAM"
require_executable_file "$BCFTOOLS_BIN"
require_executable_file "$TABIX_BIN"

safe_mkdir "$OUTPUT_DIR"

OUT_VCF_GZ="${OUTPUT_DIR}/${BARCODE_ID}.vcf.gz"
OUT_VCF_TBI="${OUT_VCF_GZ}.tbi"
TMP_VCF_GZ="${OUT_VCF_GZ}.tmp"

# restart-safe skip
if [[ -s "$OUT_VCF_GZ" && -s "$OUT_VCF_TBI" ]]; then
    log "VCF already exists for ${BARCODE_ID}, skipping SNP calling"
    exit 0
fi

rm -f "$TMP_VCF_GZ" "${TMP_VCF_GZ}.tbi"

log "Calling SNPs for ${BARCODE_ID}"

"$BCFTOOLS_BIN" mpileup \
    --threads "$BCFTOOLS_THREADS" \
    -f "$REFERENCE_FASTA" \
    -a FORMAT/AD,FORMAT/DP \
    -Ou \
    "$INPUT_BAM" \
| "$BCFTOOLS_BIN" call \
    --threads "$BCFTOOLS_THREADS" \
    -m \
    -v \
    -Oz \
    -o "$TMP_VCF_GZ"

"$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -f -t "$TMP_VCF_GZ"
mv "$TMP_VCF_GZ" "$OUT_VCF_GZ"
mv "${TMP_VCF_GZ}.tbi" "$OUT_VCF_TBI"

require_nonempty_file "$OUT_VCF_GZ"
require_nonempty_file "$OUT_VCF_TBI"

log "Finished SNP calling for ${BARCODE_ID}"
