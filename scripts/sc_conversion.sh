#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <barcode_id> <input_vcf_gz> <marker_tsv> <output_dir>" >&2
    exit 1
fi

BARCODE_ID="$1"
INPUT_VCF_GZ="$2"
MARKER_TSV="$3"
OUTPUT_DIR="$4"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "sc_conversion.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

require_nonempty_file "$INPUT_VCF_GZ"
require_nonempty_file "$MARKER_TSV"
require_executable_file "$BCFTOOLS_BIN"
require_executable_file "$AWK_BIN"
require_executable_file "$SORT_BIN"

safe_mkdir "$OUTPUT_DIR"

OUT_TXT="${OUTPUT_DIR}/input_corrected.txt"
TMP_TXT="${OUT_TXT}.tmp"

# restart-safe skip
if [[ -s "$OUT_TXT" ]]; then
    log "input_corrected.txt already exists for ${BARCODE_ID}, skipping conversion"
    exit 0
fi

rm -f "$TMP_TXT"

log "Converting marker counts for ${BARCODE_ID}"

# Query VCF as:
# CHROM  POS  REF  ALT  AD
#
# AD is expected as comma-separated allele depths, e.g. "7,5"
# We only keep biallelic sites matching the marker reference exactly.

"$BCFTOOLS_BIN" query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
    "$INPUT_VCF_GZ" \
| "$AWK_BIN" -F '\t' -v OFS='\t' '
    BEGIN {
        while ((getline line < "'"$MARKER_TSV"'") > 0) {
            split(line, f, "\t")
            if (length(f) >= 4) {
                key = f[1] "\t" f[2]
                marker_ref[key] = f[3]
                marker_alt[key] = f[4]
            }
        }
        close("'"$MARKER_TSV"'")
    }
    {
        chrom = $1
        pos   = $2
        ref   = $3
        alt   = $4
        ad    = $5

        key = chrom "\t" pos

        if (!(key in marker_ref)) next
        if (ref != marker_ref[key]) next
        if (alt != marker_alt[key]) next

        n = split(ad, arr, ",")
        if (n < 2) next

        ref_count = arr[1]
        alt_count = arr[2]

        if (ref_count == "." || alt_count == ".") next

        print chrom, pos, ref, ref_count, alt, alt_count
    }
' \
| "$SORT_BIN" -k1,1 -k2,2n > "$TMP_TXT"

move_if_success "$TMP_TXT" "$OUT_TXT"
require_nonempty_file "$OUT_TXT"

log "Finished conversion for ${BARCODE_ID}"
