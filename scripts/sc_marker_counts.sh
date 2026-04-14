#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <barcode_id> <reference_fasta> <input_bam> <marker_tsv> <output_dir>" >&2
    exit 1
fi

BARCODE_ID="$1"
REFERENCE_FASTA="$2"
INPUT_BAM="$3"
MARKER_TSV="$4"
OUTPUT_DIR="$5"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "sc_marker_counts.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

require_nonempty_file "$REFERENCE_FASTA"
require_nonempty_file "$INPUT_BAM"
require_nonempty_file "$MARKER_TSV"
require_nonempty_file "$MARKER_POS_BED"
require_executable_file "$SAMTOOLS_BIN"
require_executable_file "$AWK_BIN"
require_executable_file "$SORT_BIN"

safe_mkdir "$OUTPUT_DIR"

OUT_TXT="${OUTPUT_DIR}/input_corrected.txt"
TMP_PILEUP="${OUT_TXT}.pileup.tmp"
TMP_TXT="${OUT_TXT}.tmp"

if [[ -s "$OUT_TXT" ]]; then
    log "input_corrected.txt already exists for ${BARCODE_ID}, skipping marker counting"
    exit 0
fi

rm -f "$TMP_PILEUP" "$TMP_TXT"

log "Counting marker alleles for ${BARCODE_ID}"

# Generate mpileup only at marker positions
"$SAMTOOLS_BIN" mpileup \
    -f "$REFERENCE_FASTA" \
    -l "$MARKER_POS_BED" \
    "$INPUT_BAM" > "$TMP_PILEUP"

# Join pileup with marker definitions and count ref/alt alleles from pileup bases
"$AWK_BIN" -F '\t' -v OFS='\t' '
    BEGIN {
        while ((getline < "'"$MARKER_TSV"'") > 0) {
            key = $1 "\t" $2
            marker_ref[key] = $3
            marker_alt[key] = $4
        }
        close("'"$MARKER_TSV"'")
    }

    function count_base(bases, target,   i,c,nextc,n,clean,count,indel_len,indel_seq,j,numstr) {
        clean = bases
        gsub(/\^./, "", clean)
        gsub(/\$/, "", clean)

        n = length(clean)
        count = 0
        i = 1

        while (i <= n) {
            c = substr(clean, i, 1)

            if (c == "+" || c == "-") {
                i++
                numstr = ""
                while (i <= n && substr(clean, i, 1) ~ /[0-9]/) {
                    numstr = numstr substr(clean, i, 1)
                    i++
                }
                indel_len = numstr + 0
                i += indel_len
                continue
            }

            if (target == "REF") {
                if (c == "." || c == ",") count++
            } else {
                if (toupper(c) == toupper(target)) count++
            }

            i++
        }

        return count
    }

    {
        chrom = $1
        pos   = $2
        ref_base_from_pileup = $3
        depth = $4
        bases = $5

        key = chrom "\t" pos
        if (!(key in marker_ref)) next

        ref = marker_ref[key]
        alt = marker_alt[key]

        ref_count = count_base(bases, "REF")
        alt_count = count_base(bases, alt)

        print chrom, pos, ref, ref_count, alt, alt_count
    }
' "$TMP_PILEUP" \
| "$SORT_BIN" -k1,1 -k2,2n > "$TMP_TXT"

move_if_success "$TMP_TXT" "$OUT_TXT"
require_nonempty_file "$OUT_TXT"

rm -f "$TMP_PILEUP"

log "Finished marker counting for ${BARCODE_ID}"
