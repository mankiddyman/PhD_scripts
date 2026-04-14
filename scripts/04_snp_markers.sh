#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "04_snp_markers.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

stage_start "04_snp_markers"

safe_mkdir "$SNP_DIR"
safe_mkdir "$MARKER_DIR"
safe_mkdir "$FAILURES_DIR"

ELIGIBLE_BARCODES_TXT="${MARKER_DIR}/eligible_barcodes.txt"
MARKER_FAILURES_TSV="${FAILURES_DIR}/marker_count_failed_barcodes.tsv"
MARKER_JOBLOG="${MARKER_DIR}/parallel_marker_count.joblog"

require_nonempty_file "$DEMUX_SUMMARY_TSV"
require_nonempty_file "$MARKER_FILE"
require_nonempty_file "$HAP1_FASTA"
require_nonempty_file "$SC_MARKER_COUNTS_SCRIPT"

require_executable_file "$SAMTOOLS_BIN"
require_executable_file "$PARALLEL_BIN"
require_executable_file "$AWK_BIN"
require_executable_file "$SORT_BIN"
require_executable_file "$WC_BIN"

###############################################################################
# step 1: select barcodes eligible for marker counting
###############################################################################

if ! file_exists_and_nonempty "$ELIGIBLE_BARCODES_TXT"; then
    log "Selecting barcodes with uniq_map_reads >= ${MIN_UNIQ_MAP_READS} and uniq_ratio >= ${MIN_UNIQ_RATIO}"

    TMP_ELIGIBLE="$(make_temp_path "$ELIGIBLE_BARCODES_TXT")"

    "$AWK_BIN" -F '\t' -v min_reads="$MIN_UNIQ_MAP_READS" -v min_ratio="$MIN_UNIQ_RATIO" '
        NR == 1 { next }
        ($3 >= min_reads) && ($4 >= min_ratio) { print $1 }
    ' "$DEMUX_SUMMARY_TSV" > "$TMP_ELIGIBLE"

    eligible_count="$("$WC_BIN" -l < "$TMP_ELIGIBLE")"

    if [[ "$eligible_count" -eq 0 ]]; then
        warn "No barcodes passed stage-04 mapping filters"
        warn "Stopping stage 04 with no eligible barcodes"
        rm -f "$TMP_ELIGIBLE"
        stage_end "04_snp_markers"
        exit 0
    fi

    move_if_success "$TMP_ELIGIBLE" "$ELIGIBLE_BARCODES_TXT"
    assert_nonempty_file "$ELIGIBLE_BARCODES_TXT"

    log "Selected ${eligible_count} eligible barcodes for marker counting"
else
    log "Eligible barcode list already exists, skipping selection"
fi

###############################################################################
# helper: run marker counting for one barcode
###############################################################################

run_marker_count_one_barcode() {
    local barcode="$1"
    local bam="results/04_demultiplex/${barcode}/${barcode}.sorted.bam"
    local outdir="results/06_markers/${barcode}"

    if [[ ! -s "$bam" ]]; then
        printf '[%s]\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$barcode" "missing_demux_bam" >> "$MARKER_FAILURES_TSV"
        return 1
    fi

    bash "$SC_MARKER_COUNTS_SCRIPT" \
        "$barcode" \
        "$HAP1_FASTA" \
        "$bam" \
        "$MARKER_FILE" \
        "$outdir"
}

export SC_MARKER_COUNTS_SCRIPT
export HAP1_FASTA
export MARKER_FILE
export MARKER_FAILURES_TSV
export -f run_marker_count_one_barcode

###############################################################################
# step 2: run marker counting in parallel
###############################################################################

log "Running marker counting across eligible barcodes"
: > "$MARKER_FAILURES_TSV"

"$PARALLEL_BIN" \
    --jobs "$MARKER_COUNT_JOBS" \
    --joblog "$MARKER_JOBLOG" \
    --eta \
    run_marker_count_one_barcode :::: "$ELIGIBLE_BARCODES_TXT"

###############################################################################
# step 3: build switches.stats.tsv
###############################################################################

log "Computing marker counts and switch counts"

TMP_SWITCHES="$(make_temp_path "$SWITCH_STATS_TSV")"

{
    printf "barcode\tmarker_num\tswitch_num\tswitch_rate\n"

    while read -r barcode; do
        input_txt="${MARKER_DIR}/${barcode}/input_corrected.txt"

        if [[ ! -s "$input_txt" ]]; then
            warn "Missing input_corrected.txt for ${barcode}, skipping switches summary"
            continue
        fi

        "$AWK_BIN" -F '\t' -v bc="$barcode" '
            ($4 > 0 || $6 > 0) {
                total = $4 + $6
                if (total == 0) next

                marker_num++

                frac = $4 / total
                if (frac < 0.2) {
                    geno = 0
                } else if (frac > 0.8) {
                    geno = 1
                } else {
                    geno = 0.5
                }

                if (seen) {
                    if (geno != prev_geno) {
                        switch_num++
                    }
                }

                prev_geno = geno
                seen = 1
            }
            END {
                if (marker_num == 0) {
                    switch_num = 0
                    switch_rate = 0
                } else {
                    if (switch_num == "") switch_num = 0
                    switch_rate = switch_num / marker_num
                }

                printf "%s\t%d\t%d\t%.6f\n", bc, marker_num, switch_num, switch_rate
            }
        ' "$input_txt"

    done < "$ELIGIBLE_BARCODES_TXT"
} > "$TMP_SWITCHES"

move_if_success "$TMP_SWITCHES" "$SWITCH_STATS_TSV"
assert_nonempty_file "$SWITCH_STATS_TSV"

###############################################################################
# step 4: select barcodes for CO calling
###############################################################################

log "Selecting cells for CO calling"

TMP_SELECTED_CO="$(make_temp_path "$SELECTED_FOR_CO_TXT")"

"$AWK_BIN" -F '\t' -v min_markers="$MIN_MARKERS" -v max_switch="$MAX_SWITCH_RATE" '
    NR == 1 { next }
    ($2 >= min_markers) && ($4 <= max_switch) { print $1 }
' "$SWITCH_STATS_TSV" > "$TMP_SELECTED_CO"

move_if_success "$TMP_SELECTED_CO" "$SELECTED_FOR_CO_TXT"

selected_co_count="$("$WC_BIN" -l < "$SELECTED_FOR_CO_TXT")"
log "Selected ${selected_co_count} cells for CO calling"

log "switches.stats.tsv written to $SWITCH_STATS_TSV"
log "selected_for_co.txt written to $SELECTED_FOR_CO_TXT"
log "Marker counting failures recorded in $MARKER_FAILURES_TSV"

stage_end "04_snp_markers"
