#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "05_co_and_plots.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

stage_start "05_co_and_plots"

safe_mkdir "$CO_DIR"
safe_mkdir "${CO_DIR}/summary"
safe_mkdir "$PLOT_DIR"
safe_mkdir "$FAILURES_DIR"

require_nonempty_file "$SELECTED_FOR_CO_TXT"
require_nonempty_file "$HAPCO_SCRIPT"
require_nonempty_file "${HAP1_FASTA}.fai"
require_executable_file "$PARALLEL_BIN"

CO_FAILURES_TSV="${FAILURES_DIR}/co_failed_barcodes.tsv"
CO_JOBLOG="${CO_DIR}/parallel_co.joblog"
CO_COUNT_LIST="${CO_DIR}/summary/CO_num.list"

###############################################################################
# helper: run CO calling for one barcode
###############################################################################

run_co_one_barcode() {
    local barcode="$1"
    local prefix="${barcode%-1}"
    local input_file="${MARKER_DIR}/${barcode}/input_corrected.txt"
    local summary_file="${CO_DIR}/${prefix}_co_summary.tsv"
    local co_file="${CO_DIR}/${prefix}_allele_cnts_at_markers_sorted_co_pred.txt"

    if [[ -s "$summary_file" && -s "$co_file" ]]; then
        return 0
    fi

    if [[ ! -s "$input_file" ]]; then
        printf '[%s]\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$barcode" "missing_input_corrected" >> "$CO_FAILURES_TSV"
        return 0
    fi

    Rscript "$HAPCO_SCRIPT" \
        -i "$input_file" \
        -p "$prefix" \
        -g "${HAP1_FASTA}.fai" \
        -o "$CO_DIR" \
        -c "$HAPCO_MIN_MARKERS" \
        -s "$HAPCO_SEGMENT_SIZE" \
        -n "$HAPCO_N"
}

export PROJECT_ROOT
export HAP1_FASTA
export HAPCO_SCRIPT
export CO_DIR
export MARKER_DIR
export CO_FAILURES_TSV
export HAPCO_MIN_MARKERS
export HAPCO_SEGMENT_SIZE
export HAPCO_N
export -f run_co_one_barcode

###############################################################################
# step 1: run CO calling across selected cells
###############################################################################

log "Running hapCO identification across selected cells"
: > "$CO_FAILURES_TSV"

"$PARALLEL_BIN" \
    --jobs "$PARALLEL_JOBS" \
    --joblog "$CO_JOBLOG" \
    --eta \
    run_co_one_barcode :::: "$SELECTED_FOR_CO_TXT"

###############################################################################
# step 2: merge per-cell summaries
###############################################################################

log "Merging per-cell CO summaries"

TMP_CO_SUMMARY="$(make_temp_path "$CO_SUMMARY_TSV")"

summary_files_found=$(find "$CO_DIR" -maxdepth 1 -type f -name "*_co_summary.tsv" | wc -l)

if [[ "$summary_files_found" -eq 0 ]]; then
    die "No per-cell CO summary files found in ${CO_DIR}"
fi

{
    header_written="false"
    find "$CO_DIR" -maxdepth 1 -type f -name "*_co_summary.tsv" | sort | while read -r f; do
        if [[ "$header_written" == "false" ]]; then
            cat "$f"
            header_written="true"
        else
            awk 'NR > 1' "$f"
        fi
    done
} > "$TMP_CO_SUMMARY"

move_if_success "$TMP_CO_SUMMARY" "$CO_SUMMARY_TSV"
assert_nonempty_file "$CO_SUMMARY_TSV"

###############################################################################
# step 3: build global CO intervals BED
###############################################################################

log "Building global CO intervals BED"

TMP_CO_BED="$(make_temp_path "$CO_INTERVALS_BED")"

{
    find "$CO_DIR" -maxdepth 1 -type f -name "*_allele_cnts_at_markers_sorted_co_pred.txt" | sort | while read -r f; do
        prefix="$(basename "$f" _allele_cnts_at_markers_sorted_co_pred.txt)"
        awk -F '\t' -v OFS='\t' -v sample="$prefix" '
            NR == 1 { next }
            NF >= 3 { print $1, $2, $3, sample }
        ' "$f"
    done
} > "$TMP_CO_BED"

move_if_success "$TMP_CO_BED" "$CO_INTERVALS_BED"

###############################################################################
# step 4: build simple CO count list
###############################################################################

log "Building CO count list"

TMP_CO_COUNT_LIST="$(make_temp_path "$CO_COUNT_LIST")"

find "$CO_DIR" -maxdepth 1 -type f -name "*_allele_cnts_at_markers_sorted_co_pred.txt" | sort | while read -r f; do
    awk 'NR > 1' "$f" | wc -l
done > "$TMP_CO_COUNT_LIST"

move_if_success "$TMP_CO_COUNT_LIST" "$CO_COUNT_LIST"

###############################################################################
# step 5: report
###############################################################################

selected_n=$(wc -l < "$SELECTED_FOR_CO_TXT")
summary_n=$(awk 'NR > 1' "$CO_SUMMARY_TSV" | wc -l)

log "Selected-for-CO cells: ${selected_n}"
log "Per-cell CO summaries merged: ${summary_n}"
log "CO summary written to $CO_SUMMARY_TSV"
log "CO intervals written to $CO_INTERVALS_BED"
log "CO count list written to $CO_COUNT_LIST"
log "CO failures recorded in $CO_FAILURES_TSV"

stage_end "05_co_and_plots"
