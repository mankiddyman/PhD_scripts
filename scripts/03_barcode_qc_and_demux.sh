#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "03_barcode_qc_and_demux.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

stage_start "03_barcode_qc_and_demux"

safe_mkdir "$BARCODE_QC_DIR"
safe_mkdir "$DEMUX_DIR"
safe_mkdir "$FAILURES_DIR"

INPUT_BAM="${CELLRANGER_DIR}/outs/possorted_genome_bam.bam"
DEMUX_HEADER_SAM="${DEMUX_DIR}/_header.sam"
PARALLEL_JOBLOG="${BARCODE_QC_DIR}/parallel_demux.joblog"

require_nonempty_file "$INPUT_BAM"
require_executable_file "$SAMTOOLS_BIN"
require_executable_file "$PARALLEL_BIN"
require_executable_file "$AWK_BIN"
require_executable_file "$SORT_BIN"
require_executable_file "$WC_BIN"

###############################################################################
# step 1: make mapped-only BAM
###############################################################################

if ! file_exists_and_nonempty "$MAPPED_BAM"; then
    log "Creating mapped-only BAM"

    TMP_BAM="$(make_temp_path "$MAPPED_BAM")"

    run_cmd "$SAMTOOLS_BIN" view \
        -h -b -F 4 \
        -@ "$SAMTOOLS_THREADS" \
        "$INPUT_BAM" \
        -o "$TMP_BAM"

    move_if_success "$TMP_BAM" "$MAPPED_BAM"
    run_cmd "$SAMTOOLS_BIN" index "$MAPPED_BAM"

    assert_nonempty_file "$MAPPED_BAM"
    assert_nonempty_file "$MAPPED_BAI"
else
    log "Mapped-only BAM already exists, skipping"
fi

###############################################################################
# step 2: barcode read counts
###############################################################################

if ! file_exists_and_nonempty "$BARCODE_COUNTS_TSV"; then
    log "Counting reads per barcode"

    TMP_COUNTS="$(make_temp_path "$BARCODE_COUNTS_TSV")"

    {
        printf "barcode\ttotal_reads\n"
        "$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" "$INPUT_BAM" \
            | "$AWK_BIN" '
                {
                    for (i = 12; i <= NF; i++) {
                        if ($i ~ /^CB:Z:/) {
                            bc = $i
                            sub(/^CB:Z:/, "", bc)
                            counts[bc]++
                            break
                        }
                    }
                }
                END {
                    for (bc in counts) {
                        printf "%s\t%d\n", bc, counts[bc]
                    }
                }
            ' \
            | "$SORT_BIN" -k1,1
    } > "$TMP_COUNTS"

    move_if_success "$TMP_COUNTS" "$BARCODE_COUNTS_TSV"
    assert_nonempty_file "$BARCODE_COUNTS_TSV"
else
    log "Barcode counts TSV already exists, skipping"
fi

###############################################################################
# step 3: build barcode QC table
###############################################################################

if ! file_exists_and_nonempty "$BARCODE_READ_QC_TSV"; then
    log "Building barcode QC table"

    TMP_QC="$(make_temp_path "$BARCODE_READ_QC_TSV")"

    "$AWK_BIN" -F '\t' -v OFS='\t' '
        NR == 1 { print "barcode", "total_reads"; next }
        { print $1, $2 }
    ' "$BARCODE_COUNTS_TSV" > "$TMP_QC"

    move_if_success "$TMP_QC" "$BARCODE_READ_QC_TSV"
    assert_nonempty_file "$BARCODE_READ_QC_TSV"
else
    log "Barcode QC TSV already exists, skipping"
fi

###############################################################################
# step 4: select barcodes by total read count
###############################################################################

if [[ "$TEST_MODE" == "true" ]]; then
    BARCODE_MIN_READS="$MIN_TOTAL_READS_TEST"
else
    BARCODE_MIN_READS="$MIN_TOTAL_READS"
fi

if ! file_exists_and_nonempty "$SELECTED_BARCODES_TXT"; then
    log "Selecting barcodes with total_reads >= ${BARCODE_MIN_READS}"

    TMP_SELECTED="$(make_temp_path "$SELECTED_BARCODES_TXT")"

    "$AWK_BIN" -F '\t' -v min_reads="$BARCODE_MIN_READS" '
        NR == 1 { next }
        $2 >= min_reads { print $1 }
    ' "$BARCODE_READ_QC_TSV" > "$TMP_SELECTED"

    selected_count="$("$WC_BIN" -l < "$TMP_SELECTED")"

    if [[ "$selected_count" -eq 0 ]]; then
        warn "No barcodes passed total_reads threshold (${BARCODE_MIN_READS})"
        warn "Barcode QC completed, but there is nothing to demultiplex."
        rm -f "$TMP_SELECTED"
        stage_end "03_barcode_qc_and_demux"
        exit 0
    fi

    move_if_success "$TMP_SELECTED" "$SELECTED_BARCODES_TXT"
    assert_nonempty_file "$SELECTED_BARCODES_TXT"

    log "Selected ${selected_count} barcodes"
else
    log "Selected barcodes file already exists, skipping"
fi

###############################################################################
# step 5A: single-pass split into per-barcode SAM bodies
###############################################################################

if ! file_exists_and_nonempty "$DEMUX_HEADER_SAM"; then
    log "Extracting SAM header for demultiplexed outputs"
    "$SAMTOOLS_BIN" view -H "$MAPPED_BAM" > "$DEMUX_HEADER_SAM"
    assert_nonempty_file "$DEMUX_HEADER_SAM"
else
    log "Demux header SAM already exists, skipping"
fi

if [[ ! -f "$DEMUX_SUMMARY_TSV" ]]; then
    log "Preparing barcode directories"
    while read -r barcode; do
        mkdir -p "${DEMUX_DIR}/${barcode}"
    done < "$SELECTED_BARCODES_TXT"

    log "Single-pass split of mapped BAM into per-barcode SAM bodies"

    total_mapped_alignments="$("$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" -c "$MAPPED_BAM")"
    log "Total mapped alignments to scan: ${total_mapped_alignments}"

    "$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" "$MAPPED_BAM" \
        | "$AWK_BIN" \
            -v selected_file="$SELECTED_BARCODES_TXT" \
            -v demux_dir="$DEMUX_DIR" \
            -v progress_every="$SPLIT_PROGRESS_EVERY" \
            -v total="$total_mapped_alignments" '
            BEGIN {
                while ((getline bc < selected_file) > 0) {
                    selected[bc] = 1
                }
                close(selected_file)
            }
            {
                processed++
                barcode = ""

                for (i = 12; i <= NF; i++) {
                    if ($i ~ /^CB:Z:/) {
                        barcode = $i
                        sub(/^CB:Z:/, "", barcode)
                        break
                    }
                }

                if (barcode != "" && (barcode in selected)) {
                    outfile = demux_dir "/" barcode "/" barcode ".body.sam"
                    print $0 >> outfile
                    close(outfile)
                }

                if (processed % progress_every == 0) {
                    pct = (total > 0 ? (processed / total) * 100 : 0)
                    printf "[split] processed %d / %d alignments (%.2f%%)\n", processed, total, pct > "/dev/stderr"
                }
            }
            END {
                pct = (total > 0 ? (processed / total) * 100 : 100)
                printf "[split] processed %d / %d alignments (%.2f%%)\n", processed, total, pct > "/dev/stderr"
            }
        '
else
    log "Demux summary already exists; skipping split/sort work"
fi

###############################################################################
# helper for sorting/indexing one barcode
###############################################################################

sort_and_index_one_barcode() {
    local barcode="$1"
    local outdir="${DEMUX_DIR}/${barcode}"
    local body_sam="${outdir}/${barcode}.body.sam"
    local outbam="${outdir}/${barcode}.sorted.bam"
    local tmpbam="${outbam}.tmp"
    local statsfile="${outdir}/stats.tsv"

    if [[ -s "$outbam" && -s "$statsfile" ]]; then
        return 0
    fi

    if [[ ! -s "$body_sam" ]]; then
        printf '[%s]\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$barcode" "missing_or_empty_body_sam" >> "$DEMUX_FAILURES_TSV"
        return 0
    fi

    cat "$DEMUX_HEADER_SAM" "$body_sam" \
        | "$SAMTOOLS_BIN" sort -@ "$DEMUX_SORT_THREADS" -o "$tmpbam" -

    if [[ ! -s "$tmpbam" ]]; then
        printf '[%s]\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$barcode" "empty_sorted_bam" >> "$DEMUX_FAILURES_TSV"
        rm -f "$tmpbam"
        return 0
    fi

    mv "$tmpbam" "$outbam"
    "$SAMTOOLS_BIN" index "$outbam"

    local total_reads
    total_reads="$("$AWK_BIN" -F '\t' -v bc="$barcode" 'NR > 1 && $1 == bc { print $2 }' "$BARCODE_READ_QC_TSV")"

    local uniq_map_reads
    uniq_map_reads="$("$SAMTOOLS_BIN" view -@ "$DEMUX_SORT_THREADS" -c -q 4 "$outbam")"

    {
        printf "barcode\ttotal_reads\tuniq_map_reads\tuniq_ratio\n"
        "$AWK_BIN" -v bc="$barcode" -v total="$total_reads" -v uniq="$uniq_map_reads" 'BEGIN {
            ratio = (total > 0 ? uniq / total : 0)
            printf "%s\t%d\t%d\t%.6f\n", bc, total, uniq, ratio
        }'
    } > "$statsfile"

    rm -f "$body_sam"
    return 0
}

export SAMTOOLS_BIN
export AWK_BIN
export WC_BIN
export DEMUX_SORT_THREADS
export DEMUX_HEADER_SAM
export BARCODE_READ_QC_TSV
export DEMUX_DIR
export DEMUX_FAILURES_TSV
export -f sort_and_index_one_barcode

###############################################################################
# step 5B: sort/index each barcode BAM in parallel
###############################################################################

if [[ ! -f "$DEMUX_SUMMARY_TSV" ]]; then
    log "Sorting and indexing per-barcode BAMs with GNU Parallel"
    : > "$DEMUX_FAILURES_TSV"

    "$PARALLEL_BIN" \
        --jobs "$PARALLEL_JOBS" \
        --joblog "$PARALLEL_JOBLOG" \
        --bar \
        --eta \
        sort_and_index_one_barcode :::: "$SELECTED_BARCODES_TXT"
fi

###############################################################################
# step 6: merge per-barcode stats into summary
###############################################################################

if ! file_exists_and_nonempty "$DEMUX_SUMMARY_TSV"; then
    log "Merging per-barcode demux stats"

    TMP_SUMMARY="$(make_temp_path "$DEMUX_SUMMARY_TSV")"

    {
        printf "barcode\ttotal_reads\tuniq_map_reads\tuniq_ratio\n"
        find "$DEMUX_DIR" -mindepth 2 -maxdepth 2 -type f -name "stats.tsv" \
            | "$SORT_BIN" \
            | while read -r statsfile; do
                "$AWK_BIN" 'NR > 1' "$statsfile"
            done
    } > "$TMP_SUMMARY"

    move_if_success "$TMP_SUMMARY" "$DEMUX_SUMMARY_TSV"
    assert_nonempty_file "$DEMUX_SUMMARY_TSV"
else
    log "Demux summary already exists, skipping"
fi

log "Demux summary written to $DEMUX_SUMMARY_TSV"
log "Failed barcodes recorded in $DEMUX_FAILURES_TSV"

stage_end "03_barcode_qc_and_demux"
