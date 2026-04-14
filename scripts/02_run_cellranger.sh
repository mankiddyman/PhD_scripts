#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "02_run_cellranger.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

stage_start "02_run_cellranger"

safe_mkdir "$(dirname "$CELLRANGER_DIR")"
safe_mkdir "$TEST_FASTQ_DIR"

###############################################################################
# helpers
###############################################################################

subset_fastq_gz_by_read_count() {
    local input_fastq="$1"
    local output_fastq="$2"
    local n_reads="$3"

    local tmp_fastq
    tmp_fastq="$(make_temp_path "$output_fastq")"

    gzip -cd "$input_fastq" \
        | awk -v n="$n_reads" 'NR <= (n * 4) {print}' \
        | gzip > "$tmp_fastq"

    move_if_success "$tmp_fastq" "$output_fastq"
}

make_test_fastqs() {
    local source_dir="$1"
    local dest_dir="$2"
    local sample_name="$3"
    local n_reads="$4"

    safe_mkdir "$dest_dir"

    log "Creating test FASTQs with ${n_reads} reads per FASTQ"

    local fastq
    while IFS= read -r fastq; do
        local out_fastq
        out_fastq="${dest_dir}/$(basename "$fastq")"

        if file_exists_and_nonempty "$out_fastq" && [[ "${REBUILD_TEST_FASTQS}" != "true" ]]; then
            log "Test FASTQ already exists, skipping: ${out_fastq}"
            continue
        fi

        log "Creating test FASTQ: ${out_fastq}"
        subset_fastq_gz_by_read_count "$fastq" "$out_fastq" "$n_reads"
    done < <(find "$source_dir" -maxdepth 1 -type f -name "${sample_name}*.fastq.gz" | sort)

    local test_fastq_count
    test_fastq_count="$(find "$dest_dir" -maxdepth 1 -type f -name "${sample_name}*.fastq.gz" | wc -l)"
    [[ "$test_fastq_count" -gt 0 ]] || die "No test FASTQs were created in ${dest_dir}"
}

###############################################################################
# choose FASTQ source
###############################################################################

CELLRANGER_FASTQ_DIR="$FASTQ_DIR"

if [[ "$TEST_MODE" == "true" ]]; then
    log "TEST_MODE is true; creating/using subset FASTQs"
    make_test_fastqs "$FASTQ_DIR" "$TEST_FASTQ_DIR" "$SAMPLE_NAME" "$TEST_READ_COUNT"
    CELLRANGER_FASTQ_DIR="$TEST_FASTQ_DIR"
else
    log "TEST_MODE is false; using full FASTQ directory"
fi

print_key_path "CELLRANGER_FASTQ_DIR" "$CELLRANGER_FASTQ_DIR"
print_key_path "CELLRANGER_REF_DIR" "$CELLRANGER_REF_DIR"
print_key_path "CELLRANGER_DIR" "$CELLRANGER_DIR"

###############################################################################
# restart-safe skip
###############################################################################

if file_exists_and_nonempty "${CELLRANGER_DIR}/outs/possorted_genome_bam.bam"; then
    log "Cell Ranger BAM already exists, skipping cellranger count"
    assert_nonempty_file "${CELLRANGER_DIR}/outs/possorted_genome_bam.bam"
    stage_end "02_run_cellranger"
    exit 0
fi

###############################################################################
# run cellranger count
###############################################################################

require_dir "$CELLRANGER_REF_DIR"
require_dir "$CELLRANGER_FASTQ_DIR"

cleanup_temp_if_exists "$CELLRANGER_DIR"

(
    cd "$(dirname "$CELLRANGER_DIR")"
    run_in_env cellranger count \
        --id="$RUN_NAME" \
        --create-bam=true \
        --transcriptome="$CELLRANGER_REF_DIR" \
        --fastqs="$CELLRANGER_FASTQ_DIR" \
        --sample="$SAMPLE_NAME" \
        --include-introns="$INCLUDE_INTRONS" \
        --expect-cells="$EXPECT_CELLS" \
        --chemistry="$CHEMISTRY" \
        --localcores="$CELLRANGER_CORES" \
        --localmem="$CELLRANGER_MEMGB"
)

###############################################################################
# final checks
###############################################################################

assert_nonempty_dir "$CELLRANGER_DIR"
assert_nonempty_file "${CELLRANGER_DIR}/outs/possorted_genome_bam.bam"

log "Cell Ranger count complete"
stage_end "02_run_cellranger"
