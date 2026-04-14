#!/usr/bin/env bash
#
# run_pipeline.sh
#
# Main launcher for the CO pipeline.
#
# Responsibilities:
# - load config and shared helpers
# - run preflight validation
# - execute stages in order
# - skip stages with valid checkpoints unless force-rerun is requested
# - write stage-specific logs
#

set -euo pipefail


###############################################################################
# conda environment activation
###############################################################################

if [[ -n "${CONDA_ENV_NAME:-}" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda not found but CONDA_ENV_NAME is set" >&2
        exit 1
    fi

    # initialize conda for non-interactive shells
    eval "$(conda shell.bash hook)"

    echo "[INFO] Activating conda environment: ${CONDA_ENV_NAME}" >&2
    conda activate "${CONDA_ENV_NAME}" || {
        echo "ERROR: failed to activate conda env ${CONDA_ENV_NAME}" >&2
        exit 1
    }
fi

###############################################################################
# locate project root and source config/helpers
###############################################################################

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"
source "${LIB_DIR}/validate.sh"

###############################################################################
# startup directory preparation
###############################################################################

safe_mkdir "$LOG_DIR"
safe_mkdir "$CHECKPOINT_DIR"
safe_mkdir "$TMP_DIR"
safe_mkdir "$FAILURES_DIR"

###############################################################################
# stage runner
###############################################################################

run_stage() {
    local stage_name="$1"
    local stage_script="$2"
    local done_file="$3"
    local stage_stdout="${LOG_DIR}/${stage_name}.out"
    local stage_stderr="${LOG_DIR}/${stage_name}.err"

    require_nonempty_file "$stage_script"

    if is_done "$done_file" && ! should_force_rerun_stage "$stage_name"; then
        log "Skipping ${stage_name}: checkpoint exists at ${done_file}"
        return 0
    fi

    if should_force_rerun_stage "$stage_name"; then
        warn "Force-rerun requested for ${stage_name}"
        remove_done "$done_file"
    fi

    stage_start "$stage_name"
    log "Stage script: ${stage_script}"
    log "STDOUT log: ${stage_stdout}"
    log "STDERR log: ${stage_stderr}"

    if bash "$stage_script" >"$stage_stdout" 2>"$stage_stderr"; then
        mark_done "$done_file"
        stage_end "$stage_name"
    else
        local exit_code=$?
        error "Stage ${stage_name} failed with exit code ${exit_code}"
        error "See logs:"
        error "  STDOUT: ${stage_stdout}"
        error "  STDERR: ${stage_stderr}"
        exit "$exit_code"
    fi
}

###############################################################################
# main
###############################################################################

main() {
    stage_start "PIPELINE"

    print_key_path "PROJECT_ROOT" "$PROJECT_ROOT"
    print_key_path "RESULTS_DIR" "$RESULTS_DIR"
    print_key_path "LOG_DIR" "$LOG_DIR"
    print_key_path "CHECKPOINT_DIR" "$CHECKPOINT_DIR"
    print_key_path "TMP_DIR" "$TMP_DIR"

    log "Running preflight checks"
    run_preflight_checks

    run_stage "00_make_reference_markers" "${SCRIPTS_DIR}/00_make_reference_markers.sh" "$STAGE0_DONE"
    run_stage "01_build_reference"         "${SCRIPTS_DIR}/01_build_reference.sh"     "$STAGE1_DONE"
    run_stage "02_run_cellranger"          "${SCRIPTS_DIR}/02_run_cellranger.sh"      "$STAGE2_DONE"
    run_stage "03_barcode_qc_and_demux"    "${SCRIPTS_DIR}/03_barcode_qc_and_demux.sh" "$STAGE3_DONE"
    run_stage "04_snp_markers"             "${SCRIPTS_DIR}/04_snp_markers.sh"         "$STAGE4_DONE"
    run_stage "05_co_and_plots"            "${SCRIPTS_DIR}/05_co_and_plots.sh"        "$STAGE5_DONE"

    stage_end "PIPELINE"
    log "Pipeline completed successfully"
}

main "$@"
