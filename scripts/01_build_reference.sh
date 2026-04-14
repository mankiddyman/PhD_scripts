#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PROJECT_ROOT}/config.sh"
source "${LIB_DIR}/common.sh"

trap 'error "01_build_reference.sh failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

stage_start "01_build_reference"

safe_mkdir "$REFERENCE_DIR"

CELLRANGER_GENOME_TMP="${GENOME_NAME}_tmp"
CELLRANGER_REF_TMP="${REFERENCE_DIR}/${CELLRANGER_GENOME_TMP}"

###############################################################################
# helper: make a Cell Ranger-compatible GTF from Helixer GFF
###############################################################################

make_cellranger_gtf_from_gff() {
    local input_gff="$1"
    local output_gtf="$2"

    awk -F '\t' '
    BEGIN { OFS="\t" }

    FNR==NR {
        if ($0 ~ /^#/) next

        if ($3=="mRNA" || $3=="transcript") {
            id=""
            parent=""
            n=split($9, attrs, ";")
            for (i=1; i<=n; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", attrs[i])
                if (attrs[i] ~ /^ID=/) {
                    id=attrs[i]
                    sub(/^ID=/, "", id)
                } else if (attrs[i] ~ /^Parent=/) {
                    parent=attrs[i]
                    sub(/^Parent=/, "", parent)
                }
            }
            if (id != "" && parent != "") {
                tx2gene[id] = parent
            }
        }
        next
    }

    {
        if ($0 ~ /^#/) next

        if ($3=="mRNA" || $3=="transcript") {
            tx=""
            gene=""
            n=split($9, attrs, ";")
            for (i=1; i<=n; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", attrs[i])
                if (attrs[i] ~ /^ID=/) {
                    tx=attrs[i]
                    sub(/^ID=/, "", tx)
                } else if (attrs[i] ~ /^Parent=/) {
                    gene=attrs[i]
                    sub(/^Parent=/, "", gene)
                }
            }

            if (tx != "" && gene != "") {
                print $1, $2, "transcript", $4, $5, $6, $7, $8, \
                      "gene_id \"" gene "\"; transcript_id \"" tx "\";"
            }
        }
        else if ($3=="exon") {
            tx=""
            gene=""
            n=split($9, attrs, ";")
            for (i=1; i<=n; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", attrs[i])
                if (attrs[i] ~ /^Parent=/) {
                    tx=attrs[i]
                    sub(/^Parent=/, "", tx)
                }
            }

            if (tx != "") {
                gene = tx2gene[tx]
                if (gene == "") {
                    gene = tx
                    sub(/\.[0-9]+$/, "", gene)
                }

                print $1, $2, "exon", $4, $5, $6, $7, $8, \
                      "gene_id \"" gene "\"; transcript_id \"" tx "\";"
            }
        }
    }' "$input_gff" "$input_gff" > "$output_gtf"
}

###############################################################################
# step 1: build Cell Ranger-compatible hap1 GTF
###############################################################################

if ! file_exists_and_nonempty "$HAP1_GTF"; then
    log "Building Cell Ranger-compatible hap1 GTF from GFF"

    TMP_GTF="$(make_temp_path "$HAP1_GTF")"

    make_cellranger_gtf_from_gff "$HAP1_GFF3" "$TMP_GTF"

    awk -F '\t' '
        $3=="transcript" || $3=="exon" {
            if ($9 !~ /gene_id "/ || $9 !~ /transcript_id "/) {
                print "Bad GTF line:", NR, $0 > "/dev/stderr"
                bad=1
            }
        }
        END { exit bad }
    ' "$TMP_GTF" || die "Generated GTF failed attribute sanity check"

    move_if_success "$TMP_GTF" "$HAP1_GTF"
    assert_nonempty_file "$HAP1_GTF"
else
    log "HAP1_GTF already exists, skipping conversion"
fi

###############################################################################
# step 2: build Cell Ranger reference
###############################################################################

if ! dir_exists_and_nonempty "$CELLRANGER_REF_DIR"; then
    log "Building Cell Ranger reference"

    cleanup_temp_if_exists "$CELLRANGER_REF_TMP"
    cleanup_temp_if_exists "$CELLRANGER_REF_DIR"

    (
        cd "$REFERENCE_DIR"
        run_in_env cellranger mkref \
            --nthreads="$CELLRANGER_CORES" \
            --memgb="$CELLRANGER_MEMGB" \
            --genome="$CELLRANGER_GENOME_TMP" \
            --fasta="$HAP1_FASTA" \
            --genes="$HAP1_GTF"
    )

    [[ -d "$CELLRANGER_REF_TMP" ]] || die "Cell Ranger reference directory not found at expected path: ${CELLRANGER_REF_TMP}"

    mv "$CELLRANGER_REF_TMP" "$CELLRANGER_REF_DIR" \
        || die "Failed to rename ${CELLRANGER_REF_TMP} to ${CELLRANGER_REF_DIR}"
else
    log "Cell Ranger reference directory already exists, skipping"
fi

###############################################################################
# final checks
###############################################################################

assert_nonempty_file "$HAP1_GTF"
assert_nonempty_dir "$CELLRANGER_REF_DIR"

log "Reference build complete"
stage_end "01_build_reference"
