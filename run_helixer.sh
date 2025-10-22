#!/bin/bash
# Author: Aaryan Bhatia
# Running Helixer to get genome annotations for Drosera scorpioides in a memory-efficient way

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# ============================================================================
# Configuration
# ============================================================================
wd=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/Helixer
genome=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/Files/genome/out_JBAT.FINAL.fa
split_dir="split_genome"
results_dir="annotations"
helixer_sif=/netscratch/dep_mercier/grp_marques/Aaryan/methods/helixer_docker/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif

# State tracking files
STATE_DIR="$wd/.helixer_state"
COMPLETED_FILE="$STATE_DIR/completed.txt"
FAILED_FILE="$STATE_DIR/failed.txt"
LOCK_FILE="$STATE_DIR/running.lock"

# ============================================================================
# Functions
# ============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error_exit() {
    log "ERROR: $*"
    exit 1
}

cleanup_on_exit() {
    if [ -f "$LOCK_FILE" ]; then
        rm -f "$LOCK_FILE"
    fi
}

trap cleanup_on_exit EXIT INT TERM

mark_completed() {
    local chr_name="$1"
    echo "$chr_name" >> "$COMPLETED_FILE"
    log "✓ Completed: $chr_name"
}

mark_failed() {
    local chr_name="$1"
    echo "$chr_name" >> "$FAILED_FILE"
    log "✗ Failed: $chr_name"
}

is_completed() {
    local chr_name="$1"
    [ -f "$COMPLETED_FILE" ] && grep -Fxq "$chr_name" "$COMPLETED_FILE"
}

is_failed() {
    local chr_name="$1"
    [ -f "$FAILED_FILE" ] && grep -Fxq "$chr_name" "$FAILED_FILE"
}

check_output_valid() {
    local gff_file="$1"
    
    # Check if file exists and is not empty
    [ -f "$gff_file" ] || return 1
    [ -s "$gff_file" ] || return 1
    
    # Check if it has at least a GFF header or some content
    grep -q "^##gff-version" "$gff_file" || head -1 "$gff_file" | grep -q "^#" || return 1
    
    return 0
}

run_helixer() {
    local fasta_path="$1"
    local output_path="$2"
    local chr_name="$3"
    
    # Skip if already completed successfully
    if is_completed "$chr_name" && check_output_valid "$output_path"; then
        log "Skipping $chr_name (already completed)"
        return 0
    fi
    
    log "Processing: $chr_name"
    
    # Remove partial output if exists
    [ -f "$output_path" ] && rm -f "$output_path"
    
    # Run Helixer with error checking
    if singularity run --nv \
        -B "$PWD:$PWD" \
        "$helixer_sif" Helixer.py \
        --fasta-path "$fasta_path" \
        --lineage land_plant \
        --gff-output-path "$output_path"; then
        
        # Verify output is valid
        if check_output_valid "$output_path"; then
            mark_completed "$chr_name"
            return 0
        else
            log "WARNING: Output file appears invalid: $output_path"
            mark_failed "$chr_name"
            return 1
        fi
    else
        log "ERROR: Helixer failed for $chr_name (exit code: $?)"
        mark_failed "$chr_name"
        return 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

log "Starting Helixer annotation pipeline"

# Load necessary modules
module load singularity

# Create necessary directories FIRST
mkdir -p "$wd" && cd "$wd"
mkdir -p "$split_dir" "$results_dir" "$STATE_DIR"

# Check if another instance is running
if [ -f "$LOCK_FILE" ]; then
    error_exit "Another instance appears to be running (lock file exists: $LOCK_FILE)"
fi
touch "$LOCK_FILE"

# Verify input files
[ -f "$genome" ] || error_exit "Genome file not found: $genome"
[ -f "$helixer_sif" ] || error_exit "Helixer SIF file not found: $helixer_sif"

# Check disk space (warn if less than 50GB free)
available_space=$(df -BG "$wd" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$available_space" -lt 50 ]; then
    log "WARNING: Low disk space available: ${available_space}GB"
fi

# Split genome if not already done
if [ ! -d "$split_dir" ] || [ -z "$(ls -A "$split_dir" 2>/dev/null)" ]; then
    log "Splitting genome into chromosome and debris files..."
    
    declare -i count=0
    awk '/^>/ {
        count++
        if (count <= 16) {
            file=sprintf("%s/chr_%d.fa", "'$split_dir'", count)
        } else {
            file="'$split_dir'/debris.fa"
        }
    } 
    {print > file}' "$genome"
    
    log "Genome successfully split into $(ls -1 "$split_dir" | wc -l) files"
else
    log "Using existing split genome files"
fi

# Show status if resuming
if [ -f "$COMPLETED_FILE" ]; then
    completed_count=$(wc -l < "$COMPLETED_FILE")
    log "Resuming: $completed_count chromosome(s) already completed"
fi

if [ -f "$FAILED_FILE" ]; then
    failed_count=$(wc -l < "$FAILED_FILE")
    log "WARNING: $failed_count chromosome(s) previously failed. Will retry."
fi

# Run Helixer on each chromosome
log "Starting Helixer annotation..."

success_count=0
fail_count=0

for chr_fasta in "$split_dir"/chr_*.fa; do
    [ -f "$chr_fasta" ] || continue
    
    chr_name=$(basename "$chr_fasta" .fa)
    output_file="$results_dir/${chr_name}_helixer.gff3"
    
    if run_helixer "$chr_fasta" "$output_file" "$chr_name"; then
        ((success_count++))
    else
        ((fail_count++))
        log "Continuing with remaining chromosomes..."
    fi
done

# Process debris sequences
debris_fasta="$split_dir/debris.fa"
if [ -f "$debris_fasta" ]; then
    debris_output="$results_dir/debris_helixer.gff3"
    if run_helixer "$debris_fasta" "$debris_output" "debris"; then
        ((success_count++))
    else
        ((fail_count++))
    fi
fi

# Summary
log "Processing complete: $success_count succeeded, $fail_count failed"

if [ $fail_count -gt 0 ]; then
    log "WARNING: Some chromosomes failed. Check $FAILED_FILE"
    log "Failed chromosomes:"
    cat "$FAILED_FILE" | sed 's/^/  - /'
fi

# Merge GFF3 files only if we have successful outputs
if [ $success_count -gt 0 ]; then
    log "Merging GFF3 files..."
    
    merged_output="$results_dir/D_scorpioides_helixer_merged.gff3"
    echo "##gff-version 3" > "$merged_output"
    
    # Merge only completed files
    for chr_name in $(cat "$COMPLETED_FILE" 2>/dev/null | sort -V); do
        if [ "$chr_name" = "debris" ]; then
            gff_file="$results_dir/debris_helixer.gff3"
        else
            gff_file="$results_dir/${chr_name}_helixer.gff3"
        fi
        
        if [ -f "$gff_file" ]; then
            grep -v "^##gff-version" "$gff_file" >> "$merged_output"
        fi
    done
    
    log "Merged annotation available: $merged_output"
else
    log "ERROR: No successful annotations to merge"
    exit 1
fi

log "Pipeline complete!"

if [ $fail_count -eq 0 ]; then
    log "All chromosomes processed successfully!"
else
    log "Completed with $fail_count failures. Rerun script to retry failed chromosomes."
    exit 1
fi
