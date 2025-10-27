##!/bin/bash

# HiFi Read Decontamination Pipeline
# Filters new library reads based on k-mer overlap with trusted old library

set -euo pipefail

# Check arguments
new_hifi="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/new_hifi_reads.fastq.gz"
old_hifi="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/old_hifi.fastq.gz"
outdir="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/decontamination/kmer"
mkdir -p $outdir
mkdir -p $outdir/temp

OLD_LIB=$old_hifi
NEW_LIB=$new_hifi
OUTPUT_PREFIX=${outdir}/decontaminated_reads
TEMP_DIR=${outdir}/temp
KMER_SIZE=21
# Check dependencies
command -v kmc >/dev/null 2>&1 || { echo "Error: kmc not found in PATH"; exit 1; }
command -v kmc_tools >/dev/null 2>&1 || { echo "Error: kmc_tools not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "Error: python3 not found in PATH"; exit 1; }

# Check input files
[ -f "$OLD_LIB" ] || { echo "Error: Old library file not found: $OLD_LIB"; exit 1; }
[ -f "$NEW_LIB" ] || { echo "Error: New library file not found: $NEW_LIB"; exit 1; }

# Create temp directory
mkdir -p "$TEMP_DIR"

echo "=========================================="
echo "HiFi Read Decontamination Pipeline"
echo "=========================================="
echo "Old library: $OLD_LIB"
echo "New library: $NEW_LIB"
echo "K-mer size: $KMER_SIZE"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Temp directory: $TEMP_DIR"
echo ""

# Step 1: Build k-mer database from old library
echo "[Step 1/3] Building k-mer database from old library..."
STEP1_START=$(date +%s)
KMC_DB="${TEMP_DIR}/old_lib_kmers"

kmc -k${KMER_SIZE} -m64 -t8 -ci1 -fq "$OLD_LIB" "$KMC_DB" "$TEMP_DIR"

STEP1_END=$(date +%s)
STEP1_TIME=$((STEP1_END - STEP1_START))
echo "K-mer database built successfully (Time: ${STEP1_TIME}s)"
echo ""

# Step 2: Dump k-mers to text file for Python processing
echo "[Step 2/3] Dumping k-mer database to text file..."
STEP2_START=$(date +%s)
KMER_FILE="${TEMP_DIR}/old_lib_kmers.txt"

kmc_tools transform "$KMC_DB" dump "$KMER_FILE"

STEP2_END=$(date +%s)
STEP2_TIME=$((STEP2_END - STEP2_START))
echo "K-mer dump complete (Time: ${STEP2_TIME}s)"
echo ""

# Step 3: Filter new library using Python helper
echo "[Step 3/3] Filtering new library reads..."
STEP3_START=$(date +%s)
echo ""

# Create Python helper script
PYTHON_SCRIPT="${TEMP_DIR}/filter_reads.py"

cat > "$PYTHON_SCRIPT" << 'PYTHON_EOF'
#!/usr/bin/env python3

import gzip
import sys
from collections import defaultdict

def load_kmers(kmer_file):
    """Load k-mers from KMC dump file into a set"""
    print("Loading k-mers from old library...", file=sys.stderr)
    kmers = set()
    with open(kmer_file, 'r') as f:
        for line in f:
            kmer = line.strip().split()[0]  # KMC format: kmer\tcount
            kmers.add(kmer)
    print(f"Loaded {len(kmers):,} unique k-mers", file=sys.stderr)
    return kmers

def reverse_complement(seq):
    """Return reverse complement of sequence"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(base, 'N') for base in reversed(seq))

def extract_kmers(seq, k):
    """Extract all k-mers from a sequence (canonical form)"""
    kmers = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:  # Skip k-mers with N
            rc_kmer = reverse_complement(kmer)
            canonical = min(kmer, rc_kmer)  # KMC uses canonical k-mers
            kmers.append(canonical)
    return kmers

def process_fastq(input_fq, output_min1, output_all, trusted_kmers, k):
    """Filter reads based on k-mer overlap"""
    print("Processing reads...", file=sys.stderr)
    
    total_reads = 0
    min1_reads = 0
    all_reads = 0
    
    with gzip.open(input_fq, 'rt') as infile, \
         gzip.open(output_min1, 'wt') as out_min1, \
         gzip.open(output_all, 'wt') as out_all:
        
        while True:
            # Read 4 lines (one fastq record)
            header = infile.readline()
            if not header:
                break
            
            seq = infile.readline().strip()
            plus = infile.readline()
            qual = infile.readline().strip()
            
            total_reads += 1
            
            # Extract k-mers from read
            read_kmers = extract_kmers(seq, k)
            
            if len(read_kmers) == 0:
                # Skip reads with no valid k-mers
                continue
            
            # Count matching k-mers
            matching = sum(1 for kmer in read_kmers if kmer in trusted_kmers)
            
            # Strategy 1: At least 1 matching k-mer
            if matching >= 1:
                out_min1.write(header)
                out_min1.write(seq + '\n')
                out_min1.write(plus)
                out_min1.write(qual + '\n')
                min1_reads += 1
            
            # Strategy 2: ALL k-mers match
            if matching == len(read_kmers):
                out_all.write(header)
                out_all.write(seq + '\n')
                out_all.write(plus)
                out_all.write(qual + '\n')
                all_reads += 1
            
            # Progress update
            if total_reads % 10000 == 0:
                print(f"Processed {total_reads:,} reads... "
                      f"(min1: {min1_reads:,}, all: {all_reads:,})", 
                      file=sys.stderr)
    
    return total_reads, min1_reads, all_reads

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: filter_reads.py <kmer_file> <input_fq> <output_min1> <output_all> <k>")
        sys.exit(1)
    
    kmer_file = sys.argv[1]
    input_fq = sys.argv[2]
    output_min1 = sys.argv[3]
    output_all = sys.argv[4]
    k = int(sys.argv[5])
    
    # Load trusted k-mers
    trusted_kmers = load_kmers(kmer_file)
    print("", file=sys.stderr)
    
    # Process reads
    total, min1, all_kmers = process_fastq(input_fq, output_min1, output_all, trusted_kmers, k)
    
    print("", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
    print("FILTERING RESULTS", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
    print(f"Total reads processed: {total:,}", file=sys.stderr)
    print(f"", file=sys.stderr)
    print(f"Strategy 1 (≥1 matching k-mer):", file=sys.stderr)
    print(f"  Reads retained: {min1:,} ({100*min1/total:.2f}%)", file=sys.stderr)
    print(f"", file=sys.stderr)
    print(f"Strategy 2 (ALL k-mers match):", file=sys.stderr)
    print(f"  Reads retained: {all_kmers:,} ({100*all_kmers/total:.2f}%)", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
PYTHON_EOF

# Run Python filtering
python3 "$PYTHON_SCRIPT" \
    "$KMER_FILE" \
    "$NEW_LIB" \
    "${OUTPUT_PREFIX}_min1.fastq.gz" \
    "${OUTPUT_PREFIX}_all.fastq.gz" \
    "$KMER_SIZE"

STEP3_END=$(date +%s)
STEP3_TIME=$((STEP3_END - STEP3_START))
TOTAL_TIME=$((STEP1_TIME + STEP2_TIME + STEP3_TIME))

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo "Timing summary:"
echo "  Step 1 (KMC database): ${STEP1_TIME}s"
echo "  Step 2 (KMC dump):     ${STEP2_TIME}s"
echo "  Step 3 (Filtering):    ${STEP3_TIME}s"
echo "  Total time:            ${TOTAL_TIME}s"
echo ""
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_min1.fastq.gz  - Reads with ≥1 matching k-mer (less conservative)"
echo "  ${OUTPUT_PREFIX}_all.fastq.gz   - Reads where ALL k-mers match (super conservative)"
echo ""
echo "Temporary files in: $TEMP_DIR"
echo "You can remove temp directory when done: rm -rf $TEMP_DIR"
!/bin/bash

# HiFi Read Decontamination Pipeline
# Filters new library reads based on k-mer overlap with trusted old library

set -euo pipefail

# Check arguments
new_hifi="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/new_hifi_reads.fastq.gz"
old_hifi="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/old_hifi.fastq.gz"
outdir="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/decontamination/kmer"
mkdir -p $outdir
mkdir -p $outdir/temp

OLD_LIB=$old_hifi
NEW_LIB=$new_hifi
OUTPUT_PREFIX=${outdir}/decontaminated_reads
TEMP_DIR=${outdir}/temp
KMER_SIZE=21


source ~/.bashrc
micromamba activate python3.10
# Check dependencies
command -v kmc >/dev/null 2>&1 || { echo "Error: kmc not found in PATH"; exit 1; }
command -v kmc_tools >/dev/null 2>&1 || { echo "Error: kmc_tools not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "Error: python3 not found in PATH"; exit 1; }

# Check input files
[ -f "$OLD_LIB" ] || { echo "Error: Old library file not found: $OLD_LIB"; exit 1; }
[ -f "$NEW_LIB" ] || { echo "Error: New library file not found: $NEW_LIB"; exit 1; }

# Create temp directory
mkdir -p "$TEMP_DIR"

echo "=========================================="
echo "HiFi Read Decontamination Pipeline"
echo "=========================================="
echo "Old library: $OLD_LIB"
echo "New library: $NEW_LIB"
echo "K-mer size: $KMER_SIZE"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Temp directory: $TEMP_DIR"
echo ""

# Step 1: Build k-mer database from old library
echo "[Step 1/3] Building k-mer database from old library..."
STEP1_START=$(date +%s)
KMC_DB="${TEMP_DIR}/old_lib_kmers"

kmc -k${KMER_SIZE} -m64 -t8 -ci1 -fq "$OLD_LIB" "$KMC_DB" "$TEMP_DIR"

STEP1_END=$(date +%s)
STEP1_TIME=$((STEP1_END - STEP1_START))
echo "K-mer database built successfully (Time: ${STEP1_TIME}s)"
echo ""

# Step 2: Dump k-mers to text file for Python processing
echo "[Step 2/3] Dumping k-mer database to text file..."
STEP2_START=$(date +%s)
KMER_FILE="${TEMP_DIR}/old_lib_kmers.txt"

kmc_tools transform "$KMC_DB" dump "$KMER_FILE"

STEP2_END=$(date +%s)
STEP2_TIME=$((STEP2_END - STEP2_START))
echo "K-mer dump complete (Time: ${STEP2_TIME}s)"
echo ""

# Step 3: Filter new library using Python helper
echo "[Step 3/3] Filtering new library reads..."
STEP3_START=$(date +%s)
echo ""

# Create Python helper script
PYTHON_SCRIPT="${TEMP_DIR}/filter_reads.py"

cat > "$PYTHON_SCRIPT" << 'PYTHON_EOF'
#!/usr/bin/env python3

import gzip
import sys
from collections import defaultdict

def load_kmers(kmer_file):
    """Load k-mers from KMC dump file into a set"""
    print("Loading k-mers from old library...", file=sys.stderr)
    kmers = set()
    with open(kmer_file, 'r') as f:
        for line in f:
            kmer = line.strip().split()[0]  # KMC format: kmer\tcount
            kmers.add(kmer)
    print(f"Loaded {len(kmers):,} unique k-mers", file=sys.stderr)
    return kmers

def reverse_complement(seq):
    """Return reverse complement of sequence"""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(base, 'N') for base in reversed(seq))

def extract_kmers(seq, k):
    """Extract all k-mers from a sequence (canonical form)"""
    kmers = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' not in kmer:  # Skip k-mers with N
            rc_kmer = reverse_complement(kmer)
            canonical = min(kmer, rc_kmer)  # KMC uses canonical k-mers
            kmers.append(canonical)
    return kmers

def process_fastq(input_fq, output_min1, output_all, trusted_kmers, k):
    """Filter reads based on k-mer overlap"""
    print("Processing reads...", file=sys.stderr)
    
    total_reads = 0
    min1_reads = 0
    all_reads = 0
    
    with gzip.open(input_fq, 'rt') as infile, \
         gzip.open(output_min1, 'wt') as out_min1, \
         gzip.open(output_all, 'wt') as out_all:
        
        while True:
            # Read 4 lines (one fastq record)
            header = infile.readline()
            if not header:
                break
            
            seq = infile.readline().strip()
            plus = infile.readline()
            qual = infile.readline().strip()
            
            total_reads += 1
            
            # Extract k-mers from read
            read_kmers = extract_kmers(seq, k)
            
            if len(read_kmers) == 0:
                # Skip reads with no valid k-mers
                continue
            
            # Count matching k-mers
            matching = sum(1 for kmer in read_kmers if kmer in trusted_kmers)
            
            # Strategy 1: At least 1 matching k-mer
            if matching >= 1:
                out_min1.write(header)
                out_min1.write(seq + '\n')
                out_min1.write(plus)
                out_min1.write(qual + '\n')
                min1_reads += 1
            
            # Strategy 2: ALL k-mers match
            if matching == len(read_kmers):
                out_all.write(header)
                out_all.write(seq + '\n')
                out_all.write(plus)
                out_all.write(qual + '\n')
                all_reads += 1
            
            # Progress update
            if total_reads % 10000 == 0:
                print(f"Processed {total_reads:,} reads... "
                      f"(min1: {min1_reads:,}, all: {all_reads:,})", 
                      file=sys.stderr)
    
    return total_reads, min1_reads, all_reads

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: filter_reads.py <kmer_file> <input_fq> <output_min1> <output_all> <k>")
        sys.exit(1)
    
    kmer_file = sys.argv[1]
    input_fq = sys.argv[2]
    output_min1 = sys.argv[3]
    output_all = sys.argv[4]
    k = int(sys.argv[5])
    
    # Load trusted k-mers
    trusted_kmers = load_kmers(kmer_file)
    print("", file=sys.stderr)
    
    # Process reads
    total, min1, all_kmers = process_fastq(input_fq, output_min1, output_all, trusted_kmers, k)
    
    print("", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
    print("FILTERING RESULTS", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
    print(f"Total reads processed: {total:,}", file=sys.stderr)
    print(f"", file=sys.stderr)
    print(f"Strategy 1 (≥1 matching k-mer):", file=sys.stderr)
    print(f"  Reads retained: {min1:,} ({100*min1/total:.2f}%)", file=sys.stderr)
    print(f"", file=sys.stderr)
    print(f"Strategy 2 (ALL k-mers match):", file=sys.stderr)
    print(f"  Reads retained: {all_kmers:,} ({100*all_kmers/total:.2f}%)", file=sys.stderr)
    print("=" * 50, file=sys.stderr)
PYTHON_EOF

# Run Python filtering
python3 "$PYTHON_SCRIPT" \
    "$KMER_FILE" \
    "$NEW_LIB" \
    "${OUTPUT_PREFIX}_min1.fastq.gz" \
    "${OUTPUT_PREFIX}_all.fastq.gz" \
    "$KMER_SIZE"

STEP3_END=$(date +%s)
STEP3_TIME=$((STEP3_END - STEP3_START))
TOTAL_TIME=$((STEP1_TIME + STEP2_TIME + STEP3_TIME))

echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "=========================================="
echo "Timing summary:"
echo "  Step 1 (KMC database): ${STEP1_TIME}s"
echo "  Step 2 (KMC dump):     ${STEP2_TIME}s"
echo "  Step 3 (Filtering):    ${STEP3_TIME}s"
echo "  Total time:            ${TOTAL_TIME}s"
echo ""
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_min1.fastq.gz  - Reads with ≥1 matching k-mer (less conservative)"
echo "  ${OUTPUT_PREFIX}_all.fastq.gz   - Reads where ALL k-mers match (super conservative)"
echo ""
echo "Temporary files in: $TEMP_DIR"
echo "You can remove temp directory when done: rm -rf $TEMP_DIR"
