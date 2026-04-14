#!/usr/bin/env bash

###############################################################################
# Project / run metadata
###############################################################################

RUN_NAME="C_epithymum_6303_B"
SPECIES_NAME="Cuscuta_epithymum"

###############################################################################
# Project root and standard directories
###############################################################################

PROJECT_ROOT="/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/CO_pipeline"

SCRIPTS_DIR="${PROJECT_ROOT}/scripts"
LIB_DIR="${PROJECT_ROOT}/lib"
RESOURCES_DIR="${PROJECT_ROOT}/resources"
RESULTS_DIR="${PROJECT_ROOT}/results"
LOG_DIR="${PROJECT_ROOT}/logs"
CHECKPOINT_DIR="${PROJECT_ROOT}/checkpoints"
TMP_DIR="${PROJECT_ROOT}/tmp"
FAILURES_DIR="${RESULTS_DIR}/failures"

###############################################################################
# User-supplied input files
###############################################################################

HAP1_FASTA="/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/assemblies/C_epithymum_haphic_mar2026/HiC_scaffolding/haphic/05.post_juicebox/hap1.fasta"
HAP1_GFF3="/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/Helixer/hap1/annotations/hap1_helixer_merged.gff3"
HAP2_FASTA="/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/assemblies/C_epithymum_haphic_mar2026/HiC_scaffolding/haphic/05.post_juicebox/hap2.fasta"
HAP2_GFF3="/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/Helixer/hap2/hap2_chr/annotations/hap2_chr_helixer_merged.gff3"
FASTQ_DIR="/netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/Files/scRNA/Cuscuta_epithymum_6303B/6303_B"
SAMPLE_NAME="6303_B"

###############################################################################
# Tool paths
###############################################################################

SAMTOOLS_BIN="/usr/bin/samtools"
PARALLEL_BIN="/usr/bin/parallel"
CELLRANGER_BIN="/opt/share/software/bin/cellranger"
MINIMAP2_BIN="/opt/share/software/bin/minimap2"
BCFTOOLS_BIN="/opt/share/software/bin/bcftools"
TABIX_BIN="/opt/share/software/bin/tabix"
BGZIP_BIN="/opt/share/software/bin/bgzip"
BEDTOOLS_BIN="/usr/bin/bedtools"

AWK_BIN="/usr/bin/awk"
SORT_BIN="/usr/bin/sort"
WC_BIN="/usr/bin/wc"
GZIP_BIN="/usr/bin/gzip"

# Optional tools
SEQKIT_BIN=""
PIGZ_BIN=""

###############################################################################
# Pipeline-owned scripts
###############################################################################

MAKE_MARKERS_SCRIPT="${SCRIPTS_DIR}/00_make_reference_markers.sh"
SC_SNP_CALLING_SCRIPT="${SCRIPTS_DIR}/sc_snp_calling.sh"
SC_CONVERSION_SCRIPT="${SCRIPTS_DIR}/sc_conversion.sh"
SC_MARKER_COUNTS_SCRIPT="${SCRIPTS_DIR}/sc_marker_counts.sh"
MARKER_COUNT_JOBS=8
MARKER_COUNT_THREADS=2
HAPCO_SCRIPT="${SCRIPTS_DIR}/hapCO_identification.R"

###############################################################################
# Micromamba
###############################################################################

MAMBA_EXE="micromamba"
MAMBA_ENV_PREFIX="/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/co_pipeline_env"
MAMBA_BIN_DIR="${MAMBA_ENV_PREFIX}/bin"

###############################################################################
# Stage output directories
###############################################################################

MARKER_REF_DIR="${RESULTS_DIR}/00_reference_markers"
REFERENCE_DIR="${RESULTS_DIR}/01_reference"
CELLRANGER_DIR="${RESULTS_DIR}/02_cellranger/${RUN_NAME}"
BARCODE_QC_DIR="${RESULTS_DIR}/03_barcode_qc"
DEMUX_DIR="${RESULTS_DIR}/04_demultiplex"
SNP_DIR="${RESULTS_DIR}/05_snps"
MARKER_DIR="${RESULTS_DIR}/06_markers"
CO_DIR="${RESULTS_DIR}/07_crossovers"
PLOT_DIR="${RESULTS_DIR}/08_plots"

###############################################################################
# Derived files
###############################################################################

HAP1_GTF="${REFERENCE_DIR}/hap1.gtf"
HAP2_GTF="${REFERENCE_DIR}/hap2.gtf"

GENOME_NAME="${SPECIES_NAME}_hap1"
CELLRANGER_REF_DIR="${REFERENCE_DIR}/${GENOME_NAME}"

HAP1_FASTA_FAI="${HAP1_FASTA}.fai"
HAP2_FASTA_FAI="${HAP2_FASTA}.fai"

MARKER_FILE="${MARKER_REF_DIR}/reference_markers.tsv"
MARKER_BED="${MARKER_REF_DIR}/reference_markers.bed"
MARKER_POS_BED="${MARKER_REF_DIR}/reference_markers.pos.bed"
HAP1_HAP2_WHOLEGENOME_VCF="${MARKER_REF_DIR}/hap1_vs_hap2.vcf.gz"

###############################################################################
# Marker generation settings
###############################################################################

MARKER_MIN_MQ=20
MARKER_MIN_BQ=20
MARKER_VARIANT_TYPE="snps"
MARKER_REQUIRE_BIALLELIC="true"

###############################################################################
# Cell Ranger settings
###############################################################################

CHEMISTRY="auto"
CELLRANGER_CORES=32
CELLRANGER_MEMGB=64
EXPECT_CELLS=10000
INCLUDE_INTRONS="true"

###############################################################################
# Barcode / read filtering thresholds
###############################################################################

# real run
MIN_TOTAL_READS=10000
MIN_UNIQ_MAP_READS=5000
MIN_UNIQ_RATIO="0.2"

# test mode
MIN_TOTAL_READS_TEST=200
TEST_READ_COUNT=300000

###############################################################################
# Stage 03 files
###############################################################################

MAPPED_BAM="${BARCODE_QC_DIR}/possorted_genome_bam.mapped.bam"
MAPPED_BAI="${MAPPED_BAM}.bai"
BARCODE_COUNTS_TSV="${BARCODE_QC_DIR}/barcode_read_counts.tsv"
BARCODE_READ_QC_TSV="${BARCODE_QC_DIR}/barcode_read_qc.tsv"
SELECTED_BARCODES_TXT="${BARCODE_QC_DIR}/selected_barcodes.txt"
DEMUX_SUMMARY_TSV="${BARCODE_QC_DIR}/demux_summary.tsv"
DEMUX_FAILURES_TSV="${FAILURES_DIR}/demux_failed_barcodes.tsv"

###############################################################################
# Marker / switch / CO thresholds
###############################################################################

MIN_MARKERS=500
MAX_SWITCH_RATE="0.10"

HAPCO_MIN_MARKERS=500
HAPCO_SEGMENT_SIZE=3000000
HAPCO_N=10

###############################################################################
# Parallelism / runtime behavior
########################################################################### w####

MINIMAP2_THREADS=40
SAMTOOLS_THREADS=8
BCFTOOLS_THREADS=16
PARALLEL_JOBS=8
DEMUX_SORT_THREADS=8
SPLIT_PROGRESS_EVERY=1000000
###############################################################################
# Test mode
###############################################################################

TEST_MODE="false"
TEST_FASTQ_DIR="${RESULTS_DIR}/test_fastqs"
REBUILD_TEST_FASTQS="false"

###############################################################################
# Logging / cleanup / rerun behavior
###############################################################################

KEEP_TEMP="false"
FORCE_RERUN_STAGE=""
STRICT_MODE="true"

###############################################################################
# Standard output file names
###############################################################################

SWITCH_STATS_TSV="${MARKER_DIR}/switches.stats.tsv"
SELECTED_FOR_CO_TXT="${MARKER_DIR}/selected_for_co.txt"

CO_SUMMARY_TSV="${CO_DIR}/co_summary.tsv"
CO_INTERVALS_BED="${CO_DIR}/co_intervals.bed"

###############################################################################
# Checkpoint files
###############################################################################

STAGE0_DONE="${CHECKPOINT_DIR}/00_make_reference_markers.done"
STAGE1_DONE="${CHECKPOINT_DIR}/01_build_reference.done"
STAGE2_DONE="${CHECKPOINT_DIR}/02_run_cellranger.done"
STAGE3_DONE="${CHECKPOINT_DIR}/03_barcode_qc_and_demux.done"
STAGE4_DONE="${CHECKPOINT_DIR}/04_snp_markers.done"
STAGE5_DONE="${CHECKPOINT_DIR}/05_co_and_plots.done"
