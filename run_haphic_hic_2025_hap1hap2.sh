#!/bin/bash

# Author Aaryan Bhatia
# #haphic for D scorpioides with its new hic library, should get chromosome level assembly
# strategy is using hap1 hap2 as done for binata and paradoxa
#11:16|Mo.|Sep.|15|2025


CHR_num=16

wd=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/assemblies/Drosera_scorpioides_HiC_2025_04_29_saurab_no_hom_cov/newres/HiC_scaffolding/HapHiC_hap1hap2
mkdir -p $wd
cd $wd


#HiC files
#Hic_files
hic_7039_A=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/Files/HiC_reads/Drosera_scorpioides_7039_A
hic_6181M=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/Files/HiC_reads/Drosera_scorpioides_6181M
hic_7039_A_deep=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/Files/HiC_reads/Drosera_scorpioides_7039_A_deep

# Enable extended glob patterns
shopt -s extglob

# Combine R1 and R2 files with both naming conventions
h1_files=$(printf "%s " $hic_7039_A/*R1*.fastq.gz $hic_6181M/*R1*.fastq.gz $hic_7039_A_deep/*_@(1|3).fq.gz)
h2_files=$(printf "%s " $hic_7039_A/*R2*.fastq.gz $hic_6181M/*R2*.fastq.gz $hic_7039_A_deep/*_@(2|4).fq.gz) 


hap_dir=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/assemblies/Drosera_scorpioides_HiC_2025_04_29_saurab_no_hom_cov/newres

#defining assembly file
hap1_p_ctg=$hap_dir/new.hic.hap1.p_ctg.gfa
hap2_p_ctg=$hap_dir/new.hic.hap2.p_ctg.gfa

ls -lh $hap1_p_ctg $hap2_p_ctg

source ~/.bashrc
gfa_2_fasta $hap1_p_ctg $hap_dir/hap1_p_ctg.fa
gfa_2_fasta $hap2_p_ctg $hap_dir/hap2_p_ctg.fa

cat $hap_dir/hap1_p_ctg.fa $hap_dir/hap2_p_ctg.fa > $hap_dir/all_haps.fa

ref=$hap_dir/all_haps.fa


#test if hic files can be found
echo "=== Testing Hi-C files (R1) ==="
for f in $h1_files; do
    if [[ -f $f ]]; then
        ls -lh $f
    else
        echo "Missing: $f"
    fi
done

echo "=== Testing Hi-C files (R2) ==="
for f in $h2_files; do
    if [[ -f $f ]]; then
        ls -lh $f
    else
        echo "Missing: $f"
    fi
done

#test if p_ctg.gfa files can be found
echo "=== Testing p_ctg.gfa files ==="
if [[ -f $hap1_p_ctg ]]; then
    ls -lh $hap1_p_ctg
else
    echo "Missing: $hap1_p_ctg"
fi
if [[ -f $hap2_p_ctg ]]; then
    ls -lh $hap2_p_ctg
else
    echo "Missing: $hap2_p_ctg"
fi
#test if fasta file can be found
echo "=== Testing assembly fasta file (all_haps.fa) ==="
if [[ -f $ref ]]; then
    ls -lh $ref
else
    echo "Missing: $ref"
fi




bwa index $ref
# prepare arrays of R1 and R2 files
#

HIC_r1_files=($h1_files) # Convert space-separated list to array
HIC_r2_files=($h2_files)

#initialise an array to store BAM file paths
BAM_files=()




for i in "${!HIC_r1_files[@]}"; do
    HIC_r1=${HIC_r1_files[$i]}
    HIC_r2=${HIC_r2_files[$i]}
    BAM_file="HIC_${i}.bam"
    BAM_files+=($BAM_file)

    # print file sizes
    size_r1=$(stat -c%s "$HIC_r1")
    size_r2=$(stat -c%s "$HIC_r2")
    size_total=$((size_r1 + size_r2))
    size_total_gb=$(echo "scale=2; $size_total / (1024*1024*1024)" | bc)

    echo "Processing read pair $((i+1)) out of ${#HIC_r1_files[@]}"
    echo "  R1 size: $(du -h "$HIC_r1" | cut -f1)"
    echo "  R2 size: $(du -h "$HIC_r2" | cut -f1)"
    echo "  Total: ${size_total_gb} GB"

    # start timer
    SECONDS=0
    bwa mem -t 100 -5SP $ref "$HIC_r1" "$HIC_r2" \
        | samblaster \
        | samtools view -@ 100 -S -h -b -F 3340 -o "$BAM_file"
    runtime=$SECONDS

    runtime_min=$(echo "scale=2; $runtime / 60" | bc)
    echo "BAM file $BAM_file generated in ${runtime_min} minutes"

    # time per GB of input
    if (( $(echo "$size_total_gb > 0" | bc -l) )); then
        time_per_gb=$(echo "scale=2; $runtime_min / $size_total_gb" | bc)
        echo "Time per GB: ${time_per_gb} min/GB"
    else
        echo "Skipping time-per-GB calculation (input size 0)"
    fi
done



#merge BAM files
samtools merge -@ 40 HIC.bam ${BAM_files[@]}
#now sort the merged bam file by read name
samtools sort -n -@ 40 -o HiC.sorted.bam HIC.bam
echo "BAM files merged"


source ~/.bashrc
micromamba activate haphic_py311
/netscratch/dep_mercier/grp_marques/Aaryan/methods/HapHiC/utils/filter_bam.py HiC.sorted.bam 2 --NM 3 --threads 40 | samtools view - -b -@ 40 -o HiC.filtered.bam
echo "BAM file filtered"

haphic pipeline $ref HiC.filtered.bam $CHR_num --threads 40 --processes 40 --gfa "$hap1_p_ctg,$hap2_p_ctg"
haphic pipeline $ref HiC.filtered.bam $CHR_num --threads 40 --processes 40 --gfa "$hap1_p_ctg,$hap2_p_ctg"


cd 04.build/
module load java/jdk-17.0.10
bash juicebox.sh

