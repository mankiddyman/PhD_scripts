#!/bin/bash

# running juicer for qc purposes
#
# 1. u make a JSD
JSD="/netscratch/dep_mercier/grp_marques/Aaryan/tasks/HiC_QC/6181M"

# 2. u put ur genome (contig level) in the references folder
#
REF=${JSD}/references/reference_scorpioides.fa

# 3. u make a fastq folder containing the HiC reads in question
mkdir -p ${JSD}/HiC_work/fastq

# 4. u soft link the scripts juicer needs to run
ln -s /netscratch/dep_mercier/grp_marques/Aaryan/methods/myJuicerDir/juicer/CPU ./scripts 


prefix="D_scorp"


cd ${JSD}/references

bwa index ${REF}

samtools faidx ${REF}

awk '{print $1"\t"$2}' ${REF}.fai > ${prefix}.chrom.sizes

# Run juicer
#
cd $JSD

/netscratch/dep_mercier/grp_marques/Aaryan/methods/myJuicerDir/scripts/juicer.sh \
  -g $prefix \
  -d $JSD/HiC_work \
  -p $JSD/references/${prefix}.chrom.sizes \
  -z $REF \
  -D $JSD \
  -t 32 \
  --assembly 0


# can run 3ddna from here but im not interested in it right now

