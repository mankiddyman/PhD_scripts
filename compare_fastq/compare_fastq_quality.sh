#!/bin/sh
# compare the quality distributions of two fastq.gz files and plot in R

micromamba activate r_lang

fastq_1=/biodata/dep_mercier/grp_marques/marques/Hologen/Drosera/23022LRa012_5546D_75pM-Cell1/Drosera_binata.fastq.gz
fastq_2=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/blobtools_analysis/filtered_reads/new_reads_only_keep.fastq.gz



module load seqkit

wd=$netscratch_home/Drosera/binata/results/compare_fastq_quality
mkdir -p $wd && cd $wd


fastp -i $fastq_1 -o /dev/null --json file1.json --html /dev/null --thread 32
fastp -i $fastq_2 -o /dev/null --json file2.json --html /dev/null --thread 32
