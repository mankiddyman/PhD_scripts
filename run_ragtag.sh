#!/bin/bash

# Author: Aaryan Bhatia
# making a ragtag script using the p_utg manually corrected haphic scaffolds as reference to scaffold the hap1_p_ctg.fa and hap2_p_ctg.fa concatenation as per André's instruction
# 10:37|Mi.|Sep.|17|2025


# installing ragtag , dont run but im keeping it here for reference
#micromamba create env -n ragtag -c bioconda ragtag
micromamba activate ragtag

# running ragtag
wd=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/assemblies/Drosera_scorpioides_HiC_2025_04_29_saurab_no_hom_cov/newres/HiC_scaffolding/ragtag_putg_hap1hap2

cd $wd
#ref is gonna be the manually corrected putg fasta
ref=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/assemblies/Drosera_scorpioides_HiC_2025_04_29_saurab_no_hom_cov/newres/HiC_scaffolding/HapHiC_putg/05.post_juicebox/out_JBAT.FINAL.fa


all_haps=/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/scorpioides/results/assemblies/Drosera_scorpioides_HiC_2025_04_29_saurab_no_hom_cov/newres/HiC_scaffolding/HapHiC_hap1hap2/04.build/all_haps.fa

# list all files
ls -lhaL $ref $all_haps

# run ragtag
# # scaffold a query assembly

ragtag.py scaffold $ref $all_haps -t 49
