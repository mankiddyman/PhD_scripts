#!/bin/bash
# this is for phasing a sub and dominant genome, and some bed and gff3 and protien fastas, useful for genespace when u want to put the dominant and recessive subgenomes as different species

# Create output directories
mkdir -p dominant_subgenome recessive_subgenome

# Separate genome FASTA - dominant
awk 'NR==FNR{dom[$1]=1; next} /^>/ {header=$0; flag=(substr($1,2) in dom)} flag' dom_genome_ids.txt Nepenthes_gracilis_renamed_chr_asm.fa > dominant_subgenome/dominant_genome.fa

# Separate genome FASTA - recessive  
awk 'NR==FNR{dom[$1]=1; next} /^>/ {header=$0; flag=!(substr($1,2) in dom)} flag' dom_genome_ids.txt Nepenthes_gracilis_renamed_chr_asm.fa > recessive_subgenome/recessive_genome.fa

# Separate GFF3 - dominant
awk 'NR==FNR{dom[$1]=1; next} /^#/ || dom[$1]' dom_genome_ids.txt N_gracilis_helixer.gff3 > dominant_subgenome/dominant_helixer.gff3

# Separate GFF3 - recessive
awk 'NR==FNR{dom[$1]=1; next} /^#/ || !dom[$1]' dom_genome_ids.txt N_gracilis_helixer.gff3 > recessive_subgenome/recessive_helixer.gff3

# Separate BED - dominant
awk 'NR==FNR{dom[$1]=1; next} dom[$1]' dom_genome_ids.txt N_gracilis_helixer.bed > dominant_subgenome/dominant_helixer.bed

# Separate BED - recessive
awk 'NR==FNR{dom[$1]=1; next} !dom[$1]' dom_genome_ids.txt N_gracilis_helixer.bed > recessive_subgenome/recessive_helixer.bed

# Extract protein IDs from dominant GFF3
awk '$3=="mRNA" {match($9, /ID=([^;]+)/, arr); print arr[1]}' dominant_subgenome/dominant_helixer.gff3 > dominant_proteins.tmp

# Separate protein FASTA - dominant
awk 'NR==FNR{gsub(/-/,"_",$1); prot[$1]=1; next} /^>/ {header=$0; flag=(substr($1,2) in prot)} flag' dominant_proteins.tmp N_gracilis_helixer.peptide.fa > dominant_subgenome/dominant_helixer_peptide.fa

# Extract protein IDs from recessive GFF3  
awk '$3=="mRNA" {match($9, /ID=([^;]+)/, arr); print arr[1]}' recessive_subgenome/recessive_helixer.gff3 > recessive_proteins.tmp

# Separate protein FASTA - recessive
awk 'NR==FNR{gsub(/-/,"_",$1); prot[$1]=1; next} /^>/ {header=$0; flag=(substr($1,2) in prot)} flag' recessive_proteins.tmp N_gracilis_helixer.peptide.fa > recessive_subgenome/recessive_helixer_peptide.fa

# Clean up temporary files
rm dominant_proteins.tmp recessive_proteins.tmp

echo "Separation complete!"
echo "Dominant subgenome files in: dominant_subgenome/"
echo "Recessive subgenome files in: recessive_subgenome/"
