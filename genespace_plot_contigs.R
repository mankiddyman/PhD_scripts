# time to do this genespace ting ig
# module load genespace before u run R
library(GENESPACE)
library(Biostrings)
library(GenomicRanges)
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1_0.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.1.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.2.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.3.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.4.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.5.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.6.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.7.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.8.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.9.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.10.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.11.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.12.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.13.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.14.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.15.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.16.FINAL.fa"
#asm <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.17.FINAL.fa"
#asm<-"/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.19.FINAL.fa"
#asm<-"/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.20.FINAL.fa"
#asm<-"/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.21.FINAL.fa"
asm<-"/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/results/assemblies/oct_2025/decon_blobtools_kmers/HiC_scaffolding/haphic_putg/05.post_juicebox/v1.22.FINAL.fa"
dnass<- readDNAStringSet(asm)

#plant telomere
teloKmers <- c("TTTAGGG", "CCCTAAA")


result <- find_contigsGapsTelos(
  dnass = dnass,
  teloKmers = teloKmers,
  minContigGapSize = 100,      # Minimum size of N runs to call a gap
  maxDistBtwTelo = 20,          # Max distance between neighboring telomere kmers
  minTeloSize = 200,            # Minimum size for telomere kmer cluster
  minTeloDens = 0.75,           # Minimum density (0-1) for telomere kmers
  minChrSize=20e6,            # Minimum scaffold size to include
  maxDist2end = 10000,          # Max distance to chr end for telomere call
  verbose = TRUE)
# Filter to keep only scaffold_1 through scaffold_32
# Suppose your list is called `gr_list`
# Define valid scaffolds
valid_scaffolds <- paste0("scaffold_", 1:32)
gr_list <- result
filtered_gr_list <- lapply(gr_list, function(gr) {
  gr[seqnames(gr) %in% valid_scaffolds]
})
palette=colorRampPalette(c("blue", "green"))

plot_contigs(cgt=filtered_gr_list,nColors=4,palette=palette)
#dev.off()
# Save the plot to a file
# now checking exact coords
# getting a table of telomeres over 1kb
telo <- result$telo
telo_df <- as.data.frame(telo)
telo_df[(telo_df[,1]=="scaffold_24"),]
