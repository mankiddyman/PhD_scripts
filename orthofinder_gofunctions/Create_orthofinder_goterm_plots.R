# so this is  a script that is meant to analyse:
# 1. the orthofinder copy numbers from stefan, we want to assert that copy numbers fit our expectation based on synteny, whereradicans is effectively diploid, pubera is octoploid and the other 2 pubera zeocin mutants are aneuploid, between tetra and octoploid.


wd <- "/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/"
setwd(wd)
orthofinder_file <- "/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/genespace/orthofinderv3/Results_Feb14/Orthogroups/Orthogroups.GeneCount.tsv"

orthofinder_df <- read.table(orthofinder_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# create a new df where all numbers are divided by their value in Rr
orthofinder_ratios_df <- orthofinder_df %>%
  mutate(Rp12c_ratio = Rp12c / (Rr ),  
         Rp16c_ratio = Rp16c / (Rr ),
         RpWT_ratio = RpWT / (Rr )) %>%
  dplyr::select(Orthogroup, Rp12c_ratio, Rp16c_ratio, RpWT_ratio,Rr)


# make a frequency polygon of the relative copy number to Rr
ggplot(orthofinder_ratios_df, aes(x = Rp12c_ratio)) +
  geom_freqpoly(binwidth = 0.5, color = "red", linewidth = 0.8) +
  geom_freqpoly(aes(x = Rp16c_ratio), binwidth = 0.5, color = "blue", linewidth = 0.8) +
  geom_freqpoly(aes(x = RpWT_ratio), binwidth = 0.5, color = "green", linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Frequency Polygon of Copy Number Ratios Relative to Rr",
       x = "Copy Number Ratio (Species/Rr)",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())
library(ggplot2)
library(tidyr)
library(dplyr)

# Reshape data from wide to long format for plotting
ortho_long <- orthofinder_df %>%
  dplyr::select(Orthogroup, Rp12c, Rp16c, RpWT, Rr) %>%
  pivot_longer(cols = c(Rp12c, Rp16c, RpWT, Rr),
               names_to = "Species",
               values_to = "Copy_Number")

# Create frequency polygon
ggplot(ortho_long, aes(x = Copy_Number, color = Species)) +
  geom_freqpoly(binwidth = 1, linewidth = 0.8) +
  scale_color_manual(values = c("Rp12c" = "#E41A1C", 
                                 "Rp16c" = "#377EB8", 
                                 "RpWT" = "#4DAF4A", 
                                 "Rr" = "#984EA3")) +
  labs(title = "Frequency Polygon of Orthogroup Copy Numbers",
       x = "Copy Number",
       y = "Frequency",
       color = "Species") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())


# Calculate means for each species
library(ggrepel)

# Calculate means for each species
species_means <- ortho_long %>%
  group_by(Species) %>%
  summarise(Mean_Copy = mean(Copy_Number, na.rm = TRUE))

a <- ggplot(ortho_long, aes(x = Copy_Number, color = Species)) +
  geom_freqpoly(binwidth = 1, linewidth = 0.8) +
  scale_color_manual(values = c("Rp12c" = "#0EC6D6", 
                                 "Rp16c" = "#0AFFF7", 
                                 "RpWT" = "#0055FF", 
                                 "Rr" = "#001011")) +
  coord_cartesian(xlim = c(0, 15)) +
  labs(title = "Frequency Polygon of Orthogroup Copy Numbers (0-15)",
       x = "Copy Number",
       y = "Frequency",
       color = "Species") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())


library(ggplot2)
library(tidyr)
library(dplyr)

# Calculate ratios relative to Rr
ortho_ratios <- orthofinder_df %>%
  mutate(Rp12c_ratio = Rp12c / (Rr ),  # add small constant to avoid division by zero
         Rp16c_ratio = Rp16c / (Rr),
         RpWT_ratio = RpWT / (Rr)) %>%
  select(Orthogroup, Rp12c_ratio, Rp16c_ratio, RpWT_ratio) %>%
  pivot_longer(cols = c(Rp12c_ratio, Rp16c_ratio, RpWT_ratio),
               names_to = "Species",
               values_to = "Ratio") %>%
  mutate(Species = gsub("_ratio", "", Species))

# Plot frequency polygon
ggplot(ortho_ratios, aes(x = Ratio, color = Species)) +
  geom_freqpoly(binwidth = 0.5, linewidth = 0.8) +
  scale_color_manual(values = c("Rp12c" = "#0EC6D6", 
                                 "Rp16c" = "#0AFFF7", 
                                 "RpWT" = "#0055FF")) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Frequency Polygon of Copy Number Ratios Relative to Rr",
       x = "Copy Number Ratio (Species/Rr)",
       y = "Frequency",
       color = "Species") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())



library(ggplot2)
library(tidyr)
library(dplyr)

# Calculate ratios relative to Rr
librarlibrary(dplyr)
library(ggplot2)

# ===============================
# Libraries
# ===============================
library(dplyr)
library(tidyr)
library(ggplot2)

# ===============================
# Calculate ratios relative to Rr
# ===============================
ortho_ratios <- orthofinder_df %>%
  mutate(
    Rp12c_ratio = Rp12c / Rr,
    Rp16c_ratio = Rp16c / Rr,
    RpWT_ratio  = RpWT  / Rr
  ) %>%
  select(Orthogroup, Rp12c_ratio, Rp16c_ratio, RpWT_ratio) %>%
  pivot_longer(
    cols = c(Rp12c_ratio, Rp16c_ratio, RpWT_ratio),
    names_to = "Species",
    values_to = "Ratio"
  ) %>%
  mutate(
    Species = gsub("_ratio", "", Species),
    Species = factor(Species, levels = c("Rp16c", "Rp12c", "RpWT"))  # <- ORDER SET HERE
  ) %>%
  filter(is.finite(Ratio))

# ===============================
# Plot
# ===============================
b <- ggplot(ortho_ratios, aes(x = Ratio, fill = Species)) +

  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = Inf,
           alpha = 0.15, fill = "red") +
  annotate("rect", xmin = 1, xmax = 4, ymin = 0, ymax = Inf,
           alpha = 0.10, fill = "orange") +

  geom_vline(xintercept = 0, color = "red",
             linewidth = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 1, color = "orange",
             linewidth = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 4, color = "orange",
             linewidth = 0.5, alpha = 0.5) +

  geom_histogram(
    binwidth = 0.25,
    position = position_dodge2(preserve = "total", padding = 0.05),
    alpha = 0.85,
    color = "black",
    linewidth = 0.2
  ) +

  annotate("text", x = 0.5, y = 13500,
           label = "Complete\ngene loss \n (0 copies)",
           hjust = 0.5, size = 3.5,
           color = "gray20", fontface = "bold") +
  annotate("text", x = 2.5, y = 13500,
           label = "Partial\ngene loss \n (1-4 copies)",
           hjust = 0.5, size = 3.5,
           color = "gray20", fontface = "bold") +

  scale_fill_manual(
    values = c("Rp16c" = "#0AFFF7",
               "Rp12c" = "#0EC6D6",
               "RpWT"  = "#0055FF")
  ) +

  scale_x_continuous(breaks = 0:10) +
  coord_cartesian(xlim = c(0, 10),
                  ylim = c(0, 15500)) +
  # change ylim to make room for annotations
  labs(
    title = "Copy Number Ratios Relative to Rr",
    subtitle = "Gene loss patterns in Rp12c and Rp16c compared to RpWT (binwidth = 0.1)",
    x = "Copy Number Ratio (Species/Rr)",
    y = "Frequency",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

b
ggsave("Copy_Number_Ratios_Histogram.pdf", b, width = 8, height = 5)

# make a quick and dirty histogram, but without copy number ratios, instead copy numbers of orthogroups directly
a <- ggplot(ortho_long, aes(x = Copy_Number, fill = Species)) +
  geom_histogram(binwidth = .2, position = "dodge", color = "black", alpha = 0.7) +
  # use the same colors as the ratio histogram
  scale_fill_manual(values = c("Rp12c" = "#0EC6D6", 
                               "Rp16c" = "#0AFFF7", 
                               "RpWT" = "#0055FF", 
                               "Rr" = "#001011")) +
  labs(title = "Histogram of Orthogroup Copy Numbers",
       x = "Copy Number",
       y = "Frequency",
       fill = "Species") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 10))+ 
  #add ticks at every integer
  scale_x_continuous(breaks = 0:10)
a
ggsave("Orthogroup_Copy_Number_Histogram.pdf", a, width = 8, height = 5)

# acc pdf is better
pdf("Orthogroup_histograms.pdf", width = 10, height = 5)
a
b
dev.off()
# Select orthogroups with copy number < rpWT in Rp12C or Rp16C
library(dplyr)
library(eulerr)
library(grid)

# --- Define "partial loss" across ALL orthogroups ---
df <- orthofinder_df %>%
  mutate(
    Rp12c_loss = Rp12c < RpWT,
    Rp16c_loss = Rp16c < RpWT
  )

# --- Counts ---
N  <- nrow(df)
n1 <- sum(df$Rp12c_loss)
n2 <- sum(df$Rp16c_loss)
k  <- sum(df$Rp12c_loss & df$Rp16c_loss)  # observed overlap

# Expected overlap under independence
expected_overlap <- (n1 * n2) / N

# Null distribution for overlap given margins: X ~ Hypergeometric(N, n1, n2)
# (m = n1 "successes" in population, n = N-n1 "failures", k = n2 draws)
null_lo <- qhyper(0.025, m = n1, n = N - n1, k = n2)
null_hi <- qhyper(0.975, m = n1, n = N - n1, k = n2)

# 2x2 table + Fisher
a <- k
b <- n1 - k
c <- n2 - k
d <- N - (a + b + c)

tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
              dimnames = list(
                Rp12c = c("loss", "no_loss"),
                Rp16c = c("loss", "no_loss")
              ))

f <- fisher.test(tab)  # default two-sided

# --- Venn list ---
partial_loss_list <- list(
  Rp12c = df %>% filter(Rp12c_loss) %>% pull(Orthogroup),
  Rp16c = df %>% filter(Rp16c_loss) %>% pull(Orthogroup)
)

fit <- euler(partial_loss_list)

# --- Annotation text (includes null interval for overlap + Fisher OR CI) ---
annot <- sprintf(
  paste0(
    "N = %d orthogroups\n",
    "Rp12c loss = %d | Rp16c loss = %d\n",
    "Observed overlap = %d\n",
    "Expected overlap (mean) = %.2f\n",
    "Null overlap 95%% interval = [%d, %d]\n",
    "Fisher: OR = %.2f (95%% CI %.2f–%.2f), p = %.3g"
  ),
  N, n1, n2, k,
  expected_overlap,
  null_lo, null_hi,
  unname(f$estimate), f$conf.int[1], f$conf.int[2], f$p.value
)

pdf("OG_overlap_venn_partial_loss_stats.pdf", width = 6, height = 6)

plot(fit,
     fills = list(fill = c("#0EC6D6", "#0AFFF7"), alpha = 0.3),
     labels = list(col = "black", font = 2, cex = 1.2),
     quantities = list(col = "black", font = 2, cex = 1),
     edges = list(col = c("#0EC6D6", "#0AFFF7"), lwd = 3),
     main = list(
       label = expression(atop(
         "Orthogroups with Partial Gene Loss relative to " * italic("R. pubera WT"),
         "(Copy Number < RpWT)"
       )),
       cex = 1.15
     ))

# Put stats text in bottom-left; tweak x/y/fontsize as needed
grid.text(annot,
          x = unit(0.02, "npc"), y = unit(0.02, "npc"),
          just = c("left", "bottom"),
          gp = gpar(col = "black", fontsize = 8.5))

dev.off()
# Create a list of orthogroups for each species
# ok now functional enrichment,
# first we should collect a list of orthogrups of interest , namely those where Rp16c or Rp12c are < RpWt
lost_12c <- orthofinder_df %>%
  filter(Rp12c < RpWT) %>%
  pull(Orthogroup)
lost_16c <- orthofinder_df %>%
  filter(Rp16c < RpWT) %>%
  pull(Orthogroup)
#we observe a large intersection
lost_both <- intersect(lost_12c, lost_16c)
#now identify unique losses
lost_12c_unique <- setdiff(lost_12c, lost_16c)
lost_16c_unique <- setdiff(lost_16c, lost_12c)


# tsv where orthogroups can be related to geneIDs,
orthogroup_genes_file <- "/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/genespace/orthofinderv3/Results_Feb14/Orthogroups/Orthogroups.tsv"


# ok now we need to link to function, ill put all species as a named vector
func_files <- c("/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/Rp12c_pep/outputs_20260210_133848/results.csv","/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/Rp16c_pep/outputs_20260210_155949/results.csv","/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/RpWTpubera_pep/outputs_20260210_184357/results.csv","/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/Rr_pep/outputs_20260212_142658/results.csv")
names(func_files) <- c("Rp12c", "Rp16c", "RpWT", "Rr")

# now load the orthogroup to geneID mapping
orthogroup_genes_df <- read.table(orthogroup_genes_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(orthogroup_genes_df)

# load func files into a list of dataframes
func_dfs <- lapply(func_files, function(file) {
  read.csv(file, header = TRUE, stringsAsFactors = FALSE)
})
names(func_dfs) <- names(func_files)
head(func_dfs[[1]])



# needed later
peptide_fastas_path <- "/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/genespace/peptide"


# ok beginnign with relating function and OG


library(dplyr)
library(tidyr)
library(stringr)
library(readr)

orthogroup_genes_df <- read_tsv(
  orthogroup_genes_file,
  show_col_types = FALSE
)

og_long <- orthogroup_genes_df %>%
  pivot_longer(
    cols = -Orthogroup,
    names_to = "species",
    values_to = "genes"
  ) %>%
  mutate(genes = na_if(genes, "")) %>%
  filter(!is.na(genes)) %>%
  separate_rows(genes, sep = ",\\s*") %>%     # split on comma + optional spaces
  mutate(genes = str_trim(genes)) %>%
  filter(genes != "")

# columns: Orthogroup, species, genes
head(og_long)


gene2go_by_species <- lapply(func_dfs, function(df) {
  df %>%
    select(gene_id = query_accession, go_id, category, go_description) %>%
    distinct() %>%
    filter(!is.na(go_id), go_id != "")
})

# example:
head(gene2go_by_species$Rp12c)


og_gene_go <- og_long %>%
  rename(gene_id = genes) %>%
  left_join(
    bind_rows(
      lapply(names(gene2go_by_species), function(sp) {
        gene2go_by_species[[sp]] %>% mutate(species = sp)
      })
    ),
    by = c("species", "gene_id")
  )

og_consensus_go <- og_gene_go %>%
  filter(!is.na(go_id)) %>%
  group_by(Orthogroup) %>%
  mutate(n_annot_members = n_distinct(gene_id)) %>%
  group_by(Orthogroup, go_id) %>%
  summarise(
    n_members_with_term = n_distinct(gene_id),
    n_annot_members = first(n_annot_members),
    .groups = "drop"
  ) %>%
  filter(n_members_with_term == n_annot_members) %>%
  select(Orthogroup, go_id)
print(nrow(og_gene_go))
print(nrow(og_consensus_go))


# calculate how many unique orthogroups in og_gene_go and og_consensus_go
unique_ogs_gene_go <- og_gene_go %>%
  filter(!is.na(go_id)) %>%
  distinct(Orthogroup) %>%
  nrow()
unique_ogs_consensus_go <- og_consensus_go %>%
  distinct(Orthogroup) %>%
  nrow()
total_orthogroups <- orthogroup_genes_df %>%
  distinct(Orthogroup) %>%
  nrow()
cat("Total orthogroups in dataset:", total_orthogroups, "\n")
cat("Unique orthogroups with any GO annotation:", unique_ogs_gene_go, "\n")
cat("Unique orthogroups with consensus GO annotation:", unique_ogs_consensus_go, "\n")


#print unique go ids and their counts for orthogroup OG0000000
og_gene_go %>%
  filter(Orthogroup == "OG0000000") %>%
  group_by(go_id) %>%
  summarise(count = n_distinct(gene_id)) %>%
  arrange(desc(count)) %>%
  print(n = Inf)



# im gonna try filtering instead by a majority rule
threshold <- 0.8  

og_consensus_go <- og_gene_go %>%
  filter(!is.na(go_id)) %>%
  group_by(Orthogroup) %>%
  mutate(n_annot_members = n_distinct(gene_id)) %>%
  group_by(Orthogroup, go_id) %>%
  summarise(
    n_members_with_term = n_distinct(gene_id),
    n_annot_members = first(n_annot_members),
    .groups = "drop"
  ) %>%
  filter(n_members_with_term / n_annot_members >= threshold) %>%
  select(Orthogroup, go_id)

length(unique(og_consensus_go$Orthogroup))

unique_ogs_majority_go <- og_consensus_go %>%
  distinct(Orthogroup) %>%
  nrow()
cat("Unique orthogroups with majority consensus GO annotation:", unique_ogs_majority_go, "\n")
print(unique_ogs_majority_go / total_orthogroups * 100)


# do we have skewing for small OGs?
og_sizes <- og_long %>%
  group_by(Orthogroup) %>%
  summarise(n_members = n_distinct(genes), .groups = "drop")

kept_ogs <- og_consensus_go %>% distinct(Orthogroup)

summary_all  <- summary(og_sizes$n_members)
summary_kept <- summary(og_sizes$n_members[og_sizes$Orthogroup %in% kept_ogs$Orthogroup])

summary_all
summary_kept


# do i lose ogs with heterogenous functions?
all_annot_ogs <- og_gene_go %>%
  filter(!is.na(go_id)) %>%
  distinct(Orthogroup)

kept_ogs <- og_consensus_go %>% distinct(Orthogroup)

lost_ogs <- anti_join(all_annot_ogs, kept_ogs, by = "Orthogroup")
nrow(lost_ogs)
# yeah those r the ones lost but thats okay


library(Biostrings)
library(goseq)
library(GO.db)
library(AnnotationDbi)

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)

library(ggplot2)
library(patchwork)

library(rrvgo)  # optional but recommended


peptide_fastas_path <- "/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan/genespace/peptide"

fasta_files <- c(
  Rp12c = file.path(peptide_fastas_path, "Rp12c.fa"),
  Rp16c = file.path(peptide_fastas_path, "Rp16c.fa"),
  RpWT  = file.path(peptide_fastas_path, "RpWT.fa"),
  Rr    = file.path(peptide_fastas_path, "Rr.fa")
)

# gene lengths from fasta
length_df <- imap_dfr(fasta_files, function(fp, sp) {
  aa <- readAAStringSet(fp)
  tibble(
    species = sp,
    gene_id = names(aa),
    aa_length = Biostrings::width(aa)
  )
})

length_df %>% summarise(n = n(), n_unique = n_distinct(gene_id))


og_long <- og_long %>% rename(gene_id = genes)  # if needed

og_long <- og_long %>% rename(genes = "gene_id")  # if needed
og_length <- og_long %>%
  left_join(length_df, by = c("species", "gene_id")) %>%
  group_by(Orthogroup) %>%
  summarise(
    bias_length = median(aa_length, na.rm = TRUE),
    n_len = sum(!is.na(aa_length)),
    .groups = "drop"
  ) %>%
  filter(is.finite(bias_length), n_len > 0)

bias_vec <- og_length$bias_length
names(bias_vec) <- og_length$Orthogroup

length(bias_vec)



og_lists <- list(
  lost_12c_unique = lost_12c_unique,
  lost_16c_unique = lost_16c_unique,
  lost_both       = lost_both
)

make_og2go_consensus <- function(og_gene_go,
                                 category_letter = c("P", "F", "C"),
                                 threshold = 0.8,
                                 min_support = 2) {
  category_letter <- match.arg(category_letter)

  # Only this ontology
  df <- og_gene_go %>%
    dplyr::filter(category == category_letter, !is.na(go_id))

  # Annotated members per OG (only genes with ≥1 GO in this ontology)
  og_n_annot <- df %>%
    dplyr::distinct(Orthogroup, gene_id) %>%
    dplyr::group_by(Orthogroup) %>%
    dplyr::summarise(n_annot = dplyr::n_distinct(gene_id), .groups = "drop")

  # Count how many genes per OG have each GO term
  og_go_counts <- df %>%
    dplyr::distinct(Orthogroup, gene_id, go_id) %>%
    dplyr::group_by(Orthogroup, go_id) %>%
    dplyr::summarise(n_with = dplyr::n_distinct(gene_id), .groups = "drop") %>%
    dplyr::left_join(og_n_annot, by = "Orthogroup") %>%
    dplyr::mutate(prop = n_with / n_annot)

  # Apply 80% rule + minimum support
  og2go <- og_go_counts %>%
    dplyr::filter(
      n_annot > 0,
      prop >= threshold,
      n_with >= min_support
    ) %>%
    dplyr::select(Orthogroup, go_id) %>%
    dplyr::distinct()

  og2go
}

og2go_P <- make_og2go_consensus(og_gene_go, "P", threshold = 0.8, min_support = 2)
og2go_F <- make_og2go_consensus(og_gene_go, "F", threshold = 0.8, min_support = 2)
og2go_C <- make_og2go_consensus(og_gene_go, "C", threshold = 0.8, min_support = 2)


sapply(list(P = og2go_P, F = og2go_F, C = og2go_C),
       function(x) length(unique(x$Orthogroup)))



make_gene2cat <- function(og2go_df) {
  split(og2go_df$go_id, og2go_df$Orthogroup)
}

gene2cat_P <- make_gene2cat(og2go_P)
gene2cat_F <- make_gene2cat(og2go_F)
gene2cat_C <- make_gene2cat(og2go_C)



make_universe <- function(gene2cat, bias_vec) {
  intersect(names(gene2cat), names(bias_vec))
}

universe_P <- make_universe(gene2cat_P, bias_vec)
universe_F <- make_universe(gene2cat_F, bias_vec)
universe_C <- make_universe(gene2cat_C, bias_vec)

length(universe_P)
length(universe_F)
length(universe_C)



check_overlap <- function(og_list, universe) {
  c(
    total_in_list = length(og_list),
    in_universe   = sum(og_list %in% universe)
  )
}

check_overlap(lost_12c_unique, universe_P)
check_overlap(lost_16c_unique, universe_P)
check_overlap(lost_both, universe_P)



# how many lost_both OGs have any BP GO at all (even before consensus)?
lost_both_anyP <- og_gene_go %>%
  filter(category == "P", !is.na(go_id), Orthogroup %in% lost_both) %>%
  distinct(Orthogroup) %>%
  nrow()

c(
  lost_both_total = length(lost_both),
  lost_both_anyP  = lost_both_anyP,
  lost_both_in_universeP = sum(lost_both %in% universe_P)
)



# plotting time

# =========================
# GOseq + redundancy reduction + 2-page PDF
# Page 1: 3 panels (one per og_list), ALL ontologies pooled (BP/MF/CC together)
# Page 2: 9 panels (3 lists × 3 ontologies)
# =========================

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(goseq)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)

# ---- 1) Helpers ----

go_id_to_name <- function(go_ids) {
  out <- rep(NA_character_, length(go_ids))
  names(out) <- go_ids
  ok <- go_ids %in% keys(GO.db::GOTERM)
  out[ok] <- AnnotationDbi::Term(GO.db::GOTERM[go_ids[ok]])
  out
}

get_ancestors <- function(go_ids, ontology = c("BP", "MF", "CC")) {
  ontology <- match.arg(ontology)
  anc_db <- switch(
    ontology,
    BP = GO.db::GOBPANCESTOR,
    MF = GO.db::GOMFANCESTOR,
    CC = GO.db::GOCCANCESTOR
  )
  anc <- AnnotationDbi::mget(go_ids, anc_db, ifnotfound = NA)
  lapply(anc, function(x) {
    if (length(x) == 1 && is.na(x)) character(0) else as.character(x)
  })
}

reduce_redundancy_elim_ancestors <- function(df_terms, ontology = c("BP", "MF", "CC")) {
  # Keep best terms; drop ancestors of already-kept terms (post-hoc redundancy reduction)
  ontology <- match.arg(ontology)
  if (nrow(df_terms) == 0) return(df_terms)

  df_terms <- df_terms %>% arrange(FDR, desc(GeneRatio), desc(Count))
  anc_map <- get_ancestors(df_terms$go_id, ontology = ontology)

  blocked <- character(0)
  keep <- logical(nrow(df_terms))

  for (i in seq_len(nrow(df_terms))) {
    term <- df_terms$go_id[i]
    if (term %in% blocked) next
    keep[i] <- TRUE
    blocked <- union(blocked, anc_map[[term]])
  }

  df_terms[keep, , drop = FALSE]
}

run_goseq_one_list <- function(target_ogs, universe, gene2cat, bias_vec) {
  universe <- intersect(universe, names(bias_vec))
  target <- intersect(target_ogs, universe)

  selected <- as.integer(universe %in% target)
  names(selected) <- universe
  n_selected <- sum(selected)

  pwf <- nullp(selected, bias.data = bias_vec[universe])
  res <- goseq(pwf, gene2cat = gene2cat[universe], method = "Wallenius") %>% as_tibble()

  # Expect these columns in goseq output
  stopifnot(all(c("category", "over_represented_pvalue", "numDEInCat", "numInCat") %in% names(res)))

  res %>%
    transmute(
      go_id      = category,
      PValue     = over_represented_pvalue,
      FDR        = p.adjust(PValue, method = "BH"),
      Count      = numDEInCat,
      BgCount    = numInCat,
      ListSize   = n_selected,
      UniverseN  = length(universe),
      GeneRatio  = ifelse(n_selected > 0, Count / n_selected, NA_real_)
    ) %>%
    arrange(FDR, desc(GeneRatio), desc(Count))
}

# ---- 2) Main function ----
library(viridis)

run_og_go_enrichment_pdf <- function(
  og_lists,
  bias_vec,
  gene2cat_P, universe_P,
  gene2cat_F, universe_F,
  gene2cat_C, universe_C,
  out_pdf = "OG_GO_enrichment_2page.pdf",
  fdr_cutoff = 0.05,
  top_n = 10
) {
  ontologies <- list(
    BP = list(gene2cat = gene2cat_P, universe = universe_P),
    MF = list(gene2cat = gene2cat_F, universe = universe_F),
    CC = list(gene2cat = gene2cat_C, universe = universe_C)
  )

  # A) Run goseq: list × ontology
  all_res <- imap_dfr(og_lists, function(og_vec, list_name) {
    imap_dfr(ontologies, function(ont, ont_name) {
      run_goseq_one_list(
        target_ogs = og_vec,
        universe   = ont$universe,
        gene2cat   = ont$gene2cat,
        bias_vec   = bias_vec
      ) %>%
        mutate(CategoryList = list_name, Ontology = ont_name)
    })
  })

  # B) Significant only + GO names
  sig <- all_res %>%
    filter(!is.na(FDR), FDR < fdr_cutoff, Count > 0)

  go_names <- go_id_to_name(unique(sig$go_id))
  sig <- sig %>%
    mutate(
      Term = unname(go_names[go_id]),
      Term = ifelse(is.na(Term), go_id, Term),
      # add GO IDs + ontology tag for clarity in pooled plots
      TermLabel = paste0("[", Ontology, "] ", Term, " (", go_id, ")"),
      TermLabelWrapped = str_wrap(TermLabel, width = 55)
    )

  # C) Redundancy reduction happens AFTER testing
  #    Page 2: per panel (CategoryList × Ontology) reduce -> top_n
  reduced_top_page2 <- sig %>%
    group_by(CategoryList, Ontology) %>%
    group_modify(~{
      df <- .x %>% arrange(FDR, desc(GeneRatio), desc(Count))
      df_red <- reduce_redundancy_elim_ancestors(df, ontology = unique(df$Ontology))
      df_red %>% slice_head(n = top_n)
    }) %>%
    ungroup()

  #    Page 1: pooled ontologies per CategoryList
  #    Approach: reduce redundancy WITHIN each ontology, then combine, then pick top_n overall.
  reduced_top_page1 <- sig %>%
    group_by(CategoryList, Ontology) %>%
    group_modify(~{
      df <- .x %>% arrange(FDR, desc(GeneRatio), desc(Count))
      reduce_redundancy_elim_ancestors(df, ontology = unique(df$Ontology))
    }) %>%
    ungroup() %>%
    group_by(CategoryList) %>%
    arrange(FDR, desc(GeneRatio), desc(Count), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  # D) Plotting

  # Page 1: 3 panels (one per category list), ontologies pooled
  p1 <- ggplot(
    reduced_top_page1 %>%
      mutate(CategoryList = factor(CategoryList, levels = names(og_lists))),
    aes(x = GeneRatio, y = TermLabelWrapped)
  ) +
    geom_point(aes(size = Count, color = -log10(FDR))) +
facet_wrap(~ CategoryList, ncol = 1, scales = "free_y", drop=FALSE) +
    scale_color_viridis_c(option = "plasma", end = 0.9)+
    labs(
      title = paste0("GO enrichment (goseq length-corrected; redundancy-reduced) — Top ", top_n,
                     " terms per category (BP/MF/CC pooled)"),
      x = "GeneRatio (Count / ListSize)",
      y = NULL,
      size = "Count",
      alpha = "-log10(FDR)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  # Page 2: 9 panels (3 lists × 3 ontologies), top_n per panel
  p2 <- ggplot(
    reduced_top_page2 %>%
      mutate(
        CategoryList = factor(CategoryList, levels = names(og_lists)),
        Ontology = factor(Ontology, levels = c("BP", "MF", "CC"))
      ),
    aes(x = GeneRatio, y = TermLabelWrapped)
  ) +
    geom_point(aes(size = Count, alpha = -log10(FDR))) +
    facet_grid(CategoryList ~ Ontology, scales = "free_y") +
    scale_alpha_continuous(range = c(0.4, 1)) +
    labs(
      title = paste0("GO enrichment (goseq length-corrected; redundancy-reduced) — Top ", top_n,
                     " per panel (BP/MF/CC separate)"),
      x = "GeneRatio (Count / ListSize)",
      y = NULL,
      size = "Count",
      alpha = "-log10(FDR)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  # E) Write 2-page PDF
  grDevices::pdf(out_pdf, width = 14, height = 10, onefile = TRUE)
  print(p1)
  print(p2)
  grDevices::dev.off()

  invisible(list(
    all_results = all_res,
    significant = sig,
    reduced_top_page1 = reduced_top_page1,
    reduced_top_page2 = reduced_top_page2,
    plot_page1 = p1,
    plot_page2 = p2,
    pdf_path = out_pdf
  ))
}

# ---- 3) Run it ----
# Make sure you have:
# og_lists, bias_vec, gene2cat_P/F/C, universe_P/F/C already created

res <- run_og_go_enrichment_pdf(
  og_lists = og_lists,
  bias_vec = bias_vec,
  gene2cat_P = gene2cat_P, universe_P = universe_P,
  gene2cat_F = gene2cat_F, universe_F = universe_F,
  gene2cat_C = gene2cat_C, universe_C = universe_C,
  out_pdf = "OG_GO_enrichment_2page.pdf",
  fdr_cutoff = 0.05,
  top_n = 10
)

res$pdf_path


# save other plots to pdf
pdf("OG_copy_number_distributions.pdf", width = 12, height = 6)
print(a)
print(b)
dev.off()

pdf("OG_overlap_venn.pdf", width = 6, height = 6)
print(c)
dev.off()

