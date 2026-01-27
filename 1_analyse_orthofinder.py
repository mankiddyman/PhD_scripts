# in this script  we will analyse the results from orthofinder
# objective is the find out whether those genes that are allied in drosera regia are that allied in drosera capenisis, what chromosomes those r found on, whether they are randomly distributed or not



import pandas as pd
import os
import matplotlib.pyplot as plt
from io import StringIO
import io
# setup

working_dir="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/orthofinder/Results_Dec23"

# step 1 : retrieve all orthogroups
orthogroups_file_with_gene_counts=pd.read_csv(os.path.join(working_dir,"Orthogroups","Orthogroups.GeneCount.tsv"),sep="\t") # header is there
orthogroups_file_with_gene_counts.head()


with open(os.path.join(working_dir,"Comparative_Genomics_Statistics","Statistics_Overall.tsv"), 'r') as f:
    lines = f.readlines()
# 1. Join lines and split into blocks based on empty lines (double newline)
content = "".join(lines).strip()
blocks = content.split('\n\n')

# --- DataFrame 1: General Statistics ---
# This section has no header, just Key-Value pairs
df_stats = pd.read_csv(io.StringIO(blocks[0]), sep='\t', header=None, names=['Metric', 'Value'])

# --- DataFrame 2: Average number of genes per-species in orthogroup ---
# This section has a header. Note: We strip the single quote (') often found in this file (e.g. '1, '2)
df_genes_per_species = pd.read_csv(io.StringIO(blocks[1]), sep='\t')
# Clean the first column to remove the Excel-formatting apostrophe
first_col = df_genes_per_species.columns[0]
df_genes_per_species[first_col] = df_genes_per_species[first_col].astype(str).str.replace("'", "")

# --- DataFrame 3: Number of species in orthogroup ---
# This section has a header and simple integers
df_species_in_orthogroup = pd.read_csv(io.StringIO(blocks[2]), sep='\t')

#put these dfs in a list
statistics_overall_dfs=[df_stats,df_genes_per_species,df_species_in_orthogroup]
statistics_overall_dfs


# how many rows of single copy orthogroups are there
# where  the collumns that arent the first or the last one all equal one as first collumn is orthogroup and last collumn is total

single_copy_orthogroups=orthogroups_file_with_gene_counts[(orthogroups_file_with_gene_counts.iloc[:,1:-1] == 1).all(axis=1)]
single_copy_orthogroups.shape
print(f"Number of single copy orthogroups: {single_copy_orthogroups.shape[0]}")

# make sure all single copy orthogroups have total equal to number of species
for index, row in single_copy_orthogroups.iterrows():
    total=row['Total']
    num_species=len(orthogroups_file_with_gene_counts.columns)-2 # minus orthogroup and total
    if total!=num_species:
        print(f"Orthogroup {row['Orthogroup']} has total {total} but number of species is {num_species}")

# in the paper single copy busco genes were of the number 563, almost 100 extra, they have extra species but it still doesnt make sense that more species should let u find more single copy orthologues, they identified their genes using busco directly on the genome, maybe thats why




# now we are gonna make a densitree, well prepare the data.
# make a folder called single_copy_orthogroups_trees
single_copy_orthogroups_trees_dir=os.path.join(working_dir,"Single_Copy_Orthogroups_Trees")
os.makedirs(single_copy_orthogroups_trees_dir,exist_ok=True)

# in this folder export all the single copy orthogroup names in a text file
single_copy_orthogroups_names_file=os.path.join(single_copy_orthogroups_trees_dir,"single_copy_orthogroups_names.txt")
single_copy_orthogroups['Orthogroup'].to_csv(single_copy_orthogroups_names_file,index=False,header=False)

# now getting just single copy gene trees
# make a folder called Single_Copy_Orthogroups_Gene_Trees

# the actual trees are located in wd/Gene_Trees/ as {orthogroup_name}_tree.txt
single_copy_orthogroups_gene_trees_dir=os.path.join(working_dir,"Single_Copy_Orthogroups_Trees","Single_Copy_Orthogroups_Gene_Trees")
os.makedirs(single_copy_orthogroups_gene_trees_dir,exist_ok=True)

for orthogroup in single_copy_orthogroups['Orthogroup']:
    source_tree_file=os.path.join(working_dir,"Gene_Trees",f"{orthogroup}_tree.txt")
    dest_tree_file=os.path.join(single_copy_orthogroups_gene_trees_dir,f"{orthogroup}_tree.txt")
    if os.path.exists(source_tree_file):
        os.system(f"cp {source_tree_file} {dest_tree_file}")
    else:
        print(f"Tree file for orthogroup {orthogroup} does not exist at {source_tree_file}")


# make a subfolder within this folder with sanitized names
single_copy_orthogroups_gene_trees_sanitized_dir=os.path.join(single_copy_orthogroups_gene_trees_dir,"Sanitized_Names")
os.makedirs(single_copy_orthogroups_gene_trees_sanitized_dir,exist_ok=True)

# the problem is that the trees have tipnames  that are gene names, we need to change them to species names only
# the gene names can be found in the bed folder, the name of the .bed tells us the species
# the genename within the bed tells us the gene name    

bed_dir="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/bed"
# make a dictionary mapping gene name to species name
bed_dir = "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/bed"

gene_to_species = {}
gene_to_bed = {}  # optional, for better error messages

for bed_file in os.listdir(bed_dir):
    if not bed_file.endswith(".bed"):
        continue

    species_name = bed_file.replace(".bed", "")
    bed_path = os.path.join(bed_dir, bed_file)

    bed_df = pd.read_csv(bed_path, sep=r"\s+", header=None)

    for gene_name in bed_df[3]:
        if gene_name in gene_to_species:
            raise ValueError(
                f"Duplicate gene name detected:\n"
                f"  Gene: {gene_name}\n"
                f"  Species 1: {gene_to_species[gene_name]} ({gene_to_bed[gene_name]})\n"
                f"  Species 2: {species_name} ({bed_file})"
            )

        gene_to_species[gene_name] = species_name
        gene_to_bed[gene_name] = bed_file
# now we have the dictionary, we can go through each tree file and replace the tip names
from ete3 import Tree
import os

for tree_file in os.listdir(single_copy_orthogroups_gene_trees_dir):
    if not tree_file.endswith("_tree.txt"):
        continue

    in_path = os.path.join(single_copy_orthogroups_gene_trees_dir, tree_file)
    out_path = os.path.join(single_copy_orthogroups_gene_trees_sanitized_dir, tree_file)

    tree = Tree(in_path, format=1)  # format=1 = Newick with branch lengths

    for leaf in tree.iter_leaves():
        if leaf.name not in gene_to_species:
            raise KeyError(
                f"Gene '{leaf.name}' in tree '{tree_file}' "
                f"not found in any BED file"
            )

        leaf.name = gene_to_species[leaf.name]

    tree.write(outfile=out_path)

#finish preparing data now run 2_make_densitree.R
