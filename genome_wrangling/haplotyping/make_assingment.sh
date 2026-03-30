#
#first i must make a split paradoxa genome

BASE="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/paradoxa/Files/genomes"

mkdir -p "$BASE/raw" "$BASE/metadata" "$BASE/derived/renamed" "$BASE/derived/haplotypes" "$BASE/scripts"

# only do this if the file is currently in genomes/
mv "$BASE/out_JBAT_V2.FINAL.fa" "$BASE/raw/"




cd "$BASE"

chr_number_2n=12

fasta_summary raw/out_JBAT_V2.FINAL.fa | head -n "$chr_number_2n" | \
awk 'BEGIN{OFS="\t"} {gsub(/^>/,"",$1); print $1,$2}' \
> metadata/fasta_summary.tsv



# creating scaffold_to_chr_assignment.tsv

ASSIGN_FILE="$BASE/metadata/scaffold_to_chr_assignment.tsv"

# need it to have tabs

{
    printf "old_id\tnew_id\thaplotype\n"
    printf "scaffold_1\tchr1_h1\thap1\n"
    printf "scaffold_2\tchr1_h2\thap2\n"
    printf "scaffold_3\tchr2_h1\thap1\n"
    printf "scaffold_4\tchr2_h2\thap2\n"
    printf "scaffold_5\tchr3_h1\thap1\n"
    printf "scaffold_6\tchr3_h2\thap2\n"
    printf "scaffold_7\tchr4_h1\thap1\n"
    printf "scaffold_8\tchr4_h2\thap2\n"
    printf "scaffold_9\tchr5_h1\thap1\n"
    printf "scaffold_10\tchr5_h2\thap2\n"
    printf "scaffold_11\tchr6_h1\thap1\n"
    printf "scaffold_12\tchr6_h2\thap2\n"
} > "$ASSIGN_FILE"


