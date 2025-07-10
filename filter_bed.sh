# for retaining only chromosomes 1 2 and 3 in a bed
grep -E '^scaffold(1|2|3)\s' genome_a.bed > genome_a_filtered.bed

