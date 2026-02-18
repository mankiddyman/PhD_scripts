 #!/bin/bash
# AUTHOR: AARYANBHATIA
#13:36|Tue.|Feb.|10|2026
# we instaleld fantasia-lite from here https://github.com/CBBIO/FANTASIA-Lite currently it will only work when abhatia runs it, but we can fix that later.
#ssh gpu01
source ~/.bashrc
micromamba activate fantasia-lite

fantasia_pipeline=/netscratch/dep_mercier/grp_marques/Aaryan/methods/FANTASIA-Lite/src/fantasia_pipeline.py
BASE_WD=/netscratch/dep_mercier/grp_marques/ssteckenborn/Chromosome_breaks_project/pubera/ForAaryan

cd "$BASE_WD"
export CUDA_VISIBLE_DEVICES=2

for infile in *.fa *.fasta; do
  [[ -e "$infile" ]] || continue

  prefix="${infile%.*}"
  WD="$BASE_WD/$prefix"
  mkdir -p "$WD"

  echo "========================================"
  echo "Running FANTASIA-Lite on: $infile"
  echo "Working directory: $WD"
  echo "========================================"

  (
    cd "$WD"
    $fantasia_pipeline $BASE_WD/$infile 
  )
done
