#!/usr/bin/env bash
set -euo pipefail

Q=${1:?query_fasta}
DB=${2:?blast_db_prefix}   # made with makeblastdb -dbtype prot
OUT=${3:-blast.mixed.tsv}
SHORT_MAX=${SHORT_MAX:-30}

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

short="$tmpdir/short.fa"
long="$tmpdir/long.fa"

awk -v m="$SHORT_MAX" '
  BEGIN{RS=">"; ORS=""; FS="\n"}
  NR>1{
    hdr=$1; seq="";
    for(i=2;i<=NF;i++){gsub(/[ \t\r]/,"",$i); seq=seq $i}
    if(length(seq)<=m){print ">"hdr"\n"seq"\n" >> s}
    else              {print ">"hdr"\n"seq"\n" >> l}
  }
' s="$short" l="$long" "$Q"

: > "$OUT"

if [[ -s "$short" ]]; then
  blastp -task blastp-short -query "$short" -db "$DB" \
    -evalue 1000 -word_size 2 -seg no -comp_based_stats 0 \
    -max_hsps 1 -max_target_seqs 100000 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    >> "$OUT"
fi

if [[ -s "$long" ]]; then
  blastp -query "$long" -db "$DB" \
    -evalue 1e-5 \
    -max_hsps 1 -max_target_seqs 100000 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    >> "$OUT"
fi

