#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./blastp_mixed.sh queries.fa DBPREFIX out.tsv
#
# Example:
#   makeblastdb -in D_capensis_PASA_proteins.fa -dbtype prot -out Dcap_PASA_prot
#   ./blastp_mixed.sh knl1_all_queries.fasta Dcap_PASA_prot knl1_capensis_blastp_results.tsv

Q=${1:?query_fasta}
DB=${2:?blast_db_prefix}          # output prefix from makeblastdb -dbtype prot
OUT=${3:-blast.mixed.tsv}

# Tunables (override via env vars)
SHORT_MAX=${SHORT_MAX:-30}        # length cutoff (aa) for blastp-short
EVALUE_SHORT=${EVALUE_SHORT:-1000}
EVALUE_LONG=${EVALUE_LONG:-1e-5}
MAX_TARGETS=${MAX_TARGETS:-100000}

# Output columns (match your manual pipeline)
OUTFMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

short="$tmpdir/queries.short.fa"
long="$tmpdir/queries.long.fa"
raw="$tmpdir/blast.raw.tsv"

# Split fasta by sequence length
awk -v m="$SHORT_MAX" '
  BEGIN{RS=">"; ORS=""; FS="\n"}
  NR>1{
    hdr=$1; seq="";
    for(i=2;i<=NF;i++){gsub(/[ \t\r]/,"",$i); seq=seq $i}
    if(length(seq)<=m){print ">"hdr"\n"seq"\n" >> s}
    else              {print ">"hdr"\n"seq"\n" >> l}
  }
' s="$short" l="$long" "$Q"

: > "$raw"

# Run BLASTP-short on short queries
if [[ -s "$short" ]]; then
  blastp -task blastp-short \
    -query "$short" -db "$DB" \
    -evalue "$EVALUE_SHORT" -word_size 2 -seg no -comp_based_stats 0 \
    -max_hsps 1 -max_target_seqs "$MAX_TARGETS" \
    -outfmt "$OUTFMT" \
    >> "$raw"
fi

# Run standard BLASTP on long queries
if [[ -s "$long" ]]; then
  blastp \
    -query "$long" -db "$DB" \
    -evalue "$EVALUE_LONG" \
    -max_hsps 1 -max_target_seqs "$MAX_TARGETS" \
    -outfmt "$OUTFMT" \
    >> "$raw"
fi

# Filter: keep best per sseqid by bitscore (col 12), then sort by qcovs (col 13) desc and bitscore desc,
# then add header.
awk -F'\t' '
  {
    s=$2; bs=$12+0
    if (!(s in best) || bs > score[s]) { best[s]=$0; score[s]=bs }
  }
  END { for (s in best) print best[s] }
' "$raw" \
| sort -t$'\t' -k13,13nr -k12,12nr \
| awk -F'\t' 'BEGIN{
    OFS="\t";
    print "qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs"
  }
  { print }
' > "$OUT"

