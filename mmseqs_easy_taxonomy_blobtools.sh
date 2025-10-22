#!/bin/bash
# Author: Aaryan Bhatia
# Date 16/10/2025
# so we gonna use mmseq easytaxonomy to decontminate that one hifi library of binata
#
#
#


#fILES
to_decontaminate="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/new_hifi_reads.fastq.gz"
db="/netscratch/dep_mercier/grp_marques/Aaryan/aaryan_scratch/mmseqs_db"

ls -lha $to_decontaminate $db

# download the database, NR its a run once
mmseqs databases NR "$db" tmp --threads 48
# this is 256gb
#

# Run the easy-taxonomy workflow, keep temp files if doing next steps
# The resulting report file is KRAKEN-style, can be visualized e.g. with https://fbreitwieser.shinyapps.io/pavian/
mmseqs easy-taxonomy "$fa" "$db" mm tmp --remove-tmp-files 0 

# Optional: extract the list of plant sequences (by Viridiplantae NCBI taxon ID = 33090)
# I use seqkit here, installable as a conda module
mmseqs filtertaxdb "$db" ./tmp/latest/result mm_plants --taxon-list 33090 &&
mmseqs createtsv ./tmp/latest/query mm_plants plants_scaffolds.tsv &&
rm -r ./tmp ./mm* &&
seqkit grep -f <( cut -f1 plants_scaffolds.tsv ) "$fa" > "$( basename ${fa%.*} )"_filtered.fa



# ok that took a long time to run, couldnt identify a contaminant , no bacterial or animal enrichment 

# blobtools is in path 
# run blobtools here
#
#
#
# ok actually we are running it a bit differently
#!/bin/bash
# Author: Aaryan Bhatia
# Date: 20/10/2025
# MMseqs2 decontamination + BlobTools analysis for Drosera binata HiFi libraries

# FILES
new_hifi="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/new_hifi_reads.fastq.gz"
old_hifi="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/Files/old_hifi.fastq.gz"
db="/netscratch/dep_mercier/grp_marques/Aaryan/aaryan_scratch/mmseqs_db"
outdir="/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/binata/blobtools_analysis"
export BLASTDB="/opt/share/blastdb/ncbi"
#check that files exist
if [[ ! -f ${new_hifi} ]]; then
  echo "Error: New HiFi file not found at ${new_hifi}"
  exit 1
fi
if [[ ! -f ${old_hifi} ]]; then
  echo "Error: Old HiFi file not found at ${old_hifi}"
  exit 1
fi
# Create output directory


mkdir -p ${outdir}
cd ${outdir}


echo "=== Starting analysis at $(date) ==="
SCRIPT_START=$(date +%s)

# STEP 1: Assemble BOTH libraries separately with hifiasm
echo ""
echo "=== Assembling old HiFi library ==="
START=$(date +%s)
hifiasm -o old_hifi_asm -t 48 ${old_hifi} 2>&1 | tee old_hifi_asm.log
awk '/^S/{print ">"$2;print $3}' old_hifi_asm.bp.p_ctg.gfa > old_hifi_asm.fa
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> Old HiFi assembly completed in $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"

echo ""
echo "=== Assembling new HiFi library ==="
START=$(date +%s)
hifiasm -o new_hifi_asm -t 48 ${new_hifi} 2>&1 | tee new_hifi_asm.log
awk '/^S/{print ">"$2;print $3}' new_hifi_asm.bp.p_ctg.gfa > new_hifi_asm.fa
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> New HiFi assembly completed in $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"

# STEP 2: Map reads back to their own assemblies
echo ""
echo "=== Mapping reads back to assemblies ==="

# convert assemblies to fasta
gfa_2_fasta old_hifi_asm.bp.p_ctg.gfa old_hifi_asm.fa
gfa_2_fasta new_hifi_asm.bp.p_ctg.gfa new_hifi_asm.fa




# Old library
echo "Mapping old HiFi reads..."
START=$(date +%s)
minimap2 -ax map-hifi -t 48 old_hifi_asm.fa ${old_hifi} | \
  samtools sort -@8 -o old_hifi_mapped.bam -
samtools index old_hifi_mapped.bam
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> Old HiFi mapping completed in $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"

# New library
echo "Mapping new HiFi reads..."
START=$(date +%s)
minimap2 -ax map-hifi -t 48 new_hifi_asm.fa ${new_hifi} | \
  samtools sort -@8 -o new_hifi_mapped.bam -
samtools index new_hifi_mapped.bam
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> New HiFi mapping completed in $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"

# STEP 3: BLAST for taxonomic assignment for BlobTools
# Description: Chunked parallel BLAST of two assemblies for BlobTools
module load seqkit


###############################
# CONFIGURATION
###############################
TOTAL_THREADS=80
THREADS_PER_CHUNK=10
MAX_PARALLEL=$((TOTAL_THREADS / THREADS_PER_CHUNK))
CHUNK_COUNT=100   # number of chunks to split each assembly into
BLAST_DB="nt"
EVALUE="1e-25"
OUTFMT="6 qseqid staxids bitscore std"
###############################

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

run_blast_for_assembly() {
  local asm_file="$1"
  local label="$2"

  echo ""
  echo "=== [$(timestamp)] Running BLAST on ${label} assembly ==="

  # Assembly stats
  echo "Checking assembly stats..."
  local size contigs
  size=$(grep -v ">" "${asm_file}" | wc -c)
  contigs=$(grep -c ">" "${asm_file}")
  echo "${label} assembly: ${contigs} contigs, ${size} bp"

  local start_time=$(date +%s)
  local chunk_dir="${label}_blast_chunks_$$"
  mkdir -p "${chunk_dir}"

  echo ""
  echo "[INFO] Splitting ${label} assembly into ~${CHUNK_COUNT} chunks..."
  seqkit split2 -p "${CHUNK_COUNT}" -O "${chunk_dir}" "${asm_file}" >/dev/null
  local num_chunks
  num_chunks=$(ls "${chunk_dir}"/*.fa | wc -l)
  echo "[INFO] Created ${num_chunks} chunks for ${label} assembly"
  echo "[INFO] Running ${MAX_PARALLEL} chunks in parallel with ${THREADS_PER_CHUNK} threads each"

  local status_file="${chunk_dir}/status.txt"
  > "${status_file}"

  local chunk_num=0
  for chunk in "${chunk_dir}"/*.fa; do
    chunk_num=$((chunk_num + 1))
    (
      local chunk_start=$(date +%s)
      local chunk_name
      chunk_name=$(basename "${chunk%.fa}")

      blastn -task megablast \
        -query "${chunk}" \
        -db "${BLAST_DB}" \
        -num_threads "${THREADS_PER_CHUNK}" \
        -evalue "${EVALUE}" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -outfmt "${OUTFMT}" \
        -out "${chunk%.fa}_blast.out" \
        2> "${chunk%.fa}_blast.log"

      local chunk_end=$(date +%s)
      local chunk_elapsed=$((chunk_end - chunk_start))
      echo "${chunk_name} ${chunk_elapsed}" >> "${status_file}"
    ) &

    # Limit parallel jobs
    if (( chunk_num % MAX_PARALLEL == 0 )); then
      wait
    fi
  done

  # Monitor progress
  echo ""
  echo "[INFO] Monitoring progress..."
  while [ "$(wc -l < "${status_file}")" -lt "${num_chunks}" ]; do
    sleep 30
    local completed avg_time remaining est_remaining
    completed=$(wc -l < "${status_file}")
    avg_time=$(awk '{sum+=$2; count++} END {if(count>0) print int(sum/count); else print 0}' "${status_file}")
    remaining=$((num_chunks - completed))
    est_remaining=$(((remaining + MAX_PARALLEL - 1) / MAX_PARALLEL * avg_time))

    echo "[$(timestamp)] ${label}: ${completed}/${num_chunks} chunks done | Avg: $((avg_time / 60))m $((avg_time % 60))s per chunk | Est. remaining: $((est_remaining / 3600))h $((est_remaining % 3600 / 60))m"
  done

  wait

  echo ""
  echo "[INFO] Verifying completed chunks..."
  local found_chunks
  found_chunks=$(find "${chunk_dir}" -name '*_blast.out' | wc -l)
  if [ "${found_chunks}" -lt "${num_chunks}" ]; then
    echo "[WARNING] Some chunks missing: ${found_chunks}/${num_chunks}"
  fi

  echo "[INFO] Merging results for ${label}..."
  cat "${chunk_dir}"/*_blast.out > "${label}_hifi_blast.out"

  local avg_time_label
  avg_time_label=$(awk '{sum+=$2; count++} END {if(count>0) print int(sum/count); else print 0}' "${status_file}")

  rm -rf "${chunk_dir}"

  local end_time=$(date +%s)
  local elapsed=$((end_time - start_time))
  echo ">>> [$(timestamp)] ${label} HiFi BLAST completed in $((elapsed / 3600))h $((elapsed % 3600 / 60))m $((elapsed % 60))s"
  echo "BLAST hits found: $(wc -l < "${label}_hifi_blast.out")"

  echo "${avg_time_label}"
}

###############################
# RUN OLD ASSEMBLY
###############################
avg_time_old=$(run_blast_for_assembly "old_hifi_asm.fa" "old")

###############################
# RUN NEW ASSEMBLY
###############################
echo ""
echo "=== [$(timestamp)] Running BLAST on new assembly ==="
echo "[INFO] Estimated total time based on old assembly average chunk time: $((avg_time_old / 60))m $((avg_time_old % 60))s per chunk"
run_blast_for_assembly "new_hifi_asm.fa" "new" >/dev/null

###############################
# VERIFY OUTPUTS
###############################
echo ""
echo "[INFO] Verifying final outputs..."
if [ ! -s old_hifi_blast.out ]; then
  echo "ERROR: old_hifi_blast.out is empty or missing!"
  exit 1
fi

if [ ! -s new_hifi_blast.out ]; then
  echo "ERROR: new_hifi_blast.out is empty or missing!"
  exit 1
fi

echo ""
echo "[SUCCESS] BLAST step completed successfully for both assemblies."
echo "[INFO] Results:"
echo " - old_hifi_blast.out"
echo " - new_hifi_blast.out"


# Backup originals
cp old_hifi_blast.out old_hifi_blast.out.backup
cp new_hifi_blast.out new_hifi_blast.out.backup

# Keep only the best hit per query (highest bitscore)
awk '!seen[$1]++ || $3 > max[$1] {max[$1]=$3; best[$1]=$0} END {for (q in best) print best[q]}' \
  old_hifi_blast.out.backup > old_hifi_blast.out

awk '!seen[$1]++ || $3 > max[$1] {max[$1]=$3; best[$1]=$0} END {for (q in best) print best[q]}' \
  new_hifi_blast.out.backup > new_hifi_blast.out

echo "Filtered BLAST outputs to best hit per query"
echo "Old: $(wc -l < old_hifi_blast.out.backup) -> $(wc -l < old_hifi_blast.out) hits"
echo "New: $(wc -l < new_hifi_blast.out.backup) -> $(wc -l < new_hifi_blast.out) hits"


# STEP 4: BlobTools analysis
micromamba activate blobtools
echo ""
echo "=== Creating BlobTools databases ==="
# Old library
echo "Creating old HiFi blobtools database..."
START=$(date +%s)
blobtools create -i old_hifi_asm.fa \
  -t old_hifi_blast.out \
  -b old_hifi_mapped.bam \
  --nodes /opt/share/blastdb/taxonomy/nodes.dmp \
  --names /opt/share/blastdb/taxonomy/names.dmp \
  -o old_hifi
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> Old HiFi blobtools create completed in $((ELAPSED / 60))m $((ELAPSED % 60))s"

START=$(date +%s)
blobtools view -i old_hifi.blobDB.json -o old_hifi_table
blobtools plot -i old_hifi.blobDB.json -o old_hifi_plot
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> Old HiFi blobtools view/plot completed in $((ELAPSED / 60))m $((ELAPSED % 60))s"

# New library
echo "Creating new HiFi blobtools database..."
START=$(date +%s)
blobtools create -i new_hifi_asm.fa \
  -t new_hifi_blast.out \
  -b new_hifi_mapped.bam \
  --nodes /opt/share/blastdb/taxonomy/nodes.dmp \
  --names /opt/share/blastdb/taxonomy/names.dmp \
  -o new_hifi
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> New HiFi blobtools create completed in $((ELAPSED / 60))m $((ELAPSED % 60))s"

START=$(date +%s)
blobtools view -i new_hifi.blobDB.json -o new_hifi_table
blobtools plot -i new_hifi.blobDB.json -o new_hifi_plot
END=$(date +%s)
ELAPSED=$((END - START))
echo ">>> New HiFi blobtools view/plot completed in $((ELAPSED / 60))m $((ELAPSED % 60))s"

# STEP 5: Generate comparison statistics
echo ""
echo "=== Assembly statistics ==="
for asm in old_hifi_asm.fa new_hifi_asm.fa; do
  echo "Assembly: $asm"
  echo "Number of contigs:"
  grep ">" $asm | wc -l
  echo "Total size:"
  awk '/^>/ {next} {total+=length($0)} END {print total" bp"}' $asm
  echo ""
done

SCRIPT_END=$(date +%s)
TOTAL_ELAPSED=$((SCRIPT_END - SCRIPT_START))
echo "=== Analysis complete at $(date) ==="
echo ">>> TOTAL RUNTIME: $((TOTAL_ELAPSED / 3600))h $((TOTAL_ELAPSED % 3600 / 60))m $((TOTAL_ELAPSED % 60))s"
echo ""
echo "Check the following outputs:"
echo "  - old_hifi.blobDB.json and old_hifi_plot.*"
echo "  - new_hifi.blobDB.json and new_hifi_plot.*"
echo "  - old_hifi_table.old_hifi.blobDB.table.txt"
echo "  - new_hifi_table.new_hifi.blobDB.table.txt"
echo ""
echo "Look for differences in blob patterns between old and new libraries"
echo "Contamination will appear as separate blobs with different GC/coverage"
