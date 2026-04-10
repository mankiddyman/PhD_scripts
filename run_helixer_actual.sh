bash run_helixer_batch.sh \
  -n 7 \
  -o /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/results/Helixer \
  -s /netscratch/dep_mercier/grp_marques/Aaryan/methods/helixer_docker/helixer-docker_helixer_v0.3.6_cuda_12.2.2-cudnn8.sif \
  --skip-existing \
  /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/Files/aaryan_asm_2026/hap1.fasta \
  /netscratch/dep_mercier/grp_marques/Aaryan/Cuscuta/epithymum/Files/aaryan_asm_2026/hap2.fasta
