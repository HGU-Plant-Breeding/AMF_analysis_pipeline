#!/bin/bash
##Name of Job
#SBATCH --job-name=vsearch

##maximum run time, can be left out since there is no time limit defined on the cluster
##SBATCH --time=1:00:00

##number of servers used
#SBATCH --nodes=1

##number of cpus used on each server
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=2

##Array
#SBATCH --array=1-57

##Which partition standard, gpu or all
#SBATCH --partition=bigmem

##Where to write the log file and how to call it
## %j jobnumber %x name
#SBATCH -o ../logs/%x_%j.log

WD=/mnt/data/max/20260613_AMF
cd ${WD}
i=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")

module load vsearch

zcat trim/${i}.fastq.gz | vsearch --fastq_filter - --fastq_qmax 50 --fastaout vsearch/${i}_trim.fasta

vsearch --derep_fulllength vsearch/${i}_trim.fasta --output vsearch/${i}_derep.fasta --sizeout --minuniquesize 2

vsearch --uchime_denovo vsearch/${i}_derep.fasta --nonchimeras vsearch/${i}_nochim.fasta

vsearch --cluster_size vsearch/${i}_nochim.fasta --id 0.99 --centroids vsearch/${i}_otu.fasta --relabel OTU_

vsearch --usearch_global vsearch/${i}_trim.fasta --db vsearch/${i}_otu.fasta --id 0.97 --otutabout vsearch/${i}_otu_map.txt

awk -v OFS="\t" '{sum=0; for(i=2; i<=NF; i++) sum+=$i; print $1, sum}' vsearch/${i}_otu_map.txt > vsearch/${i}_otu_map_sum.txt
