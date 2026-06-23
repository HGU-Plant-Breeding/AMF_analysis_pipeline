#!/bin/bash

##Name of Job
#SBATCH --job-name=ngspeciesid

##maximum run time, can be left out since there is no time limit defined on the cluster
##SBATCH --time=1:00:00

##number of servers used
#SBATCH --nodes=1

##number of cpus used on each server
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=2

##Which partition standard, gpu or all
#SBATCH --partition=bigmem

##Array
#SBATCH --array=1-57

##Where to write the log file and how to call it
## %j jobnumber %x name
#SBATCH -o ../logs/%x_%j.log

WD=/mnt/data/max/20260613_AMF
SEEDFILE=${WD}/samples.list
SEED=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE | cut -f1)

module load samtools

cd ${WD}

module load ngspeciesid

NGSpeciesID --t 16 --ont --consensus --medaka --abundance_ratio 0 --fastq <(zcat ${WD}/filter/${SEED}_filter.fastq.gz) --outfolder ${WD}/ngspeciesid/${SEED}

for i in $(ls ${WD}/ngspeciesid/${SEED} | grep medaka); do cat ${WD}/ngspeciesid/${SEED}/${i}/consensus.fasta; done > ${WD}/ngspeciesid/${SEED}/all_consensus.fasta
