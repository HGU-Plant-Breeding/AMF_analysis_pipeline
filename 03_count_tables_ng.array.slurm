#!/bin/bash
##Name of Job
#SBATCH --job-name=count_tables

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

python scripts/blast_to_count_table.py -d amf -b blast/${i}.blast.out -t database/amf_taxonomy.tsv --ngspeciesid -o count_tables_ng/${i}.cluster_count_taxonomy.tsv
