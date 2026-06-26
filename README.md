# AMF_analysis_pipeline
Steps for analysis from basecalled Nanopore reads to count-matrix

## 0. Concatenante basecalled reads
Reads were life basecalled on the P2 workstation with the sup model and a quality cut-off at 10 (default)
MinKNOW always generates a lot of files so we concatenate them into one file per samples. Since we are anyways doing that we can also remove low quality reads (Q<12) and filter for only reads between 700 and 1000 bp

```
for i in {01..57}
 do
    zcat /mnt/raw/nanopore/20260611_AMF_01/no_sample_id/20260611_1829_P2S_00089-A_PBK61177_d32a5bed/fastq_pass/barcode${i}/*.fastq.gz | \
    seqkit seq -M 1000 -m 700 -Q 12 | bgzip -c > ${i}.fastq.gz
 done
```

## Optional: Filter for the amplicons within Reads
Recycled a python script that I initially wanted to use to extract marker sequences from the previti Rapid-Kit data

```
for i in {01..57}
 do
    python scripts/call_markers_from_reads_both_orientations.py \
      -i input/${i}.fastq.gz \
      -p scripts/aml-primer.txt \
      --summary filter/${i}.summary \
      --remove-primers Y -o filter/${i}_filter.fastq.gz
 done
```

Or atlernatively just trim the adapters and barcodes form the reads to prepare then for NGSpeciesID

```
for i in {01..57}
 do
   dorado trim input/${i}.fastq.gz -t 12 -k SQK-NBD114.96 | \
   bgzip -c > trim/${i}.fastq.gz
 done
```

## 1. Run NGSpeciesID
NGSpecies ID clusters reads within every samples and creates a consensus that is then polished with either racon or medaka (by deafult medaka, which we use as well here). 

**Caution**: medaka does not run on kassel due to the old processor architecture, so do not use the ```all``` or ```standard``` queues, only ```bigmem```


NGSpeciesID will only use `fastq` not `fastq.gz` hence the `zcat` injection. It will also only create one folder with one consensus file per cluster in every samples, so concate all polished clusters into one fasta-file after NGSpeciesID.

```
NGSpeciesID --t 16 --ont --consensus --medaka \
 --abundance_ratio 0 \
 --fastq <(zcat ${WD}/filter/${SEED}_filter.fastq.gz) \
 --outfolder ${WD}/ngspeciesid/${SEED}

for i in $(ls ${WD}/ngspeciesid/${SEED} | grep medaka)
 do
   cat ${WD}/ngspeciesid/${SEED}/${i}/consensus.fasta
 done > ${WD}/ngspeciesid/${SEED}/all_consensus.fasta
```

Because this approach only yields about 1-2 clusters per species, I also tried NGSpeciesID on the trimmed but not amplicon-filtered reads

```
NGSpeciesID --t 16 --ont --consensus --medaka \
 --abundance_ratio 0 \
 --fastq <(zcat ${WD}/trim/${SEED}.fastq.gz) \
 --outfolder ${WD}/ngspeciesid_trim/${SEED}

for i in $(ls ${WD}/ngspeciesid_trim/${SEED} | grep medaka)
 do
   cat ${WD}/ngspeciesid_trim/${SEED}/${i}/consensus.fasta
 done > ${WD}/ngspeciesid_trim/${SEED}/all_consensus.fasta
```

## 1.a. Run vsearch instead of NGSpeciesID
Because NGSpecies ID does not have a cut-off parameter to decide till which level of identity two reads can be considered to be from the same cluster we also tried to use vsearch instead because there appears to be high degree of similarity between the reads

vsearch does not take gziped input so we pipe it through `zcat` then we convert it to fasta because there are some very high quality reads in it we have to allow for q-scores up to 50 (`--fastq_qmax 50`).

We then dereplicate the sequences to remove PCR-duplicates, while also removing singletons

Next we remove chimeric reads which are unlikely to be present here because we already removed reads larger than 1000 bp

Then we merge clusters that have an identity to each other of at least 99%

Lastly we assign the reads to the individual clusters down to an identity of 97%, which is okay here because they will always be counted to the cluster they share the highest identity with.

Because unlike with the Ilumina data we had recently not all fastq-files have the same ID (first column of fasta header) here there will be one column per fasta ID, so one for every reads. So if we want the sum of all reads belonging to an otu/cluster we need to sum up all column from the second to the last one so we have a table with two column were the first column is the OTU-ID and the second one the sum of reads

```
zcat trim/${i}.fastq.gz | \
 vsearch --fastq_filter - \
  --fastq_qmax 50 \
  --fastaout vsearch/${i}_trim.fasta

vsearch --derep_fulllength vsearch/${i}_trim.fasta \
 --output vsearch/${i}_derep.fasta \
 --sizeout --minuniquesize 2

vsearch --uchime_denovo vsearch/${i}_derep.fasta \
 --nonchimeras vsearch/${i}_nochim.fasta

vsearch --cluster_size vsearch/${i}_nochim.fasta \
 --id 0.99 \
 --centroids vsearch/${i}_otu.fasta \
 --relabel OTU_

vsearch --usearch_global vsearch/${i}_trim.fasta \
 --db vsearch/${i}_otu.fasta \
 --id 0.97 \
 --otutabout vsearch/${i}_otu_map.txt

awk -v OFS="\t" '{sum=0; for(i=2; i<=NF; i++) sum+=$i; print $1, sum}' \
 vsearch/${i}_otu_map.txt > vsearch/${i}_otu_map_sum.txt
```

## 2. Blast against AMF database
Davide gave me an AMF database file that I prepared as a blast database with a separate taxonomy file because blast does not like long sequence names.
For that we use a custom python script

```
python reformat_amf_database.py \
 -i maarjam_database_SSU.fasta \
 -o amfdb_reformat.fasta \
 -t amf_taxonomy.tsv
```

The resulting fastafile is then used to build a blast database

```
makeblastdb -in amfdb_reformat.fasta \
 -dbtype nucl \
 -out amf_db \
 -parse_seqids
```

This database can then be used to blast the clusters against it

```
blastn -query ngspeciesid_trim/${i}/all_consensus.fasta \
 -db database/amf_db \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
 -out blast/${i}.blast.out \
 -max_target_seqs 5 -num_threads 8 \
 -max_hsps 1
```

In the same way blast is used with the vsearch OTUs

```
blastn -query vsearch/${i}_otu.fasta \
 -db database/amf_db \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
 -out blast_otu/${i}.blast.out \
 -max_target_seqs 5 -num_threads 8 \
 -max_hsps 1
```

## 3. Make per sample count tables
I improved and united the various `blast_to_count_table.py` scripts. So that the user now has the option to either specify silva, unite or amf database input, which differ on the ID and taxonomy columns. Furthermore the option `--ngsspeciesid` was added as an alternative to `-c,--counts` to specify a separate counts file because NGSpeciesID will ad the read-count to the fasta name of the polished conensus cluster.

This is how the script works with vsearch input

```
python scripts/blast_to_count_table.py -d amf \
-b blast_otu/${i}.blast.out \
-t database/amf_taxonomy.tsv \
-c vsearch/${i}_otu_map_sum.txt \
-o count_tables/${i}.cluster_count_taxonomy.tsv
```

And this is how it works with NGSpeciesID input

```
python scripts/blast_to_count_table.py -d amf \
-b blast_otu/${i}.blast.out \
-t database/amf_taxonomy.tsv \
-c vsearch/${i}_otu_map_sum.txt \
-o count_tables/${i}.cluster_count_taxonomy.tsv
```

## 4. Combine per-samples tables into one count-matrix, merging fields on blast hits
Done by another custom python script `combine_sample_taxonomy_tables_from_folder.py` for either NGSpeciesID on amplicons filtered reads (`_amp`) or just trimmed reads with NGSpeciesID (`_ng`) or vsearch output (no file name extension)

```
python scripts/combine_sample_taxonomy_tables_from_folder.py -i count_tables/ -o species_matrix_AMF.tsv
python scripts/combine_sample_taxonomy_tables_from_folder.py -i count_tables_ng/ -o species_matrix_AMF_ng.tsv
python scripts/combine_sample_taxonomy_tables_from_folder.py -i count_tables_amp/ -o species_matrix_AMF_amp.tsv
```

## 5. Optional: combine all sample cluster sequences into one file
Collect all NGSpeciesID consensus files (amplicon filtered) into one fasta-file
and rename and concat in bash

```
for i in {01..57}; do cp /mnt/data/max/20260613_AMF/ngspeciesid/${i}/all_consensus.fasta ${i}_consensus.fasta; done

for file in $(ls *_consensus.fasta | sort -V); do
     sample=$(basename "$file" _consensus.fasta)
     cluster=1
      while read line; do
         if [[ "$line" =~ ^">" ]]; then
             echo ">sample${sample}_cluster${cluster}" >> all_samples.fasta
             ((cluster++))
         else
             echo "$line" >> all_samples.fasta
         fi
     done
done
```
