#! /bin/bash

# download in parallel 3 first HMP project runs
# (please remove --n_first parameter to download all HMP runs)  
./srr-download.py \
	--n_threads 16 \
	--n_first 3 \
	./input/hmp_srr.list \
	./input/samples   

# create list of R1 FASTQ files
ls ./input/samples/*_1.fastq > ./input/samples.list

# run the pipeline against uniprotkb_NOT_Eukaryota_AND_reviewed.fasta sequence database
./assemble_and_align.py \
	--positive_index_list ./input/lookup/positive_index.list \
	--positive_index_k 23 \
	--negative_index_list ./input/lookup/contaminant_index.list \
	--negative_index_k 18 \
	--out_dir ./out-align \
	./input/samples.list \
	./input/uniprotkb_NOT_Eukaryota_AND_reviewed.fasta \
	
# extract and analyze CAS-targeted hits	
./extract_cas.py \
	--out_dir out-cas \
	./input/samples.list \
	./input/uniprotkb_NOT_Eukaryota_AND_reviewed.fasta \
	./out-align/alignments.tsv \
	./out-align/compactors.fasta \
	
