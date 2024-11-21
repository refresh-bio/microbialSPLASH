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

# concatenate partial files into single one
cat ./input/uniprot/*.fa > ./input/uniprotkb_not_eukaryota.fa

# run the pipeline against non-eukariotic UniProtKbB sequences
./assemble_and_align.py \
	--positive_index_list ./input/lookup/positive_index.list \
	--positive_index_k 23 \
	--negative_index_list ./input/lookup/contaminant_index.list \
	--negative_index_k 18 \
	--out_dir ./out-align \
	./input/samples.list \
	./input/uniprotkb_not_eukaryota.fa \
	
# extract and analyze CAS-targeted hits	
./extract_cas.py \
	--out_dir out-cas \
	./input/samples.list \
	./input/uniprotkb_not_eukaryota.fa \
	./out-align/alignments.tsv \
	./out-align/compactors.fasta \
	
