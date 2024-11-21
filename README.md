# microbialSPLASH 

We present a set of tools extending SPLASH for microbial sequencing data. The package requires [SPLASH](https://github.com/refresh-bio/SPLASH) to be installed.

microbialSPLASH includes the following modules:
* *SPLASH alignment* used for CRISPR repeats identification in the Human Microbiome Project data,
* *SPLASH clustering* employed for strain classification.

## SPLASH alignment

The SPLASH alignment pipeline consists of two scripts: `./assemble_and_align.py` and `./extract_cas.py`.

#### `./assemble_and_align.py`
The script performs the following actions:
  1. Identification of anchors from specified FASTQ files.
  2. Optional positive and negative filtering of anchors using provided FASTA files.
  3. Generation of compactors seeded by anchors independently for each FASTQ.
  4. Alignment of the compactors to a specified target database in the FASTA format.
       
#### Usage:

```bash
assemble_and_align.py [-h] [--out_dir OUT_DIR]
                   [--positive_index_list POSITIVE_INDEX_LIST]
                   [--positive_index_k POSITIVE_INDEX_K]
                   [--negative_index_list NEGATIVE_INDEX_LIST]
                   [--negative_index_k NEGATIVE_INDEX_K]
                   fastq_list target_db_fasta

```

Positional arguments:
* `fastq_list` - list of input FASTQ files (each in a separate line),
* `target_db_fasta` - FASTA file with sequences to be used as a target database,

Options:
* `-h`, `--help`- show this help message and exit
* `--out_dir OUT_DIR` - output directory with all the results (default: `out-align`)
* `--positive_index_list POSITIVE_INDEX_LIST` - list of FASTA files for positive index (default: None)
* `--positive_index_k POSITIVE_INDEX_K` - k-mer length for positive index (default: 18)
* `--negative_index_list NEGATIVE_INDEX_LIST` - list of FASTA files for negative index (default: None)
* `--negative_index_k NEGATIVE_INDEX_K` - k-mer length for negative index (default: 18)

#### Output

The script produces the following files in the output directory:
* `anchors.scores.tsv` - all significant anchors identified by SPLASH,
* `anchors.filtered.tsv` - anchors after optional positive/negative filtering,
* `compactors.fasta` - compactors in FASTA format,
* `alignments.tsv` - alignments of compactors against database sequences.

#### `./extract_cas.py`

The script performs the following actions:
1. Extraction of compactor alignments targeting CAS proteins.
2. Deduplication of overlapping alignments using a greedy algorithm.
3. Extraction of corresponding compactors.

#### Usage

```bash
extract_cas.py [-h] [--out_dir OUT_DIR]
                      fastq_list target_db_fasta alignments_tsv
                      compactors_fasta
```

Positional arguments:
* `fastq_list` - list of FASTQ files to analyze
* `target_db_fasta` - FASTA file with sequences to be used as a target database
* `alignments_tsv` - input TSV table with alignments of compactors against the target database produced by `assemble_and_align.py`
* `compactors_fasta` - input FASTA with compactors sequences produced by `assemble_and_align.py`

Options:
* `-h`, `--help` - show this help message and exit
* `--out_dir OUT_DIR` -  output directory with results of analysis (default: `out-cas`)

#### Output

The script produces the following files in the output directory:
* `cas-alignments.tsv` - deduplicated alignments of compactors against CAS proteins,
* `cas-compactors.fasta` - corresponding compactors sequences.

### Example: CRISPR repeats identification in the Human Microbiome Project data
```bash
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
```

## SPLASH clustering

In our manuscript, strain classification is performed using unsupervised clustering of sequencing information from SPLASH. This sequencing information is one of:

- Anchors' sample embeddings: SPLASH "c-vectors" optimizing sample partitioning during the SPLASH statistical test.
- Anchors' sample-target fractions. 

Here, we provide a script to generate the 5 SPLASH feature matrices described in the M. tuberculosis classification benchmark using native SPLASH outputs.

#### Installation

To set up the singularity image: 
```sh
singularity build feature_selection_packages.sif feature_selection_packages.def
```

#### Running the feature extraction script
To run the feature extraction script, provide 3 arguments: 
1. ```<SPLASH output filepath>``` a path to the native SPLASH output folder. This contains the SPLASH significant anchors file (results.after_correction.scores.tsv), the subfolder containing binarized sample-anchor-target count (SATC)  files  called result_satc, and the subfolder containing sample embeddings (c-vectors) called result_Cjs.
2. ```<anchor count>``` after filtering anchors by count, prevalence across samples, and effect size, this integer limits the total anchor count for feature selection.
3. ```<satc_dump executable>``` this is a path to the satc_dump executable included in the SPLASH release version using which the user ran SPLASH.

The script can then be run as follows: 

```sh 
singularity exec feature_selection_packages.sif python3.9 splash_feature_selection.py <SPLASH output filepath> <anchor count> <satc_dump executable>
```

#### Outputs

Outputs are written to the directory given as ```<SPLASH output filepath>```. 

1. ```clustering_anchor_list.tsv``` a list of anchors used for SATC unpack. 
2. ```unpacked_SATC_filtered_anchors_top_10_targs.tsv``` for the top 10 targets per anchor per sample, a file containing sample, anchor, target, count (the count of this anchor-target in the sample) and:
a. anchor_target_sample_fraction : the anchor-target count in this sample divided by the sum of the anchor's top 10 targets' counts in the sample. 
b. global_target_rank : for each anchor, the rank of this anchor-target in terms of its sum of anchor_target_sample_fraction across samples. 
3. ```unpacked_SATC_filtered_RANK1_anchors_top_10_targs.tsv``` (2) but restricted to the rank 1 anchor per cluster in each sample.
4. ```SATC_MSA_selected_anchors_top10_targets.tsv``` (2) but restricted to the top-ranked anchor per cluster in each sample.

The following sample-by-anchor-target-fraction feature matrices are produced:
4. ```sample_by_anchortargetfraction_RANK1_anchors_dynamic_featurization.tsv``` for only the rank 1 anchor per cluster, with dynamic target selection. 
5. ```sample_by_anchortargetfraction_RANK1_anchors_global_featurization.tsv``` for only the rank 1 anchor per cluster, with global feature selection.
6. ```sample_by_anchortargetfraction_MSA_anchors_dynamic_featurization.tsv``` for the top-ranked available anchor per cluster per sample, with dynamic target selection. 
7. ```sample_by_anchortargetfraction_MSA_anchors_global_featurization.tsv```for the top-ranked available anchor per cluster per sample, with global target selection. 

For per-cluster rank 1 anchors only, the anchor's concatenated sample embeddings (c-vectors):
8. ```anchor_by_sample_Cj_RANK1_anchors.tsv```

#### Method description
##### Anchor filtering and offset-grouping
SPLASH-significant anchors are selected via the following restrictions, where any percentile-based restriction is computed for the subset resulting from the previous restriction: anchors containing homopolymers of length > 8 reexcluded, restricted to having > ½ the 50th percentile of anchors’ counts of nonzero samples, > ½ the 50th percentile of M (anchor count), and > ½ the 50th percentile of effect size. Per user speciications, the top N anchors are  selected by descending effect size. Offset anchors are grouped at a maximum shift distance of 3. Within offset groups, anchors are ranked by descending effect size and ties broken using counts of nonzero samples. To produce anchor-target count aggregation matrices, resulting anchors are SATC-unpacked and their top 10 targets per sample retained.

##### Fixed and dynamic anchor selection
In the resulting feature matrices, each sample reports columns for the top 10 (or fewer per availability) targets of 1 anchor per offset group. Columns report anchor-target fractions computed with respect to the anchor’s top 10 anchor-targets’ sum count in the sample. In fixed anchor selection, the highest-ranked anchor is selected per offset group, even if the anchor has 0 counts in a given sample. In dynamic anchor selection, for each sample and for each offset group, the highest-ranked anchor with nonzero counts is selected. 

##### Fixed and dynamic target selection 
In fixed target selection, each anchors’ target-sample fractions are summed across samples and ranked by descending sum fraction, with lexicographic tiebreaking. In this scheme, for two samples having the same top-ranked anchor from a given offset group, the rank 1 and rank 2 targets, for instance, will also be fixed. In dynamic anchor selection, targets are ranked in each sample by their descending target fraction, with lexicographic tiebreaking. 

#### Running the clustering script
To run the clustering script, provide 3 arguments: 
1. ```<feature table>``` a path to a file, i.e. any of outputs 4-8 from the previous section, on which sample clustering might be performed.
2. ```<omit columns>``` columns in (1) which are not features for PCA, i.e. a sample identifier or metadata category. This is passed as a comma-separated string (i.e. "sample_name,lineage,sublineage").
3. ```<scale before PCA>``` a 1 or 0 corresponding to whether the user wishes to standardize columns prior to running PCA.
4. ```<PCA dimensions>``` the number of dimensions to which to reduce via PCA (this is also the number of dimensions input to UMAP).
5. ```<run UMAP>``` a 1 or 0 corresponding to whether the user wishes to run UMAP. 
6. ```<scale before UMAP>``` a 1 or 0 corresponding to whether the user wishes to standardize the columns (now PCs) prior to running UMAP.

The script can then be run as follows: 

```sh 
singularity exec feature_selection_packages_with_umap.sif python3.9 splash_feature_selection.py <feature table> <omit columns> <scale before PCA> <PCA dimensions> <run UMAP> <scale before UMAP>
```

#### Clustering outputs

Clustering outputs are written to the directory in which the clustering script was run. 

1. ```PCA_result_<PCA dimension>_PCs.tsv``` a table of principal components, with the non-feature columns (those in ```<omit columns>```) added. 
2. ```UMAP_result.tsv``` UMAP 0 and 1 of the PCA result in (1), with non-feature columns  (those in ```<omit columns>```) added.

## Citing

George Henderson, Adam Gudys, Tavor Baharav, Punit Sundaramurthy, Marek Kokot, Peter L. Wang, Sebastian Deorowicz, Allison F. Carey, Julia Salzman.
[Ultra-efficient, unified discovery from microbial sequencing with SPLASH and precise statistical assembly](https://www.biorxiv.org/content/10.1101/2024.01.18.576133v1.full)
bioRxiv 2024.01.18.576133 (2024)
