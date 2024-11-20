import numpy as np
import pandas as pd
import os
import sys
import networkx as nx
import glob
import scipy
import multiprocessing as mp
import time

#### Define inputs. 
output_directory_name = sys.argv[1]
top_N_anchors = int(sys.argv[2])

print(output_directory_name,flush=True)

#### Move to output directory. 
os.chdir(output_directory_name)

#### Load results file. 
result = pd.read_csv('result.after_correction.scores.tsv',sep='\t',
                    usecols=['anchor','most_freq_target_1',
                             'number_nonzero_samples','M',
                             'effect_size_bin'])

#### Select anchors for which there is not a homopolymer of length >= 8 in anchor-target 1. 
result = result[(~(result['anchor'] + result['most_freq_target_1']).str.contains('AAAAAAAA')) &
(~(result['anchor'] + result['most_freq_target_1']).str.contains('TTTTTTTT')) &
(~(result['anchor'] + result['most_freq_target_1']).str.contains('GGGGGGGG')) &
(~(result['anchor'] + result['most_freq_target_1']).str.contains('CCCCCCCC'))]

#### Select anchors occurring in > the 50% of median number of samples seen for a single anchor. 
result = result[result['number_nonzero_samples'] > (result['number_nonzero_samples'].median() * 0.50)]

#### Select anchors having > 50% of the median of the maximum M. 
result = result[result['M'] > (result['M'].median() * 0.50)]

#### Require effect_size_bin > 50% of the median after these filters were applied. 
result = result[result['effect_size_bin'] > (np.percentile(result['effect_size_bin'],50) * 0.50)]

#### Take the top N anchors by effect size.  
result = result.sort_values(by='effect_size_bin',ascending=False).head(top_N_anchors).reset_index(drop=True)

#### Deduplicate anchors. 
def bpToInt(x):
    if x=='A':
        return 0
    if x=='T':
        return 1
    if x=='C':
        return 2
    if x=='G':
        return 3
    return 4 #### for nan

#### i,j-th entry is 1 if i-th row of A equals j-th row of B
def compare_rows_hash(A, B):
    n, k = A.shape
    hash_A = np.array([hash(tuple(row)) for row in A])
    hash_B = np.array([hash(tuple(row)) for row in B])
    return scipy.sparse.csr_matrix(hash_A[:, None] == hash_B)

#### Function to perform shift-distance-based anchor clustering. 
def clusterAnchors(anchLst,maxShiftDist=5):
    start = time.time()
    bpArr = np.array([[bpToInt(x) for x in s] for s in anchLst], dtype=np.uint8)

    n,k=bpArr.shape
    assert maxShiftDist<=k

    simMat = scipy.sparse.csr_matrix((n, n), dtype=bool)
    for shift in range(1, maxShiftDist + 1):
        simMatUpdate = compare_rows_hash(bpArr[:, shift:], bpArr[:, :-shift])
        simMat = simMat + simMatUpdate
    print("Time until adjacency mat constructed", time.time()-start)

    ### from this similarity matrix, generate clusters
    G = nx.from_numpy_array(simMat)
    assemblies = list(nx.connected_components(G))
    print("Time until networkx done", time.time()-start)

    ### order clusters by size
    assemblies.sort(key=len,reverse=True)

    return [[anchLst[i] for i in list(cc)] for cc in assemblies] ### output sequences

#### Shift-distance-based anchor clustering result. 
anchors = clusterAnchors(result['anchor'].tolist(), 3)

#### Produce a dataframe ranking anchors based on their effect size and count of nonzero samples. 
cluster_df = pd.DataFrame({'anchor':anchors[0]})
cluster_df['cluster_id'] = 0
for i in range(1,len(anchors)):
    cluster_df_i = pd.DataFrame({'anchor':anchors[i]})
    cluster_df_i['cluster_id'] = i
    cluster_df = pd.concat([cluster_df,cluster_df_i])

#### Add effect size and number of nonzero samples to that dataframe. 
cluster_df = cluster_df.merge(result[['anchor','number_nonzero_samples','effect_size_bin']])

#### Sort by effect size and tiebreak using # nonzero samples. 
cluster_df = cluster_df.sort_values(by=['effect_size_bin','number_nonzero_samples'],ascending=[False,False])

#### Assign a rank to each anchor within its cluster. 
cluster_df['cluster_rank'] = cluster_df.groupby(['cluster_id'])['anchor'].rank(method='first',ascending=True)

#### Write the anchor list to the disk. 
anchor_file = 'clustering_anchor_list.tsv'
cluster_df[['anchor']].to_csv(anchor_file,index=None,sep='\t')

#### Run SATC unpack.
root_filename = str(np.random.rand())[2:10]

#### Produce an intermediate directory where SATCs will be unpacked.
try:
    os.mkdir('unpack_satcs_'+root_filename)
except FileExistsError:
    pass

#### A helper function to unpack each SATC using the provided anchor list. 
def unpack(satcname, root_filename, anchor_file, keep_n_targets, satc_unpack_executable):

    #### Run SATC unpack. 
    output_satcname = 'unpack_satcs_'+root_filename + '/' + '__'.join(satcname.split('/'))
    os.system(satc_unpack_executable+' --anchor_list '+anchor_file+' '+satcname+' '+output_satcname)
    satcpath = output_satcname

    #### Load the resulting SATC. Except error if no anchors are found in the SATC. 
    try:
        satc_i = pd.read_csv(satcpath, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame({'sample':[],'anchor':[],'target':[],'count':[]})

    #### Handling for 10x or non-10x SATCs.
    if len(satc_i.columns) == 5:
        satc_i[0] = satc_i[0].astype(str) + '_' + satc_i[1]
        satc_i = satc_i.rename(columns={0:'sample',2:'anchor',3:'target',4:'count'})
    if len(satc_i.columns) == 4: 
        satc_i = satc_i.rename(columns={0:'sample',1:'anchor',2:'target',3:'count'})
    satc_i = satc_i.sort_values(by=['count','target'],ascending=[False,True])
    if keep_n_targets > 1: 
        satc_i = satc_i.groupby(['anchor','sample'])[['sample','anchor','target','count']].head(keep_n_targets)
    else: 
        satc_i = satc_i.drop_duplicates(subset=['sample','anchor'])
    return satc_i 

#### Keep this number of targets per anchor. 
keep_n_targets = 10

#### Unpack SATCs. 
satc_input = glob.glob('result_satc/*')
inputiter = [(path, root_filename, anchor_file, keep_n_targets, sys.argv[3]) for path in satc_input ]
workers = int(os.environ['SLURM_JOB_CPUS_PER_NODE']) 
print(workers,flush=True)
if __name__ == "__main__":
    with mp.Pool(workers) as p:
        outs = p.starmap(unpack, inputiter)
outs = pd.concat(outs)
satc = outs

#### Compute total sample, anchor count. 
satc = satc.merge(satc.groupby(['anchor','sample'])['count'].sum().reset_index().rename(columns={'count':'sample_count_sum'}))

#### Compute the fraction of a sample counts assigned to each target.
satc['anchor_target_sample_fraction'] = satc['count'] / satc['sample_count_sum']

#### Rank targets in samples according to their fraction. 
satc = satc.sort_values(by=['anchor_target_sample_fraction','target'],ascending=[False,True]).reset_index(drop=True)
satc['anchor_sample_rank'] = satc.groupby(['anchor','sample'])['anchor_target_sample_fraction'].rank(method='first',ascending=False)
anchor_rank = cluster_df

#### Introduce the cluster rank for each anchor in each sample. 
satc = satc.merge(anchor_rank)

#### Introduce the global rank of anchor-targets by their sample fractions. 
atc = satc.groupby(['anchor','target'])['anchor_target_sample_fraction'].sum().reset_index().sort_values(by=['anchor_target_sample_fraction','target'],ascending=[False,True]).reset_index(drop=True)
atc['global_target_rank'] = atc.groupby(['anchor'])['anchor_target_sample_fraction'].rank(method='first',ascending=False).astype(int)
atc = atc[['anchor','target','global_target_rank']]
satc = satc.merge(atc)

#### Write the SATC of all filtered and clustered anchors. 
satc.to_csv('unpacked_SATC_filtered_anchors_top_10_targs.tsv',sep='\t',index=None)

#### Restrict to the best rank per cluster, per sample.
satc = satc.merge(satc.groupby(['sample','cluster_id'])['cluster_rank'].min().reset_index())

### Write the SATC of only rank 1 anchors. 
satc[satc['cluster_rank']==1].to_csv('unpacked_SATC_filtered_RANK1_anchors_top_10_targs.tsv',sep='\t',index=None)

#### Define the features we'll use for each sample (a feature being each cluster's best available anchor + targets).
satc['cluster_target_id'] = satc['cluster_id'].astype(str) + '_' + satc['anchor_sample_rank'].astype(int).astype(str)

#### Define the features we'll use using the top N targets as fixed features. 
satc['cluster_global_target_id'] = satc['cluster_id'].astype(str) + '_' + satc['global_target_rank'].astype(int).astype(str)

#### Produce a sample x anchor-target matrix of target-sample fractions for rank 1 anchors. 
satc2 = satc[satc['cluster_rank']==1]

#### Get a list of rank 1 anchors for the c-vector aggregation. 
final_anchs = satc2['anchor'].unique()

#### Write a pivot table of dynamic-features target fractions for rank 1 anchors. 
satc22 = pd.pivot_table(data=satc2,index='sample',columns='cluster_target_id',values='anchor_target_sample_fraction').fillna(0)
satc22.to_csv('sample_by_anchortargetfraction_RANK1_anchors_dynamic_featurization.tsv',sep='\t')
satc22 = 0

#### Write a pivot table of fixed-features target fractions for rank 1 anchors. 
satc2 = satc2[satc2['global_target_rank']<=10]
satc2 = pd.pivot_table(data=satc2,index='sample',columns='cluster_global_target_id',values='anchor_target_sample_fraction').fillna(0)
satc2.to_csv('sample_by_anchortargetfraction_RANK1_anchors_global_featurization.tsv',sep='\t')
satc2 = 0

#### Write this SATC for internal records! 
satc.to_csv('SATC_MSA_selected_anchors_top10_targets.tsv',sep='\t',index=None)

#### Redefine the SATC. 
satc = satc[['sample','anchor_target_sample_fraction','cluster_target_id','cluster_global_target_id','global_target_rank']]

#### Produce feature matrices using MSA. 
satc2 = pd.pivot_table(data=satc,index='sample',columns='cluster_target_id',values='anchor_target_sample_fraction').fillna(0)
satc2.to_csv('sample_by_anchortargetfraction_MSA_anchors_dynamic_featurization.tsv',sep='\t')
satc2 = pd.pivot_table(data=satc[satc['global_target_rank']<=10],index='sample',columns='cluster_global_target_id',values='anchor_target_sample_fraction').fillna(0)
satc2.to_csv('sample_by_anchortargetfraction_MSA_anchors_global_featurization.tsv',sep='\t')
satc, satc2 = 0, 0

#### Process c-vectors. 
cjPaths = 'result_Cjs/*bin*.cjs'
dfArr = []
for fname in glob.glob(cjPaths):
    start = time.time()
    dfArr.append(pd.concat(chunk[chunk.anchor.isin(final_anchs)] for chunk in pd.read_csv(fname,sep='\t',chunksize=int(1E7))))
    print(fname)
    print(time.time() - start)
dfcj = pd.concat(dfArr)
dfcj['sampleName']=dfcj['sample']
dfcj[['anchor','sample','Cj']].to_csv('unpacked_Cjs_RANK1_anchors.tsv',sep='\t',index=None)
pivotedDf = dfcj.pivot(index=['anchor'], columns='sampleName', values='Cj').fillna(0)
pivotedDf.to_csv('anchor_by_sample_Cj_RANK1_anchors.tsv',sep='\t')