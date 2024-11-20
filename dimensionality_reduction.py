import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import sys

input_dt = sys.argv[1]
omit_columns = sys.argv[2]
scale_before_pca = int(sys.argv[3])
pca_dimensions = int(sys.argv[4])
run_umap = int(sys.argv[5])
scale_before_umap = int(sys.argv[6])

input_dt = pd.read_csv(input_dt, sep='\t')
values = input_dt[[i for i in input_dt.columns if i not in omit_columns]].values

#### Scale the data if needed.
if scale_before_pca: 
    values = StandardScaler().fit_transform(values)

#### Run PCA using the desired number of dimensions. 
pca = PCA(n_components=pca_dimensions)
values = pca.fit_transform(values)

#### Report the result of PCA. 
values = pd.DataFrame(values)
values.columns = ['PC_'+str(i) for i in values.columns]
values[omit_columns] = input_dt[omit_columns]
values.to_csv('PCA_result_'+str(pca_dimensions)+'_PCs.tsv',sep='\t',index=None)

#### If we want also to run UMAP, do so.
if run_umap:
    values = values[[i for i in values.columns if i not in omit_columns]].values
    
    #### Again, scale the data if needed.
    if scale_before_umap: 
        values = StandardScaler().fit_transform(values)

    #### Run UMAP. 
    reducer = umap.UMAP(low_memory=True, verbose=True, random_state=0)
    values = reducer.fit_transform(values)

    #### Write the UMAP embeddings. 
    values = pd.DataFrame(values)
    values.columns = ['UMAP_'+str(i) for i in values.columns]
    values[omit_columns] = input_dt[omit_columns]
    values.to_csv('UMAP_result.tsv',sep='\t',index=None)