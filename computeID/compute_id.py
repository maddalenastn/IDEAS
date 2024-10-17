import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from tqdm import tqdm
import warnings

from id_methods import id_2nn, id_local_2nn, id_pca, id_local_pca

def compute_ID(adata:ad.AnnData, group:str='all', method:str='2nn', variance_e:float=0.9, sample_size:int=0, n_neighbors:int = 0,
               n_samples:int=10, normalization:bool=True, norm_sum:int=1e6, full_output:bool=False) -> list[float]:

    # Make obs names unique
    adata.obs_names_make_unique()
    
    # Normalization
    if normalization == True:
        sc.pp.normalize_total(adata, target_sum=norm_sum)

    # Convert to array matrix
    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray()

    
    # Add a obs column to compute ID of all sample if no group is specified
    if group == 'all':
        adata.obs['all'] = ['ID']*len(adata)
    
    # Get obs names and counts
    labels, counts = np.unique(adata.obs[group], return_counts = True)

    # Set default sample_size and n_neighbors
    if (sample_size == 0)|(sample_size>counts.min()):
        sample_size = int(0.8*counts.min())
        warnings.warn(f'Warning : Sample size is null or too large. The sample size is thus set as 80% of the minimum number of cells among categories : sample size = {sample_size} ')

    if (n_neighbors == 0)|(n_neighbors>counts.min()):
        n_neighbors = int(0.1*counts.min())
        warnings.warn(f'Warning : Number of neighbors (n_neighbors) is null or too large. n_neighbors is thus set as 10% of the minimum number of cells among categories : n_neighbors = {n_neighbors} ')

    # Create a list to store the ID
    id_list = []
    
    # Loop over the unique labels
    for label in tqdm(labels, desc='Computing Intrinsic Dimension', total=len(labels)):

        print(f'Computing ID for {group} {label}:')
        
        # Filter data per labels
        data_label = adata[adata.obs[group]==label]

        
        # Compute the ID
        # select method
        if method.lower() == '2nn':
            id_ = id_2nn(data_label, sample_size = sample_size, n_samples = n_samples)
        
        elif method.lower() == 'local_2nn':
            id_ = id_local_2nn(data_label, n_neighbors = n_neighbors, n_samples = n_samples)
            
        elif method.lower() == 'pca':
            id_ = id_pca(data_label,  sample_size = sample_size, n_samples = n_samples, threshold = variance_e)

        elif method.lower() == 'local_pca':
            id_ = id_local_pca(data_label, n_neighbors = n_neighbors, n_samples = n_samples)


        else:
            raise ValueError('Method not implemented')
            
        # Append the ID
        id_list.append(id_)

    # Reshape id_list
    id_array = np.asarray(id_list).reshape(len(labels),n_samples)
        
     # Format the output
    if full_output == False:

        id_mean = np.mean(id_array, axis = 1)
        id_sigma = np.std(id_array, axis = 1)
        id_comb = [id_mean, id_sigma]
        id = pd.DataFrame(id_comb)
        id.columns = labels
        id.index = ['Mean Intrinsic Dimension', 'Standard Deviation']
        
    else:
        
        id = pd.DataFrame(id_array.T, columns = labels)

    #Remove added temporary obs column
    if 'all' in adata.obs.columns :
        adata.obs.drop(columns = ['all'], inplace = True)
    
    return id