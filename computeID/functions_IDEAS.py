
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import random
from dadapy.data import Data
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
import warnings

def compute_ID(adata:ad.AnnData, group:str='all', method:str='2nn', roots:list=None, variance_e:float=0.9, sample_size:int=0, n_neighbors:int = 0,
               n_samples:int=30, normalization:bool=True, norm_sum:int=1e6, full_output:bool=True, id_score:bool=False) -> list[float]:

    # Make obs names unique
    adata = adata.copy()
    adata.obs.index = adata.obs.index.astype(str)
    adata.obs_names_make_unique() 
    
    # Add a obs column to compute ID of all sample if no group is specified
    if group == 'all':
        adata.obs['all'] = ['ID']*len(adata)
    
    # Get obs names and counts
    labels, counts = np.unique(adata.obs[group], return_counts = True)
    min_count = counts.min() # smaller cluster size
    n_obs = adata.shape[0] #total number of cells

    # Warnings
    if method in ['2nn', 'pca']:
        if roots != None or n_neighbors != 0 :
            warnings.warn(f"input parameters 'roots' and 'n_neighbors' will be ignored in the computation of non local ID")
        if sample_size == 0:
            sample_size = int(0.8 * min_count)
            warnings.warn(f"'sample_size' not specified: automatically set to {sample_size} (80% of smallest group size: {n_obs})")

    if method in ['local_2nn', 'local_pca']:
        if group != 'all' or sample_size != 0 or n_samples != 30 or full_output != True : 
            warnings.warn(f"input parameters 'group','sample_size','n_samples' and 'full_output' will be ignored in the computation of local ID")
        if n_neighbors == 0:
            n_neighbors = max(1, int(0.1 * n_obs))
            warnings.warn(f"'n_neighbors' not specified: automatically set to {n_neighbors} (10% of total cells: {n_obs})")

    if method not in ['pca', 'local_pca']:
        if variance_e != 0.9:
            warnings.warn(f"input parameter 'variance_e' will be ignored in the computation of ID with non PCA-based methods")
    
    # Normalization
    if normalization == True:
        sc.pp.normalize_total(adata, target_sum=norm_sum)

    # Convert to array matrix
    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray()


    # select all cells as roots by default
    if roots == None:
        roots = list(adata.obs_names)
    
    # Create a list to store the ID
    id_list = []
    index_list = []
    
    # Loop over the unique labels
    for label in tqdm(labels, desc='Computing Intrinsic Dimension', total=len(labels)):

        if group == 'all':
            print(f'Computing ID:')
        else:
            print(f'Computing ID for {group} {label}:')
        
        # Filter data per labels
        data_label = adata[adata.obs[group]==label]

         # select method to compute ID
        if method.lower() == '2nn':
            id_ = id_2nn(data_label, sample_size = sample_size, n_samples = n_samples)
        
        elif method.lower() == 'local_2nn':
            id_, index_ = id_local_2nn(data_label, n_neighbors = n_neighbors, roots = roots)
            
        elif method.lower() == 'pca':
            id_ = id_pca(data_label,  sample_size = sample_size, n_samples = n_samples, threshold = variance_e)

        elif method.lower() == 'local_pca':
            id_, index_ = id_local_pca(data_label, n_neighbors = n_neighbors, roots = roots)
        else:
            raise ValueError('Method not implemented')
            
        # Append the ID
        id_list.append(id_)

        if (method.lower() == 'local_2nn')|(method.lower() == 'local_pca'):
            index_list.append(index_)

    # Reshape id_list
    if 'local' in method:
        id_array = np.asarray(id_list).reshape(len(labels),len(roots))
    else:
        id_array = np.asarray(id_list).reshape(len(labels),n_samples)


    if id_score:
            range_vals = id_array.max() - id_array.min()
            id_array = (id_array - id_array.min())/range_vals
        
     # Format the output
    if full_output == False and 'local' not in method:

        id_mean = np.mean(id_array, axis = 1)
        id_sigma = np.std(id_array, axis = 1)       

        id_comb = [id_mean, id_sigma]
        id = pd.DataFrame(id_comb)
        id.columns = labels
        id.index = ['Mean Intrinsic Dimension', 'Standard Deviation']
        
    else:
        
        id = pd.DataFrame(id_array.T, columns = labels)

        if method in ['local_2nn', 'local_pca']:
            id.index = index_list
        else:
            id.index = [f"{i+1}" for i in range(n_samples)]

    #Remove temporary obs column
    if 'all' in adata.obs.columns :
        adata.obs.drop(columns = ['all'], inplace = True)
    
    return id




def id_local_2nn(adata, n_neighbors = 0, roots=None):
              
    # Create empty list of id and names
    id_samples = []

     # Create nearest neighbors network
    neighbors = NearestNeighbors(n_neighbors=n_neighbors)

    # Find the nearest neighbors for every cell
    data_matrix = adata.X
    neighbors.fit(data_matrix)
    indices = neighbors.kneighbors(data_matrix, return_distance=False)
    indices_df = pd.DataFrame(indices, index=adata.obs_names, columns=None)

    #estimate ID on neighborhood of root cells
    res = [process_neighbor(adata,indices_df,x) for x in tqdm(roots, position = 0, leave = True)]
    id_samples,_ = zip(*res)

    return id_samples, roots




def id_local_pca(adata, n_neighbors=0, roots=None):

     # Create empty list of id
    id_samples = []

   # Create nearest neighbors network
    neighbors = NearestNeighbors(n_neighbors=n_neighbors)

    # Find the nearest neighbors for every cell
    data_matrix = adata.X
    neighbors.fit(data_matrix)
    indices = neighbors.kneighbors(data_matrix, return_distance=False)
    indices_df = pd.DataFrame(indices, index=adata.obs_names, columns=None)

    #estimate ID on neighborhood of root cells
    res = [process_neighbor_pca(adata, indices_df, x) for x in tqdm(roots, position=0, leave=True)]
    id_samples, _ = zip(*res)

    return id_samples, roots




def process_neighbor(adata,indices_df,name):
    
    # Find the nearest neighbors for the specific cell rapresented by 'cell_data'
    indices = indices_df.loc[name]
    #indices = indices.flatten()
    
    # Extract the local neighborhood data
    local_data = adata[indices]
    nnz_genes = np.count_nonzero(local_data.X, axis = 0)>0
    local_data = local_data[:, nnz_genes]
    
     # initialise the "Data" class with the set of coordinates
    data = Data(np.asarray(local_data.X).astype(float)) #astype float to avoid errors in dadapy
    
    # compute the intrinsic dimension using 2nn estimator
    local_id,error,_ = data.compute_id_2NN()

    return local_id,error




def process_neighbor_pca(adata, indices_df, name):
    # Find the nearest neighbors for the specific cell rapresented by 'cell_data'
    indices = indices_df.loc[name]

     # Extract the local neighborhood data
    local_data = adata[indices]
    nnz_genes = np.count_nonzero(local_data.X, axis=0) > 0
    local_data = local_data[:, nnz_genes]

    names = list(local_data.obs_names)
    local_id = d_PCA(local_data, names)

    return local_id, names




def id_2nn(adata, sample_size , n_samples=10):

    # Create empty list of id
    id_samples = []

    # select samples and calculate ID 
    for s in  tqdm(range(n_samples), position = 0, leave = True):
        indexes = np.array(adata.obs.index)
        random.shuffle(indexes)
        names=indexes[:sample_size]

        # initialise the "Data" class with the set of coordinates
        data = Data(np.asarray(adata[names,:].X).astype(float)) #astype float to avoid errors in dadapy
        
        # compute the intrinsic dimension using 2nn estimator
        id,error,_ = data.compute_id_2NN()
        id_samples.append(id)
 
    return id_samples




def id_pca(adata, sample_size, n_samples=10,  threshold=.9):
    
    #adata = adata.copy()
    #adata.obs.reset_index(inplace=True,drop=True)

     # Create empty list of id
    id_samples = []

    # select samples and calculate ID 
    for s in tqdm(range(n_samples), position = 0, leave = True):
        indexes = np.array(adata.obs.index)
        random.shuffle(indexes)
        names=indexes[:sample_size]

        # compute the intrinsic dimension using PCA estimator
        id = d_PCA(adata, names, threshold = threshold)
        id_samples.append(id)

    return id_samples




def d_PCA(Adata, Names, threshold=.9):
    
    X=np.asarray(Adata[Names,:].X)
    
    comp = min(np.shape(X))  
    pca_tot = PCA(n_components=comp)
    pca_tot.fit(X)
    explained=pca_tot.explained_variance_ratio_
    
    for i in range(comp):
        if sum(explained[:i])>threshold:
            break 
    return i


