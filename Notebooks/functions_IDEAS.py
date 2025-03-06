import scanpy as sc
import anndata as ad
import pandas as pd
import random
import numpy as np
from dadapy.data import Data
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

def compute_ID(adata:ad.AnnData, group:str='all', method:str='2nn', variance_e:float=0.9, sample_size:int=0, n_neighbors:int = 0,
               n_samples:int=30, normalization:bool=True, norm_sum:int=1e6, full_output:bool=False) -> list[float]:

    # Make obs names unique
    adata = adata.copy()
    adata.obs.index = adata.obs.index.astype(str)
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
    
    # Create a list to store the ID
    id_list = []
    index_list = []
    
    # Loop over the unique labels
    for label in tqdm(labels, desc='Computing Intrinsic Dimension', total=len(labels)):

        print(f'Computing ID for {group} {label}:')
        
        # Filter data per labels
        data_label = adata[adata.obs[group]==label]

         # select method to compute ID
        if method.lower() == '2nn':
            id_ = id_2nn(data_label, sample_size = sample_size, n_samples = n_samples)
        
        elif method.lower() == 'local_2nn':
            id_, index_ = id_local_2nn(data_label, n_neighbors = n_neighbors, n_samples = n_samples)
            
        elif method.lower() == 'pca':
            id_ = id_pca(data_label,  sample_size = sample_size, n_samples = n_samples, threshold = variance_e)

        elif method.lower() == 'local_pca':
            id_, index_ = id_local_pca(data_label, n_neighbors = n_neighbors, n_samples = n_samples)


        else:
            raise ValueError('Method not implemented')
            
        # Append the ID
        id_list.append(id_)

        if (method.lower() == 'local_2nn')|(method.lower() == 'local_pca'):
            index_list.append(index_)

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
        id.index = index_list

    #Remove added temporary obs column
    if 'all' in adata.obs.columns :
        adata.obs.drop(columns = ['all'], inplace = True)
    
    return id

def id_local_2nn(adata, n_neighbors = 0, n_samples = 100):
              
    # Create empty list of id and names
    id_samples = []

    #extract n_local_measures names of cells
    indexes = np.array(adata.obs.index)
    random.shuffle(indexes)
    cells=indexes[:n_samples]

     # Create nearest neighbors network
    neighbors = NearestNeighbors(n_neighbors=n_neighbors)

    # Find the nearest neighbors for every cell
    data_matrix = adata.X
    neighbors.fit(data_matrix)
    
    res = [process_neighbor(adata,neighbors,x) for x in tqdm(cells, position = 0, leave = True)]
    
    id_samples,_ = zip(*res)

    return id_samples, cells

def process_neighbor(adata,neighbors,name):
    #filter adata 
    root_cell = adata[adata.obs_names == name].X

    # Find the nearest neighbors for the specific cell rapresented by 'cell_data'
    indices = neighbors.kneighbors(root_cell, return_distance=False)
    indices = indices.flatten()
    
    # Extract the local neighborhood data
    local_data = adata[indices]
    nnz_genes = np.count_nonzero(local_data.X, axis = 0)>0
    local_data = local_data[:, nnz_genes]
    
     # initialise the "Data" class with the set of coordinates
    data = Data(np.asarray(local_data.X).astype(float)) #astype float to avoid errors in dadapy
    
    # compute the intrinsic dimension using 2nn estimator
    local_id,error,_ = data.compute_id_2NN()

    return local_id,error

def id_local_pca(adata, n_neighbors = 0, n_samples = 100):
              
    # Create empty list of id
    id_samples = []
    names_list = []

    #extract n_local_measures names of cells
    indexes = np.array(adata.obs.index)
    random.shuffle(indexes)
    cells=indexes[:n_samples]

     # Create nearest neighbors network
    neighbors = NearestNeighbors(n_neighbors=n_neighbors)

    # Find the nearest neighbors for every cell
    data_matrix = adata.X
    neighbors.fit(data_matrix)

    id_samples,names_neighbor = zip(*[process_neighbor_pca(adata,neighbors,x) for x in tqdm(cells, position = 0, leave = True)])
    
    return id_samples, (cells, names_neighbor)
    

def process_neighbor_pca(adata,neighbors,name):
    #filter adata 
    root_cell = adata[adata.obs_names == name].X

    # Find the nearest neighbors for the specific cell rapresented by 'cell_data'
    indices = neighbors.kneighbors(root_cell, return_distance=False)
    indices = indices.flatten()
    
    # Extract the local neighborhood data
    local_data = adata[indices]
    nnz_genes = np.count_nonzero(local_data.X, axis = 0)>0
    local_data = local_data[:, nnz_genes]

    names = [list(adata.obs_names)[index] for index in indices]
    
    local_id = d_PCA(local_data, names)

    return local_id,names

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