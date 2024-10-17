import numpy as np
import random
from dadapy.data import Data
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA

def id_2nn(adata, sample_size, n_samples=10):
    id_samples = []
    for s in range(n_samples):
        indexes = np.array(adata.obs.index)
        random.shuffle(indexes)
        names = indexes[:sample_size]
        data = Data(np.asarray(adata[names, :].X).astype(float))
        id, _, _ = data.compute_id_2NN()
        id_samples.append(id)
    return id_samples

def id_local_2nn(adata, n_neighbors=0, n_samples=100):
    id_samples = []
    indexes = np.array(adata.obs.index)
    random.shuffle(indexes)
    cells = indexes[:n_samples]
    neighbors = NearestNeighbors(n_neighbors=n_neighbors)
    data_matrix = adata.X
    neighbors.fit(data_matrix)
    for name in cells:
        root_cell = adata[adata.obs_names == name].X
        indices = neighbors.kneighbors(root_cell, return_distance=False)
        indices = indices.flatten()
        local_data = adata[indices]
        nnz_genes = np.count_nonzero(local_data.X, axis=0) > 0
        local_data = local_data[:, nnz_genes]
        data = Data(np.asarray(local_data.X).astype(float))
        local_id, _, _ = data.compute_id_2NN()
        id_samples.append(local_id)
    return id_samples

def d_PCA(adata, names, threshold=0.9):
    X = np.asarray(adata[names, :].X)
    comp = min(np.shape(X))  
    pca_tot = PCA(n_components=comp)
    pca_tot.fit(X)
    explained = pca_tot.explained_variance_ratio_
    for i in range(comp):
        if sum(explained[:i]) > threshold:
            break 
    return i

def id_pca(adata, sample_size, n_samples=10, threshold=0.9):
    id_samples = []
    for s in range(n_samples):
        indexes = np.array(adata.obs.index)
        random.shuffle(indexes)
        names = indexes[:sample_size]
        id = d_PCA(adata, names, threshold=threshold)
        id_samples.append(id)
    return id_samples

def id_local_pca(adata, n_neighbors = 0, n_samples = 100):
              
    # Create empty list of id
    id_samples = []

    #extract n_local_measures names of cells
    indexes = np.array(adata.obs.index)
    shuffle(indexes)
    cells=indexes[:n_samples]

     # Create nearest neighbors network
    neighbors = NearestNeighbors(n_neighbors=n_neighbors)

    # Find the nearest neighbors for every cell
    data_matrix = adata.X
    neighbors.fit(data_matrix)
    
    for i,name in enumerate(tqdm(cells)):

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
        id_samples.append(local_id)

    return id_samples
    