
import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import dadapy
from dadapy.data import Data
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
import warnings
from scipy.stats import pearsonr
import pickle
import shutil

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


def local_id_across_scales(adata, pca_n_comps, k_grid='None', save = False, filename_tag = ''):


    # Computes Local Intrinsic Dimension (LID) for each cell across a range of
    # neighborhood sizes k, in order to identify an "optimal" scale for LID
    # estimation. For each cell, the function builds a local neighborhood of
    # size k, computes the LID using the 2NN estimator (via dadapy), and repeats
    # this for all k in the provided k-grid.
    #
    # Parameters
    # ----------
    # adata : AnnData
    #     The single-cell dataset. LID will be computed for each cell.
    #
    # pca_n_comps : int
    #     Number of PCA components to use. If 0, the function uses the raw
    #     expression matrix. Otherwise, PCA is computed and LID estimation is
    #     performed in the reduced space.
    #
    # k_grid : array-like or str
    #     A list/array of neighbor sizes k to test. If set to a string, the function 
    #     automatically builds a log-spaced range of k values (from ~3% of the dataset
    #     size up to a maximum of 2000 neighbors) to explore LID across small to large 
    #     neighborhood scales.
    #
    # save : bool
    #     Whether to save the resulting k-grid and LID matrix as .npy files.
    #
    # filename_tag : str
    #     Optional suffix to append to saved filenames, useful when running
    #     multiple experiments.
    #
    # Returns
    # -------
    # k_grid : ndarray
    #     The neighborhood sizes used in the computation.
    #
    # ids : ndarray (n_cells × n_k)
    #     Matrix of computed LID values. ids[i, j] is the LID of cell i using
    #     the j-th value of k in k_grid.
    # -------------------------------------------------------------------------


    if(pca_n_comps == 0):
        mtx = adata.X
    else:
        print(f"PCA dimensionality reduction: from {len(adata.var)} genes to {pca_n_comps} principal components.")
        sc.pp.pca(adata, n_comps = pca_n_comps)
        mtx = adata.obsm['X_pca']

    if(isinstance(k_grid, str)):
        n_neighbors = int(np.min([len(adata), 2000]))
        nk = 10
        print(f'Finding {n_neighbors} nearest neighbors...')
        neigh = NearestNeighbors(n_neighbors=int(n_neighbors)-1)
        neigh.fit(mtx)
        k_grid = np.logspace(np.log10(n_neighbors*0.03), 
                            np.log10(n_neighbors), 
                            num=nk).astype(int)
        
    else:
        print(f'Finding {int(np.max(k_grid))} nearest neighbors...')
        neigh = NearestNeighbors(n_neighbors=int(np.max(k_grid)-1))
        neigh.fit(mtx)
    
    print('Nearest-neighbors grid used to investigate ID across scales: ', k_grid)

    ids = []

    _, neigh_indexes = neigh.kneighbors()

    for sc_i, sc_neigh_indexes in enumerate(tqdm(neigh_indexes)):

        sc_ids = []

        for k in k_grid:
            cell_mask = np.concatenate([[sc_i], sc_neigh_indexes[:k]])
            neigh_mtx = mtx[cell_mask, :]
            data = dadapy.Data(neigh_mtx)
            id, _, _ = data.compute_id_2NN()
            sc_ids.append(id)
        
        ids.append(sc_ids)

    if(save):
        save_npy(k_grid, filename = f'optimal_scale_kgrid_out{filename_tag}')
        save_npy(ids, filename = f'optimal_scale_ids_out{filename_tag}')
    
    return k_grid, np.array(ids)

def optimal_scale(k_grid, ids, fontsize = 15):

    # Plots the mean Local Intrinsic Dimension (LID) across cells as a function
    # of neighborhood size k, together with the 25–75% quantile bands. The
    # function also computes and plots the numerical derivative of the mean LID
    # to highlight the scale at which LID becomes stable. The optimal k is
    # estimated as the point where the derivative crosses zero, indicating 
    # the transition to a stable LID regime.
    #
    # Parameters
    # ----------
    # k_grid : array-like
    #     Values of k (neighborhood sizes) used in LID estimation.
    #
    # ids : ndarray
    #     Matrix of LID values with shape (n_cells × n_k), where each column
    #     corresponds to a k in k_grid.
    #
    # fontsize : int, optional
    #     Font size used for plotting labels and tick marks.
    #
    # Returns
    # -------
    # vline_x : float
    #     The optimal k
    #
    # -------------------------------------------------------------------------

    
    ids = np.array(ids)                    
 
    mask = ~np.isnan(ids).all(axis=0)
    ids = ids[:, mask]

    y      = np.nanmean(ids, axis = 0)
    y_inf  = np.nanquantile(ids, 0.25, axis=0)
    y_sup  = np.nanquantile(ids, 0.75, axis=0)

    fig, ax = plt.subplots(2, 1, figsize = (10, 5), sharex = True)

    ax[0].plot(k_grid, y, '-', color = 'black', linewidth = 3)
    ax[0].plot(k_grid, y, 'o', markersize = 4, color = 'black')
    ax[0].plot(k_grid, y_inf, '--', color = 'black', linewidth = 2)
    ax[0].plot(k_grid, y_sup, '--', color = 'black', linewidth = 2)
    ax[0].set_ylabel('Mean LID \nacross cells', fontsize = fontsize)
    ax[0].grid(True, 'both', linestyle = '--')
    ax[0].tick_params(labelsize = fontsize)

    increments = np.diff(y)/np.diff(k_grid)
    x_increment = k_grid[:-1] + np.diff(k_grid)/2
    mask = ~np.isnan(increments)
    increments = increments[mask]
    x_increment = x_increment[mask]
    ax[1].plot(x_increment, increments, '-', color = 'black', linewidth = 3)
    ax[1].plot(x_increment, increments, 'o', markersize = 4, color = 'black')
    ax[1].set_xlabel('Neighborhood size k', fontsize=fontsize)
    ax[1].set_ylabel('Derivative \nof Mean LID', fontsize = fontsize)
    ax[1].grid(True, 'both', linestyle = '--')
    ax[1].tick_params(labelsize = fontsize)
    ax[1].hlines(y = 0, xmin = np.min(x_increment), 
                xmax = np.max(x_increment), linestyles = '-',
                linewidth = 2, color = 'blue')

    vline_x = np.array([increments[i+1]*increments[i] for i in range(len(increments)-1)])
    vline_x = x_increment[np.where(vline_x<0)[0][0]]
    
    ax[1].vlines(x = vline_x, ymin = np.min(increments), 
                ymax = np.max(increments), linestyles = '-',
                linewidth = 2, color = 'green', label =f"x={vline_x}")
    ax[1].legend(fontsize = 0.8*fontsize, loc = 'best')

    fig.subplots_adjust(hspace = 0)
    plt.show()

    return vline_x

def save_npy(vals, filename = 'optimal_scale_out'):
    temp_file_name_pickle = f'temp.npy'
    final_file_name_pickle = f'{filename}.npy'
    try:
        with open(temp_file_name_pickle, 'wb') as f:
            pickle.dump(vals, f)
        shutil.move(temp_file_name_pickle, final_file_name_pickle)
    except Exception as e:
        print(f"Error Pickle: {e}")


def load_npy(filename = 'optimal_scale_out'):

    final_file_name_pickle = f'{filename}.npy'
    vals = None
    try:
        with open(final_file_name_pickle, 'rb') as f:
            vals = pickle.load(f)
        print(f"Successfully loaded from: {final_file_name_pickle}")
    except Exception as e:
        print(f"ERROR loading: {e}")
    
    return vals