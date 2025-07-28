import pandas as pd
import numpy as np

def add_local_ID_to_adata(Adata, ID_dataframe, n_neigh=0, column_name = None):

    if n_neigh==0 : 
        n_neigh = 0.1 * len(Adata)
    
    if column_name is None:
        column_name = f'local_ID_{n_neigh}'

    if column_name not in Adata.obs.columns:
        Adata.obs[column_name] = pd.NA 

    index = [item[0] for item in ID_dataframe.index]
    for barcode,id_val in zip(index, ID_dataframe.ID.values):
        Adata.obs.loc[Adata.obs_names == barcode, column_name] = id_val

    Adata.obs[column_name] = Adata.obs[column_name].replace({pd.NA: np.nan}).astype(float)

    return Adata