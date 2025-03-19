# IDEAS

A package to estimate the Intrinsic Dimension of scRNA-seq data. 

## Notebooks  content:

- **Clustering_Preprocesing.ipynb** : analysis of the first dataset (Tran et al.) starting from the raw data (Tran.h5ad). The pipeline includes clustering, pluripotency score estimation, cluster ordering through pluripotency, pseudotime analsyis. It ends by saving two separate files (FBS_IDEAS.h5ad, A2S_IDEAS.h5ad) that contain in the obs and obsm the annotations resulting of the analysis, but have the raw counts in the adata.X.
___

- **FIGURE2.ipynb** : Produces figure 2 from the adata obtained in the "Clustering_Preprocesing.ipynb" notebook.

*Changes needed* : 
    - normalization in the ID trend plot?

___

- **Local_ID.ipynb** : local ID estimation analysis of the first dataset (Tran et al.), the third dataset (Nair et al.) and the Xenopus dataset, which is not in the article (yet). The local ID (for a certain number of neighbors chosen) is added as a obs column (es .adata.obs[local_ID_100] for 100 neighbors), having for each cell the local ID estimation obtained by using that cell as root cell. 
    - Tran : takes in input the preproccesed data  (FBS_IDEAS.h5ad, A2S_IDEAS.h5ad) and retrieves the same data with the local ID annotation in the obs (FBS_to_plot.h5ad, A2S_to_plot.h5ad). 
    - Nair : Adds the local ID annotaion in the obs (Nair_local_ID.h5ad). In addition to the local ID analysis, the pseudotime is calculated on the three different trajectories within the data (fibroblast-like, partial reprogramming, reprogramming). The final data with annotiation on both local ID and pseudotime is saved separately for the three trajectories (Nair_FIBRO_to_plot.h5ad, Nair_PARTIAL_to_plot.h5ad, Nair_REP_to_plot.h5ad).
    - Xenopus: The local ID analysis is performes separately on the two experimental conditions : in vitro fertilization (IVF) and nuclear transfer (NT). Additionally, the local ID is compared between the inner and outer layer both in the NT and in the IVF condition.

___

- **FIGURE_Local_ID.ipynb** : produces the figures from the adata obtained in "Local_ID.ipynb" for both the Tran et al. dataset and the Nair et al. dataset. It includes umap visualization of local ID, pseudotime, pluripotency score and scatter plots of the trends of local ID vs pseudotime and local ID vs pluripotency. 

## Installation

To install the package, we recommend that you first create a new environment

```bash 
conda env create  --file IDEAS_environment.yml 
```

Or if you use pip you can use 

```bash 
pip install -r requirements.txt
```

Note that the current version of python used for this repository is 3.11.6

Then you need to activate the environment

```bash 
conda activate ideas-env
```

And clone this repository by 

```bash
git clone https://github.com/maddalenastn/IDEAS.git
```

Then go to the main folder

```bash
cd IDEAS/
```

And  you can finally build the package by:

```bash
pip install -e .
```
