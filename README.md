# IDEAS

A package to estimate the Intrinsic Dimension of scRNA-seq data. 

Notebooks content: 

- "Clustering_Preprocesing": analysis of the first dataset (Tran et al.) starting from the raw data (Tran.h5ad). The pipeline includes clustering, pluripotency score estimation, cluster ordering through pluripotency, pseudotime analsyis. It ends by saving two separate files (FBS_IDEAS.h5ad, A2S_IDEAS.h5ad) that contain in the obs and obsm the result of the analysis, but have the raw counts in the adata.X.
Changes needed: check/change pluripotency score definition; it containes all teh function to compute ID written in the notebbok, it would instead import the functions_IDEAS.py

- "Local_ID" : local ID estimation analysis of the first dataset (Tran et al.), the third dataset (Nair et al.) and the Xenopus dataset, which is not in the article (yet). THe local ID is added as a obs column, having for each cell the local ID estimation obtained bu using that cell as root cell.
    - Tran : takes in input the preproccesed data  (FBS_IDEAS.h5ad, A2S_IDEAS.h5ad) and retrieves the same data with the local ID annotation in the obs (FBS_to_plot.h5ad, A2S_to_plot.h5ad). 
    - Nair : Adds the local ID annotaion in the obs (Nair_local_ID.h5ad). In addition to the local ID analysis, the pseudotime is calculated on the three different trajectories within the data (fibroblast-like, partial reprogramming, reprogramming). The final data with annotiation on both local ID and pseudotime is saved separately for the three trajectories (Nair_FIBRO_to_plot.h5ad, Nair_PARTIAL_to_plot.h5ad, /Nair_REP_to_plot.h5ad).
    - Xenopus: still needs to be changed, ignore it for now

- "FIGURE_Local_ID" : produces the figures from the adata obtained in "Local_ID" for both the Tran et al. dataset and the Nair et al. dataset. It includes umap visualization of local ID, pseudotime, pluripotency score and scatter plots of the trends of local ID vs pseudotime and local ID vs pluripotency. 

- "FIGURE2" : Produces figure 2 from the adata obtained in the "Clustering_Preprocesing" notebook.

## Installazione

To install the package, we recommend that you first create a new environment

```bash 
conda create -n compute_id python 
```

Then you need to activate the environment


```bash 
conda activate compute_id
```

Then you need to clone this repository 

```bash
git clone https://github.com/maddalenastn/IDEAS.git
```

Then going inside the main folder

```bash
cd IDEAS/
```

And finally you can build the package:

```bash
pip install -e .
```
