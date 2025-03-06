# IDEAS

A package to estimate the Intrinsic Dimension of scRNA-seq data. 

Notebooks content: 

-"Clustering_Preprocesing": analysis of the first dataset (Tran et al.) starting from the raw data (Tran.h5ad). The pipeline includes clustering, pluripotency score estimation, cluster ordering through pluripotency, pseudotime analsyis. It ends by saving two separate files (FBS_IDEAS.h5ad, A2S_IDEAS.h5ad) that contain in the obs and obsm the result of the analysis, but have the raw counts in the adata.X.

-""

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
