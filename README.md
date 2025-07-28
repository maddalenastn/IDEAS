# IDEAS

A package to estimate the Intrinsic Dimension of scRNA-seq data. 
This repository provides code, data, and tutorials for reproducing the main analyses presented in our article [Staiano et al.](https://doi.org/10.1101/2025.07.21.665922). The core of the repository is **IDEAS** (Intrinsic Dimensionality Estimation Analysis of Single-cell RNA sequencing data), which enables the study of cell potency in single-cell transcriptomic data by estimating their intrinsic dimension (ID). 


## Repository Structure:

- **`computeID/`**  
  Source code of the IDEAS pipeline.

- **`Tutorial/`**  
  A Jupyter Notebook introducing how to use IDEAS on a simple dataset (Dataset 1 – FBS condition). 

- **`Notebooks/`**  
Jupyter Notebooks for reproducing the main figures of the article.  
Each notebook is named according to the corresponding figure. `Local_ID.ipynb` computes local intrinsic dimensionality used for **Figure 5**. All notebooks require an `anndata` object created as described in the [Instructions for reproducing figures from datasets](#-datasets--reproduction-instructions) section.

- **`Metadata/`**  
  Metadata extracted from open access datasets. These correspond to the `.obs` attribute of the `anndata` objects and include clustering results, pluripotency scores, and pseudotime annotations.

---

## Input and output

| Input Name     | Type             | Default                  | Description |
|----------------|------------------|--------------------------|-------------|
| `adata`        | AnnData object   | —                        | scRNA-seq data matrix with annotations. |
| `group`        | str              | `'all'`                  | Name of `.obs` column to group cells for ID computation. If `'all'`, all cells are treated as one group. |
| `method`       | str              | `'2nn'`                  | Method for computing intrinsic dimension. Options: `'pca'`, `'local_pca'`, `'2nn'`, `'local_2nn'`. |
| `roots`        | list             | `None`                   | List of cell names used as roots for local methods. (`'local_2nn'`, `'local_pca'`). If `None`, all cells are used. Ignored for non-local methods (`'2nn'`, `'pca'`). |
| `variance_e`   | float            | `0.9`                    | Percentage of explained variance for PCA-based methods. Ignored if method is not `'pca'` or `'local_pca'`. |
| `sample_size`  | int              | 80% of the smallest group| Number of cells per random sub-sampling. Ignored for local methods. |
| `n_samples`    | int              | `30`                     | Number of random sub-samplings used to estimate ID. Ignored in local methods. |
| `n_neighbors`  | int              | 10% of total number of cells  | Number of neighbors for local sub-sampling. Ignored for non-local methods (`'2nn'`, `'pca'`).  |
| `normalization`| bool             | `True`                   | If `True`, performs total-count normalization. |
| `norm_sum`     | int              | `1e6`                    | Target sum for normalization if `normalization=True`. |
| `full_output`  | bool             | `True`                   | Determines the output formatting. Ignored for local methods. |
| `id_score`     | bool             | `False`                  | If `True`, rescales ID values between 0 and 1. |


| Method Type | `full_output` | Return Type       | Description |
|-------------|---------------|-------------------|-------------|
| `'2nn'`, `'pca'` | `False`        | `pandas.DataFrame` (2 rows × n groups) | Returns mean and standard deviation of intrinsic dimension for each group. |
| `'2nn'`, `'pca'` | `True`         | `pandas.DataFrame` (n_samples × n groups) | Returns full matrix of IDs: each row is a random sub-sample, each column is a group, values are ID estimates. |
| `'local_2nn'`, `'local_pca'` | *Ignored*      | `pandas.DataFrame` (n roots × 1) | Returns a single-column DataFrame: rows are root cell barcodes, values are local ID estimates. |

---

## Instructions for reproducing figures from datasets

To reproduce the figures from the article, follow the instructions below for each dataset:

1. Download raw data from the corresponding **GEO accession**.
2. Filter cells (`obs_names` in the anndata object) using `CellId` in the metadata CSV files.
3. Add the metadata to the obs of the anndata object.
4. Use the relevant notebook in the `Notebooks/` folder to generate the figure.

---

### Dataset 1: [Tran et al.](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)30529-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124719305297%3Fshowall%3Dtrue)
- **GEO Accession**: GSE108222
- **Metadata Files**:
  - `Metadata_FBS_fig2.csv`
  - `Metadata_A2S_fig2.csv`  
- **Resulting AnnData Files**:
  - `FBS_IDEAS.h5ad`
  - `A2S_IDEAS.h5ad`

`FBS_IDEAS.h5ad` is used in the **Tutorial notebook**.  
For **Figure 5**, use the `Local_ID.ipynb` notebook to compute local intrinsic dimension and add the results to the obs of the `anndata` object.  
Save the updated object and use it in the corresponding figure notebook.

---

### Dataset 2: [Francesconi et al.](https://elifesciences.org/articles/41627) 
- **GEO Accession**: GSE112004
- **Metadata File**: `Metadata_fig3.csv`  
- **Resulting AnnData File**: `Francesconi.h5ad`

---

### Dataset 3: [Nair et al.](https://pubmed.ncbi.nlm.nih.gov/37873116/)
- **GEO Accession**: GSE242424
- **Metadata File**: `Metadata_fig4.csv`  
- **Resulting AnnData File**: `Nair.h5ad`

---

### Dataset 4: [Zikmund et al.](https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(25)00051-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2213671125000517%3Fshowall%3Dtrue)
- **GEO Accession**: GSE269252 
- **Metadata Files**:
  - `Metadata_fig5_xenopusNT.csv`
  - `Metadata_fig5_xenopusIVF.csv`  
- **Resulting AnnData Files**:
  - `Xenopus_NT_to_plot.h5ad`
  - `Xenopus_IVF_to_plot.h5ad`

---


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
--- 

## Contact and Citation

If you use IDEAS in your research, please cite:   
🔗 [Staiano et al.](https://doi.org/10.1101/2025.07.21.665922)