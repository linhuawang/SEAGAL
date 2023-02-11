# Spatial Enrichment Analysis of Gene Association using L-index

## Install

SEAGAL works for Python3.9, it is advised to create a conda environment with Python3.9 to avoid conflicts of dependencies. Conda installation guide can be found at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html. After installing conda, run the following command to create a `seagal` environment:

```bash
conda create -n seagal python=3.9
```

Then install seagal through:

```bash
conda activate seagal
conda install pip
pip install seagal
```

Please check if `echo $PATH` includes the path to python3.9 site-packages if `Error: module seagal not found` occurs when importing.

## Input data format
1. For 10X Visium, Space Ranger `Folder` with the following contents:
  - [Folder]/spatial/tissue_positions_list.csv
  - [Folder]/filtered_feature_bc_matrix.h5
2. CSV-format files: 
  - count.csv - gene expression data frame in Pandas.DataFrame format.
  - meta.csv - spot meta data frame, with x, y columns denoting coordinates.

## Running SEAGAL with Visium Data

Please read [Tutorial1](./test/Tutorial1%20SEAGAL%20with%20Visium%20(SpaceRanger%20Output).ipynb) for instructions.

## Running SEAGAL with user-defined marker genes

Please read [Tutorial2](./test/Tutorial3%20SEAGAL%20with%20ST%20(csv%20files).ipynb) for details.

## Running SEAGAL with CSV-format Data

Please read [Tutorial3](./test/Tutorial2%20SEAGAL%20with%20Customized%20Gene%20Sets.ipynb) for details.

<!-- ## Parameters 
  ### I/O parameters
  * -f: path to the input raw count matrix (csv).
  * -o: path to save the imputed data sets.

  ### Model Parameters
  * -r: radius in Euclidean distance to consider as adjacent spots.
  * -s: whether to select thresholding parameter epsilon automatically or not. 0: no selection, use fixed. 1: select automatically.
  * -e: edge filtering parameter epsilon, range from 0 to 1. Only useful when -s was set to 0.
  * -l: normalization method. Must be one of "cpm", "logCPM", "logMed", "none". Default is "cpm".

  ### Other parameters
  * -n: number of processors to be used for parallel computing. 1-10. Default is 1. 

## Run example experiments
  
  The following code will impute the test data with 4 processors, save the imputed cpm data, raw data to the designated folder. Also, the component information will be saved to the same folder.
  
    python3 MIST.py -f test_data/raw.csv -o test_data/imputed.csv -l cpm -n 4

  After running the above code, the following files will be generated:

    1. test_data/imputed.csv -- imputed, normalized, gene filtered expression.
    2. test_data/imputed_complete.csv -- imputed, normalized, gene expression.
    3. test_data/imputed_rawCount.csv -- imputed, raw gene counts.
    4. imputed_cluster_info.csv -- region assignment of every spot.

  ### Visualize major tissue components
  
  The following code will take component information returned by the imputation pipeline and visualize the component information.
  
    python3 visualize_components.py test_data/imputed_cluster_info.csv test_data/cluster.png
  
  The above code will visualize the detected regions by giving a figure like:

  [Cluster Visualization](test_data/output/cluster.png) -->
