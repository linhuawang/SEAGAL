# import scanpy as sc
# import anndata as ad
# import squidpy as sq
# import pandas as pd
# import seaborn as sns
# import numpy as np
import warnings
from matplotlib import pyplot as plt
warnings.filterwarnings("ignore")
from ._utils import *
from ._gmod import *
from ._plotting import *

class SEAGAL(object):
    """ Spatial Enrichment Analysis of Gene Association using L-index.
    """
    def __init__(self, 
                 count_path = "", meta_path = "",
                 visium_path="/mnt/atlas_local/linhua/data/SEAGAL/test_data/"):
        


        data = load_raw(self, count_path, meta_path, visium_path)
            
        ## preprocess, and QC
        data = process_st(self, data, min_counts=150, min_cells=10)
        ## spatial process
        self.adata = spatial_process(data)
        self.co_expression = None
        self.co_expression_grouped = None
        self.hotspot = None
        self.hotspot_names = None
        # self.plot =  plot()






