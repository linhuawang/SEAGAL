import os
import scanpy as sc
import sys
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import squidpy as sq
from sklearn.mixture import GaussianMixture 
from .spatial_measure import Global_L, Local_L
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from pathlib import Path

def load_raw(self, count_path, meta_path , visium_path):
    if os.path.exists(visium_path):
            self.visium_path = visium_path
            data = sc.read_visium(visium_path)
            self.data_type = 'visium'
    elif os.path.exists(count_path) and os.path.exists(meta_path):
        self.data_type = 'st'
        count = pd.read_csv(count_path, index_col=0)
        meta = pd.read_csv(meta_path, index_col=0)
        
        var_df = pd.DataFrame({"gene": count.columns},
                                    index = count.columns)
        
        if 'x' not in meta.columns or 'y' not in meta.columns:
            sys.exit("x and y should be in the columns of the meta data frame.")
            
        if (meta.index != count.index).any():
            sys.exit("count and meta data frames should have the same index and order.")
            
        data = ad.AnnData(X=csr_matrix(count.values),
                                var=var_df, obs=meta)
        data.obsm = {'spatial': meta[['x', 'y']].values.copy()}
    else:
        sys.exit("Please either provide (1) a visium path or (2) paths to count matrix and meta file.")

    return data

def process_st(self,data, min_counts=150, min_cells=10):
    data.var_names_make_unique()
    data.var_names = data.var_names.to_series().apply(lambda x: str(x).upper())
    sc.pp.filter_cells(data, min_counts=min_counts)
    sc.pp.filter_genes(data, min_cells=min_cells)
    self.raw_adata = data.copy()
    sc.pp.normalize_total(data, target_sum=10**6)
    sc.pp.log1p(data)
    return data

def spatial_process(data):
    sq.gr.spatial_neighbors(data)
    sq.gr.spatial_autocorr(data,
            mode="moran",
            n_perms=1,
            n_jobs=1)
    data.var['moranI'] = data.uns['moranI']['I'].copy()
    return data

def spatial_pattern_genes(self, I=None, topK=None):
    adata = self.adata.copy()
    if I is None and topK is None:
        GM = GaussianMixture(n_components=2, covariance_type='spherical')
        GM.fit(np.reshape(adata.var.moranI.to_numpy(), (adata.shape[1], 1)))
        vals = np.arange(0, 1, 0.02)
        vals = np.reshape(vals, (len(vals), 1))
        pred_labs = GM.predict(vals)
        l0 = pred_labs[0]
        I = vals[np.where(pred_labs != l0)[0][0]][0]
    elif topK is not None:
        var_df = adata.var.copy()
        var_df = var_df.sort_values("moranI", ascending=False)
        I = var_df.moranI.tolist()[topK-1]

    adata.var['high_pattern_genes'] = False        
    adata.var.loc[(adata.var.moranI >= I), 'high_pattern_genes'] = True
    print(f"#Extracted genes with spatial patterns = {adata.var.high_pattern_genes.sum()}")
    f, ax = plt.subplots(1,1, figsize=(4,4))
    sns.histplot(data=adata.var, x='moranI',
                    hue='high_pattern_genes',
                    palette={True: "red", False: "lightgray"},ax=ax)
    plt.close()
    self.histI = f
    self.adata = adata
    return f

def spatial_association(self, grouped_only=True, use_pattern_genes=True,
                        genes=None, n_permutation=99, permute_ratio=0.2,
                        FDR_cutoff = 0.05, L_cutoff = 0.1):
    """ Calculate spatial associations
    
        Prameters:
        ----------
        grouped_only: bool, if True, only use it for gene sets previously defined and other parameters are irrelevant.
        use_pattern_genees: bool, if True, only calculate L for genes that have spatial patterns (based on Moran's I).
        genes: list of gene names, if not None, calculate L for this list of genes only.
        n_permutation: int, default 99, number of permutations for L's p value calculation.
        FDR_cutoff: float, 0 to 1, default 0.05, a FDR cutoff to define significance for L values.
        L_cutoff: float, a effect-size cutoff to define significance for L values.
        
        Return:
        ----------
        The spatial correlation data frame will be saved to self.co_expression.   
    """
    if grouped_only:
        grouped_data = self.grouped_adata.copy()
        grouped_data = Global_L(grouped_data, permutations=n_permutation, percent=permute_ratio, max_RAM=32)
        df = grouped_data.uns['co_expression'].copy()
        df['pair'] = df.gene_1 + "&" + df.gene_2
        vals  = df['L.FDR'].copy()
        vals = vals[vals > 0]
        if len(vals) > 0:
            nonzeromin = min(vals)
            df['-log10(FDR)'] = -np.log10(df['L.FDR'] + nonzeromin)  
        else:
            df['-log10(FDR)'] = -2
        df["Association"] = "NS"
        df.loc[(df['L.FDR'] <= FDR_cutoff) & (df['L'] >= L_cutoff), "Association"] = "SigPos"
        df.loc[(df['L.FDR'] <= FDR_cutoff) & (df['L'] <= -L_cutoff), "Association"] = "SigNeg"
        self.co_expression_grouped = df
        return
    
    if use_pattern_genes:
        adata = self.adata[:, self.adata.var.high_pattern_genes].copy()         
    elif genes is not None:
        adata = self.adata[:, genes].copy()
    else:
        adata = self.adata.copy()

    try:
        adata.X = adata.X.toarray()
    except:
        pass
    varN = len(adata.var_names)
    print(f"#Genes to be used is {varN}. #Pairs = {int(varN * (varN-1) /2)}, Estimated elapsed time: {varN/1000} hours.")
    adata = Global_L(adata, permutations=n_permutation, max_RAM=32)
    df = adata.uns['co_expression'].copy()
    df['pair'] = df.gene_1 + "&" + df.gene_2

    vals  = df['L.FDR'].copy()
    vals = vals[vals > 0]
    if len(vals) > 0:
        nonzeromin = min(vals)
        df['-log10(FDR)'] = -np.log10(df['L.FDR'] + nonzeromin)  
    else:
        df['-log10(FDR)'] = -2
        
    df["Association"] = "NS"
    df.loc[(df['L.FDR'] <= FDR_cutoff) & (df['L'] >= L_cutoff), "Association"] = "SigPos"
    df.loc[(df['L.FDR'] <= FDR_cutoff) & (df['L'] <= -L_cutoff), "Association"] = "SigNeg"
    self.co_expression = df

def group_adata_by_genes(self, ct_dict = None, inplace=True):
    ## Overwrite self.co_expression_grouped object if existed.
    if self.co_expression_grouped is not None:
        self.co_expression_grouped = None

    data = self.raw_adata.copy()
    sc.pp.normalize_total(data, target_sum=10**6)
    pool_genes = data.var_names.tolist()

    if ct_dict is None:
        ct_dict = {'B cells': ['IGKC', 'EBF1', 'BANK1', 'CD79A',
                             'IGHD', 'IGHM', 'PAX5', 'IGLC2',
                             'FCHSD2', 'CD79B', 'RALGPS2', 'LY6D',
                              'MS4A1', 'MEF2C', 'IGLC3', 'BACH2', 'MAN1A',
                               'CD55', 'BTLA', 'AFF3'],
                    'DC': ['CD209A', 'CRIP1', 'CST3', 'CBFA2T3', 
                            'TCF4', 'VIM', 'PID1', 'IRF8', 'FLT3',
                             'TAGLN2', 'AHNAK', 'S100A4', 'IFI30',
                              'LGALS1', 'CCND1', 'SYNGR2', 'H2AFY', 
                              'BST2', 'LY6C2', 'CCDC88A'],
                    'Macrophages': ['APOE', 'C1QB', 'C1QC', 'C1QA', 
                                    'CCL4', 'RGS1', 'CD81', 'MS4A7',
                                     'CTSB', 'LGMN', 'CTSS', 'SELENOP', 
                                     'CD63', 'CD72', 'CCL3', 'CTSC', 
                                     'TMEM176B', 'FTL1', 'ACP5', 'AIF1'], 
                    'Monocytes': ['PLAC8', 'SLC8A1', 'FN1', 'S100A4', 
                                    'LYZ2', 'LY6C2', 'F13A1', 'IFITM3', 
                                    'CCR2', 'ADGRE4', 'CYBB', 'CRIP1',
                                     'LDLRAD3', 'PID1', 'GPR141', 'LRP1',
                                      'CCL9', 'ZEB2', 'AHNAK', 'IFITM6'], 
                    'Neutrophils': ['S100A8', 'S100A9', 'RETNLG', 'CSF3R',
                                     'IFITM1', 'SLPI', 'HDC', 'MMP9', 
                                     'CXCR2', 'PBX1', 'WFDC21', 'MXD1', 
                                     'CLEC4D', 'GDA', 'SLFN4', 'GSR',
                                      'S100A11', 'IL1B', 'GCNT2', 'MSRB1'],
                    'NK cells': ['GZMA', 'GZMB', 'CD7', 'CCL5', 
                                    'NKG7', 'AW112010', 'IL2RB', 'KLRB1C', 
                                    'KLRE1', 'KLRK1', 'CTSW', 'KLRA7', 
                                    'KLRD1', 'NCR1', 'PRF1', 'ARSB', 
                                    'LGALS1', 'ID2', 'CST7', 'CAR2'],
                    'T cells': ['GM2682', 'LEF1', 'SKAP1', 'THEMIS',
                                 'IL7R', 'ITK', 'PRKCQ', 'BCL11B',
                                  'MS4A4B', 'CD247', 'TCF7', 'CD3E',
                                   'CD3D', 'TRBC2', 'INPP4B', 'CD3G',
                                    'GRAP2', 'CAMK4', 'SATB1', 'SIDT1']}
    cts = list(ct_dict.keys())
    if len(cts) < 2:
        print("Please include at list 2 groups to compare against each other and call the function again.")
        return
    for key in ct_dict.keys():
        vals = [v.upper() for v in ct_dict[key] if v.upper() in pool_genes]
        ct_dict[key] = vals
        print(f"{key}: {len(vals)} markers: {','.join(vals)}.")
    
    X = np.zeros((data.shape[0], len(cts)))
    
    for i in range(len(cts)):
        ct = cts[i]
        ct_genes = ct_dict[ct]

        ct_expr = np.sum(data[:, ct_genes].X, axis=1)
        ct_expr_log1p = np.log1p(ct_expr)

        # data.obs.insert(loc=0, column=ct, value=ct_expr_log1p) 
        X[:,i:i+1] = ct_expr_log1p
        
    data = ad.AnnData(X=X,
                    obs=data.obs.copy(),
                    var=pd.DataFrame({'ct':cts}, index=cts),
                    uns = data.uns.copy(),
                    obsm = data.obsm.copy())
    sq.gr.spatial_neighbors(data)
    sq.gr.spatial_autocorr(data,
        mode="moran",
        n_perms=100,
        n_jobs=1)

    if inplace:
        self.ct_dict = ct_dict
        self.grouped_adata=data
        spatial_association(self, grouped_only=True)
    else:
        return data

def save_modules(self, path):
    self.module_dict['module_df'].to_csv(path)

def save_coexpr(self, path, use_grouped=False):
    if use_grouped:
        self.co_expression_grouped.to_csv(path)
    else:
        self.co_expression.to_csv(path)

def module_pattern(self, ncols=4, alpha_img=0.5, cmap='Reds'):
    mod2gene = self.module_dict['mod2gene'].copy()
    modules = list(mod2gene.keys())
    adata = group_adata_by_genes(self, ct_dict=mod2gene, inplace=False)

    adata.obs.rename(columns={"ct": "module"}, inplace=True)
    if self.data_type == 'visium':
        sc.pl.spatial(adata, color=modules, alpha_img=alpha_img, ncols=ncols, cmap=cmap)
    else:
        import math
        nrows = math.ceil(len(modules)/ncols)
        f, axs = plt.subplots(nrows, ncols,
         figsize=(ncols * 5, nrows * 4))
        for i in range(len(modules)):
            ax = axs.flatten()[i]
            im = ax.scatter(x=adata.obs.x.copy(),
                            y=adata.obs.y.copy(),
                            c=np.ravel(adata[:, modules[i]].X.toarray()), cmap=cmap)
            clb = f.colorbar(im, ax=ax)
            clb.ax.set_title('Gene Expr',fontsize=8)

            ax.axis("off")
            ax.set_title(modules[i])

        for ax in axs.flatten():
            ax.axis("off")



def module_hotspot(self, dropout_rm=True, alpha_img=0.5, cmap='Blues_r', ncols=4):
    mod2gene = self.module_dict['mod2gene'].copy()
    modules = list(mod2gene.keys())
    n_modules = len(modules)
    adata = group_adata_by_genes(self, ct_dict=mod2gene, inplace=False)
    adata.obs.rename(columns={"ct": "module"}, inplace=True)

    # f, axs = plt.subplots(n_modules,n_modules, figsize=(4*n_modules,4*n_modules))
    pairs = []
    for i in range(0, n_modules - 1):
        for j in range(i+1, n_modules):
            mod1, mod2 = modules[i], modules[j]
            adata1 = Local_L(adata, mod1, mod2, dropout_rm=dropout_rm, max_RAM=32)
            adata.obs[f'{mod1} & {mod2}'] = adata1.uns['Local_L'].ravel() 
            pairs.append(f'{mod1} & {mod2}') 

    if self.data_type == 'visium':
        sc.pl.spatial(adata, 
                    color=pairs, 
                    cmap=cmap,
                    alpha_img=alpha_img, 
                    ncols=ncols)
    else:
        import math
        f, axs = plt.subplots(math.ceil(len(pairs)/ncols), ncols,
                         figsize=(ncols * 5, math.ceil(len(pairs)/ncols)* 4))
        for i in range(len(pairs)):
            ax = axs.flatten()[i]
            im = ax.scatter(x=adata.obs.x.copy(),
                            y=adata.obs.y.copy(),
                            c=np.ravel(adata.obs[pairs[i]].to_numpy()), cmap=cmap)
            ax.axis("off")
            ax.set_title(pairs[i])
            clb = f.colorbar(im, ax=ax)
            clb.ax.set_title('Local L-index',fontsize=8)

        for ax in axs.flatten():
            ax.axis("off")