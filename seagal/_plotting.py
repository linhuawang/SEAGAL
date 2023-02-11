from .spatial_measure import Local_L
import scanpy as sc
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def hotspot(self, var_1, var_2, dropout_rm=True, use_grouped=False, cmap='viridis', alpha_img=0.5):
    """ Calculate local L-index for var pairs (var_1 and var_2). 
    """

    if use_grouped:
        adata = self.grouped_adata.copy()
    else:
        adata = self.adata.copy()
    # calculate local L-index for specific var_name in adata
    ### 
    try:
        adata.X = adata.X.toarray()
    except:
        pass

    if var_1 not in adata.var_names or var_2 not in adata.var_names:
        print("One of the queried features is not in the adata. Check if the correct `use_grouped` parameter is used, or typos in the names.")
        return

    adata1 = Local_L(adata, var_1, var_2, dropout_rm=dropout_rm, max_RAM=32)
    adata.obs[f'{var_1} & {var_2}'] = adata1.uns['Local_L'].ravel()
    f, ax = plt.subplots(1,1, figsize=(5,4))
    if self.data_type == 'visium':
        sc.pl.spatial(adata, color=f'{var_1} & {var_2}', ax=ax, cmap=cmap, alpha_img=alpha_img)
    else:        
        im=ax.scatter(x=adata.obs.x.copy(),
                    y=adata.obs.y.copy(),
                    c=np.ravel(adata.obs[f'{var_1} & {var_2}'].to_numpy()), cmap=cmap)
        ax.axis("off")
        ax.set_title(f'{var_1} & {var_2}')
        clb = f.colorbar(im, ax=ax)
        clb.ax.set_title('Local L-index',fontsize=8)

    if use_grouped:
        self.grouped_adata = adata.copy()
    else:
        self.adata = adata.copy()
    return f, adata.obs[f'{var_1} & {var_2}']

def volcano(self, FDR_cutoff=0.05, L_cutoff=0.1, use_grouped=False):
    if use_grouped:
        df = self.co_expression_grouped.copy()
    else:
        df = self.co_expression.copy()
    
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

    f, ax = plt.subplots(1,1, figsize=(4,3))
    sns.scatterplot(data=df, x='L', y='-log10(FDR)', 
                    hue='Association',
                        palette={'NS': 'lightgray', 'SigPos': 'red', 'SigNeg': 'blue'},
                        s=20, ax=ax)
    sns.despine(trim=True, ax=ax)
    return f

def clustermap(self, figsize=(6,6), use_grouped=False):
    if use_grouped:
        df1 = self.co_expression_grouped.copy()
    else:
        df1 = self.co_expression.copy()

    df2 = df1.copy()
    df2.rename(columns={'gene_1': 'gene_2', 'gene_2': 'gene_1'}, inplace=True)

    coexpr = pd.concat([df1, df2])
    coexpr = coexpr.pivot(index='gene_1', columns='gene_2', values='L')
    coexpr.index.name=None
    coexpr.columns.name=None
    np.fill_diagonal(coexpr.values, 0)
    if not use_grouped:
        sns.clustermap(coexpr, figsize=figsize, cmap='bwr', annot=False,
               vmin=-0.8, vmax=0.8)
    else:
        annot = coexpr.round(2).astype(str)
        for (gene1, gene2, p) in zip(df1.gene_1.tolist(), df1.gene_2.tolist(), df1['L.FDR'].tolist()):
            if p <= 0.01:
                if "*" not in annot.loc[gene1, gene2]:
                    annot.loc[gene1, gene2] = f"{annot.loc[gene1, gene2]}**"
                if "*" not in annot.loc[gene2, gene1]:
                    annot.loc[gene2, gene1] = f"{annot.loc[gene2, gene1]}**"
            elif p <= 0.05:
                if "*" not in annot.loc[gene1, gene2]:
                    annot.loc[gene1, gene2] = f"{annot.loc[gene1, gene2]}*"
                if "*" not in annot.loc[gene2, gene1]:
                    annot.loc[gene2, gene1] = f"{annot.loc[gene2, gene1]}*"
            else:
                pass
        sns.clustermap(coexpr, figsize=figsize, cmap='bwr', 
                annot=annot, fmt='s',
               vmin=-0.8, vmax=0.8)
    return