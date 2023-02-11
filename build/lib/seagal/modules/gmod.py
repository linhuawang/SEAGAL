from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import pyplot as plt

def genemodules(sg, nmax=6, use_grouped=False, n_modules = None):
    if use_grouped:
        df1 = sg.co_expression_grouped.copy()
    else:
        df1 = sg.co_expression.copy()

    df2 = df1.copy()
    df2.rename(columns={'gene_1': 'gene_2', 'gene_2': 'gene_1'}, inplace=True)

    coexpr = pd.concat([df1, df2])
    coexpr = coexpr.pivot(index='gene_1', columns='gene_2', values='L')
    coexpr.index.name=None
    coexpr.columns.name=None
    np.fill_diagonal(coexpr.values, 0)
    X = coexpr.values.copy()

    if n_modules is None:
        range_n_clusters = list(range(2, nmax+1))
        sil_scores = []

        for n_clusters in range_n_clusters:
            hc = AgglomerativeClustering(n_clusters=n_clusters)
            cluster_labels = hc.fit_predict(X)
            silhouette_avg = silhouette_score(X, cluster_labels)
            sil_scores.append(silhouette_avg)
            print(
                "For n_clusters =",
                n_clusters,
                "The average silhouette_score is :",
                silhouette_avg,
            )
        sil_scores_r = sil_scores[::-1]
        i = len(sil_scores_r) - np.argmax(sil_scores_r) - 1
        nclust_opt = range_n_clusters[i]
    else:
        nclust_opt = n_modules
    
    print("Number of selected clusters = ", nclust_opt)

    hc_opt = AgglomerativeClustering(n_clusters=nclust_opt)
    cluster_labels_opt = hc_opt.fit_predict(X)
    silhouette = silhouette_score(X, cluster_labels_opt)
    cluster_labels_opt = [f'm{c}' for c in cluster_labels_opt]
    module2Genes = {}
    gene2module = {}

    for gene, lab in zip(coexpr.index.tolist(), np.ravel(cluster_labels_opt).tolist()):
        if lab not in module2Genes.keys():
            module2Genes[lab] = [gene]
        else:
            module2Genes[lab].append(gene)
        
        gene2module[gene] = lab

    module_df = pd.DataFrame({"label": cluster_labels_opt, "silhouette": silhouette}, 
                                    index=coexpr.index.copy())
    sg.gene_modules = module_df
    
    colors = {'m0': "green", 'm1': "magenta", 'm2': "darkorange",
                 'm3': "gray", 'm4': "beige", 'm5': "pink",
                  'm6': "darkgray", 'm7': "yellow",
                  'm8' : "red", 'm9': "blue"}

    g = sns.clustermap(coexpr, figsize=(6,6), cmap='bwr',
                vmin=-0.8, vmax=0.8,yticklabels=1, xticklabels=1,
                row_colors=[colors[l] for l in cluster_labels_opt])
    
    uniq_modules = list(set(cluster_labels_opt))
    handles = [Patch(facecolor=colors[l]) for l in uniq_modules]
    plt.legend(handles, uniq_modules, title='Gene Modules',
            bbox_to_anchor=(1.2, 1), bbox_transform=plt.gcf().transFigure,
             loc='upper right')

    mod_dict = {'gene2mod': gene2module, 'mod2gene': module2Genes, 'module_df': module_df, 'module_clustmap': g}
    sg.module_dict = mod_dict