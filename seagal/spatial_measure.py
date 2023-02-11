import numpy as np
import pandas as pd
import scanpy as sc
import random
import anndata as ad
import time
from .lee_vec import Spatial_Pearson, Spatial_Pearson_Local
from libpysal.weights import W


from scipy import stats
from statsmodels.stats.multitest import fdrcorrection


#=========================================================================================
# Graph/Embedding related measuring functions
#=========================================================================================


def Global_L(anndat_sp, features=None, permutations=0, percent=0.1, seed=1, max_RAM=16):
    random.seed(seed)
    np.random.seed(seed)

    print("Calculating global graph-based correlation value...")
    cong_mtx = anndat_sp.obsp['spatial_connectivities'].toarray()

    eSP = Spatial_Pearson(cong_mtx, permutations=permutations)
    index_df = pd.DataFrame({'index':np.arange(anndat_sp.n_vars)})
    index_df.index = anndat_sp.var.index
    
    if features is None:
        features = anndat_sp.var_names.to_numpy()
    features_1 = []
    features_2 = []
    for i in range(len(features)):
        for j in range(i+1, len(features)):
            features_1.append(features[i])
            features_2.append(features[j])

    features_1 = np.array(features_1)
    features_2 = np.array(features_2)

    coexp_df = pd.DataFrame({'gene_1': features_1,
                             'gene_2': features_2})

    # L for all pairs
    
    features_1_index = index_df.loc[features_1,'index'].to_numpy()
    features_2_index = index_df.loc[features_2,'index'].to_numpy()

    features_1_X = anndat_sp.X[:,features_1_index]
    features_2_X = anndat_sp.X[:,features_2_index]

    

    eSP = eSP.fit(features_1_X, features_2_X, percent=percent, seed=seed, max_RAM=max_RAM)
    coexp_df['L'] = eSP.association_
    coexp_df['L.p_value'] = eSP.significance_ if permutations else 1.0

    if permutations:
        _,coexp_df['L.FDR'] = fdrcorrection(coexp_df['L.p_value'],
                                            alpha=0.05, method='indep')
    else:
        coexp_df['L.FDR'] = 1.0

    #peaks_nearby['gene.pct'] = anndat_sp.var.loc[anndat_sp.var_names[gene_index],'Frac.all'].to_numpy()
    #peaks_nearby['peak.pct'] = anndat_sp.var.loc[anndat_sp.var_names[peak_index],'Frac.all'].to_numpy()
    anndat_sp.uns['co_expression'] = coexp_df

    # print("Following changes made to the AnnData object:")
    # print("\tGlobal L results updated in uns['co_expression'].")

    #global_L_df = pd.DataFrame({'L':eSP.association_,
    #                            'L.p_value': eSP.significance_ if permutations else 1.0})
    #global_L_df.index = peaks_nearby['genes'] + '_' + peaks_nearby['peaks']
    
    return anndat_sp




def Local_L(anndat_sp,
            var_1, var_2,
            dropout_rm=True,
            permutations=0, seed=1, max_RAM=32):
    random.seed(seed)
    np.random.seed(seed)
    start = time.time()

    print("Calculating local spatial correlation value...")
    
        
    L_mtx, L_mtx_names = _calculate_LL(anndat_sp, var_1, var_2,
                                       permutations=permutations, seed=seed, max_RAM=max_RAM)

    if dropout_rm:
        # print("Set L to 0 for cells with no expression on either feature of a certain pair...")
        GP_names = L_mtx_names
        Dropout_mtx = dropout_filter(anndat_sp, GP_names)
        L_mtx = L_mtx * Dropout_mtx

    anndat_sp.uns['Local_L'] = L_mtx
    anndat_sp.uns['Local_L_names'] = L_mtx_names

    #anndat_sp.uns['peaks_nearby'] = peaks_nearby_orig.copy()
    # print("Following changes made to the AnnData object:")
    # print("\tSpatial Pearson Correlation results saved in uns['Local_L']")
    # print("\tuns['co_expression']['Local_L'] added indicates feature selected for local L calculation or not.")

    print("Local L elapsed time: %.3fs"%(time.time()-start))
    return anndat_sp



def _calculate_LL(anndat_sp,
                  var_1, var_2,
                  permutations=0, seed=1, max_RAM=16):
    random.seed(seed)
    np.random.seed(seed)

    #anndat_sp.var['index'] = np.arange(anndat_sp.n_vars)
    
    index_df = pd.DataFrame({'index':np.arange(anndat_sp.n_vars)})
    index_df.index = anndat_sp.var.index
    
    #if features is None:
    #    features = anndat_sp.var_names.to_numpy()
    features_1 = [var_1]
    features_2 = [var_2]
    pair_names = [str(var_1) + '_' + str(var_2)]
    '''
    for i in range(len(features)):
        for j in range(i+1, len(features)):
            features_1.append(str(features[i]))
            features_2.append(str(features[j]))
            pair_names.append(str(features[i]) + '_' + str(features[j]))
    '''

    features_1 = np.array(features_1)
    features_2 = np.array(features_2)

    features_1_index = index_df.loc[features_1,'index'].to_numpy()
    features_2_index = index_df.loc[features_2,'index'].to_numpy()

    features_1_X = anndat_sp.X[:,features_1_index]
    features_2_X = anndat_sp.X[:,features_2_index]

    cong_mtx = anndat_sp.obsp['spatial_connectivities'].toarray()

    eSP2 = Spatial_Pearson_Local(cong_mtx, permutations=permutations)
    eSP2 = eSP2.fit(features_1_X, features_2_X, seed=seed, max_RAM=max_RAM)
    local_L_df = pd.DataFrame(eSP2.associations_)
    local_L_df.columns = pair_names

    if permutations:
        local_p_df = pd.DataFrame(eSP2.significance_)
        local_p_df.columns = pair_names

    L_mtx = local_L_df.to_numpy()

    return L_mtx, local_L_df.columns.to_numpy()



def dropout_filter(anndat_sp, GP_names):
    index_df = pd.DataFrame({'index':np.arange(len(anndat_sp.var_names))})
    index_df.index = anndat_sp.var_names
    GP_G = [gp.split('_')[0] for gp in GP_names]
    GP_P = [gp.split('_')[1] for gp in GP_names]
    GP_genes_index = index_df.loc[GP_G,:]['index'].to_numpy()
    GP_peaks_index = index_df.loc[GP_P,:]['index'].to_numpy()
    
    E_mtx = anndat_sp.X
    #E_mtx_dropout_value = E_mtx.min(axis=0)
    #E_mtx_dropout_value = E_mtx_dropout_value[np.newaxis,:][np.zeros(anndat_sp.n_obs).astype(int)] #.shape
    E_mtx_dropout_value = np.zeros(E_mtx.shape)
    #Dropout_mtx = (E_mtx_dropout_value != E_mtx).astype(int)
    Dropout_mtx = (~np.isclose(E_mtx_dropout_value,E_mtx, 1e-3)).astype(int)
    #print(Dropout_mtx)
    Dropout_mtx_G = Dropout_mtx[:,GP_genes_index]
    Dropout_mtx_P = Dropout_mtx[:,GP_peaks_index]
    
    Dropout_mtx = Dropout_mtx_G * Dropout_mtx_P
    
    return Dropout_mtx
