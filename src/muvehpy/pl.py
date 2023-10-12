from scipy import sparse as sp
from scipy import io
import scanpy as sc
import pandas as pd
import numpy as np

def muveh_plot(adata, variants = [], plot_type = 'genotype', incl_oth = False, genomic_names = False, **kwargs):
  variants_sel = adata.uns['muveh_lut'][adata.uns['muveh_lut'].isin(variants).any(axis=1)].drop('variant_name', axis = 1).drop_duplicates()

  variants_sel = variants_sel[variants_sel['variant'].isin(adata.uns['muveh_variants'].index)]
  
  variant_indices = adata.uns['muveh_variants'].index.get_indexer(variants_sel['variant'])
  
  if genomic_names is True:
    col_names = variants_sel['variant']
  else:
    col_names = variants_sel['plot_name']
  
  if plot_type == 'genotype':
    sc.pl.umap(sc.AnnData(obs = pd.concat([adata.obs, pd.DataFrame.sparse.from_spmatrix(pp.muveh_score_mut(adata)[:,variant_indices], \
                          index = adata.obs_names, columns = col_names).astype(str)], axis = 1), obsm = adata.obsm), \
               color = col_names, palette = {'0' : 'gray', '1' : 'yellow', '2' : 'red'}, **kwargs)
  elif plot_type == 'genotype3':
    sc.pl.umap(sc.AnnData(obs = pd.concat([adata.obs, pd.DataFrame.sparse.from_spmatrix(pp.muveh_score(adata)[:,variant_indices], \
                          index = adata.obs_names, columns = col_names).astype(str)], axis = 1), obsm = adata.obsm), 
               color = col_names, palette = {'0' : 'gray', '1' : 'yellow', '2' : 'red', '3' : 'purple'}, **kwargs)
  elif plot_type == 'depth':
    sc.pl.umap(sc.AnnData(obs = pd.concat([adata.obs, pd.DataFrame.sparse.from_spmatrix(pp.muveh_depth(adata, \
                                                                                        incl_oth = incl_oth)[:,variant_indices], \
                          index = adata.obs_names, columns = col_names)], axis = 1), obsm = adata.obsm), \
               color = col_names, **kwargs)
  else:
    print("plot_type must be one of genotype, geontype3, or depth")
