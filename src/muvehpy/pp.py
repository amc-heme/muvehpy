from scipy import sparse as sp
from scipy import io
import scanpy as sc
import pandas as pd
import numpy as np

def read_muveh_mat(ad_matrix, dp_matrix, oth_matrix, colnames, rownames, cell_prefix = '', cell_suffix = ''):
  # read column and rownames
  muveh_cols = pd.read_csv(colnames, header=None, names=['variant'])
  muveh_rows = pd.read_csv(rownames, header=None, names=['cell'])
  
  # append capture to cell barcodes (rownames)
  muveh_rows.cell = cell_prefix+muveh_rows.cell+cell_suffix
  
  # read mm matrices
  muveh_ad = sp.csc_matrix(io.mmread(ad_matrix)).transpose()
  muveh_dp = sp.csc_matrix(io.mmread(dp_matrix)).transpose()
  muveh_oth = sp.csc_matrix(io.mmread(oth_matrix)).transpose()
  
  # suffix duplicate variant names
  dup_num =  muveh_cols.groupby('variant').cumcount().add(1).astype(str)
  muveh_cols['variant'] += ('.' + dup_num).replace('.1', '')
  
  return muveh_cols, muveh_rows, muveh_ad, muveh_dp, muveh_oth

def read_muveh_mat_config(config):
  # read config csv
  config_df = pd.read_csv(config, dtype='str')
  
  if len(config_df.index) == 1:
    muveh_cols, muveh_rows, muveh_ad, muveh_dp, muveh_oth = read_muveh_mat(config_df['ad_matrix'][0], \
                                                                           config_df['dp_matrix'][0], \
                                                                           config_df['oth_matrix'][0], \
                                                                           config_df['colnames'][0], \
                                                                           config_df['rownames'][0], \
                                                                           cell_prefix = config_df['cell_prefix'][0].strip('\''), \
                                                                           cell_suffix = config_df['cell_suffix'][0].strip('\''))
    
  else:
    # initialize variables as lists
    muveh_cols, muveh_rows, muveh_ad, muveh_dp, muveh_oth = ([] for i in range(5))
    
    # iterate through config samples
    for i in range(len(config_df.index)):
      # read column and rownames
      muveh_cols.append(pd.read_csv(config_df['colnames'][i], header=None, names=['variant']))
      muveh_rows.append(pd.read_csv(config_df['rownames'][i], header=None, names=['cell']))
      # append capture to cell barcodes (rownames)
      muveh_rows[i].cell = config_df['cell_prefix'][i].strip('\'')+muveh_rows[i].cell+config_df['cell_suffix'][i].strip('\'')
      # read mm matrices
      muveh_ad.append(sp.csc_matrix(io.mmread(config_df['ad_matrix'][i])).transpose())
      muveh_dp.append(sp.csc_matrix(io.mmread(config_df['dp_matrix'][i])).transpose())
      muveh_oth.append(sp.csc_matrix(io.mmread(config_df['oth_matrix'][i])).transpose())
      # suffix duplicate variant names
      dup_num =  muveh_cols[i].groupby('variant').cumcount().add(1).astype(str)
      muveh_cols[i]['variant'] += ('.' + dup_num).replace('.1', '')
      
    muveh_cols_complete = pd.concat(muveh_cols).drop_duplicates()
    
    for i in range(len(muveh_cols)):
      # construct dummy matrix for extra cells and bind to muveh_*
      pad_matrix = sp.csr_matrix((muveh_ad[i].get_shape()[0], len(muveh_cols_complete.variant) - muveh_ad[i].get_shape()[1]))
      muveh_ad[i] = sp.hstack([muveh_ad[i], pad_matrix])
      muveh_dp[i] = sp.hstack([muveh_dp[i], pad_matrix])
      muveh_oth[i] = sp.hstack([muveh_oth[i], pad_matrix])
      del pad_matrix
      
      # construct rownames for new matrices
      new_variants = muveh_cols_complete.variant[~muveh_cols_complete.variant.isin(muveh_cols[i].variant)]
      muveh_cols_i_complete = pd.concat([muveh_cols[i],new_variants.to_frame(name = "variant")])
      
      # reorder matrices and add to anndata
      muveh_ad[i] = muveh_ad[i][:,pd.Index(muveh_cols_complete.variant).get_indexer(muveh_cols_i_complete.variant)]
      muveh_dp[i] = muveh_dp[i][:,pd.Index(muveh_cols_complete.variant).get_indexer(muveh_cols_i_complete.variant)]
      muveh_oth[i] = muveh_oth[i][:,pd.Index(muveh_cols_complete.variant).get_indexer(muveh_cols_i_complete.variant)]
      
    muveh_ad = sp.vstack(muveh_ad)
    muveh_dp = sp.vstack(muveh_dp)
    muveh_oth = sp.vstack(muveh_oth)
    muveh_rows = pd.concat(muveh_rows)
    muveh_cols = muveh_cols_complete

  return muveh_cols, muveh_rows, muveh_ad, muveh_dp, muveh_oth

def add_muveh(h5ad, \
             config = None, \
             ad_matrix = None, \
             dp_matrix = None, \
             oth_matrix = None, \
             colnames = None, \
             rownames = None, \
             lut = None, \
             cell_prefix = '', \
             cell_suffix = ''):
  
  # read matrices
  if config is None:
    muveh_cols, muveh_rows, muveh_ad, muveh_dp, muveh_oth = read_muveh_mat(ad_matrix, \
                                                                           dp_matrix, \
                                                                           oth_matrix, \
                                                                           colnames, \
                                                                           rownames, \
                                                                           cell_prefix = cell_prefix, \
                                                                           cell_suffix = cell_suffix)
  else:
    muveh_cols, muveh_rows, muveh_ad, muveh_dp, muveh_oth = read_muveh_mat_config(config)
  
  # read anndata object
  adata = sc.read_h5ad(h5ad)
  
  # subset muveh matrices on adata
  muveh_rows_adata = muveh_rows[muveh_rows.cell.isin(adata.obs_names)]
  muveh_ad = muveh_ad[np.where(muveh_rows.cell.isin(adata.obs_names))]
  muveh_dp = muveh_dp[np.where(muveh_rows.cell.isin(adata.obs_names))]
  muveh_oth = muveh_oth[np.where(muveh_rows.cell.isin(adata.obs_names))]
  
  # construct dummy matrix for extra cells and bind to muveh_*
  pad_matrix = sp.csr_matrix((len(adata.obs_names) - muveh_ad.get_shape()[0], muveh_ad.get_shape()[1]))
  muveh_ad = sp.vstack([muveh_ad, pad_matrix])
  muveh_dp = sp.vstack([muveh_dp, pad_matrix])
  muveh_oth = sp.vstack([muveh_oth, pad_matrix])
  del pad_matrix
  
  # construct rownames for new matrices
  new_cells = adata.obs_names[~adata.obs_names.isin(muveh_rows_adata.cell)]
  muveh_rows_adata_complete = pd.concat([muveh_rows_adata,new_cells.to_frame(name = "cell")])
  
  # reorder matrices and add to anndata
  adata.obsm['X_muveh_ad'] = muveh_ad[pd.Index(adata.obs_names).get_indexer(muveh_rows_adata_complete.cell),:]
  adata.obsm['X_muveh_dp'] = muveh_dp[pd.Index(adata.obs_names).get_indexer(muveh_rows_adata_complete.cell),:]
  adata.obsm['X_muveh_oth'] = muveh_oth[pd.Index(adata.obs_names).get_indexer(muveh_rows_adata_complete.cell),:]
  
  # read LUT
  muveh_lut = pd.read_csv(lut)
  
  # add variants and lut to anndata object
  adata.uns['muveh_variants'] = muveh_cols.set_index('variant')
  adata.uns['muveh_lut'] = muveh_lut
  
  return adata

def muveh_score_mut(adata):
  return (adata.obsm['X_muveh_ad'] > 0).astype(int) + \
         (adata.obsm['X_muveh_dp'] > 0).astype(int)
   
def muveh_score(adata):
  return ((adata.obsm['X_muveh_ad'] > 0).astype(int) * 2) + \
        ((adata.obsm['X_muveh_dp'] - adata.obsm['X_muveh_ad']) > 0).astype(int)
        
def muveh_depth(adata, incl_oth = False):
  if incl_oth is True:
    return adata.obsm['X_muveh_dp'] + adata.obsm['X_muveh_oth']
  else:
    return adata.obsm['X_muveh_dp']
   