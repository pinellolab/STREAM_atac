import numpy as np
import pandas as pd
import anndata as ad
import multiprocessing
import os
import gzip
from sklearn import preprocessing
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

os.environ['KMP_DUPLICATE_LIB_OK']='True'

def read(file_name,file_path='',file_format='tsv',delimiter='\t',experiment='rna-seq', workdir=None,**kwargs):
    """Read gene expression matrix into anndata object.
    
    Parameters
    ----------
    file_name: `str`
        File name. For atac-seq data, it's the z-score matrix file name.
    file_path: `str`, optional (default: '')
        File path. By default it's empty
    file_format: `str`, optional (default: 'tsv')
        File format. currently supported file formats: 'tsv','txt','tab','data','csv','mtx','h5ad','pklz','pkl'
    delimiter: `str`, optional (default: '\t')
        Delimiter to use.
    experiment: `str`, optional (default: 'rna-seq')
        Choose from {{'rna-seq','atac-seq'}}       
    workdir: `float`, optional (default: None)
        Working directory. If it's not specified, a folder named 'stream_result' will be created under the current directory
    **kwargs: additional arguments to `Anndata` reading functions
   
    Returns
    -------
    AnnData object
    """       
    _fp = lambda f:  os.path.join(file_path,f)
    if(file_format == 'pkl'):
        if file_name.split('.')[-1]=='gz':
            f = gzip.open(_fp(file_name), 'rb')
            adata = pickle.load(f)
            f.close() 
        else:
            f = open(_fp(file_name), 'rb')
            adata = pickle.load(f)
            f.close()    
    elif(file_format == 'pklz'):
        f = gzip.open(_fp(file_name), 'rb')
        adata = pickle.load(f)
        f.close()
    else:
        if(experiment not in ['rna-seq','atac-seq']):
            print('The experiment '+experiment +' is not supported')
            return         
        if(file_format in ['tsv','txt','tab','data']):
            adata = ad.read_text(_fp(file_name),delimiter=delimiter,**kwargs).T
            adata.raw = adata        
        elif(file_format == 'csv'):
            adata = ad.read_csv(_fp(file_name),delimiter=delimiter,**kwargs).T
            adata.raw = adata
        elif(file_format == 'mtx'):
            adata = ad.read_mtx(_fp(file_name),**kwargs).T 
            adata.X = np.array(adata.X.todense())
            print(_fp(os.path.join(os.path.dirname(file_name),'genes.tsv')))
            genes = pd.read_csv(_fp(os.path.join(os.path.dirname(file_name),'genes.tsv')), header=None, sep='\t')
            adata.var_names = genes[1]
            adata.var['gene_ids'] = genes[0].values
            print(_fp(os.path.join(os.path.dirname(file_name),'barcodes.tsv')))
            adata.obs_names = pd.read_csv(_fp(os.path.join(os.path.dirname(file_name),'barcodes.tsv')), header=None)[0]
            adata.raw = adata
        elif(file_format == 'h5ad'):
            adata = ad.read_h5ad(_fp(file_name),**kwargs)
        else:
            print('file format ' + file_format + ' is not supported')
            return
        adata.uns['experiment'] = experiment        
    if(workdir==None):
        workdir = os.path.join(os.getcwd(), 'stream_result')
        print("Using default working directory.")
    if(not os.path.exists(workdir)):
        os.makedirs(workdir)
    adata.uns['workdir'] = workdir
    print('Saving results in: %s' % workdir)
    return adata