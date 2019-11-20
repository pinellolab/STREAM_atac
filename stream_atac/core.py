# Authors: Huidong Chen
# Contact information: huidong.chen@mgh.harvard.edu

import numpy as np
import pandas as pd
import anndata as ad
import multiprocessing
import os
import subprocess
import sys
import gzip
from sklearn import preprocessing
from rpy2.robjects.packages import importr
import conda.cli

os.environ['KMP_DUPLICATE_LIB_OK']='True'

def preprocess_atac(file_count,file_region,file_sample, genome = 'hg19',
                    feature = 'kmer',k=7,resize_peak=False,peak_width=450,n_jobs = 1,
                    file_format='tsv',file_path='',workdir=None,**kwargs):
    """Preprocess single cell atac-seq data and genearate a scaled z-score matrix.
    
    Parameters
    ----------
    file_count: `str`
        Count file name. A compressed matrix in sparse format (column-oriented).
    file_region: `str`
        Region file name. It has one column. Each row is a cell name (No header should be included). 
    file_sample: `str`
        Sample file name. It has three columns, i.e.,chromosome names, the start position of the region, the end position of the region (No header should be included).
    genome: `str`, optional (default: 'hg19')
        Reference genome. Choose from {{'mm9', 'mm10', 'hg38', 'hg19'}} 
    feature: `str`, optional (default: 'kmers')
        Features used to have the analysis. Choose from {{'kmer', 'motif'}} 
    k: `int`, optional (default: 7)
        kmer length.  
    resize_peak: `bool`, optional (default: False)
        Resize peaks when peaks have different widths.
    peak_width: `int`, optional (default: 450)
        Specify the width of peak when resizing them. Only valid when resize_peak is True.
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run
    file_format: `str`, optional (default: 'tsv')
        File format of file_count. Currently supported file formats: 'tsv','txt','csv','mtx'.
    file_path: `str`, optional (default: '')
        File path. By default it's empty
    experiment: `str`, optional (default: 'rna-seq')
        Choose from {{'rna-seq','atac-seq'}}       
    workdir: `float`, optional (default: None)
        Working directory. If it's not specified, a folder named 'stream_result' will be created under the current directory
    **kwargs: additional arguments to `Anndata` reading functions
   
    Returns
    -------
    AnnData object with the following fields:
    
    X : `numpy.ndarray` (`adata.X`)
        A #observations × #k-mers scaled z-score matrix.
    z_score: `numpy.ndarray` (`adata.layers['z_score']`)
        A #observations × #k-mers z-score matrix.
    atac-seq: `dict` (`adata.uns['atac-seq']`)   
        A dictionary containing the following keys:
        'count': (`adata.uns['atac-seq']['count']`), dataframe in sparse format, 
                the first column specifies the rows indices (the regions) for non-zero entry. 
                the second column specifies the columns indices (the sample) for non-zero entry. 
                the last column contains the number of reads in a given region for a particular cell.
        'region': (`adata.uns['atac-seq']['region']`), dataframe
                the first column specifies chromosome names.
                the second column specifies the start position of the region.
                the third column specifies the end position of the region.
        'sample': (`adata.uns['atac-seq']['sample']`), dataframe, the name of samples
    """
    
    if(genome not in ['hg19','hg38','mm9','mm10']):
        raise ValueError("Not supported reference genome: '%s'" % genome)
    if(feature not in ['kmer','motif']):
        raise ValueError("Not supported feature: '%s'" % feature)             
    dict_genome = {'hg38':'BSgenome.Hsapiens.UCSC.hg38',
                   'hg19':'BSgenome.Hsapiens.UCSC.hg19',
                   'mm9':'BSgenome.Mmusculus.UCSC.mm9',
                   'mm10':'BSgenome.Mmusculus.UCSC.mm10'}
    if(genome in ['hg19','hg38']):
        species = "Homo sapiens"
    if(genome in ['mm9','mm10']):
        species = "Mus musculus"

    print('Checking if required packages are installed ...')
    rbase = importr('base')
    if(not rbase.requireNamespace(dict_genome[genome], quietly = True)[0]):
        pkg = 'bioconductor-'+dict_genome[genome].lower()
        print("Installing pacakge '%s' ..." % pkg)
        conda.cli.main('conda', 'install',  '-y', pkg)
    if(feature == 'motif'):
        if(not rbase.requireNamespace('motifmatchr', quietly = True)[0]):
            pkg = 'bioconductor-motifmatchr'
            print("Installing pacakge '%s' ..." % pkg)
            conda.cli.main('conda', 'install',  '-y', pkg)
        if(not rbase.requireNamespace('JASPAR2016', quietly = True)[0]):
            pkg = 'bioconductor-jaspar2016'
            print("Installing pacakge '%s' ..." % pkg)
            conda.cli.main('conda', 'install',  '-y', pkg)
    
    if(workdir==None):
        workdir = os.path.join(os.getcwd(), 'stream_result')
        print("Using default working directory.")
    if(not os.path.exists(workdir)):
        os.makedirs(workdir)
    print('Saving results in: %s' % workdir) 
    
    _fp = lambda f:  os.path.join(file_path,f)
    print('Running chromVAR pipeline ...')
    _root = os.path.abspath(os.path.dirname(__file__))
    if(resize_peak):
        preprocess_cmd='Rscript ' + os.path.join(_root,'run_preprocess.R') + ' -c {0} -r {1} -s {2} -g {3} -f {4} -k {5:d} --species "{6}" --n_jobs {7:d} -o {8} --file_format {9} --peak_width {10} --resize_peak'.format(_fp(file_count),_fp(file_region),_fp(file_sample),dict_genome[genome],feature,k,species,n_jobs,workdir,file_format,peak_width)
    else:
        preprocess_cmd='Rscript ' + os.path.join(_root,'run_preprocess.R') + ' -c {0} -r {1} -s {2} -g {3} -f {4} -k {5:d} --species "{6}" --n_jobs {7:d} -o {8} --file_format {9}'.format(_fp(file_count),_fp(file_region),_fp(file_sample),dict_genome[genome],feature,k,species,n_jobs,workdir,file_format)
    # print(preprocess_cmd)
    code = subprocess.call(preprocess_cmd,shell=True)
    if(code!=0):
        print('preprocessing failed')
        sys.exit(1)    
    df_zscores = pd.read_csv(os.path.join(workdir, 'zscores.tsv.gz'),sep='\t',index_col=0)
    df_zscores_scaled = preprocessing.scale(df_zscores,axis=1)
    df_zscores_scaled = pd.DataFrame(df_zscores_scaled,index=df_zscores.index,columns=df_zscores.columns)
    df_zscores_scaled.to_csv(os.path.join(workdir,'zscores_scaled.tsv.gz'),sep = '\t',compression='gzip')
    adata = ad.AnnData(X=df_zscores_scaled.values.T, obs={'obs_names':df_zscores_scaled.columns},var={'var_names':df_zscores_scaled.index})
    adata.raw = adata
    adata.uns['workdir'] = workdir
    adata.uns['experiment'] = 'atac-seq'
    adata.layers["zscores"] = df_zscores.values.T
    adata.layers["zscores_scaled"] = df_zscores_scaled.values.T
    return adata