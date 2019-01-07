import numpy as np
import pandas as pd
import anndata as ad
import multiprocessing
import os
import pickle
import gzip
from sklearn import preprocessing
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri



def read(file_name,file_name_sample=None,file_name_region=None,file_path='./',file_format='tsv',delimiter='\t',experiment='atac-seq', workdir=None,**kwargs):
    """Read gene expression matrix into anndata object.
    
    Parameters
    ----------
    file_name: `str`
        File name. For atac-seq data, it's the count file name.
    file_name_sample: `str`
        Sample file name. Only valid when atac_seq = True.
    file_name_region: `str`
        Region file name. Only valid when atac_seq = True.
    file_path: `str`, optional (default: './')
        File path. By default it's the current directory
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
    if(experiment == 'atac-seq'):
        if(file_format == 'pklz'):
            f = gzip.open(file_path+file_name, 'rb')
            adata = pickle.load(f)
            f.close()  
        elif(file_format == 'pkl'):
            f = open(file_path+file_name, 'rb')
            adata = pickle.load(f)
            f.close()
        else:            
            if(file_name_sample is None):
                print('sample file must be provided')
            if(file_name_region is None):
                print('region file must be provided')
            df_counts = pd.read_csv(file_name,sep='\t',header=None,names=['i','j','x'],compression= 'gzip' if file_name.split('.')[-1]=='gz' else None)
            df_regions = pd.read_csv(file_name_region,sep='\t',header=None,compression= 'gzip' if file_name_region.split('.')[-1]=='gz' else None)
            df_regions = df_regions.iloc[:,:3]
            df_regions.columns = ['seqnames','start','end']
            df_samples = pd.read_csv(file_name_sample,sep='\t',header=None,names=['cell_id'],compression= 'gzip' if file_name_sample.split('.')[-1]=='gz' else None)
            adata = ad.AnnData()
            adata.uns['atac-seq'] = dict()
            adata.uns['atac-seq']['count'] = df_counts
            adata.uns['atac-seq']['region'] = df_regions
            adata.uns['atac-seq']['sample'] = df_samples
    else:
        print('The experiment '+experiment +' is not supported')
        return        
    adata.uns['experiment'] = experiment
    if(workdir==None):
        if('workdir' in adata.uns_keys()):
            workdir = adata.uns['workdir']
        else:
            workdir = os.getcwd() + '/stream_result_atac/'
    if(not os.path.exists(workdir)):
        os.makedirs(workdir)
    adata.uns['workdir'] = workdir
    return adata


def counts_to_kmers(adata,k=7,n_jobs = multiprocessing.cpu_count()):
    """Covert counts files to kmer files.
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    k: `int`, optional (default: 7)
        k mer.  
    n_jobs: `int`, optional (default: all available cpus)
        The number of parallel jobs to run
        
    Returns
    -------
    updates `adata` with the following fields.
    
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
    chromVAR = importr('chromVAR')
    GenomicRanges = importr('GenomicRanges')
    SummarizedExperiment = importr('SummarizedExperiment')
    BSgenome_Hsapiens_UCSC_hg19 = importr('BSgenome.Hsapiens.UCSC.hg19')
    r_Matrix = importr('Matrix')
    BiocParallel = importr('BiocParallel')
    BiocParallel.register(BiocParallel.MulticoreParam(n_jobs))
    pandas2ri.activate()
    df_regions = adata.uns['atac-seq']['region']
    r_regions_dataframe = pandas2ri.py2ri(df_regions)
    regions = GenomicRanges.makeGRangesFromDataFrame(r_regions_dataframe)
    
    df_counts = adata.uns['atac-seq']['count']
    counts = r_Matrix.sparseMatrix(i = df_counts['i'], j = df_counts['j'], x=df_counts['x'])
    
    df_samples = adata.uns['atac-seq']['sample']
    samples = pandas2ri.py2ri(df_samples)
    samples.rownames = df_samples['cell_id']
    
    SE = SummarizedExperiment.SummarizedExperiment(rowRanges = regions,colData = samples,assays = robjects.ListVector({'counts':counts}))
    SE = chromVAR.addGCBias(SE, genome = BSgenome_Hsapiens_UCSC_hg19.BSgenome_Hsapiens_UCSC_hg19)
    
    # compute kmer deviations
    KmerMatch = chromVAR.matchKmers(k, SE, BSgenome_Hsapiens_UCSC_hg19.BSgenome_Hsapiens_UCSC_hg19)
    BiocParallel.register(BiocParallel.SerialParam())
    Kmerdev = chromVAR.computeDeviations(SE, KmerMatch)
    KmerdevTable = SummarizedExperiment.assays(Kmerdev)
    cn = pandas2ri.ri2py((Kmerdev.do_slot('colData')).do_slot('listData').rx2('cell_id'))
    rn = pandas2ri.ri2py(Kmerdev.do_slot('NAMES'))
    scores = pandas2ri.ri2py(KmerdevTable.do_slot('listData').rx2('deviations'))    

    df_zscores = pd.DataFrame(scores,index=rn,columns=cn)
    df_zscores_scaled = preprocessing.scale(df_zscores,axis=1)
    df_zscores_scaled = pd.DataFrame(df_zscores_scaled,index=df_zscores.index,columns=df_zscores.columns)
    adata_new = ad.AnnData(X=df_zscores_scaled.values.T, obs={'obs_names':df_zscores_scaled.columns},var={'var_names':df_zscores_scaled.index})
    adata_new.raw = adata_new
    adata_new.uns['workdir'] = adata.uns['workdir']
    adata_new.uns['experiment'] = adata.uns['experiment']
    adata_new.uns['atac-seq'] = dict()
    adata_new.uns['atac-seq']['count'] = df_counts
    adata_new.uns['atac-seq']['region'] = df_regions
    adata_new.uns['atac-seq']['sample'] = df_samples
    adata_new.layers["z_score"] = df_zscores.values.T
    return adata_new


def write(adata,file_name=None,file_path=None,file_format='pkl'):
    """Write Anndate object to file
    
    Parameters
    ----------
    adata: AnnData
        Annotated data matrix. 
    file_name: `str`, optional (default: None)
        File name. If it's not specified, a file named 'stream_result' with the specified file format will be created 
        under the working directory
    file_path: `str`, optional (default: None)
        File path. If it's not specified, it's set to working directory
    file_format: `str`, optional (default: 'pkl')
        File format. By default it's compressed pickle file. Currently two file formats are supported:
        'pklz': compressed pickle file
        'pkl': pickle file
    """

    if(file_name is None):
        file_name = 'stream_result_atac.'+file_format
    if(file_path is None):
        file_path = adata.uns['workdir']
    if(file_format == 'pklz'):
        f = gzip.open(file_path+file_name, 'wb')
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()  
    elif(file_format == 'pkl'):
        f = open(file_path+file_name, 'wb')
        pickle.dump(adata, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.close()            
    else:
        print('file format ' + file_format + ' is not supported')
        return