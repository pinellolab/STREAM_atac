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
import conda.cli

os.environ['KMP_DUPLICATE_LIB_OK']='True'

def preprocess_atac(file_count,file_region,file_sample, genome = 'hg19', motif_species = 'Homo sapiens',
                    feature = 'kmer',file_format='tsv',k=7,n_jobs = 1,file_path='',delimiter='\t',workdir=None,**kwargs):
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
    motif_species: `str`, optional (default: None)
        Species of motifs in the JASPAR database. Choose from {{'Homo sapiens','Mus musculus'}} 
    feature: `str`, optional (default: 'kmers')
        Features used to have the analysis. Choose from {{'kmer', 'motif'}} 
    file_format: `str`, optional (default: 'tsv')
        File format of file_count. Currently supported file formats: 'tsv','txt','csv','mtx'.
    k: `int`, optional (default: 7)
        k mer.  
    n_jobs: `int`, optional (default: 1)
        The number of parallel jobs to run
    delimiter: `str`, optional (default: '\t')
        Delimiter to use.
    file_path: `str`, optional (default: '')
        File path. By default it's empty
    delimiter: `str`, optional (default: '\t')
        Delimiter to use.
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
    
    if(genome not in ['mm9','mm10','hg38','hg19']):
        raise ValueError("Not supported reference genome: '%s'" % genome)   
    if(motif_species not in ['Homo sapiens','Mus musculus']):
        raise ValueError("Not supported species: '%s'" % motif_species)  
    if(feature not in ['kmer','motif']):
        raise ValueError("Not supported feature: '%s'" % feature)             
    dict_genome = {'mm9':'BSgenome.Mmusculus.UCSC.mm9',
                   'mm10':'BSgenome.Mmusculus.UCSC.mm10',
                   'hg38':'BSgenome.Hsapiens.UCSC.hg38',
                   'hg19':'BSgenome.Hsapiens.UCSC.hg19'}

    print('Importing packages...')
    rbase = importr('base')
    if(rbase.requireNamespace(dict_genome[genome], quietly = True)[0]):
        BSgenome = importr(dict_genome[genome])
    else:
        pkg = 'bioconductor-'+dict_genome[genome].lower()
        print("Installing pacakge '%s' ..." % pkg)
        conda.cli.main('conda', 'install',  '-y', pkg)
        BSgenome = importr(dict_genome[genome])
    if(feature == 'motif'):
        if(rbase.requireNamespace('motifmatchr', quietly = True)[0]):
            motifmatchr = importr('motifmatchr')
        else:
            pkg = 'bioconductor-motifmatchr'
            print("Installing pacakge '%s' ..." % pkg)
            conda.cli.main('conda', 'install',  '-y', pkg)
            motifmatchr = importr('motifmatchr')
        if(rbase.requireNamespace('JASPAR2016', quietly = True)[0]):
            JASPAR2016 = importr('JASPAR2016')
        else:
            pkg = 'bioconductor-jaspar2016'
            print("Installing pacakge '%s' ..." % pkg)
            conda.cli.main('conda', 'install',  '-y', pkg)
            JASPAR2016 = importr('JASPAR2016')            
    chromVAR = importr('chromVAR')
    GenomicRanges = importr('GenomicRanges')
    SummarizedExperiment = importr('SummarizedExperiment')
    r_Matrix = importr('Matrix')
    BiocParallel = importr('BiocParallel')
    print(str(n_jobs)+' cores are being used ...')
    BiocParallel.register(BiocParallel.MulticoreParam(n_jobs))
    pandas2ri.activate()
    rbase.set_seed(2019)
    
    if(workdir==None):
        workdir = os.path.join(os.getcwd(), 'stream_result')
        print("Using default working directory.")
    if(not os.path.exists(workdir)):
        os.makedirs(workdir)
    print('Saving results in: %s' % workdir) 
    
    _fp = lambda f:  os.path.join(file_path,f)
    print('Read in files ...')
    if(file_format in ['tsv','txt','csv']):
        df_counts = pd.read_csv(file_count,sep=delimiter,header=None,names=['i','j','x'],compression= 'gzip' if file_count.split('.')[-1]=='gz' else None)
    elif(file_format == 'mtx'):
        adata_atac = ad.read_mtx(_fp(file_count),**kwargs).T
        mat = adata_atac.X.tocoo()
        df_counts = pd.DataFrame(np.column_stack((mat.col+1,mat.row+1,mat.data)),columns=['i','j','x'],dtype=int)
    else:
        print('file format ' + file_format + ' is not supported')
        return        
    df_regions = pd.read_csv(file_region,sep=delimiter,header=None,compression= 'gzip' if file_region.split('.')[-1]=='gz' else None)
    df_regions = df_regions.iloc[:,:3]
    df_regions.columns = ['seqnames','start','end']
    df_samples = pd.read_csv(file_sample,sep=delimiter,header=None,names=['cell_id'],compression= 'gzip' if file_sample.split('.')[-1]=='gz' else None)                           

    
    r_regions_dataframe = pandas2ri.py2ri(df_regions)
    regions = GenomicRanges.makeGRangesFromDataFrame(r_regions_dataframe)
    counts = r_Matrix.sparseMatrix(i = df_counts['i'], j = df_counts['j'], x=df_counts['x'])
    samples = pandas2ri.py2ri(df_samples)

    SE = SummarizedExperiment.SummarizedExperiment(rowRanges = regions,colData = samples,assays = robjects.ListVector({'counts':counts}))
    SE = chromVAR.addGCBias(SE, genome = eval('BSgenome.'+dict_genome[genome].replace('.','_')))
    SE = chromVAR.filterPeaks(SE, non_overlapping = True)

    bg = chromVAR.getBackgroundPeaks(SE)
    if(feature == 'kmer'):
        # compute kmer deviations
        print('Computing k-mer deviations...')
        KmerMatch = chromVAR.matchKmers(k, SE, genome = eval('BSgenome.'+dict_genome[genome].replace('.','_')))
        dev = chromVAR.computeDeviations(SE, annotations = KmerMatch,background_peaks = bg)            
    if(feature == 'motif'):
        # compute motif deviations
        print('Computing motif deviations...')
        motifs = chromVAR.getJasparMotifs(species = motif_species)
        motifMatch = motifmatchr.matchMotifs(motifs, SE, genome = eval('BSgenome.'+dict_genome[genome].replace('.','_')))
        dev = chromVAR.computeDeviations(SE, annotations = motifMatch,background_peaks = bg)
    devTable = SummarizedExperiment.assays(dev)
    cn = pandas2ri.ri2py((dev.do_slot('colData')).do_slot('listData').rx2('cell_id'))
    rn = pandas2ri.ri2py(dev.do_slot('NAMES'))
    scores = pandas2ri.ri2py(devTable.do_slot('listData').rx2('deviations'))    

    df_zscores = pd.DataFrame(scores,index=rn,columns=cn)
    df_zscores_scaled = preprocessing.scale(df_zscores,axis=1)
    df_zscores_scaled = pd.DataFrame(df_zscores_scaled,index=df_zscores.index,columns=df_zscores.columns)
    df_zscores_scaled.to_csv(os.path.join(workdir,'zscore.tsv.gz'),sep = '\t',compression='gzip')
    adata = ad.AnnData(X=df_zscores_scaled.values.T, obs={'obs_names':df_zscores_scaled.columns},var={'var_names':df_zscores_scaled.index})
    adata.raw = adata
    adata.uns['workdir'] = workdir
    adata.uns['experiment'] = 'atac-seq'
    adata.layers["z_score"] = df_zscores.values.T
    adata.uns['atac-seq'] = dict()
    adata.uns['atac-seq']['count'] = df_counts
    adata.uns['atac-seq']['region'] = df_regions
    adata.uns['atac-seq']['sample'] = df_samples
    return adata