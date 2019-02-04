def preprocess_atac(file_count=None,file_region=None,file_sample=None,k=7,n_jobs = multiprocessing.cpu_count(),file_path='',delimiter='\t',workdir=None):
    """Preprocess single cell atac-seq data and genearate a scaled z-score matrix.
    
    Parameters
    ----------
    file_count: `str`
        Count file name. A compressed matrix in sparse format (column-oriented).
    file_region: `str`
        Region file name. It has one column. Each row is a cell name (No header should be included). 
    file_sample: `str`
        Sample file name. It has three columns, i.e.,chromosome names, the start position of the region, the end position of the region (No header should be included).
    k: `int`, optional (default: 7)
        k mer.  
    n_jobs: `int`, optional (default: all available cpus)
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
    
    if(sum(list(map(lambda x: x is not None,[file_count,file_region,file_sample])))<3):
        print('count file, region file, and sample file should all be provided.')
        return
    else:
        print('Importing packages...')
        chromVAR = importr('chromVAR')
        GenomicRanges = importr('GenomicRanges')
        SummarizedExperiment = importr('SummarizedExperiment')
        BSgenome_Hsapiens_UCSC_hg19 = importr('BSgenome.Hsapiens.UCSC.hg19')
        r_Matrix = importr('Matrix')
        BiocParallel = importr('BiocParallel')
        BiocParallel.register(BiocParallel.MulticoreParam(n_jobs))
        pandas2ri.activate()
        
        if(workdir==None):
            workdir = os.path.join(os.getcwd(), 'stream_result')
            print("Using default working directory.")
        if(not os.path.exists(workdir)):
            os.makedirs(workdir)
        print('Saving results in: %s' % workdir) 
        
        _fp = lambda f:  os.path.join(file_path,f)
        print('Read in files ...')
        df_counts = pd.read_csv(file_count,sep=delimiter,header=None,names=['i','j','x'],compression= 'gzip' if file_count.split('.')[-1]=='gz' else None)
        df_regions = pd.read_csv(file_region,sep=delimiter,header=None,compression= 'gzip' if file_region.split('.')[-1]=='gz' else None)
        df_regions = df_regions.iloc[:,:3]
        df_regions.columns = ['seqnames','start','end']
        df_samples = pd.read_csv(file_sample,sep=delimiter,header=None,names=['cell_id'],compression= 'gzip' if file_sample.split('.')[-1]=='gz' else None)                           

        r_regions_dataframe = pandas2ri.py2ri(df_regions)
        regions = GenomicRanges.makeGRangesFromDataFrame(r_regions_dataframe)
        counts = r_Matrix.sparseMatrix(i = df_counts['i'], j = df_counts['j'], x=df_counts['x'])
        samples = pandas2ri.py2ri(df_samples)
        samples.rownames = df_samples['cell_id']

        SE = SummarizedExperiment.SummarizedExperiment(rowRanges = regions,colData = samples,assays = robjects.ListVector({'counts':counts}))
        SE = chromVAR.addGCBias(SE, genome = BSgenome_Hsapiens_UCSC_hg19.BSgenome_Hsapiens_UCSC_hg19)

        # compute kmer deviations
        print('Computing k-mer deviations...')
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