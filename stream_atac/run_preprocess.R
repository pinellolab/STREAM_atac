suppressMessages(library(optparse,quietly = TRUE))

main <- function(){
  option_list = list(
    make_option(c("-c", "--count"), type="character", default=NULL, 
                help="scATAC-seq counts file name in .tsv or .tsv.gz format", metavar="character"),
    make_option(c("-r", "--region"), type="character", default=NULL, 
                help="scATAC-seq regions file name in .bed or .bed.gz format", metavar="character"),  
    make_option(c("-s", "--sample"), type="character", default=NULL, 
                help="scATAC-seq samples file name in .tsv or tsv.gz format", metavar="character"),  
    make_option(c("--file_format"), type="character", default='tsv', 
                help="File format of file_count. Currently supported file formats: 'tsv','txt','csv','mtx'", metavar="character"), 
    make_option(c("-g", "--genome"), type="character", default=NULL, 
                help="Reference genome. Choose from 'BSgenome.Hsapiens.UCSC.hg19','BSgenome.Hsapiens.UCSC.hg38','BSgenome.Mmusculus.UCSC.mm9','BSgenome.Mmusculus.UCSC.mm10' ", metavar="character"),  
    make_option(c("-f", "--feature"), type="character", default="kmer", 
                help="Features used to have the analysis [default = %default] . Choose from 'kmer','motif' ", metavar="character"),               
    make_option(c("-k","--k_kmer"), type="integer", default=7, 
                help="k-mer length for scATAC-seq analysis [default = %default]. Only valid when kmer is used", metavar="integer"),       
    make_option(c("--species"), type="character", default=NULL, 
                help="Species of motifs in the JASPAR database. Choose from 'Homo sapiens','Mus musculus'. Only valid when motif is used", metavar="character"),
    make_option(c("--n_jobs"), type="integer", default=1, 
                help="The number of parallel jobs to run. [default = %default]", metavar="integer"),
    make_option(c("--resize_peak"), action = "store_true",default=FALSE,
                help="resize regions(peaks) [default]"),
    make_option(c("--peak_width"), type="integer", default=450, 
                help="The peak width resized to. [default = %default]", metavar="integer"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output folder", metavar="character")    
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  if(any(is.null(opt$count),is.null(opt$region),is.null(opt$sample),is.null(opt$genome))){
    print_help(opt_parser)
    stop("counts file,regions file,samples file,reference genome must all be supplied", call.=FALSE)
  }
  # if(!(opt$genome %in% c('hg19','hg38','mm9','mm10'))){
  #   print_help(opt_parser)
  #   stop("Reference genome must be chosen from 'hg19','hg38','mm9','mm10' ", call.=FALSE)}
  # if(!(opt$feature %in% c('kmer','motif'))){
  #   print_help(opt_parser)
  #   stop("Features must be chosen from 'kmer','motif' ", call.=FALSE)}
  # if(!(opt$species %in% c('Homo sapiens','Mus musculus'))){
  #   print_help(opt_parser)
  #   stop("Species must be chosen from 'Homo sapiens','Mus musculus' ", call.=FALSE)}
  
  file.count = opt$count
  file.region = opt$region
  file.sample = opt$sample
  file.format = opt$file_format
  genome = opt$genome
  feature = opt$feature
  k = opt$k_kmer
  species = opt$species
  n_jobs = opt$n_jobs
  resize_peak = opt$resize_peak
  peak_width = opt$peak_width
  dir.output = opt$output
  
  suppressMessages(library(chromVAR))
  suppressMessages(library(Matrix))
  suppressMessages(library(SummarizedExperiment))
  suppressMessages(library(BiocParallel))
  suppressMessages(library(data.table))
  
  set.seed(2019)
  register(MulticoreParam(n_jobs))
  
  # Import Peaks
  print('Read in regions ...')
  peaks <- makeGRangesFromDataFrame(data.frame(fread(file.region,col.names=c('seqnames','start','end'))))
  if(resize_peak){
    print('resize regions/peaks ...')
    peaks <- resize(peaks, width = peak_width, fix = "center")
  }
  
  # Import Samples
  print('Read in samples ...')
  samples <- data.frame(fread(file.sample,header=FALSE))

  # Import counts
  print('Read in counts ...')
  if(file.format=='mtx'){
      counts = readMM(file.count)
      colnames(counts) <- samples[,1]
  }else{
      m <- data.matrix(data.frame(fread(cmd=paste0("zcat < ", file.count))))
      counts <- sparseMatrix(i = m[,1], j = m[,2], x=m[,3])
      colnames(counts) <- samples[,1]
  }
  
  # Make RangedSummarizedExperiment
  SE <- SummarizedExperiment(
    rowRanges = peaks,
    colData = samples,
    assays = list(counts = counts)
  )
  
  # Add GC bias
  print('Add GC bias ...')
  SE <- addGCBias(SE, genome = genome)
  
  print('Filter peaks ...')
  SE <- filterPeaks(SE, non_overlapping = TRUE)
  
  print('Get background peaks ...')
  bg <- getBackgroundPeaks(SE)
  
  if(feature == 'kmer'){
    # compute kmer deviations
    print('k-mer counting...')   
    kmer_ix <- matchKmers(k, SE, genome = genome)
    print('Computing k-mer deviations...') 
    dev <- computeDeviations(object = SE, annotations = kmer_ix,background_peaks = bg)
  }
  if(feature == 'motif'){
    # compute motif deviations
    suppressMessages(library('JASPAR2016')) 
    suppressMessages(library(motifmatchr))
    print('motif matching...')
    motifs = getJasparMotifs(species = species)    
    motif_ix <- matchMotifs(motifs, SE, genome = genome)
    print('Computing motif deviations...')
    dev <- computeDeviations(object = SE, annotations = motif_ix, background_peaks = bg)
  }
  
  devTable <- assays(dev)[["deviations"]]
  
  print('Saving zscores...')
  write.table(devTable, file = gzfile(file.path(dir.output,"zscores.tsv.gz")), sep = "\t", quote = FALSE, row.names = TRUE, col.names=NA)
}

main()