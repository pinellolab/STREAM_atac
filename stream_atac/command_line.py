#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Authors: Huidong Chen
# Contact information: huidong.chen@mgh.harvard.edu

import warnings
warnings.filterwarnings('ignore')

__tool_name__='stream_atac'

import stream_atac
import argparse
import multiprocessing
import os

os.environ['KMP_DUPLICATE_LIB_OK']='True'


print('- STREAM single cell ATAC-seq preprocessing steps -',flush=True)
print('Version %s\n' % stream_atac.__version__,flush=True)


def main():
    parser = argparse.ArgumentParser(description='%s Parameters' % __tool_name__ ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--file_count", dest="file_count",default = None,required=True,
                        help="scATAC-seq counts file name in .tsv or .tsv.gz format", metavar="FILE")
    parser.add_argument("-r", "--file_region",dest="file_region", default=None,required=True,
                        help="scATAC-seq regions file name in .bed or .bed.gz format")
    parser.add_argument("-s","--file_sample",dest="file_sample", default=None,required=True,
                        help="scATAC-seq samples file name in .tsv or tsv.gz format")
    parser.add_argument("-g","--genome",dest="genome", default='hg19',
                        help="Reference genome. Choose from {{'mm9', 'mm10', 'hg38', 'hg19'}} ")    
    parser.add_argument("-f","--feature",dest="feature", default='kmer',
                        help="Features used to have the analysis. Choose from {{'kmer', 'motif'}}")    
    parser.add_argument("-k",dest="k",type=int,default=7,
                        help="k-mer length for scATAC-seq analysis")
    parser.add_argument("--resize_peak",dest="resize_peak",action="store_true",
                        help="Resize peaks when peaks have different widths.")
    parser.add_argument("--peak_width",dest = "peak_width",type=int, default=450,
                        help="Specify the width of peak when resizing them. Only valid when resize_peak is True.")
    parser.add_argument("--n_jobs",dest = "n_jobs",type=int, default=1,
                        help="The number of parallel jobs to run. (default,1)")
    parser.add_argument("--file_format",dest="file_format", default='tsv',
                        help=" File format of file_count. Currently supported file formats: 'tsv','txt','csv','mtx'.") 
    parser.add_argument("-o","--output_folder",dest="output_folder", default=None,
                        help="Output folder")

    args = parser.parse_args()

    file_count = args.file_count
    file_region = args.file_region
    file_sample = args.file_sample
    genome = args.genome
    feature = args.feature
    k = args.k
    resize_peak = args.resize_peak
    peak_width = args.peak_width
    n_jobs = args.n_jobs
    file_format = args.file_format
    output_folder = args.output_folder #work directory
    try:
        adata = stream_atac.preprocess_atac(file_count=file_count,file_region=file_region,file_sample=file_sample, genome = genome,
                                            feature=feature, k=k, resize_peak=resize_peak, peak_width=peak_width, file_format=file_format, n_jobs=n_jobs,workdir=output_folder)
    except:
        print("An exception occurred.")
    else:
        print('Finished computation.')

if __name__ == "__main__":
    main()