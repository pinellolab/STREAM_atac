#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    parser.add_argument("-k",dest="k",type=int,default=7,
                        help="k-mer length for scATAC-seq analysis")
    parser.add_argument("--n_processes",dest = "n_processes",type=int, default=multiprocessing.cpu_count(),
                        help="The number of processes. (default, all the available cores)")
    parser.add_argument("-o","--output_folder",dest="output_folder", default=None,
                        help="Output folder")

    args = parser.parse_args()

    file_count = args.file_count
    file_region = args.file_region
    file_sample = args.file_sample
    k = args.k
    n_processes = args.n_processes
    output_folder = args.output_folder #work directory
    try:
        adata = stream_atac.preprocess_atac(file_count=file_count,file_region=file_region,file_sample=file_sample,k=k,n_jobs=n_processes,workdir=output_folder)
    except:
        print("An exception occurred.")
    else:
        print('Finished computation.')

if __name__ == "__main__":
    main()