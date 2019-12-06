[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/stream_atac/README.html)

# STREAM_atac
 STREAM Preprocessing steps for single cell atac-seq data

Installation with Bioconda
--------------------------

1)	If Anaconda (or miniconda) is already installed with **Python 3**, skip to 2) otherwise please download and install Python3 Anaconda from here: https://www.anaconda.com/download/

2)	Open a terminal and add the Bioconda channel with the following commands:

```sh
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

3)	Create an environment named `myenv` , install **stream_atac**, and activate it with the following commands:

```sh
$ conda create -n myenv python stream_atac
$ conda activate myenv
```

Usage
-----

To run `stream_atac` at the command-line interface:

* start a terminal session;

* enter ```stream_atac --help [options]```

Users can specify the following options:
```
-c, --file_count  
scATAC-seq counts file name
-r, --file_region  
scATAC-seq regions file name in .bed or .bed.gz format  
-s, --file_sample  
scATAC-seq samples file name
-g, --genome
Reference genome. Choose from {{'mm9', 'mm10', 'hg38', 'hg19'}}
-f, --feature
Features used to have the analysis. Choose from {{'kmer', 'motif'}} 
-k
k-mer length for scATAC-seq analysis. (default, 7)
--resize_peak
Resize peaks when peaks have different widths.
--peak_width  
Specify the width of peak when resizing them. Only valid when resize_peak is True.
--n_jobs  
The number of parallel jobs to run. (default, 1)
--file_format   
File format of file_count. Currently supported file formats: 'tsv','txt','csv','mtx'.
-o, --output_folder  
Output folder (default: None)
```
File Format
-----------

To perform STREAM_atac preprocess, the main input includes three files: **count file**, **region file** and **sample file**. 

**count file**, .tsv or .tsv.gz format. A tab-delimited triplet file. It contains three columns. The first column specifies the rows indices (the regions) for non-zero entry. The second column specifies the columns indices (the sample) for non-zero entry. The last column contains the number of reads in a given region for a particular cell. No header is necessary:

|        |     |  |
|--------|-----|--|
| 3735   | 96  | 1|
| 432739 | 171 | 2|
| 133126 | 292 | 1|
| 219297 | 359 | 1|
| 284936 | 1222| 1|
| 442588 | 1580| 2|

**region file**, .bed or .bed.gz format. A tab-delimited .bed file with three columns. The first column specifies chromosome names. The second column specifies the start position of the region. The third column specifies the end position of the region. The order of regions should be consistent with the regions indices in the count file. No header is necessary:

|      |       |      |
|------|-------|------|
| chr1 | 10279 | 10779|
| chr1 | 13252 | 13752|
| chr1 | 16019 | 16519|
| chr1 | 29026 | 29526|
| chr1 | 96364 | 96864|

**sample file**, .tsv or .tsv.gz format. It has one column. Each row is a cell name.  The order of the cells should be consistent with the sample indices in count file. No header is necessary:

|                                    |
|------------------------------------|
| singles-BM0828-HSC-fresh-151027-1  | 
| singles-BM0828-HSC-fresh-151027-2  | 
| singles-BM0828-HSC-fresh-151027-3  |
| singles-BM0828-HSC-fresh-151027-4  |
| singles-BM0828-HSC-fresh-151027-5  |


Tutorial
--------

* #### Example dataset can be found: [Buenrostro_2018](https://www.dropbox.com/sh/zv6z7f3kzrafwmq/AACAlU8akbO_a-JOeJkiWT1za?dl=0)

> Using *k-mers* to generate zscore matrix:  

```sh
$ stream_atac -c count_file.tsv.gz -r region_file.bed.gz -s sample_file.tsv.gz -g hg19 -f kmer -k 7 --n_jobs 3 -o stream_output
```

> Using *motifs* to generate zscore matrix:  

```sh
$ stream_atac -c count_file.tsv.gz -r region_file.bed.gz -s sample_file.tsv.gz -g hg19 -f motif --n_jobs 3 -o stream_output
```

* #### For 10X CellRanger output:

> Using *k-mers* to generate zscore matrix:  
```sh
$ stream_atac -c ./filtered_peak_bc_matrix/matrix.mtx -r ./filtered_peak_bc_matrix/peaks.bed -s ./filtered_peak_bc_matrix/barcodes.tsv --file_format mtx -g hg19 -f kmer -k 7 --n_jobs 3 -o stream_output
```

> Using *motifs* to generate zscore matrix:  
```sh
$ stream_atac -c ./filtered_peak_bc_matrix/matrix.mtx -r ./filtered_peak_bc_matrix/peaks.bed -s ./filtered_peak_bc_matrix/barcodes.tsv --file_format mtx -g hg19 -f motif --n_jobs 3 -o stream_output
```

* #### Final Output

After running stream_atac, three files will be generated, including `zscores.tsv.gz`, `zscores_scaled.tsv.gz`, and `adata.h5ad`.

For the scaled z-score file, each row represents a k-mer DNA sequence/motif and each column represents one cell. Each entry is a scaled z-score of the accessibility of each k-mer/motif across cells. E.g. the scaled z-score using k-mers looks like:

|        | singles-BM0828-HSC-fresh-151027-1 | singles-BM0828-HSC-fresh-151027-2 | singles-BM0828-HSC-fresh-151027-3 |
|--------|-----------------------------------|-----------------------------------|-----------------------------------|
| AAAAAAA|-0.15973157637808505               | 0.18950966450007853               | 0.07713107176524692               | 
| AAAAAAG|-1.3630723054479532                | -0.04770034004421244              | 0.6387323857481045                |
| AAACACG|-0.2065161126378667                | -1.3375384076872765               | 0.2660278729402342                |
| AGCGTTA|-0.496859947462221                 | 0.7181918229050274                | 0.19603357892921522               |
| ATACTCA|-1.2127919166377426                | 0.7938414496478844                | -1.2665513250104594               |

**`zscores_scaled.tsv.gz` and `adata.h5ad` can be directly used for downstream analysis.**

```python
import stream as st
adata = st.read(file_name='./zscores_scaled.tsv.gz',experiment='atac-seq')
```

or

```python
import stream as st
adata = st.read(file_name='./adata.h5ad',experiment='atac-seq',file_format='h5ad')
```

* #### **stream_atac** package can be also directly imported within python:

```python
import stream_atac
adata = stream_atac.preprocess_atac(file_count='count_file.tsv.gz',
                        file_region='region_file.bed.gz',
                        file_sample='sample_file.tsv.gz',
                        genome='hg19',feature='kmer',k=7,n_jobs=10,
                        workdir='./stream_output')
```

```python
import stream_atac
adata = stream_atac.preprocess_atac(file_count='./filtered_peak_bc_matrix/matrix.mtx',
                                    file_region='./filtered_peak_bc_matrix/peaks.bed',
                                    file_sample='./filtered_peak_bc_matrix/barcodes.tsv',
                                    genome='hg19',feature='kmer',k=7,n_jobs=10,
                                    file_format='mtx',workdir='./stream_output')
```

#### **More downstream analyses with STREAM**:

* Example for scATAC-seq(using k-mers): [STREAM_scATAC-seq_k-mers.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/4.STREAM_scATAC-seq_k-mers.ipynb?flush_cache=true)

* Example for scATAC-seq(using motifs): [STREAM_scATAC-seq_motifs.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/5.STREAM_scATAC-seq_motifs.ipynb?flush_cache=true)

