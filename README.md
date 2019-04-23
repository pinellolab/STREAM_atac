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
$ conda create -n myenv python=3.6 stream_atac
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
scATAC-seq counts file name in .tsv or .tsv.gz format  
-r, --file_region  
scATAC-seq regions file name in .bed or .bed.gz format  
-s, --file_sample  
scATAC-seq samples file name in .tsv or tsv.gz format 
-g, --genome
Reference genome. Choose from {{'mm9', 'mm10', 'hg38', 'hg19'}}
-f, --feature
Features used to have the analysis. Choose from {{'kmer', 'motif'}} 
-k
k-mer length for scATAC-seq analysis  
--ms, motif_species
Species of motifs in the JASPAR database. Choose from {{'Homo sapiens','Mus musculus'}}
--n_jobs  
The number of parallel jobs to run. (default, all the available cores)
-o, --output_folder  
Output folder (default: None)
```

Tutorial
--------

Example dataset can be found:

[Buenrostro_2018](https://www.dropbox.com/sh/zv6z7f3kzrafwmq/AACAlU8akbO_a-JOeJkiWT1za?dl=0)

Using *k-mers* to generate feature matrix:  

```sh
$ stream_atac -c count_file.tsv.gz -s sample_file.tsv.gz -r region_file.bed.gz -g hg19
```

Using *motifs* to generate feature matrix:  

```sh
$ stream_atac -c count_file.tsv.gz -s sample_file.tsv.gz -r region_file.bed.gz -g hg19 -f motif
```

**More downstream analyses with STREAM**:

* Example for scATAC-seq with k-mers: [STREAM_scATAC-seq.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/4.STREAM_scATAC-seq.ipynb?flush_cache=true)


