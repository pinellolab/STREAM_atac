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
-k
k-mer length for scATAC-seq analysis  
--n_processes  
Specify the number of processes (default, all the available cores).
-o, --output_folder  
Output folder (default: None)
```

Example dataset can be found:

[Buenrostro_2018](https://www.dropbox.com/sh/zv6z7f3kzrafwmq/AACAlU8akbO_a-JOeJkiWT1za?dl=0)


```sh
$ stream_atac -c count_file.tsv.gz -s sample_file.tsv.gz -r region_file.bed.gz
```
