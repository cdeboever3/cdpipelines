cdpipelines
=========

This repository holds various bioinformatics pipelines.

## Dependencies

Many of the pipeline dependencies can be obtained using the `prepare` submodule
(see below). Additionally, a working Python environment is needed along with
some of the common scientific python packages. I recommend using [Anaconda
python](https://store.continuum.io/cshop/anaconda/) since it includes most of
the needed packages. If you are using Anaconda, I'd recommend making new
environments for different projects. Besides the default Anaconda packages, you
will need the following (available through `conda` or `pip`):

* HTSeq
* pandas
* pysam (this is available through conda but currently it's an old version so
  you have to get it using `pip`)
* PyVCF

### `rpy2`

Installing `rpy2` can be tricky. Different versions of `R` and `rpy2` don't work
well together, so it's recommended to make a local installation of `R` using the
`prepare` submodule and compile `rpy2` against this installation. You can
install `R` using `prepare.download_r` and install `rpy2` using
`prepare.download_install_rpy2`. `prepare.download_install_rpy2` will prompt you
to set your PATH, LDFLAGS, and LD_LIBRARY_PATH to correctly install `rpy2`.
After installing `rpy2`, you need to set your PATH and LD_LIBRARY_PATH using
these commands for every bash session where you want to use this `rpy2`. I'd
recommend putting the commands in a file that you source every time you load
the project's Anaconda environment.

## Submodules

### `general`

`general` contains methods used in multiple pipelines. Some pipelines use
similar but different versions of some methods, so the pipelines will have
their own versions of those methods. Sometimes it may make sense to add options
to a particular method that is used in multiple pipelines (where each pipeline
has slightly different versions) and add the method into `general`.

### `prepare`

The `prepare` module contains functions for downloading various software and
reference files needed for the different pipelines. 

### `rnaseq`

This pipeline currently starts from fastq files and has two steps. For detailed
information on each step, so the docstrings for each method. The first step is
`align_and_sort` which (optionally) removes duplicates, aligns the reads, and
makes coverage bigwig files for use with the UCSC genome browser or IGV. The
read alignments are output in both genomic and transcriptomic coordinates. The
second step is `get_counts` which counts reads overlapping genes for gene
differential expression and exonic bins for use with DEXSeq.
