Frazer Lab Pipelines
=========

This repository holds various computational pipelines.

## Dependencies

Many of the pipeline dependencies can be obtained using the `prepare` submodule
(see below). Additionally, a working Python environment is needed along with
some of the common scientific python packages. I recommend using [Anaconda
python](https://store.continuum.io/cshop/anaconda/) since it includes most of
the needed packages. If you are using Anaconda, I'd recommend making a new
environment. Besides the default Anaconda packages, you will need the following
(available through `pip`):

* HTSeq
* pysam (you can this through conda but currently it's an old version)

## Prepare

The prepare module contains functions for downloading various software and
reference files needed for the different pipelines. Some dependencies 

## RNA-seq

This pipeline currently starts from fastq files and has two steps. For detailed
information on each step, so the docstrings for each method. The first step is
`align_and_sort` which (optionally) removes duplicates, aligns the reads, and
makes coverage bigwig files for use with the UCSC genome browser or IGV. The
read alignments are output in both genomic and transcriptomic coordinates. The
second step is `get_counts` which counts reads overlapping genes for gene
differential expression and exonic bins for use with DEXSeq.
