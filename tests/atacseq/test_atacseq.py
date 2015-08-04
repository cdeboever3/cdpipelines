import os

import pytest

import pipelines as ps

WEB_PATH_FILE = os.path.join(ps._root, 'tests', 'rnaseq', 'web_path_file.txt')

class TestAlign:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = WEB_PATH_FILE
        star_path = 'path/to/star'
        picard_path = 'path/to/picard'
        bedtools_path = 'path/to/bedtools'
        bedgraph_to_bigwig_path = 'path/to/bedgraph_to_bigwig'
        fastqc_path = 'path/to/fastqc'
        samtools_path = 'path/to/samtools'
        homer_path = 'path/to/homer'
        blacklist_bed = 'path/to/blacklist'
        conda_env = 'conda_env'
        environment = 'environment.sh'
        shell=False
        fn = ps.atacseq.align(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            star_path,
            picard_path,
            bedtools_path,
            bedgraph_to_bigwig_path,
            fastqc_path,
            samtools_path,
            blacklist_bed,
            conda_env,
            environment,
            shell=False
        )
        os.remove(fn)
        os.remove(tracklines_file)
        os.rmdir('s1_alignment')
    
    def test_run_multiple(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = ['r1_1.fastq.gz', 'r1_2.fastq.gz']
        r2_fastqs = ['r2_1.fastq.gz', 'r2_2.fastq.gz']
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = WEB_PATH_FILE
        star_path = 'path/to/star'
        picard_path = 'path/to/picard'
        bedtools_path = 'path/to/bedtools'
        bedgraph_to_bigwig_path = 'path/to/bedgraph_to_bigwig'
        fastqc_path = 'path/to/fastqc'
        samtools_path = 'path/to/samtools'
        homer_path = 'path/to/homer'
        blacklist_bed = 'path/to/blacklist'
        conda_env = 'conda_env'
        environment = 'environment.sh'
        shell=False
        fn = ps.atacseq.align(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            star_path,
            picard_path,
            bedtools_path,
            bedgraph_to_bigwig_path,
            fastqc_path,
            samtools_path,
            homer_path,
            blacklist_bed,
            conda_env,
            environment,
            shell=False
        )
        os.remove(fn)
        os.remove(tracklines_file)
        os.rmdir('s1_alignment')
    
    def test_run_trim(self):
        """Test to make sure the function at least runs with trimming"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = WEB_PATH_FILE
        star_path = 'path/to/star'
        picard_path = 'path/to/picard'
        bedtools_path = 'path/to/bedtools'
        bedgraph_to_bigwig_path = 'path/to/bedgraph_to_bigwig'
        fastqc_path = 'path/to/fastqc'
        samtools_path = 'path/to/samtools'
        homer_path = 'path/to/homer'
        blacklist_bed = 'path/to/blacklist'
        conda_env = 'conda_env'
        environment = 'environment.sh'
        shell = False
        trim = 25
        fn = ps.atacseq.align(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            star_path,
            picard_path,
            bedtools_path,
            bedgraph_to_bigwig_path,
            fastqc_path,
            samtools_path,
            homer_path,
            blacklist_bed,
            conda_env,
            environment,
            shell=shell,
            trim=trim
        )
        os.remove(fn)
        os.remove(tracklines_file)
        os.rmdir('s1_alignment')
    
    def test_run_multiple_trim(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = ['r1_1.fastq.gz', 'r1_2.fastq.gz']
        r2_fastqs = ['r2_1.fastq.gz', 'r2_2.fastq.gz']
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = WEB_PATH_FILE
        star_path = 'path/to/star'
        picard_path = 'path/to/picard'
        bedtools_path = 'path/to/bedtools'
        bedgraph_to_bigwig_path = 'path/to/bedgraph_to_bigwig'
        fastqc_path = 'path/to/fastqc'
        samtools_path = 'path/to/samtools'
        homer_path = 'path/to/homer'
        blacklist_bed = 'path/to/blacklist'
        conda_env = 'conda_env'
        environment = 'environment.sh'
        shell=False
        trim=-25
        fn = ps.atacseq.align(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            star_path,
            picard_path,
            bedtools_path,
            bedgraph_to_bigwig_path,
            fastqc_path,
            samtools_path,
            homer_path,
            blacklist_bed,
            conda_env,
            environment,
            shell=False,
            trim=trim
        )
        os.remove(fn)
        os.remove(tracklines_file)
        os.rmdir('s1_alignment')
