import os

import pytest

import pipelines as ps

class TestAlignAndSort:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        star_path = 'path/to/star'
        picard_path = 'path/to/picard'
        bedtools_path = 'path/to/bedtools'
        bedgraph_to_bigwig_path = 'path/to/bedgraph_to_bigwig'
        fastqc_path = 'path/to/fastqc'
        samtools_path = 'path/to/samtools'
        macs2_path = 'path/to/macs2'
        shell=False
        fn = ps.atacseq.align_and_call_peaks(
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
            macs2_path,
            shell=False
        )
        os.remove(fn)
        os.remove(tracklines_file)
    
    def test_run_multiple(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = ['r1_1.fastq.gz', 'r1_2.fastq.gz']
        r2_fastqs = ['r2_1.fastq.gz', 'r2_2.fastq.gz']
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        star_path = 'path/to/star'
        picard_path = 'path/to/picard'
        bedtools_path = 'path/to/bedtools'
        bedgraph_to_bigwig_path = 'path/to/bedgraph_to_bigwig'
        fastqc_path = 'path/to/fastqc'
        samtools_path = 'path/to/samtools'
        macs2_path = 'path/to/macs2'
        shell=False
        fn = ps.atacseq.align_and_call_peaks(
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
            macs2_path,
            shell=False
        )
        os.remove(fn)
        # os.remove(tracklines_file)
