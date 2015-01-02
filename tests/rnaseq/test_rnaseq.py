import os

import pytest

import pipelines as ps

class TestPbsHeader:
    def test_run(self):
        """Test to make sure the function at least runs"""
        out = 'run.out'
        err = 'run.err'
        name = 'run'
        threads = 30
        lines = ps.rnaseq._pbs_header(out, err, name, threads)

class TestCbarrettPairedDupRemoval:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        r1_nodup = 'r1_nodup.fastq.gz'
        r2_nodup = 'r2_nodup.fastq.gz'
        temp_dir = 'temp_dir'
        lines = ps.rnaseq._cbarrett_paired_dup_removal(r1_fastqs, r2_fastqs,
                                                        r1_nodup, r2_nodup,
                                                        temp_dir)
    
    def test_list_run(self):
        """Test to make sure the function at least runs with lists of files as
        input"""
        r1_fastqs = ['r1_1.fastq.gz', 'r1_2.fastq.gz']
        r2_fastqs = ['r2_1.fastq.gz', 'r2_2.fastq.gz']
        r1_nodup = 'r1_nodup.fastq.gz'
        r2_nodup = 'r2_nodup.fastq.gz'
        temp_dir = 'temp_dir'
        lines = ps.rnaseq._cbarrett_paired_dup_removal(r1_fastqs, r2_fastqs,
                                                        r1_nodup, r2_nodup,
                                                        temp_dir)

class TestStarAlign:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        sample = 's1'
        rgpl = 'ILLUMINA'
        rgpu = 'flowcell_barcode'
        star_index = 'path/to/star/index'
        star_path = 'path/to/star/executable'
        threads = 30
        lines = ps.rnaseq._star_align(r1_fastqs, r2_fastqs, sample, rgpl, rgpu,
                                       star_index, star_path, threads)

class TestPicardCoordSort:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_bam = 'test.sorted.bam'
        bam_index = 'test.sorted.bam.bai'
        picard_path = 'picard'
        picard_memory = '58G'
        temp_dir = 'temp_dir'
        lines = ps.rnaseq._picard_coord_sort(in_bam, out_bam, bam_index,
                                             picard_path, picard_memory,
                                             temp_dir)

class TestPicardIndex:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        index = 'test.bam.bai'
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        temp_dir = 'path/to/temp/dir'
        lines = ps.rnaseq._picard_index(in_bam, index, picard_memory,
                                         picard_path, temp_dir)

class TestBedgraphToBigwig:
    def test_run(self):
        """Test to make sure the function at least runs"""
        bedgraph = 'test.bg'
        bigwig = 'out.bw'
        bedgraph_to_bigwig_path = 'path/to/BedGraphToBigWig'
        bedtools_path = 'path/to/bedtools'
        lines = ps.rnaseq._bedgraph_to_bigwig(bedgraph, bigwig,
                                               bedgraph_to_bigwig_path,
                                               bedtools_path)

class TestCoverageBedgraph:
    def test_run(self):
        """Test to make sure the function at least runs"""
        bam = 'test.bam'
        bedgraph = 'test.bg'
        bedtools_path = 'path/to/bedtools'
        sample_name = 's1'
        lines = ps.rnaseq._coverage_bedgraph(bam, bedgraph, bedtools_path,
                                              sample_name)

    def test_run_plus(self):
        """Test to make sure the function at least runs for plus strand"""
        bam = 'test.bam'
        bedgraph = 'test.bg'
        bedtools_path = 'path/to/bedtools'
        sample_name = 's1'
        strand = '+'
        lines = ps.rnaseq._coverage_bedgraph(bam, bedgraph, bedtools_path,
                                              sample_name, strand)

    def test_run_minus(self):
        """Test to make sure the function at least runs for minus strand"""
        bam = 'test.bam'
        bedgraph = 'test.bg'
        bedtools_path = 'path/to/bedtools'
        sample_name = 's1'
        strand = '-'
        lines = ps.rnaseq._coverage_bedgraph(bam, bedgraph, bedtools_path,
                                              sample_name, strand)

class TestBigwigFiles:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_bigwig = 'test.bw'
        sample_name = 's1'
        bedgraph_to_bigwig_path = '/path/to/BedGraphToBigWig'
        bedtools_path = 'path/to/bedtools'
        lines = ps.rnaseq._bigwig_files(in_bam, out_bigwig, sample_name,
                                         bedgraph_to_bigwig_path, bedtools_path)

    def test_run_stranded(self):
        """Test to make sure the function at least runs for stranded output"""
        in_bam = 'test.bam'
        out_bigwig = 'plus.bw'
        sample_name = 's1'
        bedgraph_to_bigwig_path = '/path/to/BedGraphToBigWig'
        bedtools_path = 'path/to/bedtools'
        out_bigwig_minus = 'minus.bw'
        lines = ps.rnaseq._bigwig_files(in_bam, out_bigwig, sample_name,
                                         bedgraph_to_bigwig_path, bedtools_path,
                                         out_bigwig_minus=out_bigwig_minus)

class TestGenomeBrowserFiles:
    def test_run(self):
        """Test to make sure the function at least runs"""
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        coord_sorted_bam = 'test.bam'
        bam_index = 'test.bam.bai'
        bigwig = 'test.bw'
        sample_name = 's1'
        out_dir = '.'
        lines = ps.rnaseq._genome_browser_files(tracklines_file, link_dir,
                                                web_path_file, coord_sorted_bam,
                                                bam_index, bigwig, sample_name,
                                                out_dir, bigwig_minus='')
        os.remove(tracklines_file)

    def test_run(self):
        """Test to make sure the function at least runs for standed input"""
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        coord_sorted_bam = 'test.bam'
        bam_index = 'test.bam.bai'
        bigwig = 'plus.bw'
        sample_name = 's1'
        bigwig_minus = 'minus.bw'
        out_dir = '.'
        lines = ps.rnaseq._genome_browser_files(tracklines_file, link_dir,
                                                web_path_file, coord_sorted_bam,
                                                bam_index, bigwig, sample_name,
                                                out_dir,
                                                bigwig_minus=bigwig_minus)
        out_dir = '.'
        os.remove(tracklines_file)

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
        remove_dup=True, 
        strand_specific=False, 
        shell=False
        fn = ps.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            remove_dup=True, 
            strand_specific=False, 
            shell=False
        )
        os.remove(fn)
        os.remove(tracklines_file)
    
    def test_run_no_remove_dup(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        remove_dup = False,
        fn = ps.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            remove_dup=remove_dup
        )
        os.remove(fn)
        os.remove(tracklines_file)
    
    def test_run_no_strand_specific(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        strand_specific = True,
        fn = ps.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            strand_specific=strand_specific
        )
        os.remove(fn)
        os.remove(tracklines_file)
    
    def test_run_shell(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = '.'
        sample_name = 's1'
        star_index = 'path/to/index'
        tracklines_file = 'tracklines.txt'
        link_dir = '.'
        web_path_file = 'web_path_file.txt'
        shell = True 
        fn = ps.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            tracklines_file,
            link_dir,
            web_path_file,
            shell=True
        )
        os.remove(fn)
        os.remove(tracklines_file)

# class TestDexseqCount:
#     def test_run(self):
#         """Test to make sure the function at least runs"""
#         bam = 'test.bam'
#         counts_file = 'counts.tsv'
#         dexseq_annotation = 'annot.gtf'
#         samtools_path = 'path/to/samtools'
#         lines = ps.rnaseq._dexseq_count(bam, counts_file, dexseq_annotation,
#                                         samtools_path=samtools_path)
# 
#     def test_run_stranded(self):
#         """Test to make sure the function at least runs"""
#         bam = 'test.bam'
#         counts_file = 'counts.tsv'
#         dexseq_annotation = 'annot.gtf'
#         stranded = True
#         samtools_path = 'path/to/samtools'
#         lines = ps.rnaseq._dexseq_count(bam, counts_file, dexseq_annotation,
#                                         stranded=stranded,
#                                         samtools_path=samtools_path)

class TestHtseqCount:
    def test_run(self):
        """Test to make sure the function at least runs"""
        bam = 'test.bam'
        counts_file = 'counts.tsv'
        stats_file = 'stats.tsv'
        gtf = 'annot.gtf'
        samtools_path = 'path/to/samtools'
        lines = ps.rnaseq._htseq_count(bam, counts_file, stats_file, gtf,
                                       samtools_path=samtools_path)

    def test_run_stranded(self):
        """Test to make sure the function at least runs"""
        bam = 'test.bam'
        counts_file = 'counts.tsv'
        stats_file = 'stats.tsv'
        gtf = 'annot.gtf'
        stranded = True
        samtools_path = 'path/to/samtools'
        lines = ps.rnaseq._htseq_count(bam, counts_file, stats_file, gtf,
                                       stranded=stranded,
                                       samtools_path=samtools_path)
