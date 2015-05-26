import os

from general import _bedgraph_to_bigwig
from general import _bigwig_files
from general import _coverage_bedgraph
from general import _fastqc
from general import JobScript
from general import _make_softlink
from general import _pbs_header
from general import _picard_bam_index_stats
from general import _picard_coord_sort
from general import _picard_collect_multiple_metrics
from general import _picard_collect_rna_seq_metrics
from general import _picard_gc_bias_metrics
from general import _picard_index
from general import _picard_mark_duplicates
from general import _process_fastqs

def _star_align(r1_fastq, r2_fastq, sample, rgpl, rgpu, star_index, star_path,
                threads):
    """
    Align paired fastq files with STAR.

    Parameters
    ----------
    r1_fastqs : list 
        List of paths to gzipped R1 fastq file(s). 

    r2_fastqs : list 
        List of paths to gzipped R2 fastq file(s). 

    sample : str
        Sample name.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    """
    # r1_fastqs.sort()
    # r2_fastqs.sort()
    # r1_fastqs = ','.join(r1_fastqs)
    # r2_fastqs = ','.join(r2_fastqs)
    # I use threads - 2 for STAR so there are open processors for reading and
    # writing.
    line = (' \\\n'.join([star_path, 
                          '\t--runThreadN {}'.format(threads - 2),
                          '\t--genomeDir {}'.format(star_index), 
                          '\t--genomeLoad LoadAndRemove', 
                          '\t--readFilesCommand zcat',
                          '\t--readFilesIn {} {}'.format(r1_fastq, r2_fastq),
                          '\t--outSAMtype BAM Unsorted', 
                          '\t--outSAMattributes All', 
                          '\t--outSAMunmapped Within',
                          ('\t--outSAMattrRGline ID:1 ' + 
                           'PL:{} '.format(rgpl) + 
                           'PU:{} '.format(rgpu) + 
                           'LB:{0} SM:{0}'.format(sample)), 
                          '\t--outFilterMultimapNmax 20', 
                          '\t--outFilterMismatchNmax 999',
                          '\t--outFilterMismatchNoverLmax 0.04',
                          ('\t--outFilterIntronMotifs '
                           'RemoveNoncanonicalUnannotated'),
                           '\t--outSJfilterOverhangMin 6 6 6 6',
                           '\t--seedSearchStartLmax 20',
                           '\t--alignSJDBoverhangMin 1', 
                          '\t--quantMode TranscriptomeSAM']) + '\n\n') 
    return line

# deprecated
# def _genome_browser_files(tracklines_file, link_dir, web_path_file,
#                           coord_sorted_bam, bam_index, bigwig, sample_name,
#                           outdir, bigwig_minus=''):
#     """
#     Make files and softlinks for displaying results on UCSC genome browser.
# 
#     Parameters
#     ----------
#     tracklines_file : str
#         Path to file for writing tracklines. The tracklines will be added to the
#         file; the contents of the file will not be overwritten. These tracklines
#         can be pasted into the genome browser upload for custom data.
# 
#     link_dir : str
#         Path to directory where softlink should be made.
# 
#     web_path_file : str
#         File whose first line is the URL that points to link_dir. For example,
#         if we make a link to the file s1_coord_sorted.bam in link_dir and
#         web_path_file has http://site.com/files on its first line, then
#         http://site.com/files/s1_coord_sorted.bam should be available on the
#         web. If the web directory is password protected (it probably should be),
#         then the URL should look like http://username:password@site.com/files.
#         This is a file so you don't have to make the username/password combo
#         public (although I'd recommend not using a sensitive password). You can
#         just put the web_path_file in a directory that isn't tracked by git, 
#         figshare, etc.
# 
#     coord_sorted_bam : str
#         Path to coordinate sorted bam file.
# 
#     bam_index : str
#         Path to index file for coordinate sorted bam file.
# 
#     bigwig : str
#         Path to bigwig file. If bigwig_minus is provided, bigwig has the plus
#         strand coverage.
# 
#     sample_name : str
#         Sample name used for naming files.
# 
#     bigwig_minus : str
#         Path to bigwig file for minus strand. If bigwig_minus is not provided,
#         bigwig is assumed to have coverage for both plus and minus stand reads.
# 
#     Returns
#     -------
#     lines : str
#         Lines to be printed to shell/PBS script.
# 
#     """
#     lines = ''
#     link_dir = os.path.join(link_dir, 'rna')
# 
#     with open(web_path_file) as wpf:
#         web_path = wpf.readline().strip()
#     web_path = web_path + '/rna'
# 
#     # File with UCSC tracklines.
#     if os.path.exists(tracklines_file):
#         with open(tracklines_file) as f:
#             tf_lines = f.read()
#     else:
#         tf_lines = ''
#     
#     # Bam file and index.
#     temp_link_dir = os.path.join(link_dir, 'bam')
#     temp_web_path = web_path + '/bam'
#     try:
#         os.makedirs(temp_link_dir)
#     except OSError:
#         pass
#     fn = os.path.join(outdir, os.path.split(coord_sorted_bam)[1])
#     new_lines, bam_name = _make_softlink(fn, sample_name, temp_link_dir)
#     lines += new_lines
# 
#     fn = os.path.join(outdir, os.path.split(bam_index)[1])
#     new_lines, index_name = _make_softlink(fn, sample_name, temp_link_dir)
#     lines += new_lines
# 
#     tf_lines += ' '.join(['track', 'type=bam',
#                           'name="{}_bam"'.format(sample_name),
#                           'description="RNAseq for {}"'.format(sample_name),
#                           'bigDataUrl={}/{}\n'.format(temp_web_path, bam_name)])
#     
#     # Bigwig file(s).
#     temp_link_dir = os.path.join(link_dir, 'bw')
#     temp_web_path = web_path + '/bw'
#     try:
#         os.makedirs(temp_link_dir)
#     except OSError:
#         pass
#     if bigwig_minus != '':
#         fn = os.path.join(outdir, os.path.split(bigwig)[1])
#         new_lines, plus_name = _make_softlink(fn, sample_name, temp_link_dir)
#         lines += new_lines
#         
#         fn = os.path.join(outdir, os.path.split(bigwig_minus)[1])
#         new_lines, minus_name = _make_softlink(fn, sample_name, temp_link_dir)
#         lines += new_lines
#  
#         tf_lines += ' '.join(['track', 'type=bigWig',
#                               'name="{}_plus_cov"'.format(sample_name),
#                               ('description="RNAseq plus strand coverage for '
#                                '{}"'.format(sample_name)),
#                               'visibility=0',
#                               'db=hg19',
#                               'bigDataUrl={}/{}\n'.format(temp_web_path,
#                                                           plus_name)])
#         tf_lines += ' '.join(['track', 'type=bigWig',
#                               'name="{}_minus_cov"'.format(sample_name),
#                               ('description="RNAseq minus strand coverage for '
#                                '{}"'.format(sample_name)),
#                               'visibility=0',
#                               'db=hg19',
#                               'bigDataUrl={}/{}\n'.format(temp_web_path,
#                                                           minus_name)])
#     else:
#         fn = os.path.join(outdir, os.path.split(bigwig)[1])
#         new_lines, bigwig_name = _make_softlink(fn, sample_name, temp_link_dir)
#         lines += new_lines
# 
#         tf_lines += ' '.join(['track', 'type=bigWig',
#                               'name="{}_cov"'.format(sample_name),
#                               ('description="RNAseq coverage for '
#                                '{}"'.format(sample_name)),
#                               'visibility=0',
#                               'db=hg19',
#                               'bigDataUrl={}/{}\n'.format(temp_web_path,
#                                                           bigwig_name)])
#     
#     with open(tracklines_file, 'w') as tf:
#         tf.write(tf_lines)
#     
#     lines += '\n'
#     return lines

def align_and_sort(
    r1_fastqs, 
    r2_fastqs, 
    outdir, 
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
    ref_flat, 
    rrna_intervals,
    rgpl='ILLUMINA',
    rgpu='',
    tempdir='/scratch', 
    threads=10, 
    picard_memory=15,
    strand_specific=True, 
    shell=False
):
    """
    Make a PBS or shell script for aligning RNA-seq reads with STAR. The
    defaults are set for use on the Frazer lab's PBS scheduler on FLC.

    Parameters
    ----------
    r1_fastqs : list or str
        Either a list of paths to gzipped fastq files with R1 reads or path to a
        single gzipped fastq file with R1 reads.

    r2_fastqs : list or str
        Either a list of paths to gzipped fastq files with R2 reads or path to a
        single gzipped fastq file with R2 reads.

    outdir : str
        Directory to store PBS/shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

    tracklines_file : str
        Path to file for writing tracklines. The tracklines will be added to the
        file; the contents of the file will not be overwritten. These tracklines
        can be pasted into the genome browser upload for custom data.

    link_dir : str
        Path to directory where softlinks for genome browser should be made.

    web_path_file : str
        File whose first line is the URL that points to link_dir. For example,
        if we make a link to the file s1_coord_sorted.bam in link_dir and
        web_path_file has http://site.com/files on its first line, then
        http://site.com/files/s1_coord_sorted.bam should be available on the
        web. If the web directory is password protected (it probably should be),
        then the URL should look like http://username:password@site.com/files.
        This is a file so you don't have to make the username/password combo
        public (although I'd recommend not using a sensitive password). You can
        just put the web_path_file in a directory that isn't tracked by git, 
        figshare, etc.

    star_path : str
        Path to STAR aligner.

    picard_path : str
        Path to Picard tools.

    bedtools_path : str
        Path to bedtools.

    bedgraph_to_bigwig_path : str
        Path bedGraphToBigWig executable.

    ref_flat : str
        Path to refFlat file with non-rRNA genes. Can ge gzipped.

    rrna_intervals : str
        Pato to interval list file with rRNA intervals.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    tempdir : str
        Directory to store files as STAR runs.

    threads : int
        Number of threads to reserve using PBS scheduler. This number of threads
        minus 2 will be used by STAR, so this must be at least 3.

    picard_memory : int
        Amount of memory (in gb) to give Picard Tools.

    strand_specific : boolean
        If false, data is not strand specific.

    shell : boolean
        If true, make a shell script rather than a PBS script.

    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    assert threads >= 3
    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()
    web_path = web_path + '/rna'
    link_dir = os.path.join(link_dir, 'rna')
    job_suffix = 'alignment'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    shell=shell, queue='high', copy_input=True)

    # I'm going to handle the copying and deleting of the fastqs myself rather
    # than have the JobScript do it because I don't want to the fastqs to sit
    # around on the disk the whole time after I'm done with them.
    temp_r1_fastqs = _process_fastqs(r1_fastqs, job.tempdir)
    temp_r2_fastqs = _process_fastqs(r2_fastqs, job.tempdir)
    if type(r1_fastqs) != list:
        r1_fastqs = [r1_fastqs]
    if type(r2_fastqs) != list:
        r2_fastqs = [r2_fastqs]
    combined_r1 = os.path.join(
        job.tempdir, '{}_combined_R1.fastq.gz'.format(sample_name))
    combined_r2 = os.path.join(
        job.tempdir, '{}_combined_R2.fastq.gz'.format(sample_name))

    # Files that will be created.
    aligned_bam = os.path.join(job.tempdir, 'Aligned.out.bam')
    job.temp_files_to_delete.append(aligned_bam)
    coord_sorted_bam = os.path.join(
        job.tempdir, '{}_sorted.bam'.format(sample_name))
    job.temp_files_to_delete.append(coord_sorted_bam)
    job.temp_files_to_delete.append('_STARtmp')
    out_bam = os.path.join(
        job.tempdir, '{}_sorted_mdup.bam'.format(sample_name))
    bam_index = '{}.bai'.format(out_bam)
    job.output_files_to_copy += [out_bam, bam_index]

    # Other STAR Files to copy to output directory.
    # TODO: Eventually, I'd like all filenames to include the sample name.
    job.output_files_to_copy += [
        'Log.out', 
        'Log.final.out', 
        'Log.progress.out', 
        'SJ.out.tab', 
        'Aligned.toTranscriptome.out.bam',
    ]

    # if strand_specific:
        # TODO: I need to think about how to calculate coverage for strand
        # specific data to make sure I'm doing it correctly and what the best
        # way is to display the data. Maybe I could make a hub and color the two
        # strand differently?
        # out_bigwig_plus = os.path.join(job.tempdir,
        #                                '{}_plus_rna.bw'.format(sample_name))
        # out_bigwig_minus = os.path.join(job.tempdir,
        #                                 '{}_minus_rna.bw'.format(sample_name))
        # files_to_copy.append(out_bigwig_plus)
        # files_to_copy.append(out_bigwig_minus)
    out_bigwig = os.path.join(job.tempdir, '{}_rna.bw'.format(sample_name))
    job.output_files_to_copy.append(out_bigwig)
    
    with open(job.filename, "a") as f:
        # I'm going to copy the fastq files here rather than have the JobScript
        # class do it because I want to delete them as soon as I've combined
        # them into a single file.
        f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
            '\\\n\t'.join(r1_fastqs + r2_fastqs), job.tempdir))
        # f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
        #     '\\\n\t'.join(r1_fastqs + r2_fastqs),
        #     job.tempdir))

        # Combine fastq files.
        if len(temp_r1_fastqs) > 1:
            f.write('cat \\\n{} \\\n\t> {} &\n'.format(
                ' \\\n'.join(['\t{}'.format(x) for x in temp_r1_fastqs]),
                combined_r1))
        else:
            f.write('mv {} {}\n\n'.format(temp_r1_fastqs[0], combined_r1))

        if len(temp_r2_fastqs) > 1:
            f.write('cat \\\n{} \\\n\t> {}\n\n'.format(
                ' \\\n'.join(['\t{}'.format(x) for x in temp_r2_fastqs]),
                combined_r2))
        else:
            f.write('mv {} {}\n\n'.format(temp_r2_fastqs[0], combined_r2))
        f.write('wait\n\n')

        # Remove temp fastqs.
        if len(temp_r1_fastqs) > 1:
            f.write('rm \\\n{} &\n\n'.format(
                ' \\\n'.join(['\t{}'.format(x) for x in temp_r1_fastqs])))
        if len(temp_r2_fastqs) > 1:
            f.write('rm \\\n{}\n\n'.format(
                ' \\\n'.join(['\t{}'.format(x) for x in temp_r2_fastqs])))
            f.write('wait\n\n')

        # Run FASTQC.
        lines = _fastqc([combined_r1, combined_r2], threads, job.outdir,
                        fastqc_path)
        f.write(lines)
        f.write('wait\n\n')
        r1 = '.'.join(os.path.split(combined_r1)[1].split('.')[0:-2])
        job.add_softlink(os.path.join(job.outdir, r1), 
                         os.path.join(link_dir, 'fastqc', r1))
        r2 = '.'.join(os.path.split(combined_r2)[1].split('.')[0:-2])
        job.add_softlink(os.path.join(job.outdir, r2), 
                         os.path.join(link_dir, 'fastqc', r2))
        with open(tracklines_file, "a") as tf:
            tf_lines = ('{}/fastqc/{}/fastqc_report.html\n'.format(
                web_path, r1))
            tf.write(tf_lines)
            tf_lines = ('{}/fastqc/{}/fastqc_report.html\n'.format(
                web_path, r2))
            tf.write(tf_lines)
    
        # Align reads.
        lines = _star_align(combined_r1, combined_r2, sample_name, rgpl,
                            rgpu, star_index, star_path, threads)
        f.write(lines)
        f.write('wait\n\n')

        # Coordinate sort bam file.
        lines = _picard_coord_sort(aligned_bam, coord_sorted_bam, picard_path,
                                   picard_memory, job.tempdir)
        f.write(lines)
        f.write('wait\n\n')

        # Mark duplicates.
        duplicate_metrics = os.path.join(
            job.outdir, '{}_duplicate_metrics.txt'.format(sample_name))
        lines = _picard_mark_duplicates(coord_sorted_bam, out_bam,
                                        duplicate_metrics,
                                        picard_path=picard_path,
                                        picard_memory=picard_memory,
                                        tempdir=job.tempdir)
        f.write(lines)
        f.write('wait\n\n')
        name = os.path.split(out_bam)[1]
        job.add_softlink(os.path.join(job.outdir, name), 
                         os.path.join(link_dir, 'bam', name))
        with open(tracklines_file, "a") as tf:
            tf_lines = ('track type=bam name="{}_bam" '
                        'description="RNAseq for {}" '
                        'bigDataUrl={}/bam/{}\n'.format(
                            sample_name, sample_name, web_path, name))
            tf.write(tf_lines)

        # Index bam file.
        lines = _picard_index(out_bam, bam_index, picard_memory / 3,
                              picard_path, job.tempdir, bg=True)
        f.write(lines)
        f.write('wait\n\n')
        name = os.path.split(bam_index)[1]
        job.add_softlink(os.path.join(job.outdir, name), 
                         os.path.join(link_dir, 'bam', name))

        # Collect insert size metrics, bam index stats, GC bias, RNA seq QC.
        lines = _picard_collect_multiple_metrics(out_bam, sample_name,
                                                 picard_path, picard_memory / 3,
                                                 job.tempdir, bg=True)
        f.write(lines)
        for fn in ['{}.{}'.format(sample_name, x) for x in 
            'alignment_summary_metrics', 
            'quality_by_cycle.pdf', 
            'base_distribution_by_cycle.pdf', 
            'quality_by_cycle_metrics', 
            'base_distribution_by_cycle_metrics', 
            'quality_distribution.pdf', 
            'insert_size_histogram.pdf', 
            'quality_distribution_metrics', 
            'insert_size_metrics'
                  ]:
            job.output_files_to_copy.append(fn)

        metrics = os.path.join(job.outdir,
                               '{}_rna_seq_metrics.txt'.format(sample_name))
        chart = os.path.join(job.outdir,
                             '{}_5_3_coverage.pdf'.format(sample_name))
        lines = _picard_collect_rna_seq_metrics(
            out_bam, 
            metrics, 
            chart, 
            sample_name,
            picard_path, 
            picard_memory / 3,
            job.tempdir,
            ref_flat, 
            rrna_intervals,
            strand_specific=strand_specific, 
            bg=True)
        f.write(lines)

        # Make md5 hash for output bam file.
        f.write('md5sum {} > {}\n\n'.format(
            out_bam, os.path.join(job.outdir, '{}.md5'.format(out_bam))))
        
        # metrics = os.path.join(job.outdir,
        #                        '{}_gc_bias_metrics.txt'.format(sample_name))
        # chart = os.path.join(job.outdir, '{}_gc_bias.pdf'.format(sample_name))
        # out = os.path.join(job.outdir, '{}_gc_bias.txt'.format(sample_name))
        # lines = _picard_gc_bias_metrics(out_bam, metrics, chart, out,
        #                                 picard_path, picard_memory / 4,
        #                                 job.tempdir, bg=False)
        # f.write(lines)

        # Make bigwig files for displaying coverage.
        # TODO: update for strand specific eventually.
        lines = _bigwig_files(out_bam, out_bigwig, sample_name,
                              bedgraph_to_bigwig_path, bedtools_path)
        f.write(lines)
        f.write('wait\n\n')
        name = os.path.split(out_bigwig)[1]
        job.add_softlink(os.path.join(job.outdir, name), 
                         os.path.join(link_dir, 'bw', name))

        with open(tracklines_file, "a") as tf:
            tf_lines = ('track type=bigWig name="{}_cov" '
                        'description="RNAseq coverage for {}" '
                        'visibility=0 db=hg19 bigDataUrl={}/bw/{}\n'.format(
                            sample_name, sample_name, web_path, name))
            tf.write(tf_lines)

        index_out = os.path.join(job.outdir,
                                 '{}_index_stats.txt'.format(sample_name))
        index_err = os.path.join(job.outdir,
                                 '{}_index_stats.err'.format(sample_name))
        lines = _picard_bam_index_stats(out_bam, index_out, index_err,
                                        picard_path, picard_memory, tempdir,
                                        bg=True)
        f.write(lines)

        # Make softlinks and tracklines for genome browser.
        # TODO: update for strand specific eventually.
        # lines = _genome_browser_files(tracklines_file, link_dir, web_path_file,
        #                               out_bam, bam_index, out_bigwig,
        #                               sample_name, job.outdir)
        f.write(lines)
        f.write('wait\n\n')
    
    job.write_end()
    return job.filename

def _dexseq_count(bam, counts_file, dexseq_annotation, paired=True,
                  strand_specific=False, samtools_path='.'):
    """
    Count reads overlapping exonic bins for DEXSeq.

    Parameters
    ----------
    bam : str
        Path to coordinate sorted bam file to count reads for.

    counts_file : str
        File to write bin counts to.

    dexseq_annotation : str
        Path to DEXSeq exonic bins GFF file.

    paired : boolean
        True if the data is paired-end. False otherwise.

    strand_specific : boolean
        True if the data is strand-specific. False otherwise.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    name : str
        File name for the softlink.

    """
    import readline
    import rpy2.robjects as robjects
    robjects.r('suppressPackageStartupMessages(library(DEXSeq))')
    scripts = robjects.r('system.file("python_scripts", package="DEXSeq")')
    g = scripts.items()
    scripts_path = g.next()[1]
    script = os.path.join(scripts_path, 'dexseq_count.py')
    if paired:
        p = 'yes'
    else:
        p = 'no'
    if strand_specific:
        s = 'reverse'
    else:
        s = 'no'
    lines = (
        '{} view -h -f 2 {} | '.format(samtools_path, bam) +
        'cut -f1-17,20- | python {} '.format(script) + 
        '-p {} -s {} -a 0 -r pos -f sam '.format(p, s) + 
        '{} - {} &\n\n'.format(dexseq_annotation, counts_file)
    )
    return lines

def _htseq_count(bam, counts_file, stats_file, gtf, samtools_path,
                 strand_specific=False):
    """
    Count reads overlapping genes for use with DESeq etc.

    Parameters
    ----------
    bam : str
        Path to coordinate sorted bam file to count reads for.

    counts_file : str
        File to write counts to.

    stats_file : str
        File to write counting stats to.

    gtf : str
        Path to GTF file to count against. Optimized for use with Gencode GTF.

    strand_specific : boolean
        True if the data is strand-specific. False otherwise.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    name : str
        File name for the softlink.

    """
    import HTSeq
    if strand_specific:
        s = 'reverse'
    else:
        s = 'no'
    script = os.path.join(HTSeq.__path__[0], 'scripts', 'count.py')
    lines = ('python {} -f bam -r pos -s {} '.format(script, s) + 
             '-a 0 -t exon -i gene_id -m union ' + 
             '{} {} > temp_out.tsv\n'.format(bam, gtf))
    lines += 'tail -n 5 temp_out.tsv > {}\n'.format(stats_file)
    lines += 'lines=$(wc -l <temp_out.tsv)\n'
    lines += 'wanted=`expr $lines - 5`\n'
    lines += 'head -n $wanted temp_out.tsv > {}\n'.format(counts_file)
    lines += 'rm temp_out.tsv\n\n'

    return lines

def get_counts(bam, outdir, sample_name, tempdir, dexseq_annotation, gtf,
               samtools_path, conda_env='', r_env='', paired=True,
               strand_specific=False, shell=False):
    """
    Make a PBS or shell script for counting reads that overlap genes for DESeq2
    and exonic bins for DEXSeq.

    Parameters
    ----------
    bam : str
        Coordinate sorted bam file (genomic coordinates).

    outdir : str
        Directory to store PBS/shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    tempdir : str
        Directory to store temporary files.

    dexseq_annotation : str
        Path to DEXSeq exonic bins GFF file.

    gtf : str
        Path to GTF file to count against. Optimized for use with Gencode GTF.

    conda_env : str
        If provided, load conda environment with this name.

    r_env : str
        If provided, this file will be sourced to set the environment for rpy2.

    paired : boolean
        True if the data is paired-end. False otherwise.

    strand_specific : boolean
        True if the data is strand-specific. False otherwise.

    shell : boolean
        If true, make a shell script rather than a PBS script.

    """
    threads = 6

    if shell:
        pbs = False
    else: 
        pbs = True

    tempdir = os.path.join(tempdir, '{}_counts'.format(sample_name))
    outdir = os.path.join(outdir, '{}_counts'.format(sample_name))

    # I'm going to define some file names used later.
    temp_bam = os.path.join(tempdir, os.path.split(bam)[1])
    dexseq_counts = os.path.join(outdir, 'dexseq_counts.tsv')
    gene_counts = os.path.join(outdir, 'gene_counts.tsv')
    gene_count_stats = os.path.join(outdir, 'gene_count_stats.tsv')
    
    # Files to copy to output directory.
    files_to_copy = []
    
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = []

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_counts.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_counts.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir, '{}_counts.out'.format(sample_name))
        err = os.path.join(outdir, '{}_counts.err'.format(sample_name))
        job_name = '{}_counts'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))

    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))
    f.write('rsync -avz {} .\n\n'.format(bam))

    if conda_env != '':
        f.write('source activate {}\n\n'.format(conda_env))
    if r_env != '':
        f.write('source {}\n\n'.format(r_env))

    lines = _dexseq_count(temp_bam, dexseq_counts, dexseq_annotation,
                          paired=True, strand_specific=strand_specific,
                          samtools_path=samtools_path)
    f.write(lines)
    lines = _htseq_count(temp_bam, gene_counts, gene_count_stats, gtf,
                         strand_specific=strand_specific,
                         samtools_path=samtools_path)
    f.write(lines)
    f.write('wait\n\n')
    
    if len(files_to_copy) > 0:
        f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
            ' \\\n\t'.join(files_to_copy),
            outdir))
    if len(files_to_remove) > 0:
        f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))
    f.close()

    return fn

def _rsem_calculate_expression(bam, reference, rsem_path, sample_name,
                               threads=1, ci_mem=1024, strand_specific=False): 
    """
    Estimate expression using RSEM.

    Parameters
    ----------
    bam : str
        Transcriptome bam file.

    reference : str
        RSEM reference.

    rsem_path : str
        Path to RSEM executables.

    sample_name : str
        Sample name for RSEM to name files.

    ci_mem : int
        Amount of memory in mb to give RSEM for calculating confidence
        intervals. Passed to --ci-memory for RSEM.

    strand_specific : boolean
        True if the data is strand-specific. False otherwise. For now, this
        means that the R1 read is on the reverse strand.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    name : str
        File name for the softlink.

    """
    line = ('{}/rsem-calculate-expression --bam --paired-end --num-threads {} '
            '--no-bam-output --seed 3272015 --calc-pme --calc-ci '
            '--ci-memory {} --estimate-rspd {} {} {}'.format(
                rsem_path, threads, ci_mem, bam, reference, sample_name))
    if strand_specific:
        line += ' --forward-prob 0'
    line += '\n'
    return line

def rsem_expression(bam, outdir, sample_name, tempdir, rsem_path,
                    rsem_reference, ci_mem=1024, r_env='', threads=32,
                    strand_specific=False, shell=False):
    """
    Make a PBS or shell script for estimating expression using RSEM.

    Parameters
    ----------
    bam : str
        Coordinate sorted bam file (genomic coordinates).

    outdir : str
        Directory to store PBS/shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    tempdir : str
        Directory to store temporary files.

    rsem_path : str
        Path to directory with RSEM executables.

    rsem_reference : str
        RSEM reference.

    ci_mem : int
        Amount of memory in mb to give RSEM for calculating confidence
        intervals. Passed to --ci-memory for RSEM.

    r_env : str
        If provided, this file will be sourced to set the environment for rpy2.

    strand_specific : boolean
        True if the data is strand-specific. False otherwise.

    shell : boolean
        If true, make a shell script rather than a PBS script.

    """
    if shell:
        pbs = False
    else: 
        pbs = True

    tempdir = os.path.join(tempdir, '{}_rsem'.format(sample_name))
    outdir = os.path.join(outdir, '{}_rsem'.format(sample_name))

    # I'm going to define some file names used later.
    temp_bam = os.path.join(tempdir, os.path.split(bam)[1])
    
    # Files to copy to output directory.
    files_to_copy = ['{}.genes.results'.format(sample_name),
                     '{}.isoforms.results'.format(sample_name),
                     '{}.stat'.format(sample_name)]
    
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [temp_bam]

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_rsem.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_rsem.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir, '{}_rsem.out'.format(sample_name))
        err = os.path.join(outdir, '{}_rsem.err'.format(sample_name))
        job_name = '{}_rsem'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))

    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))
    f.write('rsync -avz {} .\n\n'.format(bam))

    lines = _rsem_calculate_expression(temp_bam, rsem_reference, rsem_path,
                                       sample_name, threads=threads,
                                       ci_mem=ci_mem,
                                       strand_specific=strand_specific)
    f.write(lines)
    f.write('wait\n\n')
    
    if len(files_to_copy) > 0:
        f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
            ' \\\n\t'.join(files_to_copy),
            outdir))
    if len(files_to_remove) > 0:
        f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if os.path.realpath(tempdir) != os.path.realpath(outdir):
        f.write('rm -r {}\n'.format(tempdir))
    f.close()

    return fn
