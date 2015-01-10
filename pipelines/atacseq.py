import os

from general import _bedgraph_to_bigwig
from general import _bigwig_files
from general import _fastqc
from general import _flagstat
from general import _make_softlink
from general import _pbs_header
from general import _picard_remove_duplicates
from general import _process_fastqs
from general import _samtools_index

def _star_align(r1_fastqs, r2_fastqs, sample, rgpl, rgpu, star_index, star_path,
                threads):
    """
    Align paired ATAC-seq fastq files with STAR.

    Parameters
    ----------
    r1_fastqs : str
        Gzipped R1 fastq file(s). If multiple files, each file should be
        separated by s space and should be ordered the same as the R2 files.

    r2_fastqs : str
        Gzipped R2 fastq file(s). If multiple files, each file should be
        separated by s space and should be ordered the same as the R1 files.

    sample : str
        Sample name.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    """
    # I use threads - 2 for STAR so there are open processors for reading and
    # writing. The parameters here are set in part based on the discussion here:
    # https://groups.google.com/forum/#!searchin/rna-star/chip$20seq/rna-star/E_mKqm9jDm0/ZpB6yRcWi60J.
    # --alignIntronMax 1 disables spliced alignments. --alignEndsType EndToEnd
    # prohibits soft clipping. I may or may not want to prohibit soft clipping.
    # Smaller-sized fragments may still have adapter sequence and won't align
    # without soft clipping. I think that since our reads are long, we should be
    # ok soft clipping since they should align pretty well. I'm only outputting 
    # unique alignments right now.
    line = (' \\\n'.join([star_path, 
                          '\t--runThreadN {}'.format(threads - 2),
                          '\t--alignIntronMax 1',
                          '\t--genomeDir {}'.format(star_index), 
                          '\t--genomeLoad NoSharedMemory', 
                          '\t--readFilesCommand zcat',
                          '\t--readFilesIn {} {}'.format(r1_fastqs, 
                                                         r2_fastqs),
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
                          '\t--seedSearchStartLmax 20']) + '\n\n')
    return line

def _picard_coord_sort_primary(in_bam, out_bam, picard_path, picard_memory,
                               samtools_path, temp_dir): 
    """
    Coordinate sort using Picard Tools while only keeping primary alignments.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out_bam : str
        Path to output bam file.

    """
    lines = (' \\\n'.join(['{} view -hu -F 256 {} | '.format(samtools_path,
                                                             in_bam),
                           '\tjava -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(temp_dir), 
                           '\t-jar {} SortSam'.format(picard_path),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tI=/dev/stdin',
                           '\tO={}'.format(out_bam),
                           '\tSO=coordinate\n']))
    return lines

def _genome_browser_files(tracklines_file, link_dir, web_path_file,
                          coord_sorted_bam, bam_index, bigwig, r1_fastqc,
                          r2_fastqc, narrow_peak, sample_name, out_dir):
    """
    Make files and softlinks for displaying results on UCSC genome browser.

    Parameters
    ----------
    tracklines_file : str
        Path to file for writing tracklines. The tracklines will be added to the
        file; the contents of the file will not be overwritten. These tracklines
        can be pasted into the genome browser upload for custom data.

    link_dir : str
        Path to directory where softlink should be made.

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

    coord_sorted_bam : str
        Path to coordinate sorted bam file.

    bam_index : str
        Path to index file for coordinate sorted bam file.

    bigwig : str
        Path to bigwig file. 

    r1_fastqc : str
        Path to fastqc_report.html for R1 reads.

    r2_fastqc : str
        Path to fastqc_report.html for R2 reads.

    narrow_peak : str
        Path to narrowPeak file from MACS2.

    sample_name : str
        Sample name used for naming files.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = ''
    link_dir = os.path.join(link_dir, 'atac')

    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()
    web_path = web_path + '/atac'

    # File with UCSC tracklines.
    if os.path.exists(tracklines_file):
        with open(tracklines_file) as f:
            tf_lines = f.read()
    else:
        tf_lines = ''

    # FastQC results.
    temp_link_dir = os.path.join(link_dir, 'fastqc')
    temp_web_path = web_path + '/fastqc'
    try:
        os.makedirs(temp_link_dir)
    except OSError:
        pass
    new_lines, r1_name = _make_softlink(r1_fastqc, sample_name + '_R1',
                                        temp_link_dir)
    lines += new_lines
    new_lines, r2_name = _make_softlink(r2_fastqc, sample_name + '_R2',
                                        temp_link_dir)
    lines += new_lines

    # Bam file and index.
    temp_link_dir = os.path.join(link_dir, 'bam')
    temp_web_path = web_path + '/bam'
    try:
        os.makedirs(temp_link_dir)
    except OSError:
        pass
    fn = os.path.join(out_dir, os.path.split(coord_sorted_bam)[1])
    new_lines, bam_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines

    fn = os.path.join(out_dir, os.path.split(bam_index)[1])
    new_lines, index_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines

    tf_lines += ' '.join(['track', 'type=bam',
                          'name="{}_bam"'.format(sample_name),
                          'description="ATAC-seq for {}"'.format(sample_name),
                          'visibility=0',
                          'db=hg19',
                          'bigDataUrl={}/{}\n'.format(temp_web_path, bam_name)])
    
    # Bigwig file(s).
    temp_link_dir = os.path.join(link_dir, 'bw')
    temp_web_path = web_path + '/bw'
    try:
        os.makedirs(temp_link_dir)
    except OSError:
        pass
    fn = os.path.join(out_dir, os.path.split(bigwig)[1])
    new_lines, bigwig_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines

    tf_lines += ' '.join(['track', 'type=bigWig',
                          'name="{}_cov"'.format(sample_name),
                          ('description="ATAC-seq coverage for '
                           '{}"'.format(sample_name)),
                          'visibility=0',
                          'db=hg19',
                          'bigDataUrl={}/{}\n'.format(temp_web_path,
                                                      bigwig_name)])
   
    # Peaks.
    temp_link_dir = os.path.join(link_dir, 'peak')
    temp_web_path = web_path + '/peak'
    try:
        os.makedirs(temp_link_dir)
    except OSError:
        pass
    fn = os.path.join(out_dir, os.path.split(bigwig)[1])
    new_lines, bigwig_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines

    with open(tracklines_file, 'w') as tf:
        tf.write(tf_lines)
    
    lines += '\n'
    return lines

def _macs2(bam, sample_name, out_dir):
    """
    Call peaks with MACS2. The macs2 executable is assumed to be in your path
    which it should be if you installed it using pip install MACS2 into your
    python environment.

    Parameters
    ----------
    bam : str:
        Path to paired-end, coordinate sorted bam file.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = ('macs2 callpeak -t {} -f BAMPE '.format(bam) + 
             '-g hs -n {} --outdir {}'.format(sample_name, out_dir) + 
             '--call-summits\n')
    
    # Add trackline to narrowPeak file from macs.
    out = os.path.join(out_dir, '{}_peaks.narrowPeak'.format(sample_name))
    temp = os.path.join(out_dir, 'temp.narrowPeak')
    track_line = ' '.join(['track', 'type=narrowPeak',
                           'name="{}_peaks"'.format(sample_name),
                           ('description="ATAC-seq peaks for '
                            '{}"'.format(sample_name)),
                           'visibility=0',
                           'db=hg19'])
    lines += 'cat <(echo "{}") {} > {}\n'.format(track_line, out, temp)
    lines += 'mv {} {}\n\n'.format(temp, out)
    return lines

def align_and_sort(
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
    conda_env='',
    rgpl='ILLUMINA',
    rgpu='',
    temp_dir='/scratch', 
    threads=32, 
    picard_memory=58, 
    shell=False
):
    """
    Make a PBS or shell script for aligning ATAC-seq reads with STAR. The
    defaults are set for use on the Frazer lab's PBS scheduler on FLC.

    Parameters
    ----------
    r1_fastqs : list or str
        Either a list of paths to gzipped fastq files with R1 reads or path to a
        single gzipped fastq file with R1 reads.

    r2_fastqs : list or str
        Either a list of paths to gzipped fastq files with R2 reads or path to a
        single gzipped fastq file with R2 reads.

    out_dir : str
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

    fastqc_path : str
        Path to FastQC executable.

    samtools_path : str
        Path to samtools executable.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of MACS2 is used.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    temp_dir : str
        Directory to store files as STAR runs.

    threads : int
        Number of threads to reserve using PBS scheduler. This number of threads
        minus 2 will be used by STAR, so this must be at least 3.

    picard_memory : int
        Amount of memory (in gb) to give Picard Tools.

    shell : boolean
        If true, make a shell script rather than a PBS script.

    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    assert threads >= 3

    if shell:
        pbs = False
    else: 
        pbs = True

    if type(r1_fastqs) == str:
        r1_fastqs = [r1_fastqs]
    if type(r2_fastqs) == str:
        r2_fastqs = [r2_fastqs]

    temp_dir = os.path.join(temp_dir, '{}_alignment'.format(sample_name))
    out_dir = os.path.join(out_dir, '{}_alignment'.format(sample_name))

    # I'm going to define some file names used later.
    temp_r1_fastqs = _process_fastqs(r1_fastqs, temp_dir)
    temp_r2_fastqs = _process_fastqs(r2_fastqs, temp_dir)
    combined_r1 = os.path.join(temp_dir, 
                               '{}_combined_R1.fastq.gz'.format(sample_name))
    combined_r2 = os.path.join(temp_dir, 
                               '{}_combined_R2.fastq.gz'.format(sample_name))
    aligned_bam = os.path.join(temp_dir, 'Aligned.out.bam')
    coord_sorted_bam = os.path.join(temp_dir, 'Aligned.out.coord.sorted.bam')
    no_dup_bam = os.path.join(temp_dir,
                              '{}_atac_no_dup.bam'.format(sample_name))
    bam_index = os.path.join(temp_dir,
                             '{}_atac_no_dup.bam.bai'.format(sample_name))
    out_bigwig = os.path.join(temp_dir, '{}_atac.bw'.format(sample_name))
    
    duplicate_metrics = os.path.join(out_dir, 'duplicate_metrics.txt')
    stats_file = os.path.join(out_dir,
                              '{}_atac_no_dup.bam.flagstat'.format(sample_name))
    chrom_counts = os.path.join(out_dir,
                                '{}_chrom_counts.txt'.format(sample_name))
    narrow_peak = os.path.join(out_dir,
                               '{}_peaks.narrowPeak'.format(sample_name))
    
    tn = os.path.split(combined_r1)[1]
    r1_fastqc = os.path.join(
        out_dir, 
        '{}_fastqc'.format('.'.join(tn.split('.')[0:-2])),
        'fastqc_report.html')
    tn = os.path.split(combined_r2)[1]
    r2_fastqc = os.path.join(
        out_dir, 
        '{}_fastqc'.format('.'.join(tn.split('.')[0:-2])),
        'fastqc_report.html')
    
    # Files to copy to output directory.
    files_to_copy = [no_dup_bam, bam_index, 'Log.out', 'Log.final.out',
                     'Log.progress.out', 'SJ.out.tab', out_bigwig]
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [combined_r1, combined_r2]

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(out_dir, '{}_alignment.sh'.format(sample_name))
    else:
        fn = os.path.join(out_dir, '{}_alignment.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(out_dir, '{}_alignment.out'.format(sample_name))
        err = os.path.join(out_dir, '{}_alignment.err'.format(sample_name))
        job_name = '{}_align'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))
    
    if conda_env != '':
        f.write('source activate {}\n'.format(conda_env))
    f.write('mkdir -p {}\n'.format(temp_dir))
    f.write('cd {}\n'.format(temp_dir))
    f.write('rsync -avz \\\n{} \\\n{} \\\n\t.\n\n'.format(
        ' \\\n'.join(['\t{}'.format(x) for x in r1_fastqs]),
        ' \\\n'.join(['\t{}'.format(x) for x in r2_fastqs])))
    
    # Combine fastq files and run FastQC.
    f.write('cat \\\n{} \\\n\t> {} &\n'.format(
        ' \\\n'.join(['\t{}'.format(x) for x in temp_r1_fastqs]),
        combined_r1))
    f.write('cat \\\n{} \\\n\t> {}\n\n'.format(
        ' \\\n'.join(['\t{}'.format(x) for x in temp_r2_fastqs]),
        combined_r2))
    f.write('wait\n\n')
    f.write('rm \\\n{} \\\n{}\n\n'.format(
        ' \\\n'.join(['\t{}'.format(x) for x in temp_r1_fastqs]),
        ' \\\n'.join(['\t{}'.format(x) for x in temp_r2_fastqs])))
    f.write('wait\n\n')
    lines = _fastqc([combined_r1, combined_r2], threads, out_dir, fastqc_path)
    f.write(lines)
    f.write('wait\n\n')

    # Align with STAR.
    lines = _star_align(combined_r1, combined_r2, sample_name, rgpl,
                        rgpu, star_index, star_path, threads)
    f.write(lines)
    f.write('wait\n\n')

    # Coordinate sort bam file.
    lines = _picard_coord_sort_primary(aligned_bam, coord_sorted_bam,
                                       picard_path, picard_memory,
                                       samtools_path, temp_dir)
    f.write(lines)
    f.write('wait\n\n')

    # Remove duplicates.
    lines = _picard_remove_duplicates(coord_sorted_bam, no_dup_bam,
                                      duplicate_metrics, picard_path,
                                      picard_memory, temp_dir)
    f.write(lines)
    f.write('wait\n\n')

    # Index bam file and collect flagstats.
    lines = _samtools_index(no_dup_bam, samtools_path)
    f.write(lines)
    lines = _flagstat(no_dup_bam, stats_file, samtools_path)
    f.write(lines)

    # Count the number of reads aligned to each chromosome.
    f.write('{} view {} | cut -f 3 | uniq -c > {} &\n\n'.format(
        samtools_path, no_dup_bam, chrom_counts))

    # Make bigwig files for displaying coverage.
    lines = _bigwig_files(no_dup_bam, out_bigwig, sample_name,
                          bedgraph_to_bigwig_path, bedtools_path)
    f.write(lines)

    # Call peaks
    lines = _macs2(no_dup_bam, sample_name, out_dir)
    f.write(lines)
    
    f.write('wait\n\n')

    # Make softlinks and tracklines for genome browser.
    lines = _genome_browser_files(tracklines_file, link_dir, web_path_file,
                                  no_dup_bam, bam_index, out_bigwig, r1_fastqc,
                                  r2_fastqc, narrow_peak, sample_name, out_dir)
    f.write(lines)
    f.write('wait\n\n')

    f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
        ' \\\n\t'.join(files_to_copy),
        out_dir))
    f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if temp_dir != out_dir:
        f.write('rm -r {}\n'.format(temp_dir))
    f.close()

    return fn
