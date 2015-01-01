import os

def _pbs_header(out, err, name, threads, queue='high'):
    """
    Write header for PBS script

    Parameters
    ----------
    out : str
        Path to file for recording standard out.

    err : str
        Path to file for recording standard error.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = ('\n'.join(['#PBS -q {}'.format(queue),
                        '#PBS -N {}'.format(name),
                        '#PBS -l nodes=1:ppn={}'.format(threads),
                        '#PBS -o {}'.format(out),
                        '#PBS -e {}\n\n'.format(err)]))
    return lines

def _cbarrett_paired_dup_removal(r1_fastqs, r2_fastqs, r1_nodup, r2_nodup,
                                 temp_dir):
    """
    Remove duplicates from paired fastq files using UNIX sort and uniq. Read 
    pairs with exactly the same sequences are removed such that every read pair
    has a different sequence. 

    Parameters
    ----------
    r1_fastqs : str
        R1 fastq file(s). If multiple files, each file should be separated by s
        space and should be ordered the same as the R2 files.

    r2_fastqs : str
        R2 fastq file(s). If multiple files, each file should be separated by s
        space and should be ordered the same as the R1 files.

    r1_nodup : str
        Path to write gzipped R1 fastq file with duplicates removed.

    r2_nodup : str
        Path to write gzipped R2 fastq file with duplicates removed.

    temp_dir : str
        Path to temporary directory where fastq files will be copied to.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = []
    lines.append('paste \\\n')
    lines.append('<(zcat {} | '.format(r1_fastqs) + 
                 'awk \'0==(NR+3)%4{ORS=" "; split($0,a," "); ' + 
                 'print substr(a[1],2)}0==(NR+2)%4{print} (NR!=1 && 0==NR%4)' + 
                 '{ORS="\\n";print}\') \\\n')
    lines.append('<(zcat *_R2_* | '.format(r2_fastqs) + 
                 'awk \'0==(NR+3)%4{ORS=" "; split($0,a," "); ' + 
                 'print substr(a[1],2)}0==(NR+2)%4{print} (NR!=1 && 0==NR%4)' + 
                 '{ORS="\\n";print}\') | \\\n')
    lines.append('awk \'{if ($2 < $5) printf "%s %s %s %s %s %s\\n",'
                 '$1,$3,$4,$6,$2,$5; else printf "%s %s %s %s %s %s\\n",'
                 '$1,$6,$4,$3,$5,$2}\' | \\\n')
    lines.append('sort -k 5,5 -k 6,6 -T {0} -S 30G --parallel=8 | '
                 'uniq -f 4 | \\\n'.format(temp_dir))
    lines.append('awk \'{printf "@%s\\n%s\\n+\\n%s\\n",$1,$5,$2 | '
                 '"gzip -c > ' + r1_nodup + 
                 '"; printf "@%s\\n%s\\n+\\n%s\\n",$3,$6,$4 | "gzip -c > ' + 
                  r2_nodup + '"}\'\n\n')
    return ''.join(lines)

def _star_align(r1_fastqs, r2_fastqs, sample, rgpl, rgpu, star_index, star_path,
                threads):
    """
    Align paired fastq files with STAR.

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
    # writing.
    line = (' \\\n'.join([star_path, 
                          '\t--runThreadN {}'.format(threads - 2),
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
                          ('\t--outFilterIntronMotifs '
                           'RemoveNoncanonicalUnannotated'),
                           '\t--outSJfilterOverhangMin 6 6 6 6',
                           '\t--seedSearchStartLmax 20',
                           '\t--alignSJDBoverhangMin 1', 
                          '\t--quantMode TranscriptomeSAM']) + '\n\n') 
    return line

def _picard_coord_sort(in_bam, out_bam, picard_path, picard_memory, temp_dir):
    """
    Coordinate sort using Picard Tools.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out_bam : str
        Path to output bam file.

    """
    line = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                          '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                          '\t-Djava.io.tmpdir={}'.format(temp_dir),
                          '\t-jar {} SortSam'.format(picard_path),
                          '\tVALIDATION_STRINGENCY=SILENT',
                          '\tI={}'.format(in_bam),
                          '\tO={}'.format(out_bam),
                          '\tSO=coordinate\n\n']))
    return line

def _picard_index(in_bam, index, picard_memory, picard_path, temp_dir):
    """
    Index bam file using Picard Tools.

    Parameters
    ----------
    in_bam : str
        Path to file input bam file.

    index : str
        Path to index file for input bam file.

    Returns
    -------
    index : str
        Path to index file for input bam file.

    """
    line = (' \\\n'.join(['java -Xmx{}g -jar'.format(picard_memory),
                          '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                          '\t-Djava.io.tmpdir={}'.format(temp_dir),
                          '\t-jar {} BuildBamIndex'.format(picard_path),
                          '\tI={}'.format(in_bam),
                          '\tO={}\n\n'.format(index)]))
    return line

def _bedgraph_to_bigwig(bedgraph, bigwig, bedgraph_to_bigwig_path,
                        bedtools_path):
    bedtools_genome_path = os.path.join(
        os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
        'human.hg19.genome')
    lines =  ' '.join(['{} {}'.format(bedgraph_to_bigwig_path, bedgraph),
                       '{}'.format(bedtools_genome_path),
                       '{} &\n'.format(bigwig)])
    return lines

def _coverage_bedgraph(bam, bedgraph, bedtools_path, sample_name, strand='.'):
    """
    Make lines that create a coverage bedgraph file.

    Parameters
    ----------
    bam : str
        Bam file to calculate coverage for.

    bedgraph : str
        Path to output bedgraph file.

    bedtools_path : str
        Path to bedtools.

    sample_name : str
        Sample name for naming files etc.

    strand : str
        If '+' or '-', calculate strand-specific coverage. Otherwise, calculate
        coverage using all reads.

    Returns
    -------
    lines : str
        Lines to be written to PBS/shell script.

    """
    if strand == '+' or strand == '-':
        if strand == '+':
            name = '{}_plus'.format(sample_name)
        else:
            name = '{}_minus'.format(sample_name)
        lines = ' \\\n'.join(['{} genomecov -ibam'.format(bedtools_path),
                              '\t{}'.format(bam),
                              '\t-g hg19.genome -split -bg ',
                              '\t-strand {} -trackline'.format(strand),
                              '\t-trackopts \'name="{}"\''.format(name),
                              '\t> {} &\n\n'.format(bedgraph)])
    else:
        name = sample_name
        lines = ' \\\n'.join(['{} genomecov -ibam'.format(bedtools_path),
                              '\t{}'.format(bam),
                              '\t-g hg19.genome -split -bg ',
                              '\t-trackline'.format(strand),
                              '\t-trackopts \'name="{}"\''.format(name),
                              '\t> {} &\n\n'.format(bedgraph)])
    return lines

def _bigwig_files(in_bam, out_bigwig, sample_name, bedgraph_to_bigwig_path,
                  bedtools_path, out_bigwig_minus=''):
    """
    Make bigwig coverage files.

    Parameters
    ----------
    in_bam : str
        Path to bam file to create bigwigs for.

    out_bigwig : str
        Path to output bigwig file. If out_bigwig_minus is provided, out_bigwig
        has the plus strand coverage.

    out_bigwig_minus : str
        Path to output bigwig file for minus strand. If out_bigwig_minus is not
        provided, the coverage is calculated using reads from both strands and
        written to out_bigwig.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = ''
    if out_bigwig_minus != '':
        lines += _coverage_bedgraph(in_bam, 'plus.bg', bedtools_path,
                                    sample_name, strand='+')
        lines += _coverage_bedgraph(in_bam, 'minus.bg', bedtools_path,
                                    sample_name, strand='-')
        lines += ('wait\n\n')
        lines += (_bedgraph_to_bigwig('plus.bg', out_bigwig,
                                      bedgraph_to_bigwig_path, bedtools_path))
        lines += (_bedgraph_to_bigwig('minus.bg', out_bigwig_minus,
                                      bedgraph_to_bigwig_path, bedtools_path))
        lines += ('\nwait\n\n')
        lines += ('rm plus.bg minus.bg\n\n')
    
    else:
        lines = _coverage_bedgraph(in_bam, 'both.bg', bedtools_path,
                                   sample_name)
        lines += (_bedgraph_to_bigwig('both.bg', out_bigwig_minus,
                                      bedgraph_to_bigwig_path, bedtools_path))
        lines += ('wait\n\n')
        lines += ('rm both.bg\n\n')
    return lines
    
def _process_fastqs(fastqs, temp_dir):
    """
    Create list of temporary fastq paths.

    Parameters
    ----------
    fastqs : list or str
        Either a list of paths to gzipped fastq files or path to a single
        gzipped fastq file.

    temp_dir : str
        Path to temporary directory where fastq files will be copied to.

    Returns
    -------
    fastqs : str
        Paths to original fastq files (concatenated with a space if multiple,
        e.g. 'path/fq1.fastq path/fq2.fastq').

    temp_fastqs : str
        Paths to temporary fastq files (concatenated with a space if multiple,
        e.g. 'tempdir/fq1.fastq tempdir/fq2.fastq').

    """
    if type(fastqs) == list:
        fns = [os.path.split(x)[1] for x in fastqs]
        temp_fastqs = [os.path.join(temp_dir, x) for x in fns]
        fastqs = ' '.join(fastqs)
    elif type(fastqs) == str:
        temp_fastqs = os.path.join(temp_dir, os.path.split(fastqs)[1])
    return fastqs, temp_fastqs

def _make_softlink(fn, sample_name, link_dir):
    """
    Make softlink for file fn in link_dir. sample_name followed by an underscore
    will be appended to the name of fn.

    Parameters
    ----------
    fn : str
        Full path to file to make link to.

    sample_name : str
        Sample name used for naming files.

    link_dir : str
        Path to directory where softlink should be made.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    name : str
        File name for the softlink.

    """
    name = '{}_'.format(os.path.split(fn)[1])
    lines = ('ln -s {} {}\n'.format(fn,
                                    os.path.join(link_dir, fn)))
    return lines, name


def _genome_browser_files(tracklines_file, link_dir, web_path_file,
                          coord_sorted_bam, bam_index, bigwig, sample_name,
                          bigwig_minus=''):
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
        Path to bigwig file. If bigwig_minus is provided, bigwig has the plus
        strand coverage.

    sample_name : str
        Sample name used for naming files.

    bigwig_minus : str
        Path to bigwig file for minus strand. If bigwig_minus is not provided,
        bigwig is assumed to have coverage for both plus and minus stand reads.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()

    # File with UCSC tracklines.
    if os.path.exists(tracklines_file):
        with open(tracklines_file) as f:
            lines = f.read()
    else:
        lines = ''
    tf = open(tracklines_file, 'w')
    
    # Bam file and index.
    new_lines, bam_name = _make_softlink(coord_sorted_bam, sample_name,
                                         link_dir)
    lines += new_lines
    new_lines, index_name = _make_softlink(bam_index, sample_name, link_dir)
    lines += new_lines
    tf.write(' '.join(['track', 'type=bam', 'name="{}_bam"'.format(sample_name),
                       'description="RNAseq for {}"'.format(sample_name),
                       'bigDataUrl={}/{}\n'.format(web_path, bam_name)]))
    
    # Bigwig file(s).
    if bigwig_minus != '':
        new_lines, plus_name = _make_softlink(bigwig, sample_name, link_dir)
        lines += new_lines
        new_lines, minus_name = _make_softlink(bigwig_minus, sample_name,
                                               link_dir)
        lines += new_lines
        tf.write(' '.join(['track', 'type=bigWig',
                           'name="{}_plus_cov"'.format(sample_name),
                           ('description="RNAseq plus strand coverage for '
                            '{}"'.format(sample_name)),
                           'bigDataUrl={}/{}\n'.format(web_path, plus_name)]))
        tf.write(' '.join(['track', 'type=bigWig',
                           'name="{}_minus_cov"'.format(sample_name),
                           ('description="RNAseq minus strand coverage for '
                            '{}"'.format(sample_name)),
                           'bigDataUrl={}/{}\n'.format(web_path, minus_name)]))
    else:
        new_lines, bigwig_name = _make_softlink(bigwig, sample_name, link_dir)
        lines += new_lines
        tf.write(' '.join(['track', 'type=bigWig',
                           'name="{}_cov"'.format(sample_name),
                           ('description="RNAseq coverage for '
                            '{}"'.format(sample_name)),
                           'bigDataUrl={}/{}\n'.format(web_path, bigwig_name)]))
    tf.close()
    lines += '\n'
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
    rgpl='ILLUMINA',
    rgpu='',
    star_path='',
    picard_path='',
    bedtools_path='',
    bedgraph_to_bigwig_path='',
    temp_dir='/scratch', 
    threads=32, 
    picard_memory=58, 
    remove_dup=True, 
    strand_specific=False, 
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

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    star_path : str
        Path to STAR aligner. If not provided, assumed to be in your path.

    picard_path : str
        Path to Picard tools. If not provided, assumed to be in your path.

    bedtools_path : str
        Path to bedtools. If not provided, assumed to be in your path.

    bedgraph_to_bigwig_path : str
        Path bedGraphToBigWig executable.

    temp_dir : str
        Directory to store files as STAR runs.

    threads : int
        Number of threads to reserve using PBS scheduler. This number of threads
        minus 2 will be used by STAR, so this must be at least 3.

    picard_memory : int
        Amount of memory (in gb) to give Picard Tools.

    remove_dup : boolean
        Whether to remove duplicate reads prior to alignment.

    strand_specific : boolean
        If true, make strand specific bigwig files. 

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

    temp_dir = os.path.join(temp_dir, '{}_alignment'.format(sample_name))
    out_dir = os.path.join(out_dir, '{}_alignment'.format(sample_name))

    # I'm going to define some file names used later.
    r1_fastqs, temp_r1_fastqs = _process_fastqs(r1_fastqs, temp_dir)
    r2_fastqs, temp_r2_fastqs = _process_fastqs(r2_fastqs, temp_dir)
    r1_nodup = os.path.join(temp_dir, 'nodup_R1.fastq.gz')
    r2_nodup = os.path.join(temp_dir, 'nodup_R2.fastq.gz')
    aligned_bam = os.path.join(temp_dir, 'Aligned.out.bam')
    coord_sorted_bam = os.path.join(temp_dir, 'Aligned.out.coord.sorted.bam')
    bam_index = os.path.join(temp_dir, 'Aligned.out.coord.sorted.bam.bai')
    
    # Files to copy to output directory.
    files_to_copy = [coord_sorted_bam, bam_index, 'Log.out', 'Log.final.out',
                     'Log.progress.out', 'SJ.out.tab',
                     'Aligned.toTranscriptome.out.bam']
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [temp_r1_fastqs, temp_r2_fastqs, r1_nodup, r2_nodup]

    if strand_specific:
        out_bigwig_plus = os.path.join(temp_dir,
                                       '{}_plus.bw'.format(sample_name))
        out_bigwig_minus = os.path.join(temp_dir,
                                        '{}_minus.bw'.format(sample_name))
        files_to_copy.append(out_bigwig_plus)
        files_to_copy.append(out_bigwig_minus)
    else:
        out_bigwig = os.path.join(temp_dir, '{}.bw'.format(sample_name))
        files_to_copy.append(out_bigwig)

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    fn = os.path.join(out_dir, '{}_alignment.pbs'.format(sample_name))
    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(out_dir, '{}_alignment.out'.format(sample_name))
        err = os.path.join(out_dir, '{}_alignment.err'.format(sample_name))
        job_name = '{}_align'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))

    f.write('mkdir -p {}\n'.format(temp_dir))
    f.write('cd {}\n'.format(temp_dir))
    f.write('rsync -avz {} {} .\n\n'.format(r1_fastqs, r2_fastqs))

    # Remove duplicates if desired and align.
    if remove_dup:
        lines = _cbarrett_paired_dup_removal(temp_r1_fastqs, temp_r2_fastqs,
                                             r1_nodup, r2_nodup, temp_dir)
        f.write(lines)
        f.write('wait\n\n')

        lines = _star_align(r1_nodup, r2_nodup, sample_name, rgpl, rgpu,
                            star_index, star_path, threads)
        f.write(lines)
        f.write('wait\n\n')
    else:
        lines = _star_align(temp_r1_fastqs, temp_r2_fastqs, sample_name, rgpl,
                            rgpu, star_index, star_path, threads)
        f.write(lines)
        f.write('wait\n\n')

    # Coordinate sort bam file.
    lines = _picard_coord_sort(aligned_bam, coord_sorted_bam, picard_path,
                               picard_memory, temp_dir)
    f.write(lines)
    f.write('wait\n\n')
    # Index coordinate sorted bam file.
    lines = _picard_index(coord_sorted_bam, bam_index, picard_memory,
                          picard_path, temp_dir)
    f.write(lines)
    f.write('wait\n\n')

    # Make bigwig files for displaying coverage.
    if strand_specific:
        lines = _bigwig_files(coord_sorted_bam, out_bigwig_plus, sample_name,
                              bedgraph_to_bigwig_path, bedtools_path,
                              out_bigwig_minus=out_bigwig_minus)
    else:
        lines = _bigwig_files(coord_sorted_bam, out_bigwig, sample_name,
                              bedgraph_to_bigwig_path, bedtools_path)
    f.write(lines)
    f.write('wait\n\n')

    # Make softlinks and tracklines for genome browser.
    if strand_specific:
        lines = _genome_browser_files(tracklines_file, link_dir, web_path_file,
                                      coord_sorted_bam, bam_index,
                                      out_bigwig_plus, sample_name,
                                      bigwig_minus=out_bigwig_minus)
    else:
        lines = _genome_browser_files(tracklines_file, link_dir, web_path_file,
                                      coord_sorted_bam, bam_index,
                                      out_bigwig, sample_name)
    f.write(lines)
    f.write('wait\n\n')

    f.write('rsync -avz {} {}'.format(' '.join(files_to_copy), out_dir))
    f.write('rm {}'.format(' '.join(files_to_remove)))

    if temp_dir != out_dir:
        f.write('rm -r {}\n'.format(temp_dir))
    f.close()

    return fn
