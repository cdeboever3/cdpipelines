import os

from general import _bedgraph_to_bigwig
from general import _bigwig_files
from general import _cutadapt_trim
from general import _fastqc
from general import _flagstat
from general import _picard_insert_size_metrics
from general import _make_softlink
from general import _pbs_header
from general import _picard_coord_sort
from general import _picard_query_sort
from general import _picard_mark_duplicates
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
                          '\t--genomeLoad LoadAndRemove', 
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

def _genome_browser_files(tracklines_file, link_dir, web_path_file,
                          coord_sorted_bam, bam_index, r1_fastqc,
                          r2_fastqc, narrow_peak, broad_peak, sample_name,
                          outdir):
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

    r1_fastqc : str
        Path to fastqc_report.html for R1 reads.

    r2_fastqc : str
        Path to fastqc_report.html for R2 reads.

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
    fn = os.path.join(outdir, os.path.split(coord_sorted_bam)[1])
    new_lines, bam_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines

    fn = os.path.join(outdir, os.path.split(bam_index)[1])
    new_lines, index_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines

    tf_lines += ' '.join(['track', 'type=bam',
                          'name="{}_bam"'.format(sample_name),
                          'description="ATAC-seq for {}"'.format(sample_name),
                          'visibility=0',
                          'db=hg19',
                          'bigDataUrl={}/{}\n'.format(temp_web_path, bam_name)])
   
    # HOMER bigwig.
    temp_link_dir = os.path.join(link_dir, 'bigwig')
    temp_web_path = web_path + '/bigwig'
    try:
        os.makedirs(temp_link_dir)
    except OSError:
        pass
    fn = os.path.join(outdir, '{}_tags'.format(sample_name),
                      '{}_tags.ucsc.bigWig'.format(sample_name))
    new_lines, bigwig_name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines
    tf_lines += ('track type=bigWig name="{0}_atac_cov" description="ATAC-seq '
                 'coverage for {0}" visibility=0 db=hg19 bigDataUrl='
                 '{1}/{0}_tags.ucsc.bigWig\n'.format(sample_name,
                                                     temp_web_path))

    # HOMER peaks.
    temp_link_dir = os.path.join(link_dir, 'peak')
    temp_web_path = web_path + '/peak'
    try:
        os.makedirs(temp_link_dir)
    except OSError:
        pass
    fn = os.path.join(outdir, '{}_tags'.format(sample_name),
                      '{}_atac_homer_peaks.bed'.format(sample_name))
    new_lines, name = _make_softlink(fn, sample_name, temp_link_dir)
    lines += new_lines
    tf_lines += '{}/{}\n'.format(temp_web_path, os.path.split(fn)[1])

    # Peaks from MACS2. Note that this file is just directly uploaded to UCSC so
    # we don't provide a trackline but rather just a URL to UCSC.
    # fn = os.path.join(outdir, os.path.split(narrow_peak)[1])
    # new_lines, name = _make_softlink(fn, sample_name, temp_link_dir)
    # lines += new_lines
    # fn = os.path.join(outdir, os.path.split(broad_peak)[1])
    # new_lines, name = _make_softlink(fn, sample_name, temp_link_dir)
    # lines += new_lines
    # tf_lines += '{}/{}\n'.format(temp_web_path, os.path.split(fn)[1])

    with open(tracklines_file, 'w') as tf:
        tf.write(tf_lines)
    
    lines += '\n'
    return lines

def _homer(bam, sample_name, temp_tagdir, final_tagdir, homer_path, link_dir,
           bedtools_path, bigwig=False):
    """
    Make tag directory and call peaks with HOMER. Optionally make bigwig file.

    Parameters
    ----------
    bam : str
        Path to paired-end, coordinate sorted bam file.

    temp_tagdir : str
        Path to temporary tag directory that HOMER will create.

    final_tagdir : str
        Path to final tag directory (what the temp_tagdir will be copied to).

    link_dir : str
        Path to directory where softlinks should be made. HOMER will put the
        bigwig file in the directory specified when setting up HOMER but a
        softlink to a bed file with the HOMER peaks will be made here.

    bigwig : bool
        If True, have HOMER make a bigwig file.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    temp_tagdir = os.path.realpath(temp_tagdir)
    final_tagdir = os.path.realpath(final_tagdir)
    name = '{}_atac_homer'.format(sample_name)
    lines = []
    lines.append('{}/makeTagDirectory {} {}'.format(homer_path, temp_tagdir,
                                                    bam))
    if bigwig:
        lines.append('{}/makeBigWig.pl {} hg19 -name {}'
                     ' -url www.fake.com/ -webdir {}'.format(
                         homer_path, os.path.split(temp_tagdir)[1], name, 
                         os.path.split(temp_tagdir)[0]))
        lines.append('mv {}/{}_tags.ucsc.bigWig {}'.format(
            os.path.split(temp_tagdir)[0], sample_name, temp_tagdir))
    lines.append('{}/findPeaks {} -style histone -size 75 -minDist 75 '
                 '-o auto'.format(homer_path, temp_tagdir))
    # softlink_lines = []
    posfile = os.path.join(temp_tagdir, 'regions.txt')
    bed = os.path.join(temp_tagdir, '{}_peaks.bed'.format(name))
    lines.append(_convert_homer_pos_to_bed(
        posfile, bed, name, homer_path, bedtools_path))
    bed = os.path.join(final_tagdir, '{}_peaks.bed'.format(name))
    # softlink, name = _make_softlink(bed, name, os.path.join(link_dir, 'atac',
    #                                                         'peak'))
    # softlink_lines.append(softlink)
    
    lines = '\n'.join(lines) + '\n\n'
    # softlink_lines = '\n'.join(softlink_lines)
    # return lines, softlink_lines
    return lines

def _combined_homer(input_tagdirs, combined_name, temp_tagdir, final_tagdir,
                    homer_path, link_dir, bedtools_path, tracklines_file,
                    web_path_file, bigwig=False):
    """
    Combine tag directories and call peaks with HOMER. Optionally make bigwig
    file.

    Parameters
    ----------
    input_tagdirs : list 
        A list of paths to the tagdirs to be combined for calling peaks.

    combined_name : str
        Used for naming output files and directories.

    temp_tagdir : str
        Path to temporary tag directory that HOMER will create.

    final_tagdir : str
        Path to final tag directory (what the temp_tagdir will be copied to).

    link_dir : str
        Path to directory where softlinks should be made. HOMER will put the
        bigwig file in the directory specified when setting up HOMER but a
        softlink to a bed file with the HOMER peaks will be made here.

    bigwig : bool
        If True, have HOMER make a bigwig file.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()
    web_path = web_path + '/atac'

    temp_tagdir = os.path.realpath(temp_tagdir)
    final_tagdir = os.path.realpath(final_tagdir)
    name = '{}_combined_peak'.format(combined_name)
    bed = os.path.join(temp_tagdir,
                       '{}_combined_homer_peaks.bed'.format(combined_name))
    lines = []
    tf_lines = []
    softlink_lines = []
    lines.append('{}/makeTagDirectory {} -d {}'.format(homer_path, temp_tagdir, 
                                                      ' '.join(input_tagdirs)))
    if bigwig:
        lines.append('{}/makeBigWig.pl {} hg19 -name {}'
                     ' -url www.fake.com/ -webdir {}'.format(
                         homer_path, os.path.split(temp_tagdir)[1], name, 
                         os.path.split(temp_tagdir)[0]))
        lines.append('mv {0}.ucsc.bigWig {0}'.format(temp_tagdir))
        temp_link_dir = os.path.join(link_dir, 'bigwig')
        bw = os.path.join(
            final_tagdir,
            '{}.ucsc.bigWig'.format(os.path.split(temp_tagdir)[1]))
        softlink, name = _make_softlink(bw, name, os.path.join(
            temp_link_dir, 'atac', 'bigwig'))
        softlink_lines.append(softlink)
        temp_web_path = web_path + '/bigwig'
        try:
            os.makedirs(temp_link_dir)
        except OSError:
            pass
        tf_lines += ('track type=bigWig name="{0}_atac_cov" '
                     'description="ATAC-seq '
                     'coverage for {0}" visibility=0 db=hg19 bigDataUrl='
                     '{1}/{0}.ucsc.bigWig\n'.format(
                         os.path.split(temp_tagdir)[1], temp_web_path))
    lines.append('{}/findPeaks {} -style super -size 75 '
                 ' -o auto'.format(
                     homer_path, temp_tagdir))
    lines.append('{}/findPeaks {} -style histone -size 75 -minDist 75 '
                 '-o auto'.format(homer_path, temp_tagdir))
   
    posfile = os.path.join(temp_tagdir, 'regions.txt')
    name = '{}_combined'.format(combined_name)
    bed = os.path.join(temp_tagdir, '{}_homer_atac_peaks.bed'.format(name))
    lines.append(_convert_homer_pos_to_bed(
        posfile, bed, name, homer_path, bedtools_path))
    bed = os.path.join(final_tagdir, '{}_homer_atac_peaks.bed'.format(name))
    softlink, name = _make_softlink(bed, name, os.path.join(link_dir, 'atac',
                                                            'peak'))
    softlink_lines.append(softlink)
    temp_web_path = web_path + '/peak'
    tf_lines += '{}/{}_homer_atac_peaks.bed\n'.format(temp_web_path, name)

    posfile = os.path.join(temp_tagdir, 'superEnhancers.txt')
    name = '{}_combined_super_enhancers'.format(combined_name)
    bed = os.path.join(temp_tagdir, '{}_homer_atac_peaks.bed'.format(name))
    lines.append(_convert_homer_pos_to_bed(
        posfile, bed, name, homer_path, bedtools_path))
    bed = os.path.join(final_tagdir, '{}_homer_atac_peaks.bed'.format(name))
    softlink, name = _make_softlink(bed, name, os.path.join(link_dir, 'atac',
                                                            'peak'))
    softlink_lines.append(softlink)
    tf_lines += '{}/{}_homer_atac_peaks.bed\n'.format(temp_web_path, name)
    lines = '\n'.join(lines) + '\n\n'
    softlink_lines = '\n'.join(softlink_lines)

    # Write tracklines and URLs.
    if os.path.exists(tracklines_file):
        with open(tracklines_file) as f:
            existing_lines = f.read()
    else:
        existing_lines = ''

    with open(tracklines_file, 'w') as tf:
        tf.write(existing_lines + ''.join(tf_lines))
    
    return lines, softlink_lines

def _convert_homer_pos_to_bed(posfile, bed, sample_name, homer_path,
                              bedtools_path):
    """
    Convert HOMER results file to bed file.

    Parameters
    ----------
    posfile : str
        Full path to HOMER results file (i.e. regions.txt or
        superEnhancers.txt).

    bed : str
        Name of output bed file.

    sample_name : str
        Used for naming output files.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = []
    tagdir = os.path.split(posfile)[0]
    lines.append('{}/pos2bed.pl {} | grep -v \# > temp.bed'.format(
        homer_path, posfile))
    track_line = ' '.join(['track', 'type=bed',
                           'name=\\"{}_homer_atac_peaks\\"'.format(
                               sample_name),
                           ('description=\\"HOMER ATAC-seq peaks for '
                            '{}\\"'.format(sample_name)),
                           'visibility=0',
                           'db=hg19'])
    lines.append('{} sort -i temp.bed > temp2.bed'.format(bedtools_path))
    lines.append('cat <(echo {}) temp2.bed > {}'.format(track_line, bed))
    lines.append('rm temp.bed temp2.bed')
    return '\n'.join(lines) + '\n'

def _add_macs2_trackline(bed, sample_name, outdir):
    """
    Add trackline to macs2 bed file. 

    Parameters
    ----------
    bed : str:
        Full path to macs2 bed file.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    path = os.path.split(bed)[0]
    track_type = os.path.splitext(bed)[1][1:]
    if track_type == 'bed':
        name = 'summits'
    if track_type == 'narrowPeak':
        name = 'narrowPeak'
    if track_type == 'broadPeak':
        name = 'broadPeak'
    if track_type == 'gappedPeak':
        name = 'gappedPeak'
    
    temp_bed = os.path.join(outdir, 'temp.{}'.format(name))
    track_line = (
        'track type={} name=\\"{}_macs2_atac_{}\\" '
        'description=\\"macs2 ATAC-seq {} for {}\\" '
        'visibility=0 db=hg19'.format(
            track_type, sample_name, name, name, sample_name)
    )
    lines = 'cat <(echo "{}") {} > {}\n'.format(track_line, bed, temp_bed)
    lines += 'mv {} {}\n\n'.format(temp_bed, bed)
    return lines

def _macs2(
    bam, 
    sample_name, 
    outdir,
    tracklines_file,
    link_dir,
    web_path_file,
    broad=False,
):
    """
    Call peaks with MACS2. The macs2 executable is assumed to be in your path
    which it should be if you installed it using pip install MACS2 into your
    python environment.

    Parameters
    ----------
    bam : str
        Path to paired-end, coordinate sorted bam file.

    outdir : str
        Path to directory where final macs2 output should be stored. Softlinks
        will be made to the output files in this directory.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    # Prepare some stuff for making softlinks and web links.
    link_dir = os.path.join(link_dir, 'atac', 'peak')
    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()
    web_path = web_path + '/atac/peak'
    if os.path.exists(tracklines_file):
        with open(tracklines_file) as f:
            tf_lines = f.read()
    else:
        tf_lines = ''
    try:
        os.makedirs(link_dir)
    except OSError:
        pass

    run_type = '--call-summits'
    if broad:
        run_type = '--broad'
    lines = ('macs2 callpeak --nomodel --nolambda --keep-dup all \\\n'
             '\t{} --slocal 10000 -f BAMPE -g hs \\\n'
             '\t-t {} \\\n'
             '\t-n {} \\\n'
             '\t--outdir {}\n\n'.format(
                 run_type, bam, sample_name, outdir))
    
    # Add tracklines to bed files and make softlinks.
    if not broad:
        for bed in [
            os.path.join(outdir, '{}_peaks.narrowPeak'.format(sample_name)),
            os.path.join(outdir, '{}_summits.bed'.format(sample_name))
        ]:
            lines += _add_macs2_trackline(bed, sample_name, outdir)
            new_lines, name = _make_softlink(bed, sample_name, link_dir)
            lines += new_lines
            tf_lines += '{}/{}\n'.format(web_path, name)
    elif broad:
        for bed in [
            os.path.join(outdir, '{}_peaks.broadPeak'.format(sample_name)),
            os.path.join(outdir, '{}_peaks.gappedPeak'.format(sample_name))
        ]:
            lines += _add_macs2_trackline(bed, sample_name, outdir)
            new_lines, name = _make_softlink(bed, sample_name, link_dir)
            lines += new_lines
            tf_lines += '{}/{}\n'.format(web_path, name)

    # Write tracklines and URLs.
    with open(tracklines_file, 'w') as tf:
        tf.write(tf_lines)

    return lines

def align_and_call_peaks(
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
    samtools_path,
    homer_path,
    blacklist_bed,
    conda_env='',
    rgpl='ILLUMINA',
    rgpu='',
    tempdir='/scratch', 
    threads=10, 
    picard_memory=18, 
    shell=False,
    trim=None,
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

    fastqc_path : str
        Path to FastQC executable.

    samtools_path : str
        Path to samtools executable.

    homer_path : str
        Path to HOMER bin.

    blacklist_bed : str
        String to bed file with blacklist regions. This should at least include
        the ENCODE blacklist regions. Read pairs where one or both reads
        intersect these regions are filtered out.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of MACS2 is used.

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

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    trim : int
        Positive or negative integer. Positive numbers remove bases at the front
        of the read and negative numbers remove bases at the end of the read.

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

    tempdir = os.path.join(tempdir, '{}_peaks'.format(sample_name))
    outdir = os.path.join(outdir, '{}_peaks'.format(sample_name))

    # I'm going to define some file names used later.
    temp_r1_fastqs = _process_fastqs(r1_fastqs, tempdir)
    temp_r2_fastqs = _process_fastqs(r2_fastqs, tempdir)
    combined_r1 = os.path.join(tempdir, 
                               '{}_combined_R1.fastq.gz'.format(sample_name))
    combined_r2 = os.path.join(tempdir, 
                               '{}_combined_R2.fastq.gz'.format(sample_name))
    aligned_bam = os.path.join(tempdir, 'Aligned.out.bam')
    filtered_bam = os.path.join(tempdir, 
                                '{}_atac_filtered.bam'.format(sample_name))
    coord_sorted_bam = os.path.join(
        tempdir, '{}_atac_filtered_coord_sorted.bam'.format(sample_name))
    no_dup_bam = os.path.join(tempdir,
                              '{}_atac_no_dup.bam'.format(sample_name))
    bam_index = os.path.join(tempdir,
                             '{}_atac_no_dup.bam.bai'.format(sample_name))
    qsorted_bam = os.path.join(
        tempdir, '{}_atac_no_dup_qsorted.bam'.format(sample_name))
    # out_bigwig = os.path.join(tempdir, '{}_atac.bw'.format(sample_name))
    insert_metrics = os.path.join(
        outdir, '{}_insert_metrics.txt'.format(sample_name))
    insert_hist = os.path.join(
        outdir, '{}_insert_histogram.pdf'.format(sample_name))
    duplicate_metrics = os.path.join(
        outdir, '{}_duplicate_metrics.txt'.format(sample_name))
    stats_file = os.path.join(outdir,
                              '{}_atac_no_dup.bam.flagstat'.format(sample_name))
    chrM_counts = os.path.join(outdir,
                               '{}_chrM_counts.txt'.format(sample_name))
    narrow_peak = os.path.join(outdir,
                               '{}_peaks.narrowPeak'.format(sample_name))
    broad_peak = os.path.join(outdir,
                               '{}_peaks.broadPeak'.format(sample_name))
    local_tagdir = '{}_tags'.format(sample_name)
    temp_tagdir = os.path.join(tempdir, local_tagdir)
    final_tagdir = os.path.join(outdir, local_tagdir)
    
    tn = os.path.split(combined_r1)[1]
    r1_fastqc = os.path.join(
        outdir, 
        '{}_fastqc'.format('.'.join(tn.split('.')[0:-2])),
        'fastqc_report.html')
    tn = os.path.split(combined_r2)[1]
    r2_fastqc = os.path.join(
        outdir, 
        '{}_fastqc'.format('.'.join(tn.split('.')[0:-2])),
        'fastqc_report.html')
    
    # Files to copy to output directory.
    files_to_copy = [aligned_bam, no_dup_bam, bam_index, qsorted_bam, 'Log.out',
                     'Log.final.out', 'Log.progress.out', 'SJ.out.tab',
                     temp_tagdir]
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [combined_r1, combined_r2, aligned_bam, filtered_bam,
                       coord_sorted_bam, '_STARtmp']

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_peaks.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_peaks.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir, '{}_peaks.out'.format(sample_name))
        err = os.path.join(outdir, '{}_peaks.err'.format(sample_name))
        job_name = '{}_peaks'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))
    
    if conda_env != '':
        f.write('source activate {}\n'.format(conda_env))
    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))
    f.write('rsync -avz \\\n{} \\\n{} \\\n\t.\n\n'.format(
        ' \\\n'.join(['\t{}'.format(x) for x in r1_fastqs]),
        ' \\\n'.join(['\t{}'.format(x) for x in r2_fastqs])))

    # Add HOMER executables to path because HOMER expects them there.
    f.write('export PATH="{}:$PATH"\n\n'.format(homer_path))
    
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
    lines = _fastqc([combined_r1, combined_r2], threads, outdir, fastqc_path)
    f.write(lines)
    f.write('wait\n\n')

    # Optionally trim.
    if trim:
        assert(type(trim) is int)
        trimmed_r1 = os.path.join(
            tempdir, '{}_trimmed_R1.fastq.gz'.format(sample_name))
        trimmed_r2 = os.path.join(
            tempdir, '{}_trimmed_R2.fastq.gz'.format(sample_name))
        files_to_remove.append(trimmed_r1)
        files_to_remove.append(trimmed_r2)
        lines = _cutadapt_trim(combined_r1, trim, trimmed_r1, bg=True)
        lines += _cutadapt_trim(combined_r2, trim, trimmed_r2, bg=True)
        to_align_r1 = trimmed_r1
        to_align_r2 = trimmed_r2
        f.write(lines)
        f.write('\nwait\n\n')
        lines = _fastqc([trimmed_r1, trimmed_r2], threads, outdir, fastqc_path)
        f.write(lines)
        f.write('wait\n\n')
    else:
        to_align_r1 = combined_r1
        to_align_r2 = combined_r2

    # Align with STAR.
    lines = _star_align(to_align_r1, to_align_r2, sample_name, rgpl,
                        rgpu, star_index, star_path, threads)
    f.write(lines)
    f.write('wait\n\n')

    # Count the number of primary alignments for each chromosome.
    f.write('{} view -q 255 {} | \\\n\tcut -f 3 | \\\n\tgrep chrM '
            '| \\\n\tuniq -c > {} &\n\n'.format(
        samtools_path, aligned_bam, chrM_counts))

    # Remove mitochondrial reads, read pairs that are not uniquely aligned, and
    # reads where one or both of the reads were in the ENCODE blacklist regions.
    lines = (
        '{} view -h -q 255 {} | \\\n'.format(samtools_path, aligned_bam) + 
        '\tawk \'{if ($3 != "chrM") {print} ' + 
        'else if (substr($1,1,1) == "@") {print}}\' | \\\n' + 
        '\t{} view -Su - | \\\n'.format(samtools_path) + 
        '\t{} intersect -v -ubam -abam stdin -b {} | \\\n'.format(
            bedtools_path, blacklist_bed) + 
        '\t{} view -h - | \\\n'.format(samtools_path) + 
        '\tawk \'{if (substr($1,1,1) == "@") {print} ' + 
        'else if ($1 == c1) {print prev; print}  prev=$0; c1=$1}\' | \\\n' + 
        '\t{} view -Sb - > {}\n\n'.format(samtools_path, filtered_bam)
    )
    f.write(lines)

    # Coordinate sort bam file.
    lines = _picard_coord_sort(filtered_bam, coord_sorted_bam, picard_path,
                               picard_memory, tempdir)
    f.write(lines)
    f.write('wait\n\n')

    # Remove duplicates.
    lines = _picard_mark_duplicates(coord_sorted_bam, no_dup_bam,
                                    duplicate_metrics, picard_path,
                                    picard_memory, tempdir, remove_dups=True)
    f.write(lines)
    f.write('wait\n\n')

    # Query sort.
    lines = _picard_query_sort(no_dup_bam, qsorted_bam, picard_path,
                               picard_memory, tempdir, bg=True)
    f.write(lines)

    # Index bam file and collect flagstats.
    lines = _samtools_index(no_dup_bam, samtools_path, bg=True)
    f.write(lines)
    lines = _flagstat(no_dup_bam, stats_file, samtools_path, bg=True)
    f.write(lines)
    f.write('wait\n\n')
    lines = _picard_insert_size_metrics(no_dup_bam, insert_metrics, insert_hist,
                                        picard_path, picard_memory, tempdir,
                                        bg=True)
    f.write(lines)

    # Call peaks with macs2.
    lines = _macs2(no_dup_bam, sample_name, outdir, tracklines_file, link_dir,
                   web_path_file, broad=True)
    f.write(lines)
    f.write('wait\n\n')

    # Run HOMER.
    lines = _homer(qsorted_bam, sample_name, temp_tagdir, final_tagdir,
                   homer_path, link_dir, bedtools_path, bigwig=True)
    f.write(lines)
    f.write('wait\n\n')

    if tempdir != outdir:
        f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
            ' \\\n\t'.join([x for x in files_to_copy if sample_name in 
                            os.path.split(x)[1]]),
            outdir))
        for y in [x for x in files_to_copy if sample_name not in 
             os.path.split(x)[1]]:
            f.write('rsync -avz {} {}_{}\n'.format(
                y, os.path.join(outdir, sample_name), os.path.split(y)[1]))
            f.write('rm {}\n'.format(y))
    else:
        for y in [x for x in files_to_copy if sample_name not in 
             os.path.split(x)[1]]:
            f.write('mv {} {}_{}\n'.format(
                y, os.path.join(outdir, sample_name), os.path.split(y)[1]))

    f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))
    
    # Make softlinks and tracklines for genome browser.
    lines = _genome_browser_files(tracklines_file, link_dir, web_path_file,
                                  no_dup_bam, bam_index, r1_fastqc,
                                  r2_fastqc, narrow_peak, broad_peak,
                                  sample_name, outdir)
    f.write(lines)

    f.close()
    return fn

def combined_homer_peaks(
    tagdirs, 
    outdir, 
    combined_name, 
    tracklines_file,
    link_dir,
    web_path_file,
    bedtools_path,
    homer_path,
    environment,
    conda_env='',
    tempdir='/scratch', 
    threads=6, 
    shell=False,
):
    """
    Make a PBS or shell script for combining together HOMER tag directories and
    calling peaks on the combined tags.

    Parameters
    ----------
    tagdirs : list 
        A list of paths to the tagdirs to be combined for calling peaks.

    outdir : str
        Directory to store PBS/shell file and results.

    combined_name : str
        Name used for naming files and directories for combined data.

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

    bedtools_path : str
        Path to bedtools.

    homer_path : str
        Path to HOMER bin.

    environment : str
        Bash file with PATH information that can be sourced. This should include
        the paths to executables HOMER will need like bedGraphToBigWig.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of MACS2 is used.

    tempdir : str
        Directory to store files as STAR runs.

    threads : int
        Number of threads to reserve using PBS scheduler. 

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

    tempdir = os.path.join(tempdir, '{}_combined_peaks'.format(combined_name))
    outdir = os.path.join(outdir, '{}_combined_peaks'.format(combined_name))

    # I'm going to define some file names used later.
    local_tagdir = '{}_combined_tags'.format(combined_name)
    temp_tagdir = os.path.join(tempdir, local_tagdir)
    final_tagdir = os.path.join(outdir, local_tagdir)
    
    # Files to copy to output directory.
    files_to_copy = [temp_tagdir]
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [os.path.split(os.path.realpath(x))[1] for x in tagdirs]

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_combined_peaks.sh'.format(combined_name))
    else:
        fn = os.path.join(outdir, '{}_combined_peaks.pbs'.format(combined_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir,
                           '{}_combined_peaks.out'.format(combined_name))
        err = os.path.join(outdir,
                           '{}_combined_peaks.err'.format(combined_name))
        job_name = '{}_combined_peaks'.format(combined_name)
        f.write(_pbs_header(out, err, job_name, threads))
    
    if conda_env != '':
        f.write('source activate {}\n'.format(conda_env))
    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))
    f.write('rsync -avz \\\n{} \\\n\t.\n\n'.format(
        ' \\\n'.join(['\t{}'.format(x) for x in tagdirs])))

    # Add executables to path because HOMER expects them there.
    f.write('source {}\n\n'.format(environment))
    
    # Run HOMER.
    td = [os.path.split(os.path.realpath(x))[1] for x in tagdirs]
    lines, softlink_lines = _combined_homer(td, combined_name, temp_tagdir,
                                            final_tagdir, homer_path, link_dir,
                                            bedtools_path, tracklines_file,
                                            web_path_file, bigwig=True)
    f.write(lines)
    f.write('wait\n\n')

    if tempdir != outdir:
        f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
            ' \\\n\t'.join(files_to_copy), outdir))

    if len(files_to_remove) > 0:
        f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))

    f.write(softlink_lines)

    f.close()
    return fn

def _nucleoatac(bam, bed, sample_name, fasta, threads):
    lines = ('nucleoatac run --bed {} \\\n'
             '\t--bam {} \\\n'
             '\t--out {} \\\n'
             '\t--fasta {} \\\n'
             '\t--cores {}\n\n'.format(bed, bam, sample_name, fasta, threads))
    return lines

def nucleoatac(
    bam, 
    bed, 
    sample_name, 
    fasta, 
    outdir,
    tracklines_file,
    link_dir,
    web_path_file,
    environment,
    conda_env='',
    tempdir='/scratch', 
    threads=4, 
    shell=False,
):
    """
    Make a PBS or shell script for estimating nucelosome occupancy using
    ATAC-seq data.

    Parameters
    ----------
    bam : str 
        Path to bam file with aligned reads to use for estimating occupancy.

    bed : str
        Path to bed file with positions to estimate occupancy for.

    sample_name : str
        Name used for naming files and directories.

    fasta : str
        Path to genome fasta. Must be indexed.

    outdir : str
        Directory to store directory containing PBS/shell file and results.

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

    environment : str
        Bash file with PATH information that can be sourced. This should include
        the paths to executables HOMER will need like bedGraphToBigWig.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of nucleoatac is used.

    tempdir : str
        Directory to store files as nucleoatac runs.

    threads : int
        Number of threads to reserve using PBS scheduler and for nucleoatac to
        use.

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    if shell:
        pbs = False
    else: 
        pbs = True

    tempdir = os.path.join(tempdir, '{}_nucleoatac'.format(sample_name))
    outdir = os.path.join(outdir, '{}_nucleoatac'.format(sample_name))

    # I'm going to define some file names used later.
    temp_bam = os.path.join(tempdir, os.path.split(bam)[1])
    
    # Files to copy to output directory.
    files_to_copy = []
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [temp_bam]

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_nucleoatac.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_nucleoatac.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir,
                           '{}_nucleoatac.out'.format(sample_name))
        err = os.path.join(outdir,
                           '{}_nucleoatac.err'.format(sample_name))
        job_name = '{}_nucleoatac'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))
    
    if conda_env != '':
        f.write('source activate {}\n'.format(conda_env))
    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))

    f.write('source {}\n\n'.format(environment))
    f.write('rsync -avz \\\n\t{} \\\n\t.\n\n'.format(bam))
    
    # Run nucleoatac.
    lines = _nucleoatac(temp_bam, bed, sample_name, fasta, threads)
    f.write(lines)
    f.write('wait\n\n')

    if tempdir != outdir:
        if len(files_to_copy) > 0:
            f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
                ' \\\n\t'.join(files_to_copy), outdir))

    if len(files_to_remove) > 0:
        f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
            f.write('rsync -avz {} {}\n\n'.format(os.path.join(tempdir, '*'),
                                                  outdir))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))

    f.close()
    return fn

def macs2_peak_calling(
    bam, 
    sample_name, 
    outdir,
    tracklines_file,
    link_dir,
    web_path_file,
    environment,
    conda_env='',
    tempdir='/scratch', 
    threads=4, 
    shell=False,
):
    """
    Make a PBS or shell script for calling broad and narrow peaks with macs2.

    Parameters
    ----------
    bam : str 
        Path to bam file with aligned reads to use for calling peaks.

    sample_name : str
        Name used for naming files and directories.

    outdir : str
        Directory to store directory containing PBS/shell file and results.

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

    environment : str
        Bash file with PATH information that can be sourced. This should include
        the paths to executables HOMER will need like bedGraphToBigWig.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of nucleoatac is used.

    tempdir : str
        Directory to store files as nucleoatac runs.

    threads : int
        Number of threads to reserve using PBS scheduler and for nucleoatac to
        use.

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    if shell:
        pbs = False
    else: 
        pbs = True

    tempdir = os.path.join(tempdir, '{}_macs2'.format(sample_name))
    outdir = os.path.join(outdir, '{}_macs2'.format(sample_name))

    # I'm going to define some file names used later.
    temp_bam = os.path.join(tempdir, os.path.split(bam)[1])
    
    # Files to copy to output directory.
    files_to_copy = []
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [temp_bam]

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_macs2.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_macs2.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir,
                           '{}_macs2.out'.format(sample_name))
        err = os.path.join(outdir,
                           '{}_macs2.err'.format(sample_name))
        job_name = '{}_macs2'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))
    
    if conda_env != '':
        f.write('source activate {}\n'.format(conda_env))
    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))

    f.write('source {}\n\n'.format(environment))
    f.write('rsync -avz \\\n\t{} \\\n\t.\n\n'.format(bam))

    # Run macs2 for narrow peaks.
    lines = _macs2(
        temp_bam, 
        sample_name, 
        outdir,
        tracklines_file,
        link_dir,
        web_path_file,
        broad=False,
    )
    f.write(lines)
    f.write('wait\n\n')

    # Run macs2 for broad peaks.
    lines = _macs2(
        temp_bam, 
        sample_name, 
        outdir,
        tracklines_file,
        link_dir,
        web_path_file,
        broad=True,
    )
    f.write(lines)
    f.write('wait\n\n')

    if tempdir != outdir:
        if len(files_to_copy) > 0:
            f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
                ' \\\n\t'.join(files_to_copy), outdir))

    if len(files_to_remove) > 0:
        f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
            f.write('rsync -avz {} {}\n\n'.format(os.path.join(tempdir, '*'),
                                                  outdir))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))

    f.close()
    return fn

def motif_analysis(
    bed, 
    sample_name, 
    outdir,
    tracklines_file,
    link_dir,
    web_path_file,
    environment,
    mask=False,
    conda_env='',
    tempdir='/scratch', 
    threads=4, 
    shell=False,
):
    """
    Make a PBS or shell script for analyzing motifs with HOMER.

    Parameters
    ----------
    bed : str 
        Bed file with positions to analyze.

    sample_name : str
        Name used for naming files and directories.

    outdir : str
        Directory to store directory containing PBS/shell file and results.

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

    environment : str
        Bash file with PATH information that can be sourced. This should include
        the paths to executables HOMER will need like bedGraphToBigWig.

    mask : bool
        Whether to pass the -mask parameter to HOMER.

    conda_env : str
        If provided, load conda environment with this name. This will control
        which version of nucleoatac is used.

    tempdir : str
        Directory to store files as nucleoatac runs.

    threads : int
        Number of threads to reserve using PBS scheduler and for HOMER to use.

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    if shell:
        pbs = False
    else: 
        pbs = True

    tempdir = os.path.join(tempdir, '{}_motif'.format(sample_name))
    outdir = os.path.join(outdir, '{}_motif'.format(sample_name))

    # I'm going to define some file names used later.
    
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
        fn = os.path.join(outdir, '{}_motif.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_motif.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir,
                           '{}_motif.out'.format(sample_name))
        err = os.path.join(outdir,
                           '{}_motif.err'.format(sample_name))
        job_name = '{}_motif'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))
    
    if conda_env != '':
        f.write('source activate {}\n'.format(conda_env))
    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))

    f.write('source {}\n\n'.format(environment))

    # Prepare some stuff for making softlinks and web links.
    link_dir = os.path.join(link_dir, 'atac', 'motif')
    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()
    web_path = web_path + '/atac/motif'
    if os.path.exists(tracklines_file):
        with open(tracklines_file) as tf:
            tf_lines = tf.read()
    else:
        tf_lines = ''
    try:
        os.makedirs(link_dir)
    except OSError:
        pass

    # Run HOMER motif analysis.
    if mask:
        lines = ('findMotifsGenome.pl {} hg19 {} -size given -mask '
                 '-p {}\n'.format(bed, outdir, threads))
    else:
        lines = 'findMotifsGenome.pl {} hg19 {} -size given -p {}\n'.format(
            bed, outdir, threads)
    f.write(lines)
    f.write('wait\n\n')

    new_lines, name = _make_softlink(outdir, sample_name, link_dir)
    f.write(new_lines)
    tf_lines += '{}/{}\n'.format(web_path, os.path.split(outdir)[1])

    # Write tracklines and URLs.
    with open(tracklines_file, 'w') as tf:
        tf.write(tf_lines)
    if tempdir != outdir:
        if len(files_to_copy) > 0:
            f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
                ' \\\n\t'.join(files_to_copy), outdir))

    if len(files_to_remove) > 0:
        f.write('rm -r \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
            f.write('rsync -avz {} {}\n\n'.format(os.path.join(tempdir, '*'),
                                                  outdir))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))

    f.close()
    return fn
