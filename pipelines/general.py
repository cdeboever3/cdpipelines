import os
import subprocess

def _git_info():
    """Get current git version"""
    # Necessary to pipe to cat because git will open the result in less
    # otherwise.
    d = os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2] +
                         ['.git'])
    command = ('git --git-dir {0} log -1 --pretty=oneline --decorate | '
               'cat'.format(d))
    res = subprocess.Popen(command, shell=True,
                           stdout=subprocess.PIPE).communicate()
    res = res[0][:-1].strip() 
    return res 

def _make_dir(d):
    """Make directory d if it doesn't exist"""
    try:
        os.makedirs(d)
    except OSError:
        pass

def _picard_collect_rna_seq_metrics(in_bam, metrics, chart, sample_name,
                                    picard_path, picard_memory, tempdir,
                                    ref_flat, rrna_intervals,
                                    strand_specific=True, bg=False):
    """
    Collect RNA-seq metrics using Picard. The input bam file is assumed to be
    sorted.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    metrics : str
        Path to output metrics file.

    chart : str
        Path to output PDF file.

    sample_name : str
        Sample name used to name output files.

    ref_flat : str
        Path to refFlat file with non-rRNA genes. Can ge gzipped.

    rrna_intervals : str
        Pato to interval list file with rRNA intervals.

    bg : boolean
        Whether to run the process in the background.

    """
    if strand_specific:
        ss = '\tSTRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND'
    else:
        ss = '\tSTRAND_SPECIFICITY=NONE'
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} CollectRnaSeqMetrics'.format(
                               picard_path),
                           '\tI={}'.format(in_bam),
                           '\tREF_FLAT={}'.format(ref_flat),
                           ss,
                           '\tRIBOSOMAL_INTERVALS={}'.format(rrna_intervals),
                           '\tASSUME_SORTED=TRUE',
                           '\tCHART_OUTPUT={}'.format(chart),
                           '\tO={}'.format(metrics)]))

    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
    return lines

def _picard_collect_multiple_metrics(in_bam, sample_name, picard_path,
                                     picard_memory, tempdir, bg=False):
    """
    Collect multiple metrics using Picard. The input bam file is assumed to be
    sorted.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    sample_name : str
        Sample name used to name output files.

    bg : boolean
        Whether to run the process in the background.

    """
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} CollectMultipleMetrics'.format(
                               picard_path),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tASSUME_SORTED=TRUE',
                           '\tI={}'.format(in_bam), 
                           '\tO={}'.format(sample_name)]))
    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
    return lines

def _picard_gc_bias_metrics(in_bam, metrics, chart, out, picard_path,
                            picard_memory, tempdir, bg=False):
    """
    Collect GC bias metrics using Picard. The input bam file is assumed to
    be sorted.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    metrics : str
        Path to output metrics file.

    chart : str
        Path to output PDF chart.

    out : str
        Path to picard output file.

    bg : boolean
        Whether to run the process in the background.

    """
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} CollectGcBiasMetrics'.format(
                               picard_path),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tI={}'.format(in_bam), 
                           '\tO={}'.format(out),
                           '\tCHART_OUTPUT={}'.format(chart),
                           '\tSUMMARY_OUTPUT={}'.format(metrics),
                           '\tASSUME_SORTED=TRUE']))
    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
    return lines

def _picard_bam_index_stats(in_bam, out, err, picard_path, picard_memory,
                            tempdir, bg=False):
    """
    Collect bam index stats with Picard.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out : str
        Path to write stdout to. This contains the index stats.

    err : str
        Path to write stderr to. this contains picard stuff.

    bg : boolean
        Whether to run the process in the background.

    """
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} BamIndexStats'.format(
                               picard_path),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tI={}'.format(in_bam),
                           '\t> {}'.format(out),
                           '\t2> {}'.format(err)]))
    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
    return lines

def _picard_insert_size_metrics(in_bam, out_metrics, out_hist, picard_path,
                                picard_memory, tempdir, bg=False):
    """
    Collect insert size metrics using Picard. The input bam file is assumed to
    be sorted.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out_metrics : str
        Path to output metrics file.

    out_hist : str
        Path to output histogram PDF.

    bg : boolean
        Whether to run the process in the background.

    """
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} CollectInsertSizeMetrics'.format(
                               picard_path),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tI={}'.format(in_bam), 
                           '\tO={}'.format(out_metrics),
                           '\tHISTOGRAM_FILE={}'.format(out_hist),
                           '\tASSUME_SORTED=TRUE']))
    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
    return lines

def _picard_query_sort(in_bam, out_bam, picard_path, picard_memory, tempdir,
                       bg=False):
    """
    Query sort using Picard Tools.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out_bam : str
        Path to output bam file.

    bg : boolean
        Whether to run the process in the background.

    """
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} SortSam'.format(picard_path),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tI={}'.format(in_bam), 
                           '\tO={}'.format(out_bam),
                           '\tSO=queryname']))
    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
    return lines

def _picard_coord_sort(in_bam, out_bam, picard_path, picard_memory,
                       tempdir, bam_index=None):
    """
    Coordinate sort using Picard Tools.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out_bam : str
        Path to output bam file.

    bam_index : str
        If provided, generate index file for input bam file.

    """
    if bam_index:
        lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                               '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                               '\t-Djava.io.tmpdir={}'.format(tempdir), 
                               '\t-jar {} SortSam'.format(picard_path),
                               '\tVALIDATION_STRINGENCY=SILENT',
                               '\tCREATE_INDEX=TRUE', 
                               '\tI={}'.format(in_bam), 
                               '\tO={}'.format(out_bam),
                               '\tSO=coordinate\n']))
        index = '.'.join(out_bam.split('.')[0:-1]) + '.bai'
        lines += 'mv {} {}\n\n'.format(index, bam_index)
    else:
        lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                               '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                               '\t-Djava.io.tmpdir={}'.format(tempdir), 
                               '\t-jar {} SortSam'.format(picard_path),
                               '\tVALIDATION_STRINGENCY=SILENT',
                               '\tI={}'.format(in_bam), 
                               '\tO={}'.format(out_bam),
                               '\tSO=coordinate\n']))

    return lines

def _cutadapt_trim(fastq, length, out, bg=False):
    """
    Cut a specified number of bases from a fastq file using cutadapt. Cutadapt
    should be installed in your python environment.

    Parameters
    ----------
    fastq : str
        Fastq or gzipped/bzipped fastq.

    length : int
        Positive or negative integer. Positive numbers remove bases at the front
        of the read and negative numbers remove bases at the end of the read.

    out : str
        Path to output (optionally gzipped/bzipped) fastq files.

    bg : boolean
        Whether to run the process in the background (i.e. include an ampersand
        at the end of the command).

    Returns
    -------
    lines : str
        Lines to be written to PBS/shell script.

    """
    line = 'cutadapt --cut {} -o {} {}'.format(length, out, fastq)
    if bg:
        line += ' &\n\n'
    else:
        line += '\n\n'
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

def _flagstat(bam, stats_file, samtools_path, bg=False):
    """
    Run flagstat for a bam file.

    Parameters
    ----------
    bam : str
        Bam file to calculate coverage for.

    stats_file : str
        File to write flagstats to.

    samtools_path : str
        Path to samtools executable.

    Returns
    -------
    lines : str
        Lines to be written to PBS/shell script.

    """
    lines = '{} flagstat {} > {}'.format(samtools_path, bam, stats_file)
    if bg:
        lines += ' &\n\n'
    else:
        lines += '\n\n'
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
        lines += ('wait\n\n')
        lines += (_bedgraph_to_bigwig('both.bg', out_bigwig,
                                      bedgraph_to_bigwig_path, bedtools_path))
        lines += ('wait\n\n')
        lines += ('rm both.bg\n\n')
    return lines
    
def _process_fastqs(fastqs, tempdir):
    """
    Create list of temporary fastq paths.

    Parameters
    ----------
    fastqs : list or str
        Either a list of paths to gzipped fastq files or path to a single
        gzipped fastq file.

    tempdir : str
        Path to temporary directory where fastq files will be copied to.

    Returns
    -------
    temp_fastqs : list
        List of paths to temporary fastq files.

    """
    if type(fastqs) == list:
        fns = [os.path.split(x)[1] for x in fastqs]
        temp_fastqs = sorted([os.path.join(tempdir, x) for x in fns])
    elif type(fastqs) == str:
        temp_fastqs = [os.path.join(tempdir, os.path.split(fastqs)[1])]
    return temp_fastqs

def _fastqc(fastqs, threads, outdir, fastqc_path):
    """
    Run FastQC

Parameters
    ----------
    fastqs : str or list
        Path to fastq file or list of paths to fastq files.

    threads : int
        Number of threads to run FastQC with.

    outdir : str
        Path to directory to store FastQC results to.

    fastqc_path : str
        Path to FastQC.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    if type(fastqs) == list:
        fastqs = ' '.join(fastqs)
    lines = ('{} --outdir {} --nogroup \\\n'.format(fastqc_path, outdir) + 
             '\t--extract --threads {} {}\n\n'.format(threads, fastqs))
    return lines

def _softlink(target, link):
    """
    Make softlink from target to link.

    Parameters
    ----------
    target : str
        Full path to file to make link to.

    link : str
        Full path to softlink.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = 'ln -s {} {}\n'.format(target, link)
    return lines

def _make_softlink(fn, sample_name, link_dir):
    """
    Make softlink for file fn in link_dir. sample_name followed by an underscore
    will be appended to the front of fn if the sample_name isn't in fn.

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
    if sample_name not in os.path.split(fn)[1]:
        name = '{}_{}'.format(sample_name, os.path.split(fn)[1])
    else:
        name = os.path.split(fn)[1]
    lines = 'ln -s {} {}\n'.format(fn, os.path.join(link_dir, name))
    return lines, name

class JobScript:
    def __init__(self, sample_name, job_suffix, outdir, threads, tempdir=None,
                 shell=False, queue='high', conda_env=None, environment=None,
                 copy_input=True):
        """
        Create PBS/shell script object.

        Parameters
        ----------
        sample_name : str
            Sample name used for naming directories, files, etc.

        job_suffix : str
            This suffix will be used for naming directories, files, etc.

        outdir: str
            Path to directory where final output files should be stored.

        threads : int
            Number of threads to request for PBS scripts and to use for
            multi-threaded software.

        tempdir : str
            Path to directory where temporary directory should be made.

        shell : bool
            True if making a shell script. False if making a PBS script.

        queue : str
            PBS queue to use if writing PBS script.

        conda_env : str
            Path to conda environment to load when job begins.

        environment : str
            Path to file that should be sourced at beginning of script.

        copy_input : bool
            Whether to copy input files to temp directory. 

        """
        self.sample_name = sample_name
        self.job_suffix = job_suffix
        self.jobname = '{}_{}'.format(sample_name, job_suffix)
        if tempdir:
            self.tempdir = os.path.realpath(os.path.join(tempdir, self.jobname))
        else:
            self.tempdir = tempdir
        self.outdir = os.path.realpath(os.path.join(outdir, self.jobname))
        _make_dir(self.outdir)
        assert type(threads) is int
        self.threads = threads
        self.shell = shell
        self.pbs = not shell
        self.queue = queue
        self.conda_env = conda_env
        self.environment = environment
       
        self.out = os.path.join(self.outdir, '{}.out'.format(self.jobname))
        self.err = os.path.join(self.outdir, '{}.err'.format(self.jobname))

        # # self.files contains JobScriptFile objects.
        # self.files = []
        self.input_files_to_copy = []
        self.output_files_to_copy = []
        self.temp_files_to_delete = []
        self.softlinks = []
        self._set_filename()
        self._write_header()

    def _set_filename(self):
        """Make PBS/shell script filename."""
        if self.pbs:
            self.filename = os.path.join(self.outdir, 
                                         '{}.pbs'.format(self.jobname))
        else:
            self.filename = os.path.join(self.outdir, 
                                         '{}.sh'.format(self.jobname))

    def _write_header(self):
        with open(self.filename, "a") as f:
            f.write('#!/bin/bash\n\n')
            if self.pbs:
                f.write('#PBS -q {}\n'.format(self.queue))
                f.write('#PBS -N {}\n'.format(self.jobname))
                f.write('#PBS -l nodes=1:ppn={}\n'.format(self.threads))
                f.write('#PBS -o {}\n'.format(self.out))
                f.write('#PBS -e {}\n\n'.format(self.err))
            f.write('# Git repository version:\n# {}\n\n'.format(_git_info()))
            if self.environment:
                f.write('source {}\n\n'.format(self.environment))
            if self.conda_env:
                f.write('source activate {}\n\n'.format(self.conda_env))
            if self.tempdir:
                f.write('mkdir -p {}\n'.format(self.tempdir))
                f.write('cd {}\n\n'.format(self.tempdir))

    # def add_file(self, fn, copy_to_tempdir=True, delete_temp=True,
    #              copy_to_outdir=False):
    #     f = JobScriptFile(fn, copy_to_tempdir=copy_to_tempdir,
    #                       delete_temp=delete_temp,
    #                       copy_to_outdir=copy_to_outdir, tempdir=self.tempdir)
    #     self.files.append(f)
    #     if copy_to_tempdir:
    #         self.input_files_to_copy.append(f.filename)
    #     if delete_temp:
    #         self.files_to_delete.append(f.temp_filename)
    #     if copy_to_outdir:
    #         self.output_files_to_copy.append(f.temp_filename)
    #     return f

    def add_softlink(self, target, link):
        """Add a target, link pair to the JobScript instance so that a softlink
        from target to link is made at the end of the jobscript."""
        self.softlinks.append([target, link])

    def add_input_file(self, fn):
        """Add input file to self.input_files_to_copy and return temp path"""
        self.input_files_to_copy.append(fn)
        return self.temp_file_path(fn)

    def temp_file_path(self, fn):
        """Return the path to the temporary version of an input file"""
        return os.path.join(self.tempdir, os.path.split(fn)[1])

    def copy_input_files(self):
        if len(self.input_files_to_copy) > 0:
            with open(self.filename, "a") as f:
                f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
                    '\\\n\t'.join(self.input_files_to_copy),
                    self.tempdir))
            self.temp_files_to_delete += [
                os.path.join(self.tempdir, os.path.split(x)[1]) for x in
                self.input_files_to_copy]

    def _copy_output_files(self):
        if len(self.output_files_to_copy) > 0 and self.tempdir != self.outdir:
            with open(self.filename, "a") as f:
                f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
                    '\\\n\t'.join(self.output_files_to_copy),
                    self.outdir))

    def _delete_temp_files(self):
        if len(self.temp_files_to_delete) > 0:
            if self.tempdir == self.outdir or self.tempdir is None:
                self.temp_files_to_delete = [
                    x for x in self.temp_files_to_delete if x not in
                    self.output_files_to_copy
                ]
            with open(self.filename, "a") as f:
                f.write('rm -r \\\n\t{}\n\n'.format(
                    ' \\\n\t'.join(self.temp_files_to_delete)))

    def _delete_tempdir(self):
        if self.tempdir and (self.tempdir != self.outdir):
            with open(self.filename, "a") as f:
                f.write('rm -r {}\n'.format(self.tempdir))

    def _make_softlinks(self):
        with open(self.filename, "a") as f:
            for p in self.softlinks:
                f.write(_softlink(p[0], p[1]))

    def write_end(self):
        self._copy_output_files()
        self._delete_temp_files()
        self._delete_tempdir()
        self._make_softlinks()

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

def _picard_index(in_bam, index, picard_memory, picard_path, tempdir, bg=False):
    """
    Index bam file using Picard Tools.

    Parameters
    ----------
    in_bam : str
        Path to file input bam file.

    index : str
        Path to index file for input bam file.

    bg : boolean
        Whether to run the process in the background.

    Returns
    -------
    line : str
        Line to print to shell/pbs script.

    """
    line = (' \\\n'.join(['java -Xmx{}g -jar'.format(picard_memory),
                          '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                          '\t-Djava.io.tmpdir={}'.format(tempdir),
                          '\t-jar {} BuildBamIndex'.format(picard_path),
                          '\tI={}'.format(in_bam),
                          '\tO={}'.format(index)]))
    if bg:
        line += ' &\n\n'
    else:
        line += '\n\n'
    return line

def _picard_merge(bams, out_bam, picard_memory, picard_path, tempdir, bg=False):
    """
    Merge bam files using Picard. Input bam files are assumed to be coordinate
    sorted.

    Parameters
    ----------
    bams : str
        Bam files to merge.

    out_bam : str
        Path to output merged bam file.

    bg : boolean
        Whether to run the process in the background.

    Returns
    -------
    line : str
        Line to print to shell/pbs script.

    """
    merge_in = ''.join(['\tI={} \\\n'.format(x) for x in bams])
    lines = ['java -Xmx{}g -jar'.format(picard_memory),
             '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
             '\t-Djava.io.tmpdir={}'.format(tempdir), 
             '\t-jar {} MergeSamFiles'.format(picard_path),
             '\tASSUME_SORTED=TRUE',
             '\tUSE_THREADING=TRUE']
    for bam in bams:
        lines.append('\tI={}'.format(bam))
    lines.append('\tO={}'.format(out_bam))
    line = (' \\\n'.join(lines))
    if bg:
        line += ' &\n\n'
    else:
        line += '\n\n'
    return line

def _samtools_index(in_bam, samtools_path, index=None, bg=False):
    """
    Index bam file using samtools.

    Parameters
    ----------
    in_bam : str
        Path to file input bam file.

    index : str
        Path to index file to be written. If not provided, the index is written
        to the samtools default {in_bam}.bai in the current working directory.

    Returns
    -------
    index : str
        Path to index file for input bam file.

    """
    if index: 
        line = '{} index {} {}'.format(samtools_path, in_bam, index)
    else:
        line = '{} index {}'.format(samtools_path, in_bam)
    if bg:
        line += ' &\n'
    else:
        line += '\n'
    return line

def _picard_mark_duplicates(in_bam, out_bam, duplicate_metrics, picard_path,
                            picard_memory, tempdir, remove_dups=False):
    """
    Mark and optionally remove duplicates using Picard Tools.

    Parameters
    ----------
    in_bam : str
        Path to input bam file.

    out_bam : str
        Path to output bam file.

    duplicate_metrics : str
        Path to index file for input bam file.

    """
    lines = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                           '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                           '\t-Djava.io.tmpdir={}'.format(tempdir), 
                           '\t-jar {} MarkDuplicates'.format(picard_path),
                           '\tMETRICS_FILE={}'.format(duplicate_metrics),
                           '\tVALIDATION_STRINGENCY=SILENT',
                           '\tASSUME_SORTED=TRUE',
                           '\tI={}'.format(in_bam), 
                           '\tO={}'.format(out_bam)]))
    if remove_dups:
        lines += ' \\\n\tREMOVE_DUPLICATES=TRUE\n\n'
    else:
        lines += '\n\n'
    return lines

def _wasp_snp_directory(vcf, directory, sample_name=None):
    """
    Convert VCF file into input files directory and files needed for WASP. Only
    bi-allelic heterozygous sites are used.

    Parameters:
    -----------
    vcf : str
        Path to VCF file.

    directory : str
        Output directory. This is the directory that will hold the files for
        WASP.

    sample_name : str
        If provided, use this sample name to get heterozygous SNPs from VCF
        file.

    """
    chrom = []
    pos = []
    ref = []
    alt = []
    vcf_reader = pyvcf.Reader(open(vcf, 'r'))
    if sample_name:
        def condition(record, sample_name):
            return sample_name in [x.sample for x in record.get_hets()]
    else:
        def condition(record, sample_name):
            return len(record.get_hets()) > 0
    for record in vcf_reader:
        if condition(record, sample_name):
            if len(record.ALT) == 1:
                chrom.append(record.CHROM)
                pos.append(record.POS)
                ref.append(record.REF)
                alt.append(record.ALT[0].sequence)
    df = pd.DataFrame([chrom, pos, ref, alt], 
                      index=['chrom', 'position', 'RefAllele', 'AltAllele']).T
    if not os.path.exists(directory):
        os.makedirs(directory)
    for c in set(df.chrom):
        tdf = df[df.chrom == c]
        if tdf.shape[0] > 0:
            f = gzip.open(os.path.join(directory, '{}.snps.txt.gz'.format(c)),
                          'wb')
            lines = (tdf.position.astype(str) + '\t' + tdf.RefAllele + '\t' +
                     tdf.AltAllele)
            f.write('\n'.join(lines) + '\n')
            f.close()

def wasp_allele_swap(bam, find_intersecting_snps_path, snp_dir, sample_name,
                     outdir, tempdir, conda_env=None, shell=False, threads=6):
    """
    Write pbs or shell script for identifying reads in a bam file that overlap
    specified variants and switching the variant allele. This is done using
    find_intersecting_snps.py from WASP.

    Parameters
    ----------
    bam : str
        Path to input bam file.

    find_intersecting_snps_path : str
        Path to find_intersecting_snps.py script.

    snp_dir : str
        Path to directory containing SNP input files for WASP.
    
    sample_name : str
        Sample name used for naming files etc.

    outdir : str
        Directory to store PBS/shell file and aligment results.

    tempdir : str
        Directory to store temporary files.

    conda_env : str
        If provided, load conda environment with this name.

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    threads : int
        Number of threads to request for PBS script.

    """
    job_suffix = 'wasp_allele_swap'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    shell=shell, conda_env=conda_env)
    
    # Input files.
    temp_bam = job.add_input_file(bam)
    job.copy_input_files()
    
    # Files to copy to output directory.
    prefix = os.path.splitext(os.path.split(bam)[1])[0]
    fns = [
        '{}.keep.bam'.format(prefix),
        '{}.remap.fq1.gz'.format(prefix),
        '{}.remap.fq2.gz'.format(prefix),
        '{}.to.remap.bam'.format(prefix),
        '{}.to.remap.num.gz'.format(prefix)
    ]
    job.output_files_to_copy += fns

    # Run WASP to swap alleles.
    with open(job.filename, "a") as f:
        f.write('python {} -p {} {}\n\n'.format(find_intersecting_snps_path,
                                                temp_bam, snp_dir))
    
    job.write_end()
    return job.filename

def wasp_alignment_compare(to_remap_bam, to_remap_num, remapped_bam,
                           filter_remapped_reads_path, sample_name, outdir,
                           tempdir, picard_path, picard_memory=58, conda_env='',
                           shell=False, threads=6):
    """
    Write pbs or shell script for checking original mapping position of reads
    against remapping after swapping alleles using WASP.

    Parameters
    ----------
    to_remap_bam : str
        Bam file from find_intersecting_snps.py that has reads that will be
        remapped (e.g. *.to.remap.bam).

    to_remap_num : str
        Gzipped text file from find_intersecting_snps.py (e.g.
        *.to.remap.num.gz).

    remapped_bam : str
        Bam file with remapped reads.

    filter_remapped_reads_path : str
        Path to filter_remapped_reads.py script.

    sample_name : str
        Sample name used for naming files etc.

    outdir : str
        Directory to store PBS/shell file and aligment results.

    tempdir : str
        Directory to store temporary files.

    conda_env : str
        If provided, load conda environment with this name.

    shell : boolean
        If true, make a shell script rather than a PBS script.

    threads : int
        Number of threads to request for PBS script.

    """
    job_suffix = 'wasp_alignment_compare'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    shell=shell, conda_env=conda_env)
    # class JobScript:
    #     def __init__(sample_name, job_suffix, tempdir, outdir, threads, shell=False,
    #                  queue='high', conda_env=None, environment=None,
    #                  copy_input=True):

    # Input files.
    temp_to_remap_bam = job.add_input_file(to_remap_bam)
    temp_to_remap_num = job.add_input_file(to_remap_num)
    temp_remapped_bam = job.add_input_file(remapped_bam)
    job.copy_input_files()

    # Files that will be created.
    temp_filtered_bam = os.path.join(
        tempdir, '{}_filtered.bam'.format(sample_name))
    job.temp_files_to_delete.append(temp_filtered_bam)
    coord_sorted_bam = os.path.join(
        tempdir, '{}_filtered_coord_sorted.bam'.format(sample_name))
    job.output_files_to_copy.append(coord_sorted_bam)
    bam_index = coord_sorted_bam + '.bai'
    job.output_files_to_copy.append(bam_index)
   
    with open(job.filename, "a") as f:
        # Run WASP alignment compare.
        f.write('python {} -p {} {} {} {}\n\n'.format(
            filter_remapped_reads_path, temp_to_remap_bam, temp_remapped_bam,
            temp_filtered_bam, temp_to_remap_num))

        # Coordinate sort and index.
        lines = _picard_coord_sort(temp_filtered_bam, coord_sorted_bam,
                                   picard_path, picard_memory, tempdir,
                                   bam_index=bam_index)
        f.write(lines)
        f.write('\nwait\n\n')
    
    job.write_end()
    return job.filename

def wasp_remap(
    r1_fastq, 
    r2_fastq, 
    outdir, 
    sample_name, 
    star_index,
    star_path,
    picard_path,
    samtools_path,
    seq_type,
    conda_env='',
    rgpl='ILLUMINA',
    rgpu='',
    tempdir='/scratch', 
    threads=32, 
    picard_memory=58, 
    shell=False,
):
    """
    Make a PBS or shell script for aligning ATAC-seq reads with STAR. The
    defaults are set for use on the Frazer lab's PBS scheduler on FLC.

    Parameters
    ----------
    r1_fastq : str
        R1 reads from find_intersecting_snps.py to be remapped.

    r2_fastq : str
        R2 reads from find_intersecting_snps.py to be remapped.

    outdir : str
        Directory to store PBS/shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

    star_path : str
        Path to STAR aligner.

    picard_path : str
        Path to Picard tools.

    samtools_path : str
        Path to samtools executable.

    seq_type : str
        Type of data. Currently supports ATAC and RNA.

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
    
    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    assert threads >= 3
    seq_types = ['ATAC', 'RNA']
    assert seq_type in seq_types, ('Only {} currently support for '
                                   'seq_type'.format(', '.join(seq_types)))
    job_suffix = 'wasp_remap'
    job = JobScript(sample_name, job_suffix, outdir, threads, tempdir=tempdir,
                    shell=False, queue='high', conda_env=conda_env,
                    environment=None, copy_input=True)
    
    # Input files.
    temp_r1 = job.add_input_file(r1_fastq)
    temp_r2 = job.add_input_file(r2_fastq)
    job.copy_input_files()

    # Files that will be created.
    aligned_bam = os.path.join(tempdir, 'Aligned.out.bam')
    job.temp_files_to_delete.append(aligned_bam)
    coord_sorted_bam = os.path.join(
        tempdir, '{}_sorted.bam'.format(sample_name))
    job.output_files_to_copy.append(coord_sorted_bam)
    job.temp_files_to_delete.append('_STARtmp')

    # Files to copy to output directory.
    job.output_files_to_copy += ['Log.out', 'Log.final.out', 'Log.progress.out',
                                 'SJ.out.tab']

    # Run WASP remapping.
    with open(job.filename, "a") as f:
        # Align with STAR and coordinate sort.
        if seq_type == 'RNA':
            from rnaseq import _star_align
            lines = _star_align([temp_r1], [temp_r2], sample_name, rgpl,
                                rgpu, star_index, star_path, threads)
            f.write(lines)
            f.write('wait\n\n')
            lines = _picard_coord_sort(aligned_bam, coord_sorted_bam,
                                       picard_path, picard_memory, tempdir)
            f.write(lines)
            f.write('wait\n\n')
        
        elif seq_type == 'ATAC':
            from atacseq import _star_align
            lines = _star_align(temp_r1, temp_r2, sample_name, rgpl,
                                rgpu, star_index, star_path, threads)
            f.write(lines)
            f.write('wait\n\n')

            # Coordinate sort bam file.
            lines = _picard_coord_sort(aligned_bam, coord_sorted_bam,
                                       picard_path, picard_memory, tempdir)
            f.write(lines)
            f.write('wait\n\n')

    job.write_end()
    return job.filename

def _mbased(infile, locus_outfile, snv_outfile, sample_name, 
            is_phased=False, num_sim=1000000, threads=1):
    """
    Make a PBS or shell script for running MBASED to determine allelic bias from
    sequencing reads.

    Parameters
    ----------
    infile : str
        Tab-separated file with following columns: chrom, pos, ref_allele,
        alt_allele, locus, name, ref_count, alt_count.

    locus_outfile : str
        Path to file to store locus-level results.

    snv_outfile : str
        Path to file to store SNV-level results.

    sample_name : str
        Sample name used for naming files etc.

    is_phased : bool
        Whether the input file is phased. If so, the reference alleles are
        assumed to be in phase. Note that this only matter locus by locus.

    num_sim : int
        Number of simulations for MBASED to perform.

    threads : int
        Number of threads for MBASED to use.
    
    Returns
    -------
    lines : str
        Lines to be printed to PBS/shell script.

    """
    from __init__ import scripts
    is_phased = str(is_phased).upper()
    script = os.path.join(scripts, 'mbased.R')
    lines = 'Rscript '
    lines += ' '.join([script, infile, locus_outfile, snv_outfile,
                      sample_name, is_phased, str(num_sim), str(threads)])
    lines += '\n'
    return lines

def run_mbased(
    infile, 
    outdir, 
    sample_name, 
    environment=None, 
    is_phased=False,
    num_sim=1000000,
    threads=6, 
    shell=False,
):
    """
    Make a PBS or shell script for running MBASED to determine allelic bias from
    sequencing reads.

    Parameters
    ----------
    infile : str
        Tab-separated file with following columns: chrom, pos, ref_allele,
        alt_allele, locus, name, ref_count, alt_count.

    outdir : str
        Directory to store PBS/shell file and MBASED results.

    sample_name : str
        Sample name used for naming files etc.

    environment : str
        This file will be sourced to set PATH variables for R.

    is_phased : bool
        Whether the input file is phased. If so, the reference alleles are
        assumed to be in phase. Note that this only matter locus by locus.

    num_sim : int
        Number of simulations for MBASED to perform.

    threads : int
        Number of threads to reserve using PBS scheduler and for MBASED to use.

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    assert threads >= 1
    job_suffix = 'mbased'
    job = JobScript(sample_name, job_suffix, outdir, threads, shell=False,
                    queue='high', environment=environment, copy_input=True)
    
    # I'm going to define some file names used later.
    locus_outfile = os.path.join(outdir, '{}_locus.tsv'.format(sample_name))
    snv_outfile = os.path.join(outdir, '{}_snv.tsv'.format(sample_name))
    
    with open(job.filename, "a") as f:
        lines = _mbased(infile, locus_outfile, snv_outfile, sample_name, 
                        is_phased=is_phased, num_sim=num_sim, threads=threads)
        f.write(lines)
        f.write('wait\n\n')
    
    job.write_end()
    return job.filename

def convert_sra_to_fastq(
    sra_files, 
    outdir, 
    sample_name, 
    sra_toolkit_path,
    remove_sra_files=False,
    max_threads=32,
    threads_per_sra=4,
    shell=False,
):
    """
    Make a PBS or shell script for converting one or more SRA files into fastq
    files. All R1 and R2 files will be concatenated and gzipped into two single
    output files.

    Parameters
    ----------
    sra_files : list
        List of SRA files to convert.

    outdir : str
        Directory to store PBS/shell file and gzipped fastq files.

    sample_name : str
        Sample name used for naming files etc.

    sra_toolkit_path : str
        Path to the SRA toolkit bin.

    remove_sra_files : bool
        Whether to remove original SRA files after conversion is complete.

    max_threads : int
        Maximum number of threads to request from PBS scheduler.

    threads_per_sra : int
        Request this many threads per SRA input file. max_threads /
        threads_per_sra SRA files will be converted at a time. For instance, if
        you provide 3 SRA files, threads_per_sra=4, and max_threads=10, then the
        first two SRA files will be converted at the same time. After they are
        done, the last file will be converted. This is done naively using
        background processes (&).

    shell : boolean
        If true, make a shell script rather than a PBS script.
    
    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    threads = min(threads_per_sra * len(sra_files), max_threads)

    if shell:
        pbs = False
    else: 
        pbs = True

    tempdir = os.path.join(tempdir, '{}_sra_fastq'.format(sample_name))
    outdir = os.path.join(outdir, '{}_sra_fastq'.format(sample_name))

    # I'm going to define some file names used later.
    r1 = os.path.join(tempdir, '{}.R1.fastq.gz'.format(sample_name))
    r2 = os.path.join(tempdir, '{}.R2.fsatq.gz'.format(sample_name))
    temp_sra_files = [os.path.join(tempdir, os.path.split(x)[1]) for x in
                      sra_files]
    
    # Files to copy to output directory.
    files_to_copy = [r1, r2]
    
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = temp_sra_files

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_sra_fastq.sh'.format(sample_name))
    else:
        fn = os.path.join(outdir, '{}_sra_fastq.pbs'.format(sample_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir, '{}_sra_fastq.out'.format(sample_name))
        err = os.path.join(outdir, '{}_sra_fastq.err'.format(sample_name))
        job_name = '{}_sra_fastq'.format(sample_name)
        f.write(_pbs_header(out, err, job_name, threads))

    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))
    f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
        ' \\\n\t'.join(sra_files), tempdir))

    def chunks(l, n):
        """Yield successive n-sized chunks from l."""
        for i in xrange(0, len(l), n):
            yield l[i:i+n]
    
    c = chunks(temp_sra_files, threads / threads_per_sra)
    while True:
        try:
            n = c.next()
            for sra in n:
                f.write('{} {} --split-files &\n'.format(
                    os.path.join(sra_toolkit_path, 'fastq-dump'), sra))
            f.write('\nwait\n\n')
        except StopIteration:
            continue
        
    # Concatenate files, pass through awk to remove unneeded stuff, gzip.
    f.write('cat *_1.fastq | awk \'{if (NR % 4 == 1) {print "@"$2} '
            'if (NR % 4 == 2 || NR % 4 == 0) {print $1} '
            'if (NR % 4 == 3) {print "+"}}\' | '
            'gzip -c > ' + r1 + ' &\n\n')
    f.write('cat *_2.fastq | awk \'{if (NR % 4 == 1) {print "@"$2} '
            'if (NR % 4 == 2 || NR % 4 == 0) {print $1} '
            'if (NR % 4 == 3) {print "+"}}\' | '
            'gzip -c > ' + r2 + '\n\n')
    f.write('wait\n\n')

    if remove_sra_files:
        f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(sra_files)))
            
    f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
        ' \\\n\t'.join(files_to_copy),
        outdir))
    f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))
    f.close()

    return fn

def merge_bams(
    bams, 
    outdir, 
    tempdir,
    merged_name, 
    picard_path,
    picard_memory,
    index=True,
    copy_bams=True,
    threads=8,
    shell=False,
):
    """
    Make a PBS or shell script for combining multiple bam files using Picard.

    Parameters
    ----------
    bams : list
        List of SRA files to convert.

    outdir : str
        Directory to store PBS/shell file and merged bam file.

    merged_name : str
        Name used for output directory, files etc.

    picard_path : str
        Path to Picard.

    index : bool
        Whether to index the merged bam file.

    copy_bams : bool
        Whether to copy the input bam files to the temp directory. Not
        necessarcy if temp directory is on the same file system as bam files.

    threads : int
        Number of threads to request from PBS scheduler.

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

    tempdir = os.path.join(tempdir, '{}_merged_bam'.format(merged_name))
    outdir = os.path.join(outdir, '{}_merged_bam'.format(merged_name))

    # I'm going to define some file names used later.
    merged_bam = os.path.join(tempdir,
                                '{}_merged.bam'.format(merged_name))
    merged_bam_index = os.path.join(tempdir,
                                '{}_merged.bam.bai'.format(merged_name))
    
    # Files to copy to output directory.
    files_to_copy = [merged_bam]
    if index:
        files_to_copy.append(merged_bam_index)
    
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = []

    try:
        os.makedirs(outdir)
    except OSError:
        pass

    if shell:
        fn = os.path.join(outdir, '{}_merged_bam.sh'.format(merged_name))
    else:
        fn = os.path.join(outdir, '{}_merged_bam.pbs'.format(merged_name))

    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(outdir, '{}_merged_bam.out'.format(merged_name))
        err = os.path.join(outdir, '{}_merged_bam.err'.format(merged_name))
        job_name = '{}_merged_bam'.format(merged_name)
        f.write(_pbs_header(out, err, job_name, threads))

    f.write('mkdir -p {}\n'.format(tempdir))
    f.write('cd {}\n'.format(tempdir))
    if copy_bams:
        f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
            ' \\\n\t'.join(bams), tempdir))
        bams = [os.path.split(x)[1] for x in bams]
        files_to_remove += bams

    lines = _picard_merge(bams, merged_bam, picard_memory, picard_path,
                          tempdir)
    f.write(lines)

    lines = _picard_index(merged_bam, merged_bam_index, picard_memory,
                          picard_path, tempdir)
    f.write(lines)

    f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format(
        ' \\\n\t'.join(files_to_copy),
        outdir))
    if len(files_to_remove) > 0:
        f.write('rm \\\n\t{}\n\n'.format(' \\\n\t'.join(files_to_remove)))

    if tempdir != outdir:
        f.write('rm -r {}\n'.format(tempdir))
    f.close()

    return fn
