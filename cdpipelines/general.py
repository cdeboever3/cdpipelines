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

class JobScript:
    def __init__(
        self, 
        sample_name, 
        job_suffix, 
        outdir, 
        threads, 
        memory,
        linkdir=None, 
        webpath=None, 
        tempdir=None, 
        queue=None,
        conda_env=None, 
        modules=None, 
        wait_for=None, 
        copy_input=False,
    ):
        """
        Create SGE/shell script object.

        Parameters
        ----------
        sample_name : str
            Sample name used for naming directories, files, etc.

        job_suffix : str
            This suffix will be used for naming directories, files, etc.

        outdir: str
            Path to directory where final output files should be stored.

        threads : int
            Number of threads to request for SGE scripts and to use for
            multi-threaded software.

        memory : int
            Amount of memory in Gb to request for SGE scripts.

        linkdir : str
            Path to directory to make softlinks for files. This is useful if you
            want to softlink a subset of output files into a specific directory
            for viewing on the web for instance.
        
        webpath : str
            URL that points to linkdir. A file test.bam that is linked into
            linkdir is assumed to be available at webpath/test.bam. 

        tempdir : str
            Path to directory where temporary directory should be made. If not
            provided, the output directory will be used as the temp directory.

        queue : str
            SGE queue to use if writing SGE script. If not provided, jobs will
            go into the default week queue.

        conda_env : str
            Path to conda environment to load when job begins.

        modules : str
            Modules (separated by commas e.g. bedtools,samtools) to load at
            beginning of script.

        wait_for : list
            A list of jobnames to wait for before starting this job. This is
            accomplished using -hold_jid.

        copy_input : bool
            Whether to copy input files to temp directory. 

        """
        # Sample name used for naming files.
        self.sample_name = sample_name
        # Brief description of job.
        self.job_suffix = job_suffix
        # Job name for SGE.
        self.jobname = 'job_{}_{}'.format(sample_name, job_suffix)
        # Directory for storing output files.
        self.outdir = outdir
        # Directory for storing temp files. The job is executed in this
        # directory.
        if tempdir:
            self.tempdir = os.path.realpath(os.path.join(tempdir, self.jobname))
        else:
            self.tempdir = self.outdir
        _make_dir(self.outdir)
        # Number of threads to request for job from SGE.
        assert type(threads) is int
        self.threads = threads
        # Amount of memory to request for job from SGE.
        assert type(memory) is int
        self.memory = memory
        # Directory to make softlinks in.
        _make_dir(linkdir)
        self.linkdir = linkdir
        # Path on the web where output files linked to self.linkdir will be
        # available.
        self.webpath = webpath
        # Queue to request for job from SGE.
        self.queue = queue
        # Anaconda Python environment to load at the beginning of the shell
        # script.
        self.conda_env = conda_env
        # Software modules to load at the beginning of the shell script.
        if modules:
            self.modules = modules.split(',')
        else:
            self.modules = None
        # Make directory to store stdout and stderr files. 
        _make_dir(os.path.join(os.path.split(outdir)[0], 'logs'))
        self.out = os.path.join(os.path.split(self.outdir)[0], 'logs',
                                '{}.out'.format(self.jobname))
        self.err = os.path.join(os.path.split(self.outdir)[0], 'logs',
                                '{}.err'.format(self.jobname))
        # List of job names to wait for when submitting to SGE.
        self.wait_for = wait_for
        # Whether to copy input files to tempdir.
        self.copy_input = copy_input
        # List of input files to copy to tempdir.
        self._input_files_to_copy = []
        # List of output files to copy from tempdir to outdir at end of job.
        self._output_files_to_copy = []
        # List of temp files to delete at the end of the job. It's good practice
        # to delete anything you don't want to copy to the outdir.
        self._temp_files_to_delete = []
        # List of [target, link_name] pairs to create softlinks for at the end
        # of the shell script.
        self._softlinks = []
        # Whether to delete shell script after we create it.
        self.delete_sh = False
        # File to write web URLs and tracklines to.
        self.links_tracklines = os.path.join(
            self.outdir, '{}_links_tracklines.txt'.format(self.sample_name))
        # Set shell script file name.
        self._set_filename()
        # Write shell/SGE header.
        self._write_header()

    def _set_filename(self):
        """Make SGE/shell script filename. If a shell script already exists, the
        current shell script will be stored as a temp file."""
        _make_dir(os.path.join(os.path.split(self.outdir)[0], 'sh'))
        self.filename = os.path.join(os.path.split(self.outdir)[0], 'sh',
                                     '{}.sh'.format(self.jobname[4:]))
        # If the shell script already exists we'll assume that the script is
        # just being created to get the output filenames etc. so we'll write to
        # a temp file.
        if os.path.exists(self.filename):
            import tempfile
            self.filename = tempfile.NamedTemporaryFile(delete=False).name
            self.delete_sh = True
    
    def _write_header(self):
        with open(self.filename, "a") as f:
            f.write('#!/bin/bash\n\n')
            if self.queue:
                f.write('#$ -l {}\n'.format(self.queue))
            f.write('#$ -N {}\n'.format(self.jobname))
            f.write('#$ -l h_vmem={}G\n'.format(
                self.memory / float(self.threads)))
            f.write('#$ -pe smp {}\n'.format(self.threads))
            f.write('#$ -S /bin/bash\n')
            f.write('#$ -o {}\n'.format(self.out))
            f.write('#$ -e {}\n\n'.format(self.err))
            if self.modules:
                for module in self.modules:
                    f.write('module load {}\n\n'.format(module))
            if self.conda_env:
                f.write('source activate {}\n\n'.format(self.conda_env))
            if self.tempdir:
                f.write('mkdir -p {}\n'.format(self.tempdir))
                f.write('cd {}\n\n'.format(self.tempdir))

    def _copy_output_files(self):
        if (len(self._output_files_to_copy) > 0 and 
            os.path.realpath(self.tempdir) != os.path.realpath(self.outdir)):
            with open(self.filename, "a") as f:
                f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
                    '\\\n\t'.join(self._output_files_to_copy),
                    self.outdir))

    def _delete_temp_files(self):
        if len(self._temp_files_to_delete) > 0:
            if (os.path.realpath(self.tempdir) == os.path.realpath(self.outdir) 
                or self.tempdir is None):
                self._temp_files_to_delete = [
                    x for x in self._temp_files_to_delete if x not in
                    self._output_files_to_copy
                ]
            with open(self.filename, "a") as f:
                f.write('rm -r \\\n\t{}\n\n'.format(
                    ' \\\n\t'.join(list(set(self._temp_files_to_delete)))))

    def _delete_tempdir(self):
        if self.tempdir and (os.path.realpath(self.tempdir) !=
                             os.path.realpath(self.outdir)):
            with open(self.filename, "a") as f:
                f.write('rm -r {}\n'.format(self.tempdir))

    def _make_softlinks(self):
        """Softlinks are made at the end of shell script when all of the files
        are in their final locations."""
        with open(self.filename, "a") as f:
            for p in self._softlinks:
                self.softlink(p[0], p[1])

    def softlink(self, target, link):
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
        link : str
            Full path to softlink.
    
        """
        lines = 'ln -s \\\n\t{} \\\n\t{}\n\n'.format(target, link)
        with open(self.filename, "a") as f:
            f.write(lines)
        return link
    
    def add_softlink(self, target, link_name=None):
        """
        Add a softlink for file target. All softlinks are made at the end of the
        job. If link_name is provided, the softlink will be made between target
        and link_name. If only target is provided, the link will be made from
        target to a file with the same name in self.linkdir. If the file name
        does not include self.sample_name, the sample name and an underscore
        will be appended to the front of the file name for linkdir.
    
        Parameters
        ----------
        fn : str
            Full path to file to make link to.
    
        link_name : str
            Full path for link.
    
        Returns
        -------
        link : str
            Path to softlink.
    
        """
        if self.sample_name not in os.path.split(target)[1]:
            name = '{}_{}'.format(self.sample_name, os.path.split(target)[1])
            link = os.path.join(self.linkdir, name)
        else:
            name = os.path.split(target)[1]
            link = os.path.join(self.linkdir, name)
        self._softlinks.append([target, link])
        return link

    def sge_submit_command(self):
        """Get command to submit script."""
        if self.wait_for:
            return 'qsub -hold_jid {} {}'.format(','.join(self.wait_for),
                                                 self.filename)
        else:
            return 'qsub {}'.format(self.filename)

    def add_input_file(self, fn, copy=None, delete_original=False):
        """
        If self.copy_input == True or copy == True, copy input file to temp
        directory. If copy == False, do not copy input file to temp directory.
    
        Parameters
        ----------
        fn : str
            Path to input file.
    
        copy : boolean
            Whether to copy the input file to the temp directory.

        delete_original : boolean
            Whether to delete THE ORIGINAL file (i.e. the file at the original
            real path) at the end of the script. There will be no copy of this
            input file when the script complete if this is True.

        Returns
        -------
        out : str
            If the file is copied, the temp path is returned. If the file is not
            copied, the real path is returned.
    
        """
        out = os.path.realpath(fn)
        if copy or (copy is None and self.copy_input):
            self._input_files_to_copy.append(fn)
            out = self.temp_file_path(fn)
        if delete_original:
            self.add_temp_file(fn)
        return out

    def add_temp_file(self, fn):
        """Add temporary file to self.temp_files_to_delete."""
        self._temp_files_to_delete.append(fn)

    def add_output_file(self, fn):
        """Add output file to be copied to self.outdir. Returns final path (in
        self.outdir)."""
        self._output_files_to_copy.append(fn)
        return self.output_file_path(fn)

    def output_file_path(self, fn):
        """Return the path to the output version of a file."""
        return os.path.join(self.outdir, os.path.split(fn)[1])

    def temp_file_path(self, fn):
        """Return the path to the temporary version of a file."""
        return os.path.join(self.tempdir, os.path.split(fn)[1])

    def copy_input_files(self):
        if len(self._input_files_to_copy) > 0:
            with open(self.filename, "a") as f:
                f.write('rsync -avz \\\n\t{} \\\n \t{}\n\n'.format( 
                    '\\\n\t'.join(self._input_files_to_copy),
                    self.tempdir))
            # We will delete any input files from the temp directory that we
            # copy over.
            self.temp_files_to_delete += [
                self.temp_file_path(x) for x in self._input_files_to_copy]

    def write_end(self):
        self._copy_output_files()
        self._delete_temp_files()
        self._delete_tempdir()
        self._make_softlinks()
        if self.delete_sh:
            os.remove(self.filename)

    def make_het_vcf(
        vcf, 
        regions=None,
        vcf_sample_name=None,
        gatk_fai=None, 
        vcf_chrom_conv=None,
        bcftools_path='bcftools',
    ):
        """
        Extract all biallelic heterozygous variants for this sample.
    
        Parameters:
        -----------
        vcf : str
            Path to VCF file.
    
        regions : str
            Path to bed file to define regions of interests (e.g. exons, peaks,
            etc.). 
    
        vcf_sample_name : str
            Use this sample name to get heterozygous SNPs from VCF file.
    
        gatk_fai : str
            Path to karyotypically sorted fasta index (fai) file that works with
            GATK.  The output VCF file will be sorted in the order of this fai
            file for compatibility with GATK. Assumed to have associated fai and
            dict files.
    
        vcf_chrom_conv : str
            File with VCF chromosomes in first column and corresponding RNA-seq
            chromosomes in second column (no header). This is needed if the VCF
            and RNA-seq data have different chromosome naming.
    
        """
        vcf_out = os.path.join(self.tempdir,
                               '{}_hets.vcf'.format(self.sample_name))
        if regions:
            r = '\\\n\t-R {} '.format(regions)
        else:
            r = ''
        lines = ('{0} view -O u -m2 -M2 {1}\\\n\t-s {2} \\\n\t{3} \\\n\t'
                 '| {0} view -g het '.format(
                     bcftools_path, regions, vcf_sample_name, vcf))
        if vcf_chrom_conv:
            lines += '-Ou \\\n\t| {} annotate --rename-chrs {} '.format(
                bcftools_path, vcf_chrom_conv)
        lines += '> {}'.format(vcf_out)
        with open(self.filename, "a") as f:
            f.write(lines)
        return vcf_out

    def make_md5sum(
        self,
        fn,
    ):
        """Make an md5sum for fn. Returns path to md5 file."""
        lines = 'md5sum {0} > {0}.md5\n\n'.format(fn)
        with open(self.filename, "a") as f:
            f.write(lines)
        return fn + '.md5'

    def merge_bed(
        self,
        bed, 
        bedtools_path='bedtools',
    ):
        """Merge a bed file using bedtools."""
        root = os.path.splitext(os.path.split(bed)[1])[0]
        out_bed = os.path.join(self.tempdir, '{}_merged.bed'.format(root))
        lines = '{0} sort -i {1} | {0} merge -i stdin > {2}'.format(
            bedtools_path, bed, out_bed)
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bed

    def convert_sra_to_fastq(
        self,
        sra_files, 
        fastq_dump_path='fastq-dump',
    ):
        """
        Converting one or more SRA files into fastq files.  All R1 and R2 files
        will be concatenated and gzipped into two single output files.
    
        Parameters
        ----------
        sra_files : list
            List of SRA files or URLs to convert. If URLs are provided, they
            must start with ftp:// or http://.
    
        max_threads : int
            Maximum number of threads to request from SGE scheduler.
    
        fastq-dump : str
            Path to fastq-dump from the SRA toolkit.

        Returns
        -------
        r1 : str
            Path to R1 fastq file.
    
        r2 : str
            Path to R1 fastq file.
    
        """
        lines = ''
        temp_files = []
        r1 = os.path.join(self.tempdir,
                          '{}.R1.fastq.gz'.format(self.sample_name))
        r2 = os.path.join(self.tempdir,
                          '{}.R2.fastq.gz'.format(self.sample_name))
        temp_files += ['{}_1.fastq'.format(
            os.path.join(self.tempdir, os.path.split(x)[1].split('.')[0])) for x
            in sra_files]
        temp_files += ['{}_2.fastq'.format(
            os.path.join(self.tempdir, os.path.split(x)[1].split('.')[0])) for x
            in sra_files]

        # If any of the input are URLs, download the SRA file. We'll delete it
        # when we're done converting it.
        for i,fn in enumerate(sra_files):
            if fn[0:6] == 'ftp://' or fn[0:7] == 'http://':
                lines += 'wget {}\n'.format(fn)
                sra_files[i] = os.path.join(self.tempdir, os.path.split(fn)[1])
                temp_files.append(sra_files[i])
        
        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in xrange(0, len(l), n):
                yield l[i:i+n]
       
        sc = chunks(sra_files, self.threads)
        while True:
            try:
                n = sc.next()
                for sra in n:
                    lines += ('{} {} --split-files &\n'.format(fastq_dump_path,
                                                               sra))
                lines += '\nwait\n\n'
            except StopIteration:
                break
            
        # Concatenate files, pass through awk to remove unneeded stuff, gzip.
        lines += ('cat *_1.fastq \\\n\t| awk \'{if (NR % 4 == 1) {print "@"$2} '
                  'if (NR % 4 == 2 || NR % 4 == 0) {print $1} '
                  'if (NR % 4 == 3) {print "+"}}\' \\\n\t| '
                  'gzip -c \\\n\t> ' + r1 + ' &\n\n')
        lines += ('cat *_2.fastq \\\n\t| awk \'{if (NR % 4 == 1) {print "@"$2} '
                  'if (NR % 4 == 2 || NR % 4 == 0) {print $1} '
                  'if (NR % 4 == 3) {print "+"}}\' \\\n\t| '
                  'gzip -c \\\n\t> ' + r2 + '\n\n')
        lines += 'wait\n\n'

        # Remove raw fastq files from fastq-dump.
        lines += 'rm \\\n\t{}\n\n'.format(' \\\n\t'.join(temp_files))
    
        with open(self.filename, "a") as f:
            f.write(lines)
        return r1, r2

    def homer_motif_analysis(
        self,
        bed,
        size='given',
        mask=True,
        web_available=True,
        write_to_outdir=True,
    ):
        """
        Perform a motif analysis using HOMER on the input bed file. The HOMER
        scripts and other necessary HOMER stuff are assumed to be in your path.

        Parameters
        ----------
        bed : str
            Path to bed file defining intervals to use for motif analysis.

        size : int or str
            Size for HOMER -size parameter. See HOMER documentation for more
            information.

        mask : bool
            Whether to pass the -mask parameter to HOMER.

        web_available : bool
            If True, write URL to self.links_tracklines, make softlink to
            self.linkdir, and set write_to_outdir = True.

        write_to_outdir : bool
            If True, write output files directly to self.outdir.
    
        Returns
        -------
        dy : str
            Path to HOMER output diretory that contains output files.

        """
        if write_to_outdir or web_available:
            dy = self.outdir
        dy = os.path.join(self.tempdir,
                          '{}_homer_motif'.format(self.sample_name))
        lines = ('findMotifsGenome.pl \\\n\t {} \\\n\thg19 \\\n\t{} \\\n\t'
                 '-size {} \\\n\t-p {}'.format(bed, dy, size, self.threads))
        if mask:
            lines += ' \\\n\t-mask'
        lines += '\n\n'
        
        with open(self.filename, "a") as f:
            f.write(lines)
        
        if web_available:
            link = self.add_softlink(dy)
            url = self.webpath + '/' + os.path.split(dy)[1]
            with open(self.links_tracklines, "a") as f:
                f.write(url + '\n')

        return dy

    def scale_bedgraph(
        self,
        bg,
        bam,
        expected_num,
    ):
        """
        Scale a bedgraph file using a STAR Log.final.out file by multipyling by
        expected_num / actual_num. This is designed to make normalized bigwig
        files that can be roughly compared between samples.
    
        Parameters
        ----------
        bg : str 
            Path to input bedgraph file.
    
        bam : str 
            Path to bam file used to make bedgraph file.
    
        expected_num : int
            Number of expected reads (pairs). This only needs to be a rough
            estimate. You should keep this number constant for all samples
            across an experiment/project (i.e. any samples you want to compare
            against each other). For example, say you are doing RNA-seq for 10
            samples on a lane, the lane yields 100,000,000 pairs of reads, and
            you expect 90% of the read pairs to align uniquely. Then you would
            want expected_num = 9,000,000. In this case, you would set
            actual_num to the actual number of uniquely mapped reads. The reason
            you want this number close to the number of expected number of
            uniquely mapped pairs is so the normalization doesn't drastically
            change the coverage for samples near the expected number of input
            read pairs. For instance, if expected_num=20M and a sample has
            exactly 20M uniquely mapped read pairs, then the coverage will not
            be changed.  If expected_num=20M and a sample has only 10M input
            read PAIRS, it will be normalized so that it looks like it had 20M
            input read pairs (i.e. all coverages will be multiplied by 2).
    
        Returns
        -------
        out_bg : str
            Path to output scaled bedgraph file.
        
        """
        root = os.path.splitext(os.path.split(bg)[1])[0]
        out_bg = os.path.join(self.tempdir, '{}_scaled.bg'.format(root))
        lines = 'scale_bedgraph \\\n\t{} \\\n\t{} \\\n\t{} \\\n\t{}\n\n'.format(
            bg, bam, out_bg, expected_num)
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bg

    def featureCounts_count(
        self,
        features,
        bam,
        sort=True,
        filetype='gtf',
        both=False,
        strand_specific=0,
        featureCounts_path='featureCounts',
        root=None,
    ):
        """
        Run featureCounts to count the number of fragments that overlap
        intervals defined in a bed file. The bed file is converted to saf format
        for use with featureCounts.
    
        Parameters
        ----------
        features : str
            Path to bed or gtf file defining regions to count for. The file
            type is inferred using the extension or using the filetype parameter
            below.
        
        bam : str
            Path to bam file to count reads for.

        sort : bool
            If True, featureCounts will sort the bam file so that reads from the
            same pair are next to each other. If False, featureCounts will
            assume reads from the same pair are next to each other.

        filetype : str
            File type of the features file. Either gtf or bed.

        both : bool
            Use featurecounts -B option: count read pairs that have both ends
            successfully aligned only.

        strand_specific : 0, 1, or 2
            Passed to featureCounts -s parameter: Perform strand-specific read
            counting. Possible values:  0 (unstranded), 1 (stranded) and 2
            (reversely stranded).
    
        root : str
            If provided, use for naming the featureCounts output file.
    
        Returns
        -------
        out : str
            Path to output counts file.
    
        out_summary : str
            Path to output counts summary file.
    
        """
        if root is None:
            root = os.path.splitext(os.path.split(bam)[1])[0]
        out = os.path.join(self.tempdir,
                           '{}_featureCounts.tsv'.format(root))
        out_summary = os.path.join(self.tempdir, '{}.summary'.format(out))
        lines = '{} -p -T {} '.format(featureCounts_path, self.threads)
        if both:
            lines += '-B '
        if strand_specific != 0:
            lines += '-s {} '.format(strand_specific)
        if sort == False:
            lines += '--donotsort '
        if os.path.splitext(features)[1] == '.bed' or filetype == 'bed':
            saf = self.convert_bed_to_saf(features)
            self.add_temp_file(saf)
            lines += ('-F SAF \\\n\t-a {} \\\n\t-o {} '
                      '\\\n\t{}\n\n'.format(saf, out, bam))
        elif os.path.splitext(features)[1] == '.gtf' or filetype == 'gtf':
            lines += '\\\n\t-a {} \\\n\t-o {} \\\n\t{}\n\n'.format(
                features, out, bam)
        else:
            assert True == False
        with open(self.filename, "a") as f:
            f.write(lines)
        return out, out_summary

    def convert_bed_to_saf(
        self,
        bed,
    ):
        """
        Convert bed file to saf file for use with featureCounts.
    
        Parameters
        ----------
        bed : str
            Path to bed file to convert.
    
        Returns
        -------
        saf : str
            Path to output saf file.
    
        """
        saf = os.path.join(self.tempdir,
                           os.path.splitext(os.path.split(bed)[1])[0] + '.saf')
        lines = 'convert_bed_to_saf \\\n\t{} \\\n\t{}\n\n'.format(bed, saf)
        with open(self.filename, "a") as f:
            f.write(lines)
        return saf

    def bedgraph_from_bam(
        self,
        bam, 
        bedtools_path='bedtools',
        sambamba_path='sambamba',
    ):
        """
        Create a coverage bedgraph file from a bam file. Only uniquely mapped
        read pairs are used (quality score greater than 255 which is uniquely
        mapped for STAR).
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.
    
        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        bedgraph : str
            Path to output bedgraph file.
    
        """
        bedgraph = os.path.join(self.tempdir, '{}.bg'.format(self.sample_name))

        if bedtools_path == 'bedtools':
            genome_file = 'human.hg19.genome'
        else:
            genome_file = os.path.join(
                os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
                'human.hg19.genome')

        lines = ('{} view -f bam -F "not (unmapped or mate_is_unmapped) and '
                 'mapping_quality >= 255" \\\n\t{} \\\n\t| '
                 '{} genomecov -ibam stdin \\\n\t-g {} -split -bg \\\n\t'
                 '-trackline -trackopts \'name="{}"\' '.format(
                     sambamba_path, bam, bedtools_path, genome_file,
                     self.sample_name))
        lines += ' \\\n\t> {}\n\n'.format(bedgraph)
        with open(self.filename, "a") as f:
            f.write(lines)
        return bedgraph
    
    def picard_collect_rna_seq_metrics(
        self,
        in_bam, 
        ref_flat, 
        rrna_intervals,
        picard_path='$picard',
        strand_specific=True, 
        bg=False,
    ):
        """
        Collect RNA-seq metrics using Picard. The input bam file is assumed to
        be coordinate sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input coordinate sorted bam file.
    
        ref_flat : str
            Path to refFlat file with non-rRNA genes. Can be gzipped.
    
        rrna_intervals : str
            Pato to interval list file with rRNA intervals.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        metrics : str
            Path to output metrics file.
    
        chart : str
            Path to output PDF file.
    
        """
        metrics = os.path.join(self.tempdir,
                               '{}_rna_seq_metrics.txt'.format(self.sample_name))
        chart = os.path.join(self.tempdir,
                             '{}_5_3_coverage.pdf'.format(self.sample_name))
        if strand_specific:
            ss = 'SECOND_READ_TRANSCRIPTION_STRAND'
        else:
            ss = 'NONE'
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} CollectRnaSeqMetrics'.format(picard_path),
            'I={}'.format(in_bam),
            'REF_FLAT={}'.format(ref_flat),
            'STRAND_SPECIFICITY={}'.format(ss),
            'RIBOSOMAL_INTERVALS={}'.format(rrna_intervals),
            'ASSUME_SORTED=TRUE',
            'CHART_OUTPUT={}'.format(chart),
            'O={}'.format(metrics)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return metrics, chart

    def picard_collect_multiple_metrics(
        self,
        in_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Collect multiple metrics using Picard. The input bam file is assumed to
        be sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input coordinate sorted bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        output : tuple
            Tuple of the paths to the following output files:
                alignment_summary_metrics quality_by_cycle.pdf
                base_distribution_by_cycle.pdf quality_by_cycle_metrics
                base_distribution_by_cycle_metrics quality_distribution.pdf
                insert_size_histogram.pdf quality_distribution_metrics
                insert_size_metrics.
    
        """
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} CollectMultipleMetrics'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'ASSUME_SORTED=TRUE',
            'I={}'.format(in_bam), 
            'O={}'.format(self.sample_name)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        output = [os.path.join(self.tempdir, '{}.{}'.format(self.sample_name, x))
                               for x in [
                                   'alignment_summary_metrics',
                                   'quality_by_cycle.pdf',
                                   'base_distribution_by_cycle.pdf',
                                   'quality_by_cycle_metrics',
                                   'base_distribution_by_cycle_metrics',
                                   'quality_distribution.pdf',
                                   'insert_size_histogram.pdf',
                                   'quality_distribution_metrics',
                                   'insert_size_metrics']]
        return tuple(output)
    
    def sambamba_index(
        self,
        in_bam, 
        sambamba_path='sambamba',
    ):
        """
        Index bam file using sambamba.
    
        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        index : str
            Path to output index file.
    
        """
        index = os.path.join(self.tempdir, os.path.split(in_bam)[1] + '.bai')
        lines = '{} index -t {} \\\n\t{} \\\n\t{}\n\n'.format(
            sambamba_path, self.threads, in_bam, index)
        with open(self.filename, "a") as f:
            f.write(lines)
        return index
    
    def sambamba_merge(
        self,
        bams, 
        sambamba_path='sambamba',
    ):
        """
        Merge bam files using sambamba.
    
        Parameters
        ----------
        bams : list
            List of paths to bam files to merge.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        out : str
            Path to output merged bam file.
    
        """
        out = os.path.join(self.tempdir, '{}.bam'.format(self.sample_name))
        lines = '{} merge -t {} \\\n\t{} \\\n\t{}\n\n'.format(
            sambamba_path, self.threads, out, ' \\\n\t'.join(bams))
        with open(self.filename, "a") as f:
            f.write(lines)
        return out
    
    def picard_index(
        self,
        in_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Index bam file using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
    
        index : str
            Path to index file for input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        index : str
            Path to output index file.
    
        """
        index = os.path.join(self.tempdir, os.path.split(in_bam)[1] + '.bai')
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar'.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} BuildBamIndex'.format(picard_path),
            'I={}'.format(in_bam),
            'O={}'.format(index)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return index
    
    def picard_merge(
        self,
        bams, 
        out_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Merge bam files using Picard. Input bam files are assumed to be
        coordinate sorted.
    
        Parameters
        ----------
        bams : str
            Bam files to merge.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        out_bam : str
            Path to output merged bam file.
    
        """
        merge_in = ''.join(['\tI={} \\\n'.format(x) for x in bams])
        lines = ['java -Xmx{}g -jar'.format(int(self.memory * 0.8)),
                 '\t-XX:ParallelGCThreads=1',
                 '\t-Djava.io.tmpdir={}'.format(self.tempdir), 
                 '\t-jar {} MergeSamFiles'.format(picard_path),
                 '\tASSUME_SORTED=TRUE',
                 '\tUSE_THREADING=TRUE']
        for bam in bams:
            lines.append('\tI={}'.format(bam))
        lines.append('\tO={}'.format(out_bam))
        lines = (' \\\n'.join(lines))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam
    
    def samtools_index(
        self,
        in_bam, 
        bg=False,
        samtools_path='samtools',
    ):
        """
        Index bam file using samtools.
    
        Parameters
        ----------
        in_bam : str
            Path to file input bam file.
    
        index : str
            Path to index file to be written. If not provided, the index is
            written to the samtools default {in_bam}.bai in the current working
            directory.
    
        Returns
        -------
        index : str
            Path to index file for input bam file.
    
        """
        index = os.path.join(self.tempdir, os.path.splitext(in_bam)[0] + '.bai')
        line = '{} index {}'.format(samtools_path, in_bam)
        if bg:
            line += ' &\n\n'
        else:
            line += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return index

    def biobambam2_mark_duplicates(
        self,
        in_bam,
        remove_duplicates=False,
        bammarkduplicates_path='bammarkduplicates',
    ):
        """
        Mark and optionally remove duplicates using biobambam2.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        bammarkduplicates_path : str
            Path to bammarkduplicates.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        dup_metrics : str
            Path to index file for input bam file.

        removed_reads : str
            Path to bam file with removed reads if remove_duplicates == True.
    
        """
        t = 'mdup'
        if remove_duplicates:
            t = 'r' + t
        mdup_bam = os.path.join(
            self.tempdir, '{}_sorted_{}.bam'.format(self.sample_name, t))
        dup_metrics = os.path.join(
            self.tempdir, '{}_duplicate_metrics.txt'.format(self.sample_name))
        lines = ('{} markthreads={} \\\n\ttmpfile={} \\\n\tI={} \\\n\t'
                 'O={} \\\n\tM={}'.format(
                     bammarkduplicates_path, self.threads, self.tempdir, in_bam,
                     mdup_bam, dup_metrics))
        if remove_duplicates:
            removed_reads = os.path.join(
                self.tempdir,
                '{}_removed_duplicates.bam'.format(self.sample_name))
            lines += ('\\\n\trmdup=1 \\\n\tD={}'.format(removed_reads))
        lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        if remove_duplicates:
            return mdup_bam, dup_metrics, removed_reads
        else:
            return mdup_bam, dup_metrics
    
    def picard_mark_duplicates(
        self,
        in_bam, 
        picard_path='$picard',
        remove_dups=False,
    ):
        """
        Mark and optionally remove duplicates using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        dup_metrics : str
            Path to index file for input bam file.
    
        """
        mdup_bam = os.path.join(
            self.tempdir, '{}_sorted_mdup.bam'.format(self.sample_name))
        dup_metrics = os.path.join(
            self.tempdir, '{}_duplicate_metrics.txt'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} MarkDuplicates'.format(picard_path),
            'METRICS_FILE={}'.format(dup_metrics),
            'VALIDATION_STRINGENCY=SILENT',
            'ASSUME_SORTED=TRUE',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam)]))
        if remove_dups:
            lines += ' \\\n\tREMOVE_DUPLICATES=TRUE\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return mdup_bam, dup_metrics

    def picard_gc_bias_metrics(
        self,
        in_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Collect GC bias metrics using Picard. The input bam file is assumed to
        be sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        metrics : str
            Path to output metrics file.
    
        chart : str
            Path to output PDF chart.
    
        out : str
            Path to picard output file.
    
        """
        metrics = os.path.join(self.tempdir,
                               '{}_gc_bias_metrics.txt'.format(self.sample_name))
        chart = os.path.join(self.tempdir,
                             '{}_gc_bias.pdf'.format(self.sample_name))
        out = os.path.join(self.tempdir,
                           '{}_gc_bias_metrics_out.txt'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} CollectGcBiasMetrics'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out),
            'CHART_OUTPUT={}'.format(chart),
            'SUMMARY_OUTPUT={}'.format(metrics),
            'ASSUME_SORTED=TRUE']))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return metrics, chart, out
    
    def picard_bam_index_stats(
        self,
        in_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Collect bam index stats with Picard.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.

        Returns
        -------
        out : str
            Path to write stdout to. This contains the index stats.
    
        err : str
            Path to write stderr to. This contains picard stuff.
    
        """
        out = os.path.join(self.outdir,
                           '{}_index_stats.txt'.format(self.sample_name))
        err = os.path.join(self.outdir,
                           '{}_index_stats.err'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} BamIndexStats'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam),
            '> {}'.format(out),
            '2> {}'.format(err)]))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out, err
    
    def picard_insert_size_metrics(
        self,
        in_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Collect insert size metrics using Picard. The input bam file is assumed
        to be sorted.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.
   
        Returns
        -------
        metrics : str
            Path to output metrics file.
    
        hist : str
            Path to output histogram PDF.
    
        """
        metrics = os.path.join(
            self.tempdir, '{}_insert_size_metrics.txt'.format(self.sample_name))
        hist = os.path.join(self.tempdir,
                            '{}_insert_size.pdf'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} CollectInsertSizeMetrics'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(metrics),
            'HISTOGRAM_FILE={}'.format(hist),
            'ASSUME_SORTED=TRUE']))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return metrics, hist
    
    def picard_query_sort(
        in_bam, 
        picard_path='$picard',
        bg=False,
    ):
        """
        Query sort using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        bg : boolean
            Whether to run the process in the background.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        """
        out_bam = os.path.join(self.tempdir,
                               '{}_qsorted.bam'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} SortSam'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam),
            'SO=queryname']))
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam
    
    def sambamba_sort(
        self,
        in_bam, 
        queryname=False,
        tempdir='.',
        root=None,
        sambamba_path='sambamba',
    ):
        """
        Coordinate sort using sambamba.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
        
        queryname : bool
            If True, sort by query name.

        tempdir : str
            Path to directory to use for temporary files. Default is current
            directory.

        root : str
            If provided, use for naming the file [root]_query_sorted.bam or
            [root]_sorted.bam.
    
        sambamba_path : str
            Path to sambaba executable.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        """
        if not root:
            root = self.sample_name
        if queryname:
            out_bam = os.path.join(
                self.tempdir, '{}_query_sorted.bam'.format(root))
        else:
            out_bam = os.path.join(
                self.tempdir, '{}_sorted.bam'.format(root))
        lines = '{} sort -m {}GB -t {}'.format(sambamba_path, self.memory - 2,
                                               self.threads)
        if queryname:
            lines += ' -n'
        lines += ' \\\n\t'
        lines += '--tmpdir {} \\\n\t'.format(tempdir)
        lines += '{} \\\n\t'.format(in_bam)
        lines += '-o {}\n\n'.format(out_bam)

        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam
    
    def picard_coord_sort(
        self,
        in_bam, 
        index=False,
        picard_path='$picard',
    ):
        """
        Coordinate sort using Picard Tools.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        out_index : str
            Path to output index file. Only returned if index == True.
    
        """
        out_bam = os.path.join(self.tempdir,
                               '{}_sorted.bam'.format(self.sample_name))
        if index:
            out_index = os.path.join(out_bam + '.bai')
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} SortSam'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam),
            'SO=coordinate']))
        if index:
            lines += ' \\\n\tCREATE_INDEX=TRUE\n\n' 
            old_index = '.'.join(out_bam.split('.')[0:-1]) + '.bai'
            lines += 'mv {} {}'.format(old_index, out_index)
        lines += '\n\n'

        with open(self.filename, "a") as f:
            f.write(lines)
        if index:
            return out_bam, out_index
        else: 
            return out_bam
    
    def picard_reorder(
        self,
        in_bam, 
        fasta,
        picard_path='$picard',
    ):
        """
        Reorder bam file according to seq_dict.
    
        Parameters
        ----------
        in_bam : str
            Path to input bam file.
    
        picard_path : str
            Path to picard jar file. Default assumes the environmental variable
            $picard contains the path to the jar file.
    
        Returns
        -------
        out_bam : str
            Path to output bam file.
    
        """
        out_bam = os.path.join(self.tempdir,
                               '{}_reordered.bam'.format(self.sample_name))
        lines = (' \\\n\t'.join([
            'java -Xmx{}g -jar '.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {} ReorderSam'.format(picard_path),
            'VALIDATION_STRINGENCY=SILENT',
            'I={}'.format(in_bam), 
            'O={}'.format(out_bam),
            'REFERENCE={}'.format(fasta)])) + '\n\n'

        with open(self.filename, "a") as f:
            f.write(lines)
        return out_bam
    
    def cutadapt_trim(
        self,
        fastq, 
        length, 
        out, 
        bg=False,
    ):
        """
        Cut a specified number of bases from a fastq file using cutadapt.
        Cutadapt should be installed in your Python environment.
    
        Parameters
        ----------
        fastq : str
            Fastq or gzipped/bzipped fastq.
    
        length : int
            Positive or negative integer. Positive numbers remove bases at the
            front of the read and negative numbers remove bases at the end of
            the read.
    
        out : str
            Path to output (optionally gzipped/bzipped) fastq files.
    
        bg : boolean
            Whether to run the process in the background (i.e. include an
            ampersand at the end of the command).
    
        Returns
        -------
        out : str
            Path to output (optionally gzipped/bzipped) fastq files.
    
        """
        line = 'cutadapt --cut {} -o {} {}'.format(length, out, fastq)
        if bg:
            line += ' &\n\n'
        else:
            line += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out
    
    def bedgraph_to_bigwig(
        self,
        bedgraph, 
        bedgraph_to_bigwig_path='bedGraphToBigWig',
        bedtools_path='bedtools',
    ):
        """
        Convert bedgraph file to bigwig file.
    
        Parameters
        ----------
        bedgraph : str
            Input bedgraph file.
    
        bedgraph_to_bigwig_path : str
            Path bedGraphToBigWig executable from UCSC.

        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        bigwig : str
            Path to output bigwig file.
    
        """
        bigwig = os.path.join(
            self.tempdir,
            '{}.bw'.format(os.path.splitext(os.path.split(bedgraph)[1])))
        # If bedtools is in the path, I'll assume the genome file is as well.
        if bedtools_path == 'bedtools':
            bedtools_genome_path = 'human.hg19.genome'
        else:
            bedtools_genome_path = os.path.join(
                os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
                'human.hg19.genome')
        lines =  ' \\\n\t'.join([
            '{} {}'.format(bedgraph_to_bigwig_path, bedgraph),
            '{}'.format(bedtools_genome_path),
            '{} &\n'.format(bigwig)])
        with open(self.filename, "a") as f:
            f.write(lines)
        return bigwig
    
    def flagstat(
        self,
        bam, 
        samtools_path='samtools',
        bg=False,
    ):
        """
        Run flagstat for a bam file.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.
    
        samtools_path : str
            Path to samtools executable.
    
        Returns
        stats_file : str
            File to write flagstats to.
    
        -------
        """
        stats_file = os.path.join(
            self.tempdir, '{}_flagstat.txt'.format(
                os.path.splitext(os.path.split(bam)[1])[0]))
        lines = '{} flagstat {} > {}'.format(samtools_path, bam, stats_file)
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return stats_file
    
    def combine_fastqs(
        self,
        fastqs,
        suffix=None,
        bg=False,
    ):
        """
        Cat the fastqs together into a single file. If fastqs only contains one
        fastq file, make a softlink to that file.

        Parameters
        ----------
        fastqs : list
            List of paths to gzipped fastq files.
    
        suffix : str
            Add this to the combined file name (for instance, R1 or R2).
    
        Returns
        -------
        out_fastq : str
            Path to single output fastq file.
    
        """
        root = '{}_combined'.format(self.sample_name)
        if suffix:
            root += '_' + suffix
        out_fastq = os.path.join(self.tempdir, root + '.fastq.gz')
        fastqs = sorted(fastqs)
        if len(fastqs) > 1:
            lines = 'cat \\\n\t{} \\\n\t> {}'.format(' \\\n\t'.join(fastqs),
                                                     out_fastq)
        else:
            lines = 'ln -s \\\n\t{} \\\n\t{}'.format(fastqs[0], out_fastq)
        if bg:
            lines += ' &\n\n'
        else:
            lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return out_fastq
    
    def fastqc(
        self,
        fastqs, 
        web_available=True,
        write_to_outdir=True,
        fastqc_path='fastqc',
    ):
        """
        Run FastQC for fastq files in fastqs. Links to linkdir are automatically
        created for the output html files if linkdir != None.
    
        Parameters
        ----------
        fastqs : list
            List of paths to fastq files.
    
        web_available : bool
            If True, write URL to self.links_tracklines, make softlink to
            self.linkdir, and set write_to_outdir = True.

        write_to_outdir : bool
            If True, write output files directly to self.outdir.
    
        fastqc_path : str
            Path to FastQC.
    
        Returns
        -------
        fastqc_html : list
            List of paths to the html files for all input fastq files.

        fastqc_zip : str or list
            List of paths to the zip files for all input fastq files.
    
        """
        if write_to_outdir or web_available:
            dy = self.outdir
        else:
            dy = self.tempdir
        fastqc_html = []
        fastqc_zip = []
        for fq in fastqs:
            fastqc_html.append(
                os.path.join(dy, '{}_fastqc.html'.format(
                    '.'.join(os.path.split(fq)[1].split('.')[0:-2]))))
            fastqc_zip.append(
                os.path.join(dy, '{}_fastqc.zip'.format(
                    '.'.join(os.path.split(fq)[1].split('.')[0:-2]))))
        fastqs = ' \\\n\t'.join(fastqs)

        lines = ('{} --outdir {} --nogroup --threads {} \\\n'
                 '\t{}\n'.format(fastqc_path, dy, self.threads, fastqs))
        
        with open(self.filename, "a") as f:
            f.write(lines)

        if web_available:
            if self.linkdir:
                for html in fastqc_html:
                    link = self.add_softlink(html)
            url = self.webpath + '/' + os.path.split(html)[1]
            with open(self.links_tracklines, "a") as f:
                f.write(url + '\n')

        return fastqc_html, fastqc_zip
    
    def wasp_allele_swap(
        self,
        bam, 
        find_intersecting_snps_path, 
        vcfs, 
        bed,
        gatk_fai=None,
        vcf_sample_name=None, 
        vcf_chrom_conv=None,
        samtools_path='samtools',
        bcftools_path='bcftools',
    ):
        """
        Swap alleles for reads that overlap heterozygous variants.
    
        Parameters
        ----------
        bam : str
            Path to input bam file.
    
        find_intersecting_snps_path : str
            Path to find_intersecting_snps.py script.

        vcf : list
            List of VCF files to extract heterozygous variants from.

        bed : str
            Path to bed file defining regions to look for variants in. This file
            will be merged to make sure there are no overlapping regions.
       
        vcf_sample_name : str
            Sample name of this sample in the VCF file (if different than
            sample_name). For instance, the sample name in the VCF file may be
            the sample name for WGS data which may differ from the RNA-seq
            sample name.

        vcf_chrom_conv : str
            File with VCF chromosomes in first column and corresponding RNA-seq
            chromosomes in second column (no header). This is needed if the VCF
            and RNA-seq data have different chromosome naming.

        gatk_fai : str
            Path to karyotypically sorted fasta index (fai) file that works with
            GATK.  The output VCF file will be sorted in the order of this fai
            file for compatibility with GATK. Assumed to have associated fai and
            dict files.

        vcf_sample_name : str
            Sample name of this sample in the VCF file (if different than
            sample_name).

        Returns
        -------
        snp_directory : str
            Path to WASP input SNP directory.
        
        vcf_out : str
            Path to VCF file with heterozygous variants for this subject.

        keep_bam : str
            Path to WASP bam file of reads that do not overlap heterozygous
            variants.

        wasp_r1_fastq : str
            Path to WASP fastq file of reads to remap.

        wasp_r2_fastq : str
            Path to WASP fastq file of reads to remap.

        to_remap_bam : str
            Path to WASP bam file of reads that need to be remapped.

        to_remap_num
            Path to WASP file that specifies how many realignments need to be
            done per read pair.
        """
        if not vcf_sample_name:
            vcf_sample_name = sample_name

        # Files that will be created.
        uniq_bam = os.path.join(
            self.tempdir, '{}_uniq.bam'.format(self.sample_name))
        prefix = '{}_uniq'.format(self.sample_name)
        keep_bam = os.path.join(self.tempdir, '{}.keep.bam'.format(prefix))
        wasp_r1_fastq = os.path.join(self.tempdir,
                                     '{}.remap.fq1.gz'.format(prefix))
        wasp_r2_fastq = os.path.join(self.tempdir,
                                     '{}.remap.fq2.gz'.format(prefix))
        to_remap_bam = os.path.join(self.tempdir,
                                    '{}.to.remap.bam'.format(prefix))
        to_remap_num = os.path.join(self.tempdir,
                                    '{}.to.remap.num.gz'.format(prefix))
        vcf_out = os.path.join(self.tempdir,
                               '{}_hets.vcf'.format(vcf_sample_name))

        # Make SNP directory needed for WASP.
        snp_directory = os.path.join(self.tempdir, 'snps')
        lines = ' \\\n\t'.join([
            'make_wasp_input',
            vcf_out,
            vcf_sample_name,
            snp_directory,
            bed,
            '-v ' + ' \\\n\t-v '.join(vcfs),
            '-g {}'.format(gatk_fai),
            '-b {}'.format(bcftools_path),
            '-t {}'.format(self.tempdir),
        ]) 
        if vcf_chrom_conv:
            lines += ' \\\n\t-c {}'.format(vcf_chrom_conv)
        lines += '\n\n'
        lines += ('{} view -b -q 255 -F 1024 \\\n\t{} '
                  '\\\n\t> {}\n\n'.format(
                    samtools_path, bam, uniq_bam))
        lines += 'wait\n\n'
        # Run WASP to swap alleles.
        lines += ('python {} -s -p \\\n\t{} \\\n\t{}\n\n'.format(
            find_intersecting_snps_path, uniq_bam, snp_directory))
        with open(self.filename, "a") as f:
            f.write(lines)
        return (snp_directory, vcf_out, keep_bam, wasp_r1_fastq, wasp_r2_fastq,
                to_remap_bam, to_remap_num)

    def wasp_alignment_compare(
        self,
        to_remap_bam, 
        to_remap_num,
        remapped_bam, 
        filter_remapped_reads_path, 
    ):
        """
        Check original mapping position of reads against remapping after
        swapping alleles using WASP, then count allele coverage for each
        variant.
    
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

        Returns
        -------
        wasp_filtered_bam : str
            Path to bam file after filtering against re-mapped reads.
    
        """
        # Files that will be created.
        wasp_filtered_bam = os.path.join(
            self.tempdir, '{}_filtered.bam'.format(self.sample_name))
        
        # Run WASP alignment compare.
        lines = ('python {} -p \\\n\t{} \\\n\t{} \\\n\t{} '
                 '\\\n\t{}\n\n'.format(
                     filter_remapped_reads_path, to_remap_bam,
                     remapped_bam, wasp_filtered_bam, to_remap_num))
        with open(self.filename, "a") as f:
            f.write(lines)
        return wasp_filtered_bam
        
    def count_allele_coverage(
        self,
        bam, 
        vcf,
        fasta, 
        gatk_path='$GATK',
    ):
        """
        Count alleles using GATK's ASEReadCounter.

        Parameters
        ----------
        bam : str
            Alignment file to count alleles from.

        vcf : str
            Path to VCF file with heterozygous SNVs.
    
        fasta : str
            Path to fasta file used to align data.

        Returns
        -------
        Path to output counts file.

        """
    
        # Count allele coverage.
        counts = os.path.join(
            self.tempdir, '{}_allele_counts.tsv'.format(self.sample_name))
        lines = ' \\\n\t'.join([
            'java -Xmx{}g -jar'.format(int(self.memory * 0.8)),
            '-XX:ParallelGCThreads=1',
            '-Djava.io.tmpdir={}'.format(self.tempdir), 
            '-jar {}'.format(gatk_path),
            '-R {}'.format(fasta),
            '-T ASEReadCounter',
            '-o {}'.format(counts),
            '-I {}'.format(bam),
            '-sites {}'.format(vcf),
            '-overlap COUNT_FRAGMENTS_REQUIRE_SAME_BASE',
            '-U ALLOW_N_CIGAR_READS',
        ])
        lines += '\n\nwait\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return counts

    def mbased(
        self,
        allele_counts, 
        feature_bed, 
        is_phased=False, 
        num_sim=1000000, 
        vcfs=None,
        vcf_sample_name=None, 
        vcf_chrom_conv=None,
        mappability=None,
        bigWigAverageOverBed_path='bigWigAverageOverBed',
    ):
        """
        Make a shell script for running MBASED to determine allelic bias from
        sequencing reads.
    
        Parameters
        ----------
        allele_counts : str
            Output file from GATK's ASEReadCounter.
    
        feature_bed : str
            Path to bed file for assigning heterozygous SNVs to features.
        
        is_phased : bool
            Whether the input file is phased. If so, the reference alleles are
            assumed to be in phase. Note that this only matters locus by locus.
    
        num_sim : int
            Number of simulations for MBASED to perform.
    
        vcfs : str
            List of paths to gzipped, indexed VCF file with all variant calls
            (not just heterozygous calls). These will be used to filter out
            heterozygous variants that are near other heterzygous or homozygous
            alternate variants that may cause mapping problems.
    
        vcf_sample_name : str
            If vcf is provided, this must be provided to specify the sample name
            of this sample in the VCF file. Required if vcf is provided.
    
        vcf_chrom_conv : str
            File with VCF chromosomes in first column and corresponding RNA-seq
            chromosomes in second column (no header). This is needed if the VCF
            and RNA-seq data have different chromosome naming.
    
        mappability : str
            Path to bigwig file with mappability scores. A score of one should
            mean uniquely mapping.
    
        bigWigAverageOverBed_path : str
            Path to bigWigAverageOverBed. Required if mappability is provided.
    
        Returns
        -------
        mbased_infile : str
            Path to save MBASED input file to.
    
        locus_outfile : str
            Path to file to store locus-level results.
    
        snv_outfile : str
            Path to file to store SNV-level results.
    
        """
        mbased_infile = os.path.join(
            self.tempdir, '{}_mbased_input.tsv'.format(self.sample_name))
        locus_outfile = os.path.join(
            self.tempdir, '{}_locus.tsv'.format(self.sample_name))
        snv_outfile = os.path.join(
            self.outdir, '{}_snv.tsv'.format(self.sample_name))
        is_phased = str(is_phased).upper()
        lines = 'make_mbased_input \\\n\t{} \\\n\t{} \\\n\t{}'.format(
            allele_counts, mbased_infile, feature_bed)
        if vcfs:
            lines += ' \\\n\t-v {} \\\n\t-s {}'.format(' \\\n\t-v '.join(vcfs), 
                                                       vcf_sample_name)
        if mappability:
            lines += ' \\\n\t-m {} \\\n\t-p {}'.format(
                mappability, bigWigAverageOverBed_path)
        if vcf_chrom_conv:
            lines += ' \\\n\t-c {}'.format(vcf_chrom_conv)
        lines += '\n\n'
        from __init__ import _scripts
        script = os.path.join(_scripts, 'mbased.R')
        lines += 'Rscript '
        lines += ' \\\n\t'.join([script, mbased_infile, locus_outfile,
                                 snv_outfile, self.sample_name, is_phased,
                                 str(num_sim), str(self.threads)])
        lines += '\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return mbased_infile, locus_outfile, snv_outfile
# 
# def merge_bams(
#     bams, 
#     outdir, 
#     tempdir,
#     merged_name, 
#     index=True,
#     bigwig=False,
#     copy_bams=True,
#     threads=8,
#     bedgraph_to_bigwig_path='bedGraphToBigWig',
#     bedtools_path='bedtools',
#     picard_path='$picard',
#     picard_memory=2,
# ):
#     """
#     Make a shell script for combining multiple bam files using Picard.
# 
#     Parameters
#     ----------
#     bams : list
#         List of SRA files to convert.
# 
#     outdir : str
#         Directory to store shell file and merged bam file.
# 
#     merged_name : str
#         Name used for output directory, files etc.
# 
#     index : bool
#         Whether to index the merged bam file.
# 
#     bigwig : bool
#         Whether to make bigwig file from merged bam file.
# 
#     copy_bams : bool
#         Whether to copy the input bam files to the temp directory. Not
#         necessary if temp directory is on the same file system as bam files.
# 
#     threads : int
#         Number of threads to request from SGE scheduler.
# 
#     Returns
#     -------
#     fn : str
#         Path to shell script.
# 
#     """
#     if bigwig:
#         index = True
#         assert bedgraph_to_bigwig_path
#         assert bedtools_path
# 
#     job_suffix = 'merged_bam'
#     job = JobScript(merged_name, job_suffix, outdir, threads, tempdir=tempdir,
#                     copy_input=False)
#     
#     # I'm going to define some file names used later.
#     merged_bam = job.add_temp_file('{}_merged.bam'.format(merged_name),
#                                    copy=True)
#     # merged_bam = os.path.join(tempdir,
#     #                             '{}_merged.bam'.format(merged_name))
#     # job.output_files_to_copy.append(merged_bam)
#     if index:
#         merged_bam_index = job.add_temp_file(
#             '{}_merged.bam.bai'.format(merged_name), copy=True)
#         # merged_bam_index = os.path.join(tempdir,
#         #                             '{}_merged.bam.bai'.format(merged_name))
#         # job.output_files_to_copy.append(merged_bam_index)
#     
#     if bigwig:
#         merged_bigwig = job.add_temp_file(
#             '{}_merged.bw'.format(merged_name), copy=True)
#     
#     temp_bams = []
#     for bam in bams:
#         temp_bams.append(job.add_input_file(bam))
#     # self.input_files_to_copy += bams
# 
#     if copy_bams:
#         job.copy_input_files()
# 
#     with open(job.filename, "a") as f:
#         lines = _picard_merge(temp_bams, merged_bam, picard_memory,
#                               picard_path=picard_path, picard_tempdir=tempdir)
#         f.write(lines)
#         
#         if index:
#             lines = _picard_index(merged_bam, merged_bam_index,
#                                   picard_path=picard_path,
#                                   picard_memory=picard_memory,
#                                   picard_tempdir=tempdir, bg=bigwig)
#             f.write(lines)
# 
#         if bigwig:
#             lines = _bigwig_files(
#                 merged_bam, merged_bigwig, merged_name,
#                 bedgraph_to_bigwig_path=bedgraph_to_bigwig_path,
#                 bedtools_path=bedtools_path)
#             f.write(lines)
# 
#     job.write_end()
#     return job.filename
