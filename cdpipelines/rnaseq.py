import os

from general import _make_dir
from general import JobScript

class RNAJobScript(JobScript):
    def star_align(
        self,
        r1_fastq, 
        r2_fastq, 
        rgpl, 
        rgpu, 
        star_index, 
        threads,
        genome_load='LoadAndRemove',
        transcriptome_align=True,
        star_path='STAR',
    ):
        """
        Align paired fastq files with STAR.
    
        Parameters
        ----------
        r1_fastq : str 
            Path to R1 fastq file.
    
        r2_fastq : str 
            Path to R2 fastq file.
    
        rgpl : str
            Read Group platform (e.g. illumina, solid). 
    
        rgpu : str
            Read Group platform unit (eg. run barcode). 

        Returns
        -------
        bam : str
            Path to output alignment bam file.
        
        log_out : str
            Path to log file.

        log_final_out : str
            Path to final log file.
        
        log_progress_out : str
        Path to progress log file.
        
        sj_out : str
            Path to output SJ.out.tab file.

        transcriptome_bam : str
            Path to output transcriptome alignment bam file. This is returned
            only if transcriptome_align == True.
    
        """
        lines = (' \\\n\t'.join([
            star_path, 
            '--runThreadN {}'.format(threads),
            '--genomeDir {}'.format(star_index), 
            '--genomeLoad {}'.format(genome_load),
            '--readFilesCommand zcat',
            '--readFilesIn {} {}'.format(r1_fastq, r2_fastq),
            '--outSAMattributes All', 
            '--outSAMunmapped Within',
            '--outSAMattrRGline ID:1 PL:{} PU:{} LB:{} SM:{}'.format(
                rgpl, rgpu, self.sample_name, self.sample_name),
            '--outFilterMultimapNmax 20', 
            '--outFilterMismatchNmax 999',
            '--alignIntronMin 20',
            '--alignIntronMax 1000000',
            '--alignMatesGapMax 1000000',
            '--outSAMtype BAM Unsorted']))
        if transcriptome_align:
            lines +=  ' \\\n\t--quantMode TranscriptomeSAM'
        lines += '\n\n'
        lines += 'if [ -d _STARtmp ] ; then rm -r _STARtmp ; fi\n\n'
        bam = os.path.join(
            self.tempdir, '{}.bam'.format(self.sample_name))
        log_out = os.path.join(
            self.tempdir, '{}_Log.out'.format(self.sample_name))
        log_final_out = os.path.join(
            self.tempdir, '{}_Log.final.out'.format(self.sample_name))
        log_progress_out = os.path.join(
            self.tempdir, '{}_Log.progress.out'.format(self.sample_name))
        sj_out = os.path.join(
            self.tempdir, '{}_SJ.out.tab'.format(self.sample_name))
        transcriptome_bam = os.path.join(
            self.tempdir, '{}_transcriptome.bam'.format(self.sample_name))
        lines += 'mv Aligned.out.bam {}\n'.format(bam)
        lines += 'mv Log.out {}\n'.format(log_out)
        lines += 'mv Log.final.out {}\n'.format(log_final_out)
        lines += 'mv Log.progress.out {}\n'.format(log_progress_out)
        lines += 'mv SJ.out.tab {}\n'.format(sj_out)
        if transcriptome_align:
            lines += 'mv Aligned.toTranscriptome.out.bam {}\n'.format(
                transcriptome_bam)
        lines += '\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        if transcriptome_align:
            return (bam, log_out, log_final_out, log_progress_out, sj_out,
                    transcriptome_bam)
        else:
            return bam, log_out, log_final_out, log_progress_out, sj_out

    def rsem_calculate_expression(
        self,
        bam, 
        reference, 
        threads=1, 
        calc_ci=False,
        ci_mem=1024, 
        strand_specific=True,
        rsem_calculate_expression_path='rsem-calculate-expression',
    ):
        """
        Estimate expression using RSEM.
    
        Parameters
        ----------
        bam : str
            Transcriptome bam file.
    
        reference : str
            RSEM reference.

        calc_ci : bool
            Whether to calculate confidence intervals.
    
        ci_mem : int
            Amount of memory in mb to give RSEM for calculating confidence
            intervals. Passed to --ci-memory for RSEM.
    
        strand_specific : boolean
            True if the data is strand-specific. False otherwise. For now, this
            means that the R1 read is on the reverse strand.
    
        Returns
        -------
        genes : str
            Path to genes output file.

        isoforms : str
            Path to isoforms output file.

        stats : str
            Path to output stats files.

        """
        genes = os.path.join(self.tempdir,
                             '{}.genes.results'.format(self.sample_name))
        isoforms = os.path.join(self.tempdir,
                                '{}.isoforms.results'.format(self.sample_name))
        stats = os.path.join(self.tempdir, '{}.stat'.format(self.sample_name))
        line = ('{} --bam --paired-end --num-threads {} '
                '--no-bam-output --seed 3272015 --calc-ci '
                '--ci-memory {} --estimate-rspd \\\n\t{} \\\n\t{} {}'.format(
                    rsem_calculate_expression_path, threads, ci_mem, bam,
                    reference, self.sample_name))
        lines = ('{} --bam --paired-end --num-threads {} \\\n\t--no-bam-output '
                 '--seed 3272015 --estimate-rspd \\\n\t'.format(
                     rsem_calculate_expression_path, threads))
        if calc_ci:
                lines += '--calc-ci --ci-memory {} \\\n\t'.format(ci_mem)
        if strand_specific:
            lines += '--forward-prob 0 \\\n\t'
        lines += '{} \\\n\t{} \\\n\t{}\n\n'.format(bam, reference,
                                                   self.sample_name)
        with open(self.filename, "a") as f:
            f.write(lines)
        return genes, isoforms, stats

    def dexseq_count(
        self,
        bam, 
        dexseq_annotation, 
        paired=True,
        strand_specific=True, 
        dexseq_count_path=None,
        samtools_path='samtools',
    ):
        """
        Count reads overlapping exonic bins for DEXSeq.
    
        Parameters
        ----------
        bam : str
            Path to coordinate sorted bam file to count reads for.
    
        dexseq_annotation : str
            Path to DEXSeq exonic bins GFF file.
    
        paired : boolean
            True if the data is paired-end. False otherwise.
    
        strand_specific : boolean
            True if the data is strand-specific. False otherwise.

        dexseq_count_path : str
            Path to dexseq_count.py script. If not provided, rpy2 will look for
            the path in R.
    
        Returns
        -------
        counts_file : str
            Path to file with bin counts.
    
        """
        counts_file = os.path.join(
            self.tempdir, '{}_dexseq_counts.tsv'.format(self.sample_name))
        if dexseq_count_path is None:
            import readline
            import rpy2.robjects as robjects
            robjects.r('suppressPackageStartupMessages(library(DEXSeq))')
            scripts = robjects.r('system.file("python_scripts", package="DEXSeq")')
            g = scripts.items()
            scripts_path = g.next()[1]
            dexseq_count_path = os.path.join(scripts_path, 'dexseq_count.py')
        if paired:
            p = 'yes'
        else:
            p = 'no'
        if strand_specific:
            s = 'reverse'
        else:
            s = 'no'
        lines = (
            '{} view -h -f 2 {} \\\n\t'.format(samtools_path, bam) +
            '| cut -f1-16,20- \\\n\t| python {} \\\n\t'.format(dexseq_count_path) + 
            '-p {} -s {} -a 0 -r pos -f sam \\\n\t'.format(p, s) + 
            '{} \\\n\t- {}\n\n'.format(dexseq_annotation, counts_file)
        )
        with open(self.filename, "a") as f:
            f.write(lines)
        return counts_file
    
    def htseq_count(
        self,
        bam, 
        gtf, 
        strand_specific=False,
        samtools_path='samtools',
    ):
        """
        Count reads overlapping genes for use with DESeq etc.
    
        Parameters
        ----------
        bam : str
            Path to coordinate sorted bam file to count reads for.
    
        gtf : str
            Path to GTF file to count against. Optimized for use with Gencode
            GTF.
    
        strand_specific : boolean
            True if the data is strand-specific. False otherwise.
    
        Returns
        -------
        lines : str
            Lines to be printed to shell script.
    
        name : str
            File name for the softlink.
   
        Returns
        -------
        counts_file : str
            Path to file with gene counts.
    
        stats_file : str
            Path to file with counting stats.
    
        """
        counts_file = os.path.join(
            self.tempdir, '{}_gene_counts.tsv'.format(self.sample_name))
        stats_file = os.path.join(
            self.tempdir, '{}_gene_count_stats.tsv'.format(self.sample_name))
        import HTSeq
        if strand_specific:
            s = 'reverse'
        else:
            s = 'no'
        script = os.path.join(HTSeq.__path__[0], 'scripts', 'count.py')
        lines = ('python {} \\\n\t-f bam -r pos -s {} '.format(script, s) + 
                 '-a 0 -t exon -i gene_id -m union \\\n\t' + 
                 '{} \\\n\t{} \\\n\t> temp_out.tsv\n'.format(bam, gtf))
        lines += 'tail -n 5 temp_out.tsv > {}\n'.format(stats_file)
        lines += 'lines=$(wc -l <temp_out.tsv)\n'
        lines += 'wanted=`expr $lines - 5`\n'
        lines += 'head -n $wanted temp_out.tsv > {}\n'.format(counts_file)
        lines += 'rm temp_out.tsv\n\n'
        with open(self.filename, "a") as f:
            f.write(lines)
        return counts_file, stats_file

    def bedgraph_from_bam(
        self,
        bam, 
        strand=None,
        scale=None,
        bedtools_path='bedtools',
        sambamba_path='sambamba',
    ):
        """
        Make lines that create a coverage bedgraph file.
    
        Parameters
        ----------
        bam : str
            Bam file to calculate coverage for.
    
        strand : str
            If '+' or '-', calculate strand-specific coverage. Otherwise,
            calculate coverage using all reads.
    
        scale : float
            Scale the bigwig by this amount.

        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        bedgraph : str
            Path to output bedgraph file.
    
        """
        fn_root = self.sample_name
        if strand == '+':
            fn_root += '_plus'
        elif strand == '-':
            fn_root += '_minus'
        if scale:
            fn_root += '_scaled'.format(scale)
        bedgraph = os.path.join(self.tempdir, '{}.bg'.format(fn_root))

        if bedtools_path == 'bedtools':
            genome_file = 'human.hg19.genome'
        else:
            genome_file = os.path.join(
                os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
                'human.hg19.genome')

        if strand == '+':
            lines = (
                '{} view -f bam -F "((first_of_pair and mate_is_reverse_strand) '
                'or (second_of_pair and reverse_strand)) and mapping_quality >= '
                '255" \\\n\t{} \\\n\t| {} '
                'genomecov -ibam stdin -g {} \\\n\t-split -bg -trackline '
                '-trackopts \'name="{}"\' '.format(
                    sambamba_path, bam, bedtools_path, genome_file, fn_root))
        elif strand == '-':
            lines = (
                '{} view -f bam -F "((second_of_pair and mate_is_reverse_strand) '
                'or (first_of_pair and reverse_strand)) and mapping_quality >= '
                '255" \\\n\t{} \\\n\t| {} genomecov '
                '-ibam stdin -g {} \\\n\t-split -bg -trackline -trackopts '
                '\'name="{}"\' '.format(
                    sambamba_path, bam, bedtools_path, genome_file, fn_root))
        else:
            lines = ('{} view -f bam -F "not (unmapped or mate_is_unmapped) and '
                     'mapping_quality >= 255" \\\n\t{} \\\n\t| '
                     '{} genomecov -ibam stdin \\\n\t-g {} -split -bg \\\n\t'
                     '-trackline -trackopts \'name="{}"\' '.format(
                         sambamba_path, bam, bedtools_path, genome_file,
                         self.sample_name))

        if scale:
            lines += ' -scale {}'.format(scale)
        
        lines += ' \\\n\t> {}\n\n'.format( bedgraph)

        with open(self.filename, "a") as f:
            f.write(lines)
        return bedgraph
    
    def bigwig_from_bedgraph(
        self,
        bedgraph,
        strand=None,
        scale=False,
        web_available=True,
        write_to_outdir=False,
        bedGraphToBigWig_path='bedGraphToBigWig',
        bedtools_path='bedtools',
    ):
        """
        Make bigwig coverage file from bam file.
    
        Parameters
        ----------
        bedgraph : str
            Path to bedgraph file to create bigwig for.
        
        strand : str
            If '+' or '-', add this information to trackline.
    
        scale : bool
            If True, add note to trackline that data is scaled.

        web_available : bool
            If True, write trackline to self.links_tracklines, make softlink to
            self.linkdir, and set write_to_outdir = True.

        write_to_outdir : bool
            If True, write output files directly to self.outdir.
    
        bedGraphToBigWig_path : str
            Path to bedGraphToBigWig executable.

        bedtools_path : str
            Path to bedtools. If bedtools_path == 'bedtools', it is assumed that
            the hg19 human.hg19.genome file from bedtools is also in your path.

        Returns
        -------
        bigwig : str
            Path to output bigwig file.
    
        """
        if write_to_outdir or web_available:
            dy = self.outdir
        else:
            dy = self.tempdir
        if bedtools_path == 'bedtools':
            genome_file = '$(which human.hg19.genome)'
        else:
            genome_file = os.path.join(
                os.path.split(os.path.split(bedtools_path)[0])[0], 'genomes',
                'human.hg19.genome')
        root = os.path.splitext(os.path.split(bedgraph)[1])[0]
        bigwig = os.path.join(dy, '{}.bw'.format(root))
        lines = '{} \\\n\t{} \\\n\t{} \\\n\t{}\n\n'.format(
            bedGraphToBigWig_path, bedgraph, genome_file, bigwig)
        with open(self.filename, "a") as f:
            f.write(lines)
        if web_available:
            name = '{}_rna'.format(self.sample_name)
            desc = 'RNAseq coverage for {}.'.format(self.sample_name)
            if strand == '+':
                name += '_plus'
                desc = desc.replace('RNAseq coverage', 
                                    'RNAseq plus strand coverage')
            elif strand == '-':
                name += '_minus'
                desc = desc.replace('RNAseq coverage', 
                                    'RNAseq minus strand coverage')
            if scale:
                name += '_scaled'
                desc = desc.replace('RNAseq', 'Scaled RNAseq')

            url = self.webpath + '/' + os.path.split(bigwig)[1]
            t_lines = (
                'track type=bigWig name="{}" '
                'description="{}" '
                'visibility=0 db=hg19 bigDataUrl={}\n'.format(
                    name, desc, url))
            with open(self.links_tracklines, "a") as f:
                f.write(t_lines)
            link = self.add_softlink(bigwig)

        return bigwig
     
    def bigwig_hub(
        self,
        plus_bw,
        minus_bw,
        scale=False,
    ):
        """
        Make a bigwig hub for the plus and minus strand bigwigs. Hub directory
        is written directly to self.outdir and a link is made to linkdir. The
        URL of the hub is written to self.links_tracklines.
    
        Parameters
        ----------
        plus_bw : str
            Path to output plus strand bigwig file.
        
        minus_bw : str
            Path to output minus strand bigwig file.
        
        scale : bool
            If True, add labels indicated scaled.

        """
        dy = os.path.join(self.outdir, '{}_rna_hub'.format(self.sample_name))
        if scale:
            dy += '_scaled'
        _make_dir(dy)
        # Make genomes.txt file.
        with open(os.path.join(dy, 'genomes.txt'), 'w') as f:
            f.write('genome hg19\ntrackDb hg19/trackDb.txt\n')
        # Make hub.txt file.
        lines = '\n'.join([
            'hub {}_rna_cov'.format(self.sample_name),
            'shortLabel RNA cov {}'.format(self.sample_name.split('-')[0]),
            ('longLabel RNA-seq coverage for plus and minus strands '
             'for {}.'.format(self.sample_name)),
            'genomesFile genomes.txt',
            'email cdeboeve@ucsd.edu',
        ]) + '\n'
        with open(os.path.join(dy, 'hub.txt'), 'w') as f:
            f.write(lines)
        # Make trackDb.txt file.
        _make_dir(os.path.join(dy, 'hg19'))
        lines = '\n'.join([
            'track {}_rna_cov'.format(self.sample_name),
            'type bigWig',
            'container multiWig',
            'shortLabel RNA cov {}'.format(self.sample_name.split('-')[0]),
            ('longLabel RNA-seq coverage for plus and minus strands '
             'for {}.'.format(self.sample_name)),
            'visibility full',
            'aggregate transparentOverlay',
            'showSubtrackColorOnUi on',
            'viewLimits 1:400',
        ]) + '\n\n'

        url = self.webpath + '/' + os.path.split(plus_bw)[1]
        lines += '\t' + '\n\t'.join([
            'track plus',
            'bigDataUrl {}'.format(url),
            'shortLabel Plus cov {}'.format(self.sample_name.split('-')[0]),
            ('longLabel Coverage for {} from fragments originating from '
             'the plus strand.'.format(self.sample_name)),
            'parent {}_rna_cov'.format(self.sample_name),
            'type bigWig',
            'color 255,0,0',
        ]) + '\n\n'

        url = self.webpath + '/' + os.path.split(minus_bw)[1]
        lines += '\t' + '\n\t'.join([
            'track minus',
            'bigDataUrl {}'.format(url),
            'shortLabel Minus cov {}'.format(self.sample_name.split('-')[0]),
            ('longLabel Coverage for {} from fragments originating from '
             'the minus strand.'.format(self.sample_name)),
            'parent {}_rna_cov'.format(self.sample_name),
            'type bigWig',
            'color 0,0,255',
        ]) + '\n'
        if scale:
            lines = lines.replace('_rna_cov', '_rna_cov_scaled')
            lines = lines.replace('RNA cov', 'RNA scov')
            lines = lines.replace('RNA-seq coverage', 'Scaled RNA-seq coverage')
            lines = lines.replace(' Coverage ', ' Scaled coverage ')
            lines = lines.replace('us cov', 'us scov')

        _make_dir(os.path.join(dy, 'hg19'))
        with open(os.path.join(dy, 'hg19', 'trackDb.txt'), 'w') as f:
            f.write(lines)

        self.add_softlink(dy)
        url = self.webpath + '/' + os.path.split(dy)[1] + '/' + 'hub.txt'
        with open(self.links_tracklines, "a") as f:
            f.write(url + '\n')

def pipeline(
    r1_fastqs, 
    r2_fastqs, 
    outdir, 
    sample_name, 
    star_index,
    ref_flat, 
    rrna_intervals,
    dexseq_annotation,
    gene_gtf,
    gene_bed,
    exon_bed,
    rsem_reference,
    sra_files=None,
    find_intersecting_snps_path=None,
    filter_remapped_reads_path=None,
    gatk_fasta=None,
    linkdir=None,
    webpath_file=None,
    vcfs=None,
    vcf_sample_name=None,
    vcf_chrom_conv=None,
    is_phased=False,
    conda_env=None,
    modules=None,
    queue='frazer',
    star_genome_load='LoadAndRemove',
    rgpl='ILLUMINA',
    rgpu='',
    strand_specific=True, 
    tempdir=None,
    mappability=None,
    expected_unique_pairs=20000000,
    dexseq_count_path=None,
    star_path='STAR',
    picard_path='$picard',
    bedtools_path='bedtools',
    bedGraphToBigWig_path='bedGraphToBigWig',
    fastqc_path='fastqc',
    samtools_path='samtools',
    sambamba_path='sambamba',
    rsem_calculate_expression_path='rsem-calculate-expression',
    gatk_path='$GATK',
    bigWigAverageOverBed_path='bigWigAverageOverBed',
    bcftools_path='bcftools',
    bammarkduplicates_path='bammarkduplicates',
    featureCounts_path='featureCounts',
    fastq_dump_path='fastq-dump',
):
    """
    Make SGE/shell scripts for running the entire RNA-seq pipeline. The defaults
    are set for use on the Frazer lab's SGE scheduler on flh1/flh2.

    Parameters
    ----------
    r1_fastqs : list or str
        Either a list of paths to gzipped fastq files with R1 reads or path to a
        single gzipped fastq file with R1 reads. If you want to process SRA
        files, pass None here.

    r2_fastqs : list or str
        Either a list of paths to gzipped fastq files with R2 reads or path to a
        single gzipped fastq file with R2 reads. If you want to process SRA
        files, pass None here.

    outdir : str
        Directory to store shell scripts, stdout/stderr logs, and output files
        and directories.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

    ref_flat : str
        Path to refFlat file with non-rRNA genes. Can ge gzipped.

    rrna_intervals : str
        Path to interval list file with rRNA intervals.

    dexseq_annotation : str
        Path to DEXSeq exonic bins GFF file.

    gene_gtf : str
        Path to GTF file with gene annotations.

    gene_bed : str
        Path to bed file with gene definitions. The gene ID is used for
        assigning variants to genes for ASE.

    exon_bed : str
        Path to bed file with exon definitions. This is used for filtering
        variants used for ASE.

    rsem_reference : str
        Directory with RSEM reference.

    sra_files : list
        List of SRA file paths or URLs. If r1_fastqs and r2_fastqs are None, the
        pipeline will run for these SRA files (concatenating all files into one
        sample).

    find_intersecting_snps_path : str
        Path to find_intersecting_snps.py from WASP.
    
    filter_remapped_reads_path : str
        Path to filter_remapped_reads.py from WASP.

    gatk_fasta : str
        Fasta file that corresponds to the fasta file used for STAR but is
        sorted karyotypically for GATK. Assumed to have associated dict and fai
        files. This is needed for ASE.

    linkdir : str
        Path to directory where softlinks should be made. Some pipeline parts
        may make softlinks output files here for display on the web.

    webpath_file : str
        File whose first line is the URL that points to linkdir. For example,
        if we make a link to the file s1_coord_sorted.bam in linkdir and
        webpath_file has http://site.com/files on its first line, then
        http://site.com/files/s1_coord_sorted.bam should be available on the
        web. If the web directory is password protected (it probably should be),
        then the URL should look like http://username:password@site.com/files.
        This is a file so you don't have to make the username/password combo
        public (although I'd recommend not using a sensitive password). You can
        just put the webpath_file in a directory that isn't tracked by git, 
        figshare, etc.

    vcfs : list
        List of VCF files containing exonic variants used for ASE. The VCF files
        will be concatenated so they shouldn't overlap in sites (e.g. they
        should be different chromosomes, genomic regions, etc.).
    
    vcf_sample_name : str
        Sample name of this sample in the VCF files (if different than
        sample_name). For instance, the sample name in the VCF file may be the
        sample name for WGS data which may differ from the RNA-seq sample name.

    vcf_chrom_conv : str
        File with VCF chromosomes in first column and corresponding RNA-seq
        chromosomes in second column (no header). This is needed if the VCF and
        RNA-seq data have different chromosome naming.

    conda_env : str
        Conda environment to load at the beginning of the script.

    modules : str
        Comma-separated list of modules to load at the beginning of the script.

    queue : str
        Name of queue to submit jobs to. The default value of "frazer" will
        submit different jobs to the different queues on the Frazer lab cluster
        in an optimal way.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    strand_specific : boolean
        If false, data is not strand specific.

    tempdir : str
        Directory to store temporary files.

    expected_unique_pairs : int
        Number of expected mapped READ PAIRS. This only needs to be a rough 
        estimate. You should keep this number constant for all samples across 
        an experiment/project (i.e. any samples you want to compare against 
        each other). For example, say you are doing RNA-seq for 10 samples on 
        a lane, the lane yields 100,000,000 pairs of reads, and you expect 
        90%% of the read pairs to align uniquely. Then you would want 
        expected_num = 9,000,000. In this case, you would set actual_num to 
        the actual number of uniquely mapped reads. The reason you want this 
        number close to the number of expected number of uniquely mapped 
        pairs is so the normalization does not drastically change the 
        coverage for samples near the expected number of input read pairs. 
        For instance, if expected_num=20M and a sample has exactly 20M 
        uniquely mapped read pairs, then the coverage will not be changed. If 
        expected_num=20M and a sample has only 10M input read PAIRS, it will 
        be normalized so that it looks like it had 20M input read pairs (i.e. 
        all coverages will be multiplied by 2).

    star_path : str
        Path to STAR aligner.

    picard_path : str
        Path to Picard tools.

    bedtools_path : str
        Path to bedtools.

    bedGraphToBigWig_path : str
        Path bedGraphToBigWig executable.

    dexseq_count_path : str
        Path to dexseq_count.py script. If not provided, rpy2 will look for the
        path in R.
    
    Returns
    -------
    fn : str
        Path to submission shell script.

    """
    with open(webpath_file) as wpf:
        webpath = wpf.readline().strip()

    # Bash commands to submit jobs. I'll collect these as I make the jobs and
    # then write them to a file at the end.
    submit_commands = []

    # Set default queue for Frazer lab settings. None will just go to the
    # default queue. For jobs that need a specific queue, I'll set the queue
    # below.
    if queue == 'frazer':
        default_queue = None
    else:
        default_queue = queue

    ##### Job 1: Combine fastqs and align with STAR. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix='alignment',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=8, 
        memory=32,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue, 
        conda_env=conda_env,
        modules=modules,
    )
    alignment_jobname = job.jobname
    
    # If needed, convert SRA files to fastq files. 
    # TODO: Eventually, I can also take bam files as input and convert them to
    # fastq if needed.
    if r1_fastqs is None and r1_fastqs is None and sra_files:
        r1,r2 = job.convert_sra_to_fastq(sra_files,
                                         fastq_dump_path=fastq_dump_path)
        r1_fastqs = [r1]
        r2_fastqs = [r2]
    
    # Input files.
    for fq in r1_fastqs + r2_fastqs:
        job.add_input_file(fq)

    # Combine R1 and R2 fastqs.
    if type(r1_fastqs) == str:
        r1_fastqs = [r1_fastqs]
    if type(r2_fastqs) == str:
        r2_fastqs = [r2_fastqs]
    r1_fastqs = [os.path.realpath(x) for x in r1_fastqs]
    r2_fastqs = [os.path.realpath(x) for x in r2_fastqs]
    combined_r1 = job.combine_fastqs(r1_fastqs, suffix='R1', bg=True)
    combined_r2 = job.combine_fastqs(r2_fastqs, suffix='R2', bg=True)
    with open(job.filename, "a") as f:
            f.write('\nwait\n\n')
    # We don't want to keep the fastqs indefinitely, but we need them for the
    # fastQC step later.
    combined_r1 = job.add_output_file(combined_r1)
    combined_r2 = job.add_output_file(combined_r2)

    # Align reads.
    (star_bam, log_out, log_final_out, log_progress_out, sj_out, 
     transcriptome_bam) = \
            job.star_align(combined_r1, combined_r2, rgpl, rgpu, star_index,
                            job.threads, genome_load=star_genome_load)
    star_bam = job.add_output_file(star_bam)
    transcriptome_bam = job.add_output_file(transcriptome_bam)
    log_final_out = job.add_output_file(log_final_out)
    [job.add_output_file(x) for x in [log_out, log_progress_out, sj_out]]
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ##### Job 2: Run fastQC. ##### 
    job = RNAJobScript(
        sample_name, 
        job_suffix='fastqc', 
        outdir=os.path.join(outdir, 'qc'), 
        threads=1, 
        memory=4,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue, 
        conda_env=conda_env,
        modules=modules, 
        wait_for=[alignment_jobname]
    )
    fastqc_jobname = job.jobname
 
    # Input files.
    job.add_input_file(combined_r1, delete_original=True)
    job.add_input_file(combined_r2, delete_original=True)

    # Run fastQC.
    fastqc_html, fastqc_zip = job.fastqc([combined_r1, combined_r2],
                                         fastqc_path)
    [job.add_output_file(x) for x in fastqc_html + fastqc_zip]
        
    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ##### Job 3: Coordinate sort, mark duplicates and index bam. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix = 'sort_mdup_index',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=4, 
        memory=6,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[alignment_jobname]
    )
    sort_mdup_index_jobname = job.jobname

    # Input files.
    star_bam = job.add_input_file(star_bam, delete_original=True)

    # Coordinate sort.
    coord_sorted_bam = job.sambamba_sort(
        star_bam, 
        tempdir=job.tempdir,
        sambamba_path=sambamba_path,
    )
    job.add_temp_file(coord_sorted_bam)

    # Index sorted bam file.
    temp_bam_index = job.sambamba_index(coord_sorted_bam, sambamba_path)
    job.add_temp_file(temp_bam_index)

    # Mark duplicates.
    mdup_bam, duplicate_metrics = job.biobambam2_mark_duplicates(
        coord_sorted_bam,
        bammarkduplicates_path=bammarkduplicates_path)
    outdir_mdup_bam = job.add_output_file(mdup_bam)
    job.add_output_file(duplicate_metrics)

    # Add softlink to bam file in outdir and write URL and trackline.
    link = job.add_softlink(outdir_mdup_bam)
    name = '{}_rna'.format(job.sample_name)
    desc = 'RNAseq alignment for {}.'.format(job.sample_name)
    url = job.webpath + '/' + os.path.split(outdir_mdup_bam)[1]
    t_lines = (
        'track type=bam name="{}" '
        'description="{}" '
        'visibility=0 db=hg19 bigDataUrl={}\n'.format(
            name, desc, url))
    with open(job.links_tracklines, "a") as f:
        f.write(t_lines)
        f.write(url + '\n')

    # Index bam file.
    bam_index = job.sambamba_index(mdup_bam, sambamba_path)
    outdir_bam_index = job.add_output_file(bam_index)
    link = job.add_softlink(outdir_bam_index)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
    
    # These files will be input for upcoming scripts.
    mdup_bam = outdir_mdup_bam
    bam_index = outdir_bam_index

    ##### Job 4: Collect Picard metrics. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix='picard_metrics',
        outdir=os.path.join(outdir, 'qc'),
        threads=1, 
        memory=7, 
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue,
        conda_env=conda_env, 
        modules=modules,
        wait_for=[sort_mdup_index_jobname],
    )
    picard_metrics_jobname = job.jobname
    
    # Input files.
    mdup_bam = job.add_input_file(mdup_bam)

    # Collect several different Picard metrics including insert size.
    metrics_files = job.picard_collect_multiple_metrics(
        mdup_bam, 
        picard_path=picard_path, 
        bg=False,
    )
    for fn in metrics_files:
        job.add_output_file(fn)

    # Collect RNA seq metrics.
    metrics, chart = job.picard_collect_rna_seq_metrics(
        mdup_bam, 
        ref_flat, 
        rrna_intervals,
        picard_path=picard_path,
        strand_specific=strand_specific, 
        bg=False,
    )
    job.add_output_file(metrics)
    job.add_output_file(chart)

    # Collect index stats.
    index_out, index_err = job.picard_bam_index_stats(
        mdup_bam, 
        picard_path=picard_path,
        bg=False,
    )
    job.add_output_file(index_out)
    job.add_output_file(index_err)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())

    ##### Job 5: Make md5 hash for final bam file. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix='md5',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=1, 
        memory=1,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[sort_mdup_index_jobname],
    )
    md5_jobname = job.jobname
    
    # Input files.
    mdup_bam = job.add_input_file(outdir_mdup_bam)

    # Make md5 hash for output bam file.
    md5sum = job.make_md5sum(mdup_bam)
    job.add_output_file(md5sum)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
       
    ##### Job 6: Make bigwig files for final bam file. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix='bigwig',
        outdir=os.path.join(outdir, 'alignment'), 
        threads=1, 
        memory=4,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[sort_mdup_index_jobname],
    )
    bigwig_jobname = job.jobname
        
    # Input files.
    mdup_bam = job.add_input_file(mdup_bam)

    # First make bigwig from both strands.
    bg = job.bedgraph_from_bam(
        mdup_bam, 
        bedtools_path=bedtools_path,
        sambamba_path=sambamba_path,
    )
    job.add_temp_file(bg)
    bw = job.bigwig_from_bedgraph(
        bg,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    job.add_output_file(bw)

    # Now for genes on the plus strand.
    plus_bg = job.bedgraph_from_bam(
        mdup_bam, 
        strand='+',
        bedtools_path=bedtools_path,
        sambamba_path=sambamba_path,
    )
    job.add_temp_file(plus_bg)
    plus_bw = job.bigwig_from_bedgraph(
        plus_bg,
        strand='+',
        scale=False,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    plus_bw = job.add_output_file(plus_bw)

    # Now for genes on the minus strand.
    minus_bg = job.bedgraph_from_bam(
        mdup_bam, 
        strand='-',
        bedtools_path=bedtools_path,
        sambamba_path=sambamba_path,
    )
    job.add_temp_file(minus_bg)
    minus_bw = job.bigwig_from_bedgraph(
        minus_bg,
        strand='-',
        scale=None,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    minus_bw = job.add_output_file(minus_bw)
    # I'll make a track hub for the plus/minus strand bigwigs.
    job.bigwig_hub(
        plus_bw,
        minus_bw,
        scale=False,
    )

    # Now we'll make scaled versions. First I'll read the star Log.final.out
    # file and find the number of uniquely mapped reads.
    # Both strands scaled.
    scaled_bg = job.scale_bedgraph(
        bg,
        mdup_bam,
        expected_unique_pairs,
    )
    job.add_temp_file(scaled_bg)
    scaled_bw = job.bigwig_from_bedgraph(
        scaled_bg,
        scale=True,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    job.add_output_file(scaled_bw)

    # Plus strand scaled. Note that I divide by two because we expected about
    # half of the reads to map to each strand.
    plus_scaled_bg = job.scale_bedgraph(
        plus_bg,
        mdup_bam,
        expected_unique_pairs / 2,
    )
    job.add_temp_file(plus_scaled_bg)
    plus_scaled_bw = job.bigwig_from_bedgraph(
        plus_scaled_bg,
        strand='+',
        scale=True,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    plus_scaled_bw = job.add_output_file(plus_scaled_bw)

    # Minus strand scaled. Note that I divide by two because we expected about
    # half of the reads to map to each strand.
    minus_scaled_bg = job.scale_bedgraph(
        minus_bg,
        mdup_bam,
        expected_unique_pairs / 2,
    )
    job.add_temp_file(minus_scaled_bg)
    minus_scaled_bw = job.bigwig_from_bedgraph(
        minus_scaled_bg,
        strand='-',
        scale=True,
        bedGraphToBigWig_path=bedGraphToBigWig_path,
        bedtools_path=bedtools_path,
    )
    minus_scaled_bw = job.add_output_file(minus_scaled_bw)
    # I'll make a track hub for the plus/minus strand bigwigs.
    job.bigwig_hub(
        plus_scaled_bw,
        minus_scaled_bw,
        scale=True,
    )

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
    
    ##### Job 7: Get featureCounts and DEXSeq counts. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix = 'counts',
        outdir=os.path.join(outdir, 'counts'), 
        threads=1, 
        memory=4,
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, 
        queue=default_queue, 
        conda_env=conda_env,
        modules=modules,
        wait_for=[sort_mdup_index_jobname],
    )
    counts_jobname = job.jobname
    
    # Input files.
    mdup_bam = job.add_input_file(mdup_bam)

    # Get gene counts.
    if strand_specific:
        ss = 2
    else:
        ss = 0
    gene_counts, gene_count_stats = job.featureCounts_count(
        gene_gtf,
        mdup_bam,
        both=True,
        strand_specific=ss,
        featureCounts_path=featureCounts_path,
    )
    job.add_output_file(gene_counts)
    job.add_output_file(gene_count_stats)

    # Get DEXSeq bin counts.
    dexseq_counts = job.dexseq_count(
        mdup_bam, 
        dexseq_annotation,
        paired=True, 
        strand_specific=strand_specific,
        dexseq_count_path=dexseq_count_path,
        samtools_path=samtools_path,
    )
    job.add_output_file(dexseq_counts)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
    
    ##### Job 8: Run RSEM. #####
    job = RNAJobScript(
        sample_name, 
        job_suffix = 'rsem',
        outdir=os.path.join(outdir, 'rsem'),
        threads=8, 
        memory=32, 
        linkdir=linkdir,
        webpath=webpath,
        tempdir=tempdir, queue=default_queue,
        conda_env=conda_env, 
        modules=modules,
        wait_for=[sort_mdup_index_jobname],
    )
    rsem_jobname = job.jobname
    
    # Input files.
    transcriptome_bam = job.add_input_file(transcriptome_bam)

    # Run RSEM.
    genes, isoforms, stats = job.rsem_calculate_expression(
        transcriptome_bam, 
        rsem_reference, 
        threads=job.threads, 
        ci_mem=1024,
        strand_specific=strand_specific,
        rsem_calculate_expression_path=rsem_calculate_expression_path,
    )
    job.add_output_file(genes)
    job.add_output_file(isoforms)
    job.add_output_file(stats)

    job.write_end()
    if not job.delete_sh:
        submit_commands.append(job.sge_submit_command())
   
    # We'll only go through the ASE steps if a VCF was provided.
    if vcfs:
        ##### Job 9: WASP first step. #####
        job = RNAJobScript(
            sample_name, 
            job_suffix = 'wasp_allele_swap',
            outdir=os.path.join(outdir, 'wasp'),
            threads=1, 
            memory=4, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir, 
            queue=default_queue,
            conda_env=conda_env, 
            modules=modules,
            wait_for=[sort_mdup_index_jobname],
        )
        wasp_allele_swap_jobname = job.jobname
           
        # Input files.
        mdup_bam = job.add_input_file(mdup_bam)
        # The VCFs might be large so we probably don't want to copy it ever.
        input_vcfs = []
        for vcf in vcfs:
            input_vcfs.append(job.add_input_file(vcf, copy=False))
        # The exon bed file is small so we don't need to copy it ever.
        exon_bed = job.add_input_file(exon_bed, copy=False)

        # Run WASP allele swap.
        if not vcf_sample_name:
            vcf_sample_name = sample_name
        (snp_directory, hets_vcf, keep_bam, wasp_r1_fastq, wasp_r2_fastq,
         to_remap_bam, to_remap_num) = job.wasp_allele_swap(
             mdup_bam, 
             find_intersecting_snps_path, 
             input_vcfs, 
             exon_bed,
             gatk_fai=gatk_fasta + '.fai',
             vcf_sample_name=vcf_sample_name, 
             vcf_chrom_conv=vcf_chrom_conv,
             samtools_path=samtools_path,
             bcftools_path=bcftools_path,
        )
        # WASP outputs a file (keep_bam) that has reads that don't overlap
        # variants. I'm going to discard that file. I'll discard the WASP SNP
        # directory as well since I have the hets in a VCF file.
        job.add_temp_file(keep_bam)
        snp_directory = job.add_temp_file(snp_directory)
        hets_vcf = job.add_output_file(hets_vcf)
        wasp_r1_fastq = job.add_output_file(wasp_r1_fastq)
        wasp_r2_fastq = job.add_output_file(wasp_r2_fastq)
        to_remap_bam = job.add_output_file(to_remap_bam)
        to_remap_num = job.add_output_file(to_remap_num)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())
        
        ##### Job 10: WASP second step. #####
        job = RNAJobScript(
            sample_name, 
            job_suffix='wasp_remap',
            outdir=os.path.join(outdir, 'wasp'),
            threads=8, 
            memory=32, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir,
            queue=default_queue, 
            conda_env=conda_env, 
            modules=modules,
            wait_for=[wasp_allele_swap_jobname],
        )
        wasp_remap_jobname = job.jobname
        
        # Input files.
        wasp_r1_fastq = job.add_input_file(wasp_r1_fastq, delete_original=True)
        wasp_r2_fastq = job.add_input_file(wasp_r2_fastq, delete_original=True)

        # Realign allele-swapped fastqs.
        remapped_bam, log_out, log_final_out, log_progress_out, sj_out = \
                job.star_align(wasp_r1_fastq, wasp_r2_fastq, rgpl, rgpu,
                               star_index, job.threads,
                               genome_load=star_genome_load,
                               transcriptome_align=False)
        job.add_output_file(remapped_bam)
        job.add_output_file(log_out)
        job.add_output_file(log_final_out)
        job.add_output_file(log_progress_out)
        job.add_temp_file(sj_out)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())
        
        ##### Job 11: WASP third step. #####
        job = RNAJobScript(
            sample_name, 
            job_suffix = 'wasp_alignment_compare',
            outdir=os.path.join(outdir, 'wasp'),
            threads=1, 
            memory=5, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir, 
            queue=default_queue,
            conda_env=conda_env, 
            modules=modules,
            wait_for=[wasp_remap_jobname],
        )
        wasp_alignment_compare_jobname = job.jobname

        # Input files.
        to_remap_bam = job.add_input_file(to_remap_bam, delete_original=True)
        to_remap_num = job.add_input_file(to_remap_num, delete_original=True)
        remapped_bam = job.add_input_file(remapped_bam, delete_original=True)

        # Compare alignments.
        temp_filtered_bam = job.wasp_alignment_compare(
            to_remap_bam, 
            to_remap_num,
            remapped_bam, 
            filter_remapped_reads_path,
        )
        job.add_temp_file(temp_filtered_bam)
            
        # Coordinate sort and index filtered bam file.
        wasp_filtered_bam, wasp_bam_index = job.picard_coord_sort(
            temp_filtered_bam, 
            index=True,
            picard_path=picard_path,
        )
        # I'll keep this bam file and its index as a record of which reads were
        # used to calculate ASE. It might be useful for visualization. The bam
        # file is pretty small anyway.
        job.add_output_file(wasp_filtered_bam)
        job.add_output_file(wasp_bam_index)

        # Reorder bam file so it will work with GATK.
        reordered_bam = job.picard_reorder(
            wasp_filtered_bam, 
            fasta=gatk_fasta,
            picard_path=picard_path,
        )
        job.add_temp_file(reordered_bam)

        # Index reordered bam.
        reordered_index = job.picard_index(
            reordered_bam, 
            picard_path=picard_path,
            bg=False,
        )
        job.add_temp_file(reordered_index)

        # Get allele counts.
        allele_counts = job.count_allele_coverage(
            reordered_bam, 
            hets_vcf,
            gatk_fasta, 
            gatk_path=gatk_path,
        )
        job.add_output_file(allele_counts)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())
        
        ##### Job 12: Run MBASED for ASE. #####
        if queue == 'frazer':
            q = 'opt'
        else:
            q = queue
        job = RNAJobScript(
            sample_name,
            job_suffix='mbased',
            outdir=os.path.join(outdir, 'mbased'),
            threads=16, 
            memory=32, 
            linkdir=linkdir,
            webpath=webpath,
            tempdir=tempdir, 
            queue=q,
            conda_env=conda_env, 
            modules=modules,
            wait_for=[wasp_alignment_compare_jobname],
        )
        mbased_jobname = job.jobname
    
        # Input files.
        allele_counts = job.add_input_file(allele_counts)

        mbased_infile, locus_outfile, snv_outfile = job.mbased(
            allele_counts, 
            gene_bed, 
            is_phased=is_phased, 
            num_sim=1000000, 
            vcfs=input_vcfs,
            vcf_sample_name=vcf_sample_name, 
            vcf_chrom_conv=vcf_chrom_conv,
            mappability=mappability,
            bigWigAverageOverBed_path=bigWigAverageOverBed_path,
        )
        job.add_output_file(mbased_infile)
        job.add_output_file(locus_outfile)
        job.add_output_file(snv_outfile)

        job.write_end()
        if not job.delete_sh:
            submit_commands.append(job.sge_submit_command())

    ##### Submission script #####
    # Now we'll make a submission script that submits the jobs with the
    # appropriate dependencies.
    import datetime as dt
    now = str(dt.datetime.now())
    now = now.replace('-', '_').replace(' ', '_').replace(':', '_').replace('.', '_')
    submit_fn = os.path.join(outdir, 'sh', '{}_submit_{}.sh'.format(
        sample_name, now))
    if len(submit_commands) > 0:
        with open(submit_fn, 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write('\n'.join(submit_commands))
        return submit_fn
    else:
        return None
