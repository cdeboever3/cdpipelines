import os

from general import _bigwig_files
from general import _combine_fastqs
from general import _coverage_bedgraph
from general import _fastqc
from general import JobScript
from general import _make_dir
from general import _make_softlink
from general import _mbased
from general import _picard_coord_sort
from general import _picard_index

def _star_align(
    r1_fastq, 
    r2_fastq, 
    sample, 
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
                          '\t--genomeLoad {}'.format(genome_load),
                          '\t--readFilesCommand zcat',
                          '\t--readFilesIn {} {}'.format(r1_fastq, r2_fastq),
                          '\t--outSAMattributes All', 
                          '\t--outSAMunmapped Within',
                          ('\t--outSAMattrRGline ID:1 ' + 
                           'PL:{} '.format(rgpl) + 
                           'PU:{} '.format(rgpu) + 
                           'LB:{0} SM:{0}'.format(sample)), 
                          '\t--outFilterMultimapNmax 20', 
                          '\t--outFilterMismatchNmax 999',
                          '\t--alignIntronMin 20',
                          '\t--alignIntronMax 1000000',
                          '\t--alignMatesGapMax 1000000',
                          '\t--outSAMtype BAM Unsorted'
                         ]))
    if transcriptome_align:
        line +=  ' \\\n\t--quantMode TranscriptomeSAM'
    line += '\n\n'
    line += 'if [ -d _STARtmp ] ; then rm -r _STARtmp ; fi\n\n'
    return line

def pipeline(
    r1_fastqs, 
    r2_fastqs, 
    outdir, 
    sample_name, 
    star_index,
    link_dir,
    web_path_file,
    ref_flat, 
    rrna_intervals,
    dexseq_annotation,
    gene_gtf,
    exon_bed,
    rsem_reference,
    find_intersecting_snps_path, 
    filter_remapped_reads_path,
    genome_fasta,
    vcf, 
    vcf_sample_name=None,
    is_phased=False,
    conda_env=None,
    modules=None,
    queue=None,
    star_genome_load='LoadAndRemove',
    rgpl='ILLUMINA',
    rgpu='',
    strand_specific=True, 
    tempdir=None,
    mappability=None,
    star_path='STAR',
    picard_path='$picard',
    bedtools_path='bedtools',
    bedgraph_to_bigwig_path='bedGraphToBigWig',
    fastqc_path='fastqc',
    samtools_path='samtools',
    rsem_calculate_expression_path='rsem-calculate-expression',
    gatk_path='$GATK',
    bigWigAverageOverBed_path='bigWigAverageOverBed',
    bcftools_path='bcftools',
):
    """
    Make a shell script for aligning RNA-seq reads with STAR. The defaults are
    set for use on the Frazer lab's SGE scheduler on flh1/flh2.

    Parameters
    ----------
    r1_fastqs : list or str
        Either a list of paths to gzipped fastq files with R1 reads or path to a
        single gzipped fastq file with R1 reads.

    r2_fastqs : list or str
        Either a list of paths to gzipped fastq files with R2 reads or path to a
        single gzipped fastq file with R2 reads.

    outdir : str
        Directory to store shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

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

    ref_flat : str
        Path to refFlat file with non-rRNA genes. Can ge gzipped.

    rrna_intervals : str
        Path to interval list file with rRNA intervals.

    dexseq_annotation : str
        Path to DEXSeq exonic bins GFF file.

    gene_gtf : str
        Path to GTF file with gene annotations.

    exon_bed : str
        Path to bed file with exon definitions. The exons should be merged so
        that no bed file entries overlap each other.

    rsem_reference : str
        Directory with RSEM reference.

    find_intersecting_snps_path : str
        Path to find_intersecting_snps.py from WASP.
    
    filter_remapped_reads_path : str
        Path to filter_remapped_reads.py from WASP.

    vcf : str
        VCF file containing exonic variants used for ASE.
    
    vcf_sample_name : str
        Sample name of this sample in the VCF file (if different than
        sample_name). For instance, the sample name in the VCF file may be the
        sample name for WGS data which may differ from the RNA-seq sample name.

    conda_env : str
        Conda environment to load at the beginning of the script.

    modules : str
        Comma-separated list of modules to load at the beginning of the script.

    rgpl : str
        Read Group platform (e.g. illumina, solid). 

    rgpu : str
        Read Group platform unit (eg. run barcode). 

    strand_specific : boolean
        If false, data is not strand specific.

    star_path : str
        Path to STAR aligner.

    picard_path : str
        Path to Picard tools.

    bedtools_path : str
        Path to bedtools.

    bedgraph_to_bigwig_path : str
        Path bedGraphToBigWig executable.

    tempdir : str
        Directory to store temporary files.

    Returns
    -------
    fn : str
        Path to shell script.

    """
    with open(web_path_file) as wpf:
        web_path = wpf.readline().strip()
    tracklines_file = os.path.join(outdir, 'tracklines.txt')
    
    # I'm going to start by naming all of the jobs and their dependencies. Note
    # that if the way that JobScript names jobs changes, this will have to
    # change to reflect that. I'll make an ordered dict whose keys are job
    # scripts and whose values are lists of job names that the script must wait
    # on (i.e. using -hold_jid).
    from collections import OrderedDict
    job_holds = OrderedDict()

    alignment_jobname = '{}_{}'.format(sample_name, 'alignment')
    job_holds[alignment_jobname] = None

    fastqc_jobname = '{}_{}'.format(sample_name, 'fastqc')
    job_holds[fastqc_jobname] = [alignment_jobname]

    sort_mdup_index_jobname = '{}_{}'.format(sample_name, 'sort_mdup_index')
    job_holds[sort_mdup_index_jobname] = [alignment_jobname]

    picard_metrics_jobname = '{}_{}'.format(sample_name, 'picard_metrics')
    job_holds[picard_metrics_jobname] = [sort_mdup_index_jobname]

    md5_jobname = '{}_{}'.format(sample_name, 'md5')
    job_holds[md5_jobname] = [sort_mdup_index_jobname]

    bigwig_jobname = '{}_{}'.format(sample_name, 'bigwig')
    job_holds[bigwig_jobname] = [sort_mdup_index_jobname]

    counts_jobname = '{}_{}'.format(sample_name, 'counts')
    job_holds[counts_jobname] = [sort_mdup_index_jobname]

    rsem_jobname = '{}_{}'.format(sample_name, 'rsem')
    job_holds[rsem_jobname] = [sort_mdup_index_jobname]

    wasp_allele_swap_jobname = '{}_{}'.format(sample_name, 'wasp_allele_swap')
    job_holds[wasp_allele_swap_jobname] = [sort_mdup_index_jobname]

    wasp_remap_jobname = '{}_{}'.format(sample_name, 'wasp_remap')
    job_holds[wasp_remap_jobname] = [wasp_allele_swap_jobname]

    wasp_alignment_compare_jobname = '{}_{}'.format(sample_name,
                                                    'wasp_alignment_compare')
    job_holds[wasp_alignment_compare_jobname] = [wasp_remap_jobname]

    mbased_jobname = '{}_{}'.format(sample_name, 'mbased')
    job_holds[mbased_jobname] = [wasp_alignment_compare_jobname]
    
    ##### Job 1: Combine fastqs and align with STAR. #####
    job_suffix = 'alignment'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(alignment_jobname))):
        job = JobScript(sample_name, job_suffix, 
                        os.path.join(outdir, 'alignment'), threads=8, memory=32,
                        tempdir=tempdir, queue=queue, conda_env=conda_env,
                        modules=modules)
        
        # Files that will be created.
        combined_r1 = os.path.join(
            job.tempdir, '{}_R1.fastq.gz'.format(sample_name))
        combined_r2 = os.path.join(
            job.tempdir, '{}_R2.fastq.gz'.format(sample_name))
        temp_bam = os.path.join(job.tempdir, 'Aligned.out.bam')
        coord_sorted_bam = os.path.join(
            job.tempdir, '{}_coord_sorted.bam'.format(sample_name))

        # Other STAR Files to copy to output directory.
        job.output_files_to_copy += [
            'Log.out', 
            'Log.final.out', 
            'Log.progress.out', 
            'SJ.out.tab', 
            'Aligned.toTranscriptome.out.bam',
        ]
        transcriptome_bam = os.path.join(job.outdir,
                                         'Aligned.toTranscriptome.out.bam')

        with open(job.filename, "a") as f:
            # If multiple fastq files, we want to cat them together.
            lines = _combine_fastqs(r1_fastqs, combined_r1, bg=True)
            f.write(lines)
            lines = _combine_fastqs(r2_fastqs, combined_r2, bg=True)
            f.write(lines)
            f.write('wait\n\n')

            # Align reads.
            lines = _star_align(combined_r1, combined_r2, sample_name, rgpl,
                                rgpu, star_index, 8,
                                genome_load=star_genome_load)
            f.write(lines)
            f.write('wait\n\n')
        job.write_end()

    ##### Job 2: Run fastQC. ##### 
    job_suffix = 'fastqc'
    if not os.path.exists(
        os.path.join(outdir, 'sh', '{}.sh'.format(sort_mdup_index_jobname))):
        job = JobScript(sample_name, job_suffix, 
                        os.path.join(outdir, 'qc'), threads=1, memory=4,
                        tempdir=tempdir, queue=queue, conda_env=conda_env,
                        modules=modules)
        
        with open(job.filename, "a") as f:
            # Run fastQC.
            job.temp_files_to_delete.append(combined_r1)
            job.temp_files_to_delete.append(combined_r2)
            lines = _fastqc([combined_r1, combined_r2], job.outdir, job.threads,
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
        

    ##### Job 3: Coordinate sort, mark duplicates and index bam. #####
    job_suffix = 'sort_mdup_index'
    if not os.path.exists(
        os.path.join(outdir, 'sh', '{}.sh'.format(sort_mdup_index_jobname))):
        job = JobScript(sample_name, job_suffix, 
                        os.path.join(outdir, 'alignment'), threads=1, memory=4,
                        tempdir=tempdir, queue=queue, conda_env=conda_env,
                        modules=modules)
        
        with open(job.filename, "a") as f:
            # Coordinate sort.
            job.temp_files_to_delete.append(temp_bam)
            lines = _picard_coord_sort(temp_bam, coord_sorted_bam,
                                       picard_path=picard_path,
                                       picard_memory=job.memory,
                                       picard_tempdir=job.tempdir)
            f.write(lines)
            f.write('wait\n\n')

            # Mark duplicates.
            mdup_bam = os.path.join(
                job.tempdir, '{}_sorted_mdup.bam'.format(sample_name))
            job.output_files_to_copy.append(mdup_bam)
            job.temp_files_to_delete.append(coord_sorted_bam)
            duplicate_metrics = os.path.join(
                job.outdir, '{}_duplicate_metrics.txt'.format(sample_name))
            lines = _picard_mark_duplicates(coord_sorted_bam, mdup_bam,
                                            duplicate_metrics,
                                            picard_path=picard_path,
                                            picard_memory=job.memory,
                                            picard_tempdir=job.tempdir)
            f.write(lines)
            f.write('wait\n\n')
            name = os.path.split(mdup_bam)[1]
            job.add_softlink(os.path.join(job.outdir, name), 
                             os.path.join(link_dir, 'bam', name))
            with open(tracklines_file, "a") as tf:
                tf_lines = ('track type=bam name="{}_rna_bam" '
                            'description="RNAseq for {}" '
                            'bigDataUrl={}/bam/{}\n'.format(
                                sample_name, sample_name, web_path, name))
                tf.write(tf_lines)

            # Index bam file.
            bam_index = '{}.bai'.format(mdup_bam)
            job.output_files_to_copy.append(bam_index)
            lines = _picard_index(mdup_bam, bam_index, picard_path=picard_path,
                                  picard_memory=job.memory,
                                  picard_tempdir=job.tempdir, bg=False)
            f.write(lines)
            f.write('wait\n\n')
            name = os.path.split(bam_index)[1]
            job.add_softlink(os.path.join(job.outdir, name), 
                             os.path.join(link_dir, 'bam', name))
        job.write_end()

    ##### Job 4: Collect Picard metrics. #####
    job_suffix = 'picard_metrics'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(picard_metrics_jobname))):
        job = JobScript(sample_name, job_suffix, os.path.join(outdir, 'qc'),
                        threads=1, memory=4, tempdir=tempdir, queue=queue,
                        conda_env=conda_env, modules=modules)
        
        with open(job.filename, "a") as f:
            # Collect insert size metrics, bam index stats, RNA seq QC.
            lines = _picard_collect_multiple_metrics(mdup_bam, sample_name,
                                                     picard_path=picard_path, 
                                                     picard_memory=job.memory,
                                                     picard_tempdir=job.tempdir,
                                                     bg=False)
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
                mdup_bam, 
                metrics, 
                chart, 
                sample_name,
                ref_flat, 
                rrna_intervals,
                picard_path=picard_path,
                picard_memory=job.memory,
                picard_tempdir=job.tempdir,
                strand_specific=strand_specific, 
                bg=False)
            f.write(lines)

            index_out = os.path.join(job.outdir,
                                     '{}_index_stats.txt'.format(sample_name))
            index_err = os.path.join(job.outdir,
                                     '{}_index_stats.err'.format(sample_name))
            lines = _picard_bam_index_stats(mdup_bam, index_out, index_err,
                                            picard_path=picard_path,
                                            picard_memory=job.memory,
                                            picard_tempdir=job.tempdir,
                                            bg=False)
            f.write(lines)
            f.write('wait\n\n')
        job.write_end()

    ##### Job 5: Make md5 has for final bam file. #####
    job_suffix = 'md5'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(md5_jobname))):
        job = JobScript(sample_name, job_suffix, 
                        os.path.join(outdir, 'alignment'), threads=1, memory=4,
                        tempdir=tempdir, queue=queue, conda_env=conda_env,
                        modules=modules)
        
        with open(job.filename, "a") as f:
            # Make md5 hash for output bam file.
            f.write('md5sum {} > {}\n'.format(
                mdup_bam, os.path.join(job.outdir, '{}.md5'.format(
                    os.path.split(mdup_bam)[1]))))
        job.write_end()
       
    ##### Job 6: Make bigwig for final bam file. #####
    job_suffix = 'bigwig'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(bigwig_jobname))):
        job = JobScript(sample_name, job_suffix, 
                        os.path.join(outdir, 'alignment'), threads=1, memory=4,
                        tempdir=tempdir, queue=queue, conda_env=conda_env,
                        modules=modules)
        
        with open(job.filename, "a") as f:
            # Make bigwig file for displaying coverage.
            out_bigwig = os.path.join(job.tempdir, '{}_rnaseq.bw'.format(sample_name))
            job.output_files_to_copy.append(out_bigwig)
        
            lines = _bigwig_files(
                mdup_bam, out_bigwig, sample_name, 
                bedgraph_to_bigwig_path=bedgraph_to_bigwig_path,
                bedtools_path=bedtools_path)
            f.write(lines)
            f.write('wait\n\n')
            name = os.path.split(out_bigwig)[1]
            job.add_softlink(os.path.join(job.outdir, name), 
                             os.path.join(link_dir, 'bw', name))

            with open(tracklines_file, "a") as tf:
                tf_lines = ('track type=bigWig name="{}_rna_cov" '
                            'description="RNAseq coverage for {}" '
                            'visibility=0 db=hg19 bigDataUrl={}/bw/{}\n'.format(
                                sample_name, sample_name, web_path, name))
                tf.write(tf_lines)
        job.write_end()
    
    ##### Job 7: Get HTSeq and DEXSeq counts. #####
    job_suffix = 'counts'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(counts_jobname))):
        job = JobScript(sample_name, job_suffix, 
                        os.path.join(outdir, 'counts'), threads=1, memory=4,
                        tempdir=tempdir, queue=queue, conda_env=conda_env,
                        modules=modules)
        
        with open(job.filename, "a") as f:
            gene_counts = os.path.join(job.outdir, 'gene_counts.tsv')
            gene_count_stats = os.path.join(job.outdir, 'gene_count_stats.tsv')
            lines = _htseq_count(
                mdup_bam, 
                gene_counts, 
                gene_count_stats, 
                gene_gtf, 
                strand_specific=strand_specific,
                samtools_path=samtools_path)
            f.write(lines)
            dexseq_counts = os.path.join(
                job.outdir, '{}_dexseq_counts.tsv'.format(job.sample_name))
            lines = _dexseq_count(mdup_bam, dexseq_counts, dexseq_annotation,
                                  paired=True, strand_specific=strand_specific,
                                  samtools_path=samtools_path)
            f.write(lines)
        job.write_end()
    
    ##### Job 8: Run RSEM. #####
    job_suffix = 'rsem'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(rsem_jobname))):
        job = JobScript(sample_name, job_suffix, os.path.join(outdir, 'rsem'),
                        threads=8, memory=4, tempdir=tempdir, queue=queue,
                        conda_env=conda_env, modules=modules)
        
        with open(job.filename, "a") as f:
            job.output_files_to_copy += [
                '{}.genes.results'.format(sample_name),
                '{}.isoforms.results'.format(sample_name),
                '{}.stat'.format(sample_name)]
            lines = _rsem_calculate_expression(
                transcriptome_bam, rsem_reference, sample_name,
                threads=job.threads, ci_mem=1024,
                strand_specific=strand_specific,
                rsem_calculate_expression_path=rsem_calculate_expression_path,
            )
            f.write(lines)
        job.write_end()
    
    ##### Job 9: WASP first step. #####
    job_suffix = 'wasp_allele_swap'
    if not os.path.exists(
        os.path.join(outdir, 'sh', '{}.sh'.format(wasp_allele_swap_jobname))):
        job = JobScript(sample_name, job_suffix, os.path.join(outdir, 'wasp'),
                        threads=1, memory=4, tempdir=tempdir, queue=queue,
                        conda_env=conda_env, modules=modules)
        
        with open(job.filename, "a") as f:
            if not vcf_sample_name:
                vcf_sample_name = sample_name

            # Files that will be created.
            temp_uniq_bam = os.path.join(
                job.tempdir, '{}_uniq.bam'.format(sample_name))
            job.temp_files_to_delete.append(temp_uniq_bam)
            
            # Files to copy to output directory.
            prefix = '{}_uniq'.format(sample_name)
            fns = [
                '{}.keep.bam'.format(prefix),
                '{}.remap.fq1.gz'.format(prefix),
                '{}.remap.fq2.gz'.format(prefix),
                '{}.to.remap.bam'.format(prefix),
                '{}.to.remap.num.gz'.format(prefix)
            ]
            job.output_files_to_copy += fns
            wasp_r1_fastq = os.path.join(job.outdir,
                                         '{}.remap.fq1.gz'.format(prefix))
            wasp_r2_fastq = os.path.join(job.outdir,
                                         '{}.remap.fq2.gz'.format(prefix))
            to_remap_bam = os.path.join(job.outdir,
                                        '{}.to.remap.bam'.format(prefix))
            to_remap_num = os.path.join(job.outdir,
                                        '{}.to.remap.num.gz'.format(prefix))

            from __init__ import scripts
            input_script = os.path.join(scripts, 'make_wasp_input.py')

            # Run WASP to swap alleles.
            with open(job.filename, "a") as f:
                snp_directory = os.path.join(job.tempdir, 'snps')
                # temp_bam = os.path.join(job.outdir, 'uniq.bam')
                f.write('python {} \\\n\t{} \\\n\t{} \\\n\t{} '
                        '\\\n\t{} \\\n\t-b {}\n\n'.format(
                            input_script, vcf, vcf_sample_name, snp_directory,
                            exon_bed, bcftools_path))
                f.write('{} view -b -q 255 -F 1024 \\\n\t{} '
                        '\\\n\t> {}\n\n'.format(
                            samtools_path, mdup_bam, temp_uniq_bam))
                f.write('wait\n\n')
                f.write('python {} -s -p \\\n\t{} \\\n\t{}\n\n'.format(
                    find_intersecting_snps_path, temp_uniq_bam, snp_directory))
        job.write_end()
    
    ##### Job 10: WASP second step. #####
    job_suffix = 'wasp_remap'
    if not os.path.exists(os.path.join(outdir, 'sh',
                                       '{}.sh'.format(wasp_remap_jobname))):
        job = JobScript(sample_name, job_suffix, os.path.join(outdir, 'wasp'),
                        threads=8, memory=10, tempdir=tempdir,
                        queue=queue, conda_env=conda_env, modules=modules)
        
        with open(job.filename, "a") as f:
            # Input files.
            temp_r1 = job.add_input_file(wasp_r1_fastq)
            temp_r2 = job.add_input_file(wasp_r2_fastq)
            job.copy_input_files()

            # Files that will be created.
            remapped_bam = os.path.join(job.tempdir, 'Aligned.out.bam')
            job.output_files_to_copy.append(remapped_bam)

            # Files to copy to output directory.
            job.output_files_to_copy += ['Log.out', 'Log.final.out',
                                         'Log.progress.out', 'SJ.out.tab']

            # Remap using STAR.
            lines = _star_align(temp_r1, temp_r2, sample_name, rgpl,
                                rgpu, star_index, threads=job.threads,
                                genome_load=star_genome_load,
                                transcriptome_align=False)
            f.write(lines)
            f.write('wait\n\n')
        job.write_end()
    
    ##### Job 11: WASP third step. #####
    job_suffix = 'wasp_alignment_compare'
    fn = os.path.join(outdir, 'sh',
                      '{}.sh'.format(wasp_alignment_compare_jobname))
    if not os.path.exists(fn):
        job = JobScript(sample_name, job_suffix, os.path.join(outdir, 'wasp'),
                        threads=1, memory=4, tempdir=tempdir, queue=queue,
                        conda_env=conda_env, modules=modules)
        
        with open(job.filename, "a") as f:
            # Files that will be created.
            temp_filtered_bam = os.path.join(
                job.tempdir, '{}_filtered.bam'.format(sample_name))
            job.temp_files_to_delete.append(temp_filtered_bam)
            wasp_filtered_bam = os.path.join(
                job.tempdir, '{}_filtered_coord_sorted.bam'.format(sample_name))
            job.output_files_to_copy.append(wasp_filtered_bam)
            bam_index = wasp_filtered_bam + '.bai'
            job.output_files_to_copy.append(bam_index)
   
            with open(job.filename, "a") as f:
                # Run WASP alignment compare.
                f.write('python {} -p \\\n\t{} \\\n\t{} \\\n\t{} '
                        '\\\n\t{}\n\n'.format(
                            filter_remapped_reads_path, to_remap_bam,
                            remapped_bam, temp_filtered_bam, to_remap_num))

                # Coordinate sort and index.
                lines = _picard_coord_sort(temp_filtered_bam, wasp_filtered_bam,
                                           bam_index=bam_index,
                                           picard_path=picard_path,
                                           picard_memory=job.memory,
                                           picard_tempdir=job.tempdir)

                f.write(lines)
                f.write('\nwait\n\n')
                
                # Count allele coverage.
                counts = os.path.join(
                    job.outdir, '{}_allele_counts.tsv'.format(sample_name))
                f.write('java -jar {} \\\n'.format(gatk_path))
                f.write('\t-R {} \\\n'.format(genome_fasta))
                f.write('\t-T ASEReadCounter \\\n')
                f.write('\t-o {} \\\n'.format(counts))
                f.write('\t-I {} \\\n'.format(wasp_filtered_bam))
                f.write('\t-sites {} \\\n'.format(vcf))
                f.write('\t-overlap COUNT_FRAGMENTS_REQUIRE_SAME_BASE \\\n')
                f.write('\t-U ALLOW_N_CIGAR_READS \n')
                f.write('\nwait\n\n')
        job.write_end()
    
    ##### Job 12: Run MBASED for ASE. #####
    job_suffix = 'mbased'
    fn = os.path.join(outdir, 'sh', '{}.sh'.format(mbased_jobname))
    if not os.path.exists(fn):
        job = JobScript(sample_name, job_suffix, os.path.join(outdir, 'mbased'),
                        threads=8, memory=16, tempdir=tempdir, queue=queue,
                        conda_env=conda_env, modules=modules)
        mbased_infile = os.path.join(job.outdir,
                                     '{}_mbased_input.tsv'.format(sample_name))
        locus_outfile = os.path.join(job.outdir,
                                     '{}_locus.tsv'.format(sample_name))
        snv_outfile = os.path.join(job.outdir, '{}_snv.tsv'.format(sample_name))
    
        with open(job.filename, "a") as f:
            lines = _mbased(wasp_filtered_bam, gene_gtf, mbased_infile,
                            locus_outfile, snv_outfile, sample_name,
                            is_phased=is_phased, threads=job.threads, vcf=vcf,
                            vcf_sample_name=vcf_sample_name,
                            mappability=mappability,
                            bigWigAverageOverBed_path=bigWigAverageOverBed_path)
            f.write(lines)
            f.write('wait\n\n')
        job.write_end()

    ##### Submission script #####
    # Now we'll make a submission script that submits the jobs with the
    # appropriate dependencies.
    submit_fn = os.path.join(outdir, 'sh', 'submit.sh')
    with open(submit_fn, 'w') as f:
        f.write('#!/bin/bash\n\n')
        while True:
            try:
                jn,holds = job_holds.popitem(False)
            except:
                break
            if holds:
                f.write('qsub -hold_jid {} {}.sh\n'.format(
                    ','.join(holds), os.path.join(outdir, jn)))
            else:
                f.write('qsub {}.sh\n'.format(os.path.join(outdir, jn)))

    return submit_fn

def _dexseq_count(
    bam, 
    counts_file, 
    dexseq_annotation, 
    paired=True,
    strand_specific=True, 
    samtools_path='samtools'):
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
        Lines to be printed to shell script.

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
        '{} view -h -f 2 {} \\\n\t'.format(samtools_path, bam) +
        '| cut -f1-16,20- \\\n\t| python {} \\\n\t'.format(script) + 
        '-p {} -s {} -a 0 -r pos -f sam \\\n\t'.format(p, s) + 
        '{} \\\n\t- {}\n\n'.format(dexseq_annotation, counts_file)
    )
    return lines

def _htseq_count(
    bam, 
    counts_file, 
    stats_file, 
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
        Lines to be printed to shell script.

    name : str
        File name for the softlink.

    """
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

    return lines

# TODO: I want to split out the DEXSeq counting. I don't need the HTSeq counting
# anymore. This needs to be refactored for JobScript as well.
def get_counts(
    bam, 
    outdir, 
    sample_name, 
    tempdir, 
    dexseq_annotation, 
    gtf,
    threads=1,
    conda_env=None,
    r_env='', # TODO: get rid of this and replace with modules
    paired=True,
    strand_specific=False,
    samtools_path='samtools',
):
    """
    Make a shell script for counting reads that overlap genes for DESeq2 and
    exonic bins for DEXSeq.

    Parameters
    ----------
    bam : str
        Coordinate sorted bam file (genomic coordinates).

    outdir : str
        Directory to store shell file and aligment results.

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

    """
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

    if conda_env:
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

def _rsem_calculate_expression(
    bam, 
    reference, 
    sample_name,
    threads=1, 
    ci_mem=1024, 
    strand_specific=False,
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
        Lines to be printed to shell script.

    name : str
        File name for the softlink.

    """
    line = ('{} --bam --paired-end --num-threads {} '
            '--no-bam-output --seed 3272015 --calc-ci '
            '--ci-memory {} --estimate-rspd \\\n\t{} \\\n\t{} {}'.format(
                rsem_calculate_expression_path, threads, ci_mem, bam, reference,
                sample_name))
    if strand_specific:
        line += '\\\n\t--forward-prob 0'
    line += '\n'
    return line

# Needs to be refactored for JobScript.
def rsem_expression(
    bam, 
    outdir, 
    sample_name, 
    rsem_reference, 
    ci_mem=1024, 
    r_env='', # TODO: replace with modules
    threads=32,
    tempdir=None,
    strand_specific=False,
    rsem_calculate_expression_path='rsem-calculate-expression',
):
    """
    Make a shell script for estimating expression using RSEM.

    Parameters
    ----------
    bam : str
        Coordinate sorted bam file (genomic coordinates).

    outdir : str
        Directory to store shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    tempdir : str
        Directory to store temporary files.

    rsem_reference : str
        RSEM reference.

    ci_mem : int
        Amount of memory in mb to give RSEM for calculating confidence
        intervals. Passed to --ci-memory for RSEM.

    r_env : str
        If provided, this file will be sourced to set the environment for rpy2.

    strand_specific : boolean
        True if the data is strand-specific. False otherwise.

    """
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

    lines = _rsem_calculate_expression(temp_bam, rsem_reference,
                                       rsem_calculate_expression_path,
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
