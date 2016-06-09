import argparse
import os
import subprocess

import cdpipelines as ps

def _wasp_snp_directory(
    vcfs, 
    directory, 
    vcf_sample_name, 
    regions, 
    vcf_out,
    gatk_fai=None, 
    vcf_chrom_conv=None,
    tempdir='.', 
    bcftools_path='bcftools',
):
    """
    Convert VCF files into input files directory and files needed for WASP. Only
    bi-allelic heterozygous sites are used. Both SNPs and indels are included.

    Parameters:
    -----------
    vcfs : list
        List of path to VCF files.

    directory : str
        Output directory. A directory snps will be output in this directory with
        the variants for WASP.

    vcf_sample_name : str
        Use this sample name to get heterozygous SNPs from VCF file.

    regions : str
        Path to bed file to define regions of interests (e.g. exons, peaks,
        etc.). 

    vcf_out : str
        Path to output vcf file that contains heterozygous variants for this
        sample. 

    gatk_fai : str
        Path to karyotypically sorted fasta index (fai) file that works with
        GATK.  The output VCF file will be sorted in the order of this fai file
        for compatibility with GATK. Assumed to have associated fai and dict
        files.

    vcf_chrom_conv : str
        File with VCF chromosomes in first column and corresponding RNA-seq
        chromosomes in second column (no header). This is needed if the VCF
        and RNA-seq data have different chromosome naming.

    tempdir : str
        Path to temporary directory. 

    """
    import glob
    import pybedtools as pbt
    import vcf as pyvcf
    
    # Collapse bed file.
    bt = pbt.BedTool(regions)
    bt = bt.merge()
    # If a vcf_chrom_conv file is provided, we should check to see which
    # chromosome naming scheme was used for this bed file. If it's not the
    # correct naming scheme for the VCF, we'll convert.
    if vcf_chrom_conv:
        import pandas as pd
        conv = pd.read_table(vcf_chrom_conv, header=None, index_col=1,
                             squeeze=True)
        df = bt.to_dataframe()
        # Check to see whether any chromosomes have RNA-seq naming convention.
        # If so, we'll assume the bed file is in the RNA-seq naming convention.
        # Any chromosomes not in our conversion index will be dropped from the
        # bed file.
        if len(set(df.chrom) & set(conv.index)) > 0:
            df = df[df.chrom.apply(lambda x: x in conv.index)]
            df['chrom'] = df.chrom.apply(lambda x: conv[x])
            df = df.astype(str)
            s = '\n'.join(df.apply(lambda x: '\t'.join(x), axis=1)) + '\n'
            bt = pbt.BedTool(s, from_string=True)

    # Extract all heterozygous variants for this sample. We'll write the files
    # needed for WASP as well as a VCF with just the hets for this sample.
    if not os.path.exists(directory):
        os.makedirs(directory)

    temp_vcfs = []
    for vcf in vcfs:
        # I'll check to see whether the sample is in the VCF. This might be
        # useful for sex chromosome calls that are split up by sex etc.
        if vcf[-3:]:
            vcf_reader = pyvcf.Reader(open(vcf), compressed=True)
        else:
            vcf_reader = pyvcf.Reader(open(vcf))
        if vcf_sample_name in vcf_reader.samples:
            root = os.path.splitext(os.path.split(vcf)[1])[0]
            tvcf = os.path.join(tempdir, '{}_hets.tsv'.format(root))
            temp_vcfs.append(tvcf)
            c = ('{} view -O u -m2 -M2 \\\n\t-R {} \\\n\t-s {} \\\n\t{} \\\n\t'
                 '| {} view -g het '.format(
                     bcftools_path, bt.fn, vcf_sample_name, vcf, bcftools_path))
            if vcf_chrom_conv:
                c += '-Ou \\\n\t| {} annotate --rename-chrs {} '.format(
                    bcftools_path, vcf_chrom_conv)
                 
            c += ('\\\n\t| tee {} \\\n\t| grep -v ^\\# \\\n\t| cut -f1,2,4,5 '
                  '\\\n\t| awk \'{{print $2"\\t"$3"\\t"$4 >> '
                  '("{}/"$1".snps.txt")}}\''.format(
                     tvcf, directory))
            subprocess.check_call(c, shell=True)
        else:
            vcfs.remove(vcf)

    # If multiple VCFs, concatenate the temp VCFs to make the final VCF.
    # Otherwise, just rename the single temp VCF.
    if len(vcfs) > 1:
        c = '{} concat {} > {}'.format(bcftools_path, ' '.join(temp_vcfs),
                                                               vcf_out)
        subprocess.check_call(c, shell=True)
        [os.remove(x) for x in temp_vcfs]
    else:
        os.rename(temp_vcfs[0], vcf_out)

    # Now we gzip the files.
    fns = glob.glob(os.path.join(directory, '*.snps.txt'))
    for fn in fns:
        subprocess.check_call('gzip {}'.format(fn), shell=True)

    # If a gatk_fai file is provided, we need to reorder the VCF to match
    # the fai.
    if gatk_fai:
        vcf_sorted = os.path.join(tempdir,
                                  '{}_sorted.vcf'.format(vcf_sample_name))
        from __init__ import _scripts
        sortByRef_path = os.path.join(_scripts, 'sortByRef.pl')
        c = 'perl {} {} {} > {}'.format(sortByRef_path, vcf_out, gatk_fai,
                                           vcf_sorted)
        subprocess.check_call(c, shell=True)
        vcf_header = os.path.join(tempdir,
                                  '{}_header.vcf'.format(vcf_sample_name))
        c = 'grep ^\\# {} > {}'.format(vcf_sorted, vcf_header)
        subprocess.check_call(c, shell=True)
        vcf_body = os.path.join(tempdir,
                                  '{}_body.vcf'.format(vcf_sample_name))
        c = 'grep -v ^\\# {} > {}'.format(vcf_sorted, vcf_body)
        subprocess.check_call(c, shell=True)
        c = 'cat {} {} > {}'.format(vcf_header, vcf_body, vcf_out)
        subprocess.check_call(c, shell=True)
        os.remove(vcf_sorted)
        os.remove(vcf_header)
        os.remove(vcf_body)

def main():
    parser = argparse.ArgumentParser(description=(
        'This script takes a VCF and sample name and makes the '
        'SNP directory needed for the mapping part of WASP.'))
    parser.add_argument('vcf_out', help=(
        'Output VCF file with heterozygous variants for this sample.'))
    parser.add_argument('sample_name', help=('Sample name in VCF file.'))
    parser.add_argument('snp_directory', help=(
        'Output WASP SNP directory. This directory contains files with'
        ' all heterozygous variants for this sample.'))
    parser.add_argument('regions', help=(
        'Path to bed file to define regions of interests (e.g. exons, peaks, '
        'etc.).'))
    parser.add_argument(
        '-v', 
        metavar='vcfs', 
        required=True, 
        action='append', 
        help=(
            'VCF files with variants. Multiple -v VCFs can be provided (-v '
            'chr1.vcf.gz -v chr2.vcf.gz) but they shouldn\'t overlap in '
            'genomic coordinates (e.g. they should be for separate chromosomes '
            'etc.).'),
    )
    parser.add_argument('-g', metavar='gatk_fasta', help=(
        'Path to karyotypically sorted fasta index (fai) file that works with '
        'GATK. The output VCF file will be sorted in the order of this fai '
        'file for compatibility with GATK. Assumed to have associated fai and '
        'dict files.'), default=None)
    parser.add_argument('-c', metavar='chrom_conv', help=(
        'File with VCF chromosomes in first column and corresponding RNA-seq '
        'chromosomes in second column (no header). This is needed if the VCF '
        'and RNA-seq data have different chromosome naming.'))
    parser.add_argument('-t', metavar='tempdir', help=(
        'Path to temporary directory. Only used when sequence_dict is '
        'provided.'), default='.')
    parser.add_argument('-b', metavar='bcftools_path', help=(
        'Path to bcftools executable. By default, assumed to be in your path.'),
        default='bcftools')
        
    args = parser.parse_args()
    vcfs = args.v
    vcf_out = args.vcf_out
    sample_name = args.sample_name
    directory = args.snp_directory
    regions = args.regions
    gatk_fai = args.g
    chrom_conv = args.c
    tempdir = args.t
    bcftools_path = args.b

    _wasp_snp_directory(
        vcfs, 
        directory, 
        sample_name, 
        regions,
        vcf_out, 
        gatk_fai=gatk_fai,
        vcf_chrom_conv=chrom_conv,
        tempdir=tempdir,
        bcftools_path=bcftools_path,
    )

if __name__ == '__main__':
    main()
