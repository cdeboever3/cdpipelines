import argparse

import cdpipelines as ps

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
    vcfs = args.vcfs
    vcf_out = args.vcf_out
    sample_name = args.sample_name
    directory = args.snp_directory
    regions = args.regions
    gatk_fai = args.g
    chrom_conv = args.c
    tempdir = args.t
    bcftools_path = args.b

    ps.general._wasp_snp_directory(
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
