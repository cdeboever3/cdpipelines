import argparse

import cdpipelines as ps

def main():
    parser = argparse.ArgumentParser(description=(
        'This script takes a VCF and sample name and makes the '
        'SNP directory needed for the mapping part of WASP.'))
    parser.add_argument('vcf', help='VCF file with variants.')
    parser.add_argument('vcf_out', help=(
        'Output VCF file with heterozygous variants for this sample.'))
    parser.add_argument('sample_name', help=('Sample name in VCF file.'))
    parser.add_argument('snp_directory', help=(
        'Output WASP SNP directory. This directory contains files with'
        ' all heterozygous variants for this sample.'))
    parser.add_argument('regions', help=(
        'Path to bed file to define regions of interests (e.g. exons, peaks, '
        'etc.).'))
    parser.add_argument('-s', metavar='sequence_dict', help=(
        'Sequence dictionary from Picard CreateSequenceDictionary. If '
        'provided, the output VCF will be sorted to match the contig order in '
        'the sequence_dict file. When sequence_dict is provided, some temporary'
        ' files are written in the current working directory while the script '
        'runs.'), default=None)
    parser.add_argument('-m', metavar='picard_memory', help=(
        'Amount of memory in Gb to give Picard. Only used when sequence_dict is'
        ' provided.'), default=4)
    parser.add_argument('-t', metavar='tempdir', help=(
        'Path to temporary directory. Only used when sequence_dict is '
        'provided.'), default='.')
    parser.add_argument('-b', metavar='bcftools_path', help=(
        'Path to bcftools executable. By default, assumed to be in your path.'),
        default='bcftools')
    parser.add_argument('-p', metavar='picard_path', help=(
        'Path to picard jarfile. By default, assumed to be $picard.'),
        default='$picard')
        
    args = parser.parse_args()
    vcf = args.vcf
    vcf_out = args.vcf_out
    sample_name = args.sample_name
    directory = args.snp_directory
    regions = args.regions
    sequence_dict = args.s
    picard_memory = args.m
    tempdir = args.t
    bcftools_path = args.b
    picard_path = args.p

    ps.general._wasp_snp_directory(
        vcf, 
        directory, 
        sample_name, 
        regions,
        vcf_out, 
        sequence_dict=sequence_dict,
        picard_memory=picard_memory, 
        tempdir=tempdir,
        bcftools_path=bcftools_path,
        picard_path=picard_path,
    )

if __name__ == '__main__':
    main()
