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
        'etc.). These regions should be non-overlapping.'))
    parser.add_argument('-f', metavar='geneome_fai', help=(
        'Path genome fasta fai file. If provided, the VCF will be sorted to '
        'match the contig order in the fai file. When fai is provided, some '
        'temporary files are written in the current working directory while '
        'the script runs.'), default=None)
    parser.add_argument('-b', metavar='bcftools_path', help=(
        'Path bcftools executable. By default, assumed to be in your path.'),
        default='bcftools')

    args = parser.parse_args()
    vcf = args.vcf
    vcf_out = args.vcf_out
    sample_name = args.sample_name
    directory = args.snp_directory
    regions = args.regions
    fai = args.f
    bcftools_path = args.b

    ps.general._wasp_snp_directory(vcf, directory, sample_name, regions,
                                   vcf_out, fai=fai,
                                   bcftools_path=bcftools_path)

if __name__ == '__main__':
    main()
