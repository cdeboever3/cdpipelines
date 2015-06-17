import argparse

import pipelines as ps

def main():
    parser = argparse.ArgumentParser(description=(
        'This script takes a VCF and optional sample name and makes the '
        'SNP directory needed for the mapping part of WASP.'))
    parser.add_argument('vcf', help='VCF file with exonic variants.')
    parser.add_argument('snp_directory', help='Output WASP SNP directory.')
    parser.add_argument('all_snps', help='Output file for combined SNP info.')
    parser.add_argument('-s', metavar='sample_name', help=(
        'Sample name if VCF file contains variant calls for more than one '
        'sample.'))
    args = parser.parse_args()
    vcf = args.vcf
    directory = args.snp_directory
    all_snps = args.all_snps
    sample_name = args.s

    ps.general._wasp_snp_directory(vcf, directory, all_snps,
                                   sample_name=sample_name)


if __name__ == '__main__':
    main()
