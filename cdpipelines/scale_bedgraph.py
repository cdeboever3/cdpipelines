import argparse

def scale_bedgraph(bg, out_bg, actual_num, expected_num):
    """
    Scale bedgraph file by multipyling by expected_num / actual_num. This is
    designed to make normalized bigwig files that can be roughly compared
    between samples.
    
    Parameters
    ----------
    bg : str
        Full path to bedgraph file.

    out_bg : str
        Path to output bedgraph file.
    
    actual_num : int
        Number of uniquely mapped read pairs, input read pairs etc. You can use
        whatever value makes sense here as long as you have a corresponding
        expected value. Ideally, this value should be a constant factor
        different than coverage and that factor should be expected to be the
        same between experiments. For instance, if you use number of input read
        pairs for ATAC-seq but some ATAC-seq experiments have 50% pairs align
        uniquely and others have 90% reads map uniquely, the number of input
        read pairs will not work well for normalizing coverage since coverage is
        dependent on the percent uniquely mapped reads.
    
    expected_num : int
        Number of expected reads (pairs). This only needs to be a rough
        estimate. You should keep this number constant for all samples across
        an experiment/project (i.e. any samples you want to compare against each
        other). For example, say you are doing RNA-seq for 10 samples on a lane,
        the lane yields 100,000,000 pairs of reads, and you expect 90% of the
        read pairs to align uniquely. Then you would want expected_num =
        9,000,000. In this case, you would set actual_num to the actual number
        of uniquely mapped reads. The reason you want this number close to the
        number of expected number of uniquely mapped pairs is so the
        normalization doesn't drastically change the coverage for samples near
        the expected number of input read pairs. For instance, if
        expected_num=20M and a sample has exactly 20M uniquely mapped read
        pairs, then the coverage will not be changed.  If expected_num=20M and a
        sample has only 10M input read PAIRS, it will be normalized so that it
        looks like it had 20M input read pairs (i.e. all coverages will be
        multiplied by 2).
    
    Returns
    -------
    link : str
        Full path to softlink.
    
    """
    import subprocess
    factor = expected_num / float(actual_num)
    c = ('grep -v ^track {} | awk \'OFS="\\t" {{$4 = $4 * {}}} '
         '{{print $0}}\' > {}'.format(bg, factor, out_bg))
    subprocess.check_call(c, shell=True)

def num_reads(fn, sambamba_path='sambamba'):
    """Get the number of read pairs in a bam file. The total number of reads is
    divided by 2 to return read pairs."""
    import subprocess
    c = ('{} view -c -F "not (unmapped or mate_is_unmapped) and '
         'mapping_quality >= 255" {}'.format(sambamba_path, fn))
    count = subprocess.check_output(c, shell=True)
    return int(count.strip()) / 2

def main():
    parser = argparse.ArgumentParser(description=(
        'This script takes a bedgraph file and scales the coverage by '
        'multiplying by a factor. The purpose of doing this is generally '
        'to be able to roughly compare different samples with different '
        'library depths.'))
    parser.add_argument('bg', help='Input bedgraph file.')
    parser.add_argument('bam', help='Bam file used to create bedgraph file.')
    parser.add_argument('out_bg', help='Output bedgraph file.')
    parser.add_argument('expected_num', type=float, help=(
        'Number of expected mapped READ PAIRS. This only needs to be a rough '
        'estimate. You should keep this number constant for all samples across '
        'an experiment/project (i.e. any samples you want to compare against '
        'each other). For example, say you are doing RNA-seq for 10 samples on '
        'a lane, the lane yields 100,000,000 pairs of reads, and you expect '
        '90%% of the read pairs to align uniquely. Then you would want '
        'expected_num = 9,000,000. In this case, you would set actual_num to '
        'the actual number of uniquely mapped reads. The reason you want this '
        'number close to the number of expected number of uniquely mapped '
        'pairs is so the normalization does not drastically change the '
        'coverage for samples near the expected number of input read pairs. '
        'For instance, if expected_num=20M and a sample has exactly 20M '
        'uniquely mapped read pairs, then the coverage will not be changed. If '
        'expected_num=20M and a sample has only 10M input read PAIRS, it will '
        'be normalized so that it looks like it had 20M input read pairs (i.e. '
        'all coverages will be multiplied by 2).'))
    parser.add_argument('-s', metavar='sambamba_path',
                        default='sambamba',
                        help=('Path to sambamba.'))

    import os

    args = parser.parse_args()
    bg = os.path.realpath(args.bg)
    bam = os.path.realpath(args.bam)
    out_bg = os.path.realpath(args.out_bg)
    expected_num = args.expected_num
    sambamba_path = args.s
    
    actual_num = num_reads(bam, sambamba_path)
    scale_bedgraph(bg, out_bg, actual_num, expected_num)

if __name__ == '__main__':
    main()
