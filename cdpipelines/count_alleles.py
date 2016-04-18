import argparse

def main():
    parser = argparse.ArgumentParser(description=(
        'This script takes a bam file and list of positions and counts the '
        'the number of each nucleotide observed for that position. Overlapping '
        'mate pairs will only be counted once.'))
    parser.add_argument('bam', help='Bam file to count for.')
    parser.add_argument('positions', help=(
        'TSV file whose first column is chromosome and second is one-based '
        'position. This file should have a header.'))
    parser.add_argument('output', help='Output file path.')
    parser.add_argument('--stranded', action='store_true', help=(
        'Use if data is stranded.'))
    args = parser.parse_args()
    bam = args.bam
    positions = args.positions
    output = args.output
    stranded = args.stranded

    import cdpybio as cpb
    import pandas as pd
    import pybedtools as pbt

    tdf = pd.read_table(positions)
    lines = tdf.ix[:, 0] + '\t' + tdf.ix[:, 1].astype(str) + '\n'
    lines = ''.join(tdf.ix[:, 0] + '\t' + (tdf.ix[:, 1] - 1).astype(str) + 
                    '\t' + tdf.ix[:, 1].astype(str) + '\n')
    bt = pbt.BedTool(lines, from_string=True)
    counts = cpb.pysamext.nt_counts(bam, bt, stranded=stranded)
    counts.to_csv(output, sep='\t')

if __name__ == '__main__':
    main()
