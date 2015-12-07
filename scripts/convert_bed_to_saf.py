import argparse

def bed_to_saf(bed, out):
    """
    Convert bed file to SAF file and write to out.
    
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
    # SAF format seems to be inclusive on both ends according to how
    # featureCounts works (chr1 1 2 has length two for instance). Note that if
    # the bed file did not have names, I used the coordinates from the bed file
    # to make the names since working with zero-based, inclusive-exclusive
    # coordinates is more typical. I'm also guessing the SAF format is one-based
    # since it takes GTF by default and GTF is one-based end-inclusive.
    import pandas as pd
    # import pybedtools as pbt

    # bt = pbt.BedTool(bed)
    # df = bt.to_dataframe()
    with open(bed) as f:
        line = f.readline()
    if line.split()[0] == 'track':
        df = pd.read_table(bed, skiprows=1, header=None)
    else:
        df = pd.read_table(bed, header=None)
    if df.shape[1] > 4:
        df.columns = (['chrom', 'start', 'end', 'name', 'strand'] +
                      range(df.shape[1] - 5))
    elif df.shape[1] == 4:
        df.columns = ['chrom', 'start', 'end', 'name']
    else:
        df.columns = ['chrom', 'start', 'end']
    if 'name' not in df.columns:
        df['name'] = (df.chrom + ':' + df.start.astype(str) + '-' +
                      df.end.astype(str))
    elif len(set(df.name)) != df.shape[0]:
        df['name'] = (df.chrom + ':' + df.start.astype(str) + '-' +
                      df.end.astype(str))
    if 'strand' not in df.columns:
        df['strand'] = '+'
    new_df = pd.DataFrame(index=range(df.shape[0]))
    new_df['GeneID'] = df.name
    new_df['Chr'] = df.chrom
    # Make one-based.
    new_df['Start'] = (df.start + 1).astype(int)
    # No need to add one because this is inclusive.
    new_df['End'] = df.end.astype(int)
    new_df['Strand'] = df.strand
    new_df.to_csv(out, index=False, sep='\t')

def main():
    parser = argparse.ArgumentParser(description=(
        'This script converts a bed file to saf format for use with '
        'featureCounts.'))
    parser.add_argument('bed', help='Bed file to convert')
    parser.add_argument('saf', help='Path to write output saf file to.')
    
    args = parser.parse_args()
    bed = args.bed
    saf = args.saf

    bed_to_saf(bed, saf)

if __name__ == '__main__':
    main()
