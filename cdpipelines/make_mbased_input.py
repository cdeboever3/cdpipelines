import argparse

import pandas as pd
pd.options.mode.chained_assignment = None

def _ref_freq(df):
    """
    Estimate the reference bias for all 12 possible substitutions and add them
    to the count matrix.

    Parameters
    ----------
    df : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    Returns
    -------
    counts : pandas.DataFrame
        ASEReadCounter results with estimated reference allele percentage added.

    """
    import copy
    counts = copy.deepcopy(df)
    counts['refCount_norm'] = counts.refCount
    counts['altCount_norm'] = counts.altCount
    counts['totalCount_norm'] = counts.totalCount
    t = counts[counts.totalCount > 30]
    t['refCount_norm'] = (t.refCount / t.totalCount) * 30
    t['altCount_norm'] = (t.altCount / t.totalCount) * 30
    counts.ix[t.index, 'refCount_norm'] = t.refCount_norm
    counts.ix[t.index, 'altCount_norm'] = t.altCount_norm
    counts.ix[t.index, 'totalCount_norm'] = 30
    counts['snp_combo'] = (counts.refAllele + counts.altAllele)
    freq = {}
    for c in set(counts.snp_combo):
        t = counts[counts.snp_combo == c]
        freq[c] = (t.refCount_norm / t.totalCount_norm).mean()
    freq = pd.Series(freq)
    df['expectedRefFreq'] = freq[df.refAllele + df.altAllele].values
    return df

def _binomial_test(counts):
    """
    Add binomial p-value for each SNV.

    Parameters
    ----------
    counts : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    Returns
    -------
    counts : pandas.DataFrame
        ASEReadCounter results with binomial p-values added.

    """
    from scipy.stats import binom_test
    pvals = []
    for i in counts.index:
        p = binom_test(counts.ix[i, 'refCount'], counts.ix[i, 'totalCount'],
                       counts.ix[i, 'expectedRefFreq'])
        pvals.append(p)
    counts['binomialPValue'] = pvals
    return counts

def _vcf_filter(
    counts, 
    vcfs, 
    sample_name, 
    chrom_conv=None, 
    min_dist=10,
):
    """
    Remove heterozygous variants that are within min_dist of another variant
    for which the sample is heterozygous or homozygous alternative. The idea
    here is that the nearby variant could affect mapping and affect ASE
    estimation.

    Parameters
    ----------
    counts : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    vcfs : list
        List of paths to gzipped, indexed VCF file with all variant calls.

    sample_name : str
        Sample name for this sample in the VCF file.

    min_dist : int
        Filter variants with other non-reference variants within this window.

    Returns
    -------
    counts : pandas.DataFrame
        Filtered ASEReadCounter results with SNVs near other non-reference
        variants removed.

    """
    import vcf as pyvcf
    if chrom_conv:
        conv = pd.read_table(chrom_conv, header=None, index_col=1, squeeze=True)
    readers = [pyvcf.Reader(open(vcf), compressed=True) for vcf in vcfs]
    remove = []
    for i in counts.index:
        for reader in readers:
            chrom, pos = counts.ix[i, ['contig', 'position']]
            start = int(pos) - 10
            end = int(pos) + 10
            if chrom_conv:
                chrom = conv[chrom]
            try:
                res = reader.fetch(chrom, start, end)
                count = 0
                while True:
                    try:
                        r = res.next()
                        hets = [x.sample for x in r.get_hets()]
                        hom_alts = [x.sample for x in r.get_hom_alts()]
                        if (r.CHROM == chrom and r.POS < end and 
                            (sample_name in hets or sample_name in hom_alts)):
                            count += 1
                        else:
                            break
                    except StopIteration:
                        break
            except ValueError:
                continue
        if count > 1:
            remove.append(i)
    counts = counts.drop(remove)
    return counts

def _mappability_filter(counts, mappability, bigWigAverageOverBed_path):
    """
    Remove heterozygous SNVs that are not in unique regions according to the
    mappability bigwig file.

    Parameters
    ----------
    counts : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    mappability : str
        Path to bigwig file with mappability scores. A score of one should mean
        uniquely mapping.

    bigWigAverageOverBed_path : str
        Path to bigWigAverageOverBed executable.

    Returns
    -------
    counts : pandas.DataFrame
        Filtered ASEReadCounter results with SNVs in low mappability regions
        removed.

    """
    import subprocess
    import tempfile
    tempfile.NamedTemporaryFile()
    temp_bed = tempfile.NamedTemporaryFile()
    tdf = counts[['contig', 'position']]
    tdf['start'] = tdf.position - 1
    tdf['name'] = tdf.contig + ':' + tdf.position.astype(str)
    tdf = tdf.sort_values(by=['contig', 'start'])
    tdf[['contig', 'start', 'position', 'name']].to_csv(temp_bed.name,
                                                        index=None, header=None,
                                                        sep='\t')
    out = tempfile.NamedTemporaryFile()
    c = ('{} {} {} {}'.format(bigWigAverageOverBed_path, mappability,
                              temp_bed.name, out.name))
    subprocess.check_call(c, shell=True)
    res = pd.read_table(out.name, index_col=0, header=None)
    counts = counts.ix[res[res[5] == 1].index]
    out.close()
    temp_bed.close()
    return counts

def _min_dist_filter(counts, min_dist=300):
    """
    Remove heterozygous SNVs with min_dist base pairs of each (genomic
    distance). The idea here is that we don't want to include SNVs that are both
    covered by the same read pairs because these are not independent
    observations. A future counting scheme could collapse such variants into a
    single site and count the number of independent reads coming from each
    haplotype.

    Parameters
    ----------
    counts : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    Returns
    -------
    counts : pandas.DataFrame
        Filtered ASEReadCounter results with nearby sites removed.

    """
    import scipy.spatial.distance as dist
    keep = set()
    for c in set(counts.contig):
        t = counts[counts.contig == c]
        d = pd.DataFrame(dist.cdist(t.position.values.reshape([t.shape[0], 1]), 
                                    t.position.values.reshape([t.shape[0], 1])),
                         index=t.index, columns=t.index)
        se = d[d != 0].min()
        while se.min() <= min_dist:
            se = se[se <= 300]
            i = 0
            remove = []
            while i < se.shape[0] - 1:
                a = se.ix[i]
                b = se.ix[i + 1]
                if a == b:
                    if t.ix[se.index[i], 'totalCount'] < t.ix[se.index[i + 1],
                                                              'totalCount']:
                        remove.append(se.index[i])
                    else:
                        remove.append(se.index[i + 1])
                    i += 2
                else:
                    i += 1
            t = t.drop(remove)
            d = pd.DataFrame(
                dist.cdist(t.position.values.reshape([t.shape[0], 1]),
                           t.position.values.reshape([t.shape[0], 1])),
                index=t.index, columns=t.index)
            se = d[d != 0].min()
        keep |= set(t.index)
    return counts.ix[keep]

def _assign_features(counts, bed):
    """
    Assign each SNV to a feature for used with MBASED. A feature could be a
    gene, a histone peak, etc. If a SNV overlaps multiple features, it will be
    removed.

    Parameters
    ----------
    counts : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    bed : str
        Path to bed file with feature information. The feature name should be
        the fourth column.

    Returns
    -------
    counts : pandas.DataFrame
        Filtered ASEReadCounter results with nearby sites removed.

    """
    import pybedtools as pbt
    features = pbt.BedTool(bed)
    features = features.sort()
    tdf = counts[['contig', 'position']]
    tdf['start'] = tdf.position - 1
    tdf['name'] = tdf.contig + ':' + tdf.position.astype(str)
    s = '\n'.join(tdf.contig + '\t' + tdf.start.astype(str) + 
                  '\t' + tdf.position.astype(str) + 
                  '\t' + tdf.name) + '\n'
    bt = pbt.BedTool(s, from_string=True)
    bt = bt.sort()
    res = bt.intersect(features, wo=True, sorted=True)
    
    has_name = len(features[0].fields) > 3
    snv_to_feature = dict()
    for r in res:
        snv = r.fields[3]
        if has_name:
            feature = r.fields[7]
        else:
            feature = r.fields[4] + ':' + r.fields[5] + '-' + r.fields[6]
        snv_to_feature[snv] = snv_to_feature.get(r, []) + [feature]
        
    se = pd.Series(snv_to_feature)
    se = se[se.apply(lambda x: len(x) == 1)]
    se = se.apply(lambda x: x[0])
    counts = counts.ix[se.index]
    counts.ix[se.index, 'feature'] = se
    return counts

def make_mbased_input(
    counts, 
    out, 
    bed, 
    vcfs=None, 
    chrom_conv=None,
    sample_name=None,
    mappability=None, 
    min_cov=8,
    min_ref_percent=0.02,
    max_ref_percent=0.98,
    min_near_dist=10,
    bigWigAverageOverBed_path='bigWigAverageOverBed',
):
    """
    Make MBASED input file based on ASEReadCounter output. Also output file like
    ASEReadCounter file but only containing filtered sites and binomial test
    results for all sites that pass filtering. The default parameters for this
    script are based on the GTEx and MBASED papers (10.1038/nature12531,
    10.1126/science.1262110, 10.1186/s13059-014-0405-3).

    Parameters
    ----------
    counts : pandas.DataFrame
        ASEReadCounter results as a dataframe.

    out : str
        Output tsv file for MBASED.

    vcfs : str
        List of paths to gzipped, indexed VCF files with all variant calls (not
        just heterozygous calls).

    sample_name : str
        If vcf is provided, this must be provided to specify the sample name of
        this sample in the VCF file. Required if vcf is provided.

    mappability : str
        Path to bigwig file with mappability scores. A score of one should mean
        uniquely mapping.

    min_cov : int
        Sites with coverage less than this amount will be filtered out.   

    min_ref_percent : float
        Sites where the percentage of reads covering the reference allele is
        less than this amount will be filtered out. 

    max_ref_percent : float
        Sites where the percentage of reads covering the reference allele is
        greater than this amount will be filtered out.

    min_near_dist : int
        SNVs that are within min_near_dist of other variants will be removed to
        avoid mapping biases.

    bigWigAverageOverBed_path : str
        Path to bigWigAverageOverBed. Required if mappability is provided.

    """
    if vcfs is not None:
        assert sample_name is not None
    if type(counts) == str:
        counts = pd.read_table(counts)
    assert type(counts) == pd.DataFrame
    
    # This method applies several filters (depending on the input provided) and
    # the makes the final input file for MBASED. I've mentioned the references
    # for the filters in comments below. The GTEx methods are primarily in
    # 10.1038/nature12531 and 10.1126/science.1262110. The MBASED paper is
    # 10.1186/s13059-014-0405-3.

    counts.index = counts.contig + ':' + counts.position.astype(str)
    # Only keep SNVs with 8 or more counts. GTEx filter.
    counts = counts[counts.totalCount >= min_cov]
    # Only keep SNVs where the percentage of reads that are reference is between
    # 2% and 98%. GTEx filter.
    counts = counts[(min_ref_percent <= counts.refCount / counts.totalCount) &
                    (max_ref_percent >= counts.refCount / counts.totalCount)]
    # Remove SNVs that are too close to each other to avoid double-counting read
    # pairs. Necessary for using MBASED since it assumes counts at different
    # variants are independent.
    counts = _min_dist_filter(counts)

    if mappability:
        # Remove SNVs with non-unique mappability. GTEx filter.
        counts = _mappability_filter(counts, mappability,
                                     bigWigAverageOverBed_path)
    if vcfs:
        # Remove SNVs that are within 10 bp of other variants to avoid mapping
        # biases. MBASED filter.
        counts = _vcf_filter(counts, vcfs, sample_name, chrom_conv=chrom_conv,
                             min_dist=min_near_dist)
    # Assign SNVs to features.
    counts = _assign_features(counts, bed)
    # Add expected reference frequency for each ref->alt substitution.
    counts = _ref_freq(counts)
    # Add binomial p-values.
    counts = _binomial_test(counts)
    # Add total number of reads for each locus.
    t = counts[['totalCount', 'feature']].groupby('feature').sum()
    counts['totalFeatureCount'] = t.ix[counts.feature].values
    counts.to_csv(out, sep='\t')

def main():
    parser = argparse.ArgumentParser(description=(
        'This script takes allele counts from GATK\'s ASEReadCounter and '
        'creates an input file for MBASED. Binomal ASE p-values are also '
        'generated for comparison with the MBASED results.'))
    parser.add_argument('counts', help='Output file from ASEReadCounter.')
    parser.add_argument('out', help='Output tsv file for MBASED.')
    h = ('Bed file for assigning SNVs to features. If a SNV overlaps more than '
         'one feature, it will be removed.')
    parser.add_argument('bed', help=h)
    parser.add_argument(
        '-v', 
        metavar='vcfs', 
        required=True, 
        action='append', 
        help=(
            'VCF files with all variant calls. Heterozygous SNVs within 10 bp '
            'of another variant will be removed to avoid mapping bias. '
            'Multiple -v VCFs can be provided (-v '
            'chr1.vcf.gz -v chr2.vcf.gz) but they shouldn\'t overlap in '
            'genomic coordinates (e.g. they should be for separate chromosomes '
            'etc.).'),
    )
    parser.add_argument('-c', metavar='chrom_conv', help=(
        'File with VCF chromosomes in first column and corresponding RNA-seq '
        'chromosomes in second column (no header). This is needed if the VCF '
        'and RNA-seq data have different chromosome naming.'))
    parser.add_argument('-s', metavar='sample_name', default=None, 
                        help=('If -v is provided, this is the sample name for '
                              'this sample in the VCF file. Required with -v.'))
    h = ('Bigwig file with mappability scores (1 means unique mappability). '
         'SNVs that do not have unique mappability will be removed.')
    parser.add_argument('-m', metavar='mappability', default=None, help=h)
    parser.add_argument('-p', metavar='bigWigAverageOverBed_path',
                        default='bigWigAverageOverBed',
                        help=('Path to bigWigAverageOverBed. Only needed if -m '
                              'provided. Default: bigWigAverageOverBed.'))
    parser.add_argument('-mc', metavar='min_cov', default=8, 
                        help=('Sites with coverage less than this amount will '
                              'be filtered out. Default: 8.'))
    parser.add_argument('-mirp', metavar='min_ref_percent', default=0.02,
                        help=('Sites where the percentage of reads covering '
                              'the reference allele is less than this '
                              'amount will be filtered out. Default: 0.02.'))
    parser.add_argument('-marp', metavar='max_ref_percent', default=0.98,
                        help=('Sites where the percentage of reads covering '
                              'the reference allele is greater than this '
                              'amount will be filtered out. Default: 0.98.'))
    parser.add_argument('-mnd', metavar='min_near_dist', default=10,
                        help=('SNVs that are within min_near_dist of other '
                              'variants will be removed to avoid mapping '
                              'biases. Default: 10.'))

    args = parser.parse_args()
    counts = args.counts
    out = args.out
    bed = args.bed
    vcfs = args.v
    chrom_conv = args.c
    sample_name = args.s
    mappability = args.m
    bigWigAverageOverBed_path = args.p
    min_cov = args.mc
    min_ref_percent = args.mirp
    max_ref_percent = args.marp
    min_near_dist = args.mnd

    make_mbased_input(
        counts, 
        out, 
        bed, 
        vcfs=vcfs,
        chrom_conv=chrom_conv,
        sample_name=sample_name, 
        mappability=mappability,
        bigWigAverageOverBed_path=bigWigAverageOverBed_path,
        min_cov=min_cov,
        min_ref_percent=min_ref_percent,
        max_ref_percent=max_ref_percent,
        min_near_dist=min_near_dist)

if __name__ == '__main__':
    main()
