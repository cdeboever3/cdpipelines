# This notebook has commands to download external data and software or parse
# external data. These shouldn't be stored on figshare because they are
# available over the web. In some instances I parse some of the external data
# into more useful formats.

import os
import shutil
import subprocess
import sys
from urllib2 import urlopen

def _download_file(url, dest):
    req = urlopen(url)
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)

def _download_and_gunzip(url, dest):
    """
    Download a gzipped file url to dest and gunzip it.

    Parameters
    ----------
    url : str
        URL for gzipped file to download.

    dest : str
        Full path to save gzipped file to. This file will be gunzipped.

    """
    try:
        os.makedirs(os.path.split(dest)[0])
    except OSError:
        pass
    req = urlopen(url)
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['gunzip', dest])

def _download_and_untar(url, dest, outdir, remove_tarball=False):
    """
    Download a tarball url to dest and decompress it in outdir.

    Parameters
    ----------
    url : str
        URL to tarball to download.

    dest : str
        Full path to save tarball to.

    outdir : str
        Directory to save tarball to and decompress to.

    remove_tarball : bool
        If True, remove tarball after decompressing.

    """
    req = urlopen(url)
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call('tar -xf {} -C {}'.format(dest, outdir), shell=True)
    if remove_tarball:
        os.remove(dest)

def make_rna_seq_metrics_files(
    outdir, 
    gencode_gtf, 
    genome_fasta, 
    picard_path='$picard',
    gtfToGenePred_path='gtfToGenePred',
):
    """Make refFlat file and rRNA interval list"""
    import HTSeq
    import itertools as it

    rrna_fn = 'rrna.bed'
    rest_fn = 'rest.gtf'
    rrna = open(rrna_fn, 'w')
    rest = open(rest_fn, 'w')
    gene_rrna = False

    gtf = it.islice(HTSeq.GFF_Reader(gencode_gtf), None)
    line = gtf.next()
    while line != '':
        if line.type == 'gene':
            if (line.attr['gene_type'] == 'rRNA' or
                line.attr['gene_type'] == 'Mt_rRNA'):
                gene_rrna = True
            else:
                gene_rrna = False
                rest.write(line.get_gff_line())
        else:
            if gene_rrna:
                if line.type == 'exon':
                    rrna.write('\t'.join([
                        line.iv.chrom, str(line.iv.start - 1), str(line.iv.end),
                        line.name, '.', line.iv.strand]) + '\n')
            else:
                rest.write(line.get_gff_line())
        try:
            line = gtf.next()
        except StopIteration:
            line = ''

    rrna.close()
    rest.close()

    command = ('{} -ignoreGroupsWithoutExons -genePredExt -geneNameAsName2 {} '
               'refFlat.tmp.txt'.format(gtfToGenePred_path, rest_fn))
    subprocess.check_call(command, shell=True)
    out_fn = os.path.join(outdir, 'gencode_no_rRNA.txt.gz')
    command = ('paste <(cut -f 12 refFlat.tmp.txt) <(cut -f 1-10 '
               'refFlat.tmp.txt) | gzip > {}'.format(out_fn))
    subprocess.check_call(command, shell=True, executable='/bin/bash')
    os.remove(rest_fn)
    os.remove('refFlat.tmp.txt')

    command = ('java -Xmx20g -jar -XX:-UseGCOverheadLimit -XX:-UseParallelGC '
               '-Djava.io.tmpdir=. -jar {} CreateSequenceDictionary '
               'REFERENCE={} OUTPUT=hg19.sam'.format(picard_path, genome_fasta))
    subprocess.check_call(command, shell=True)
    command = ('java -Xmx20g -jar -XX:-UseGCOverheadLimit -XX:-UseParallelGC '
               '-Djava.io.tmpdir=. -jar {} BedToIntervalList '
               'I={} SEQUENCE_DICTIONARY=hg19.sam OUTPUT={}'.format(
                   picard_path, rrna_fn, os.path.join(outdir, 'rrna.interval')))
    subprocess.check_call(command, shell=True)
    os.remove(rrna_fn)
    os.remove('hg19.sam')

def download_igvtools(outdir):
    src = ('http://data.broadinstitute.org/igv/projects/downloads/'
           'igvtools_2.3.55.zip')
    dest = os.path.join(outdir, 'igvtools_2.3.55.zip')
    req = urlopen(src)
    with open(dest, 'w') as f:
        shutil.copyfileobj(req, f)
    subprocess.check_call(['unzip', '-d', outdir, dest])

def download_grasp_database(outdir):
    src = ('https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip')
    dest = os.path.join(outdir, 'GraspFullDataset2.zip')
    req = urlopen(src)
    with open(dest, 'w') as f:
        shutil.copyfileobj(req, f)
    subprocess.check_call(['unzip', '-d', outdir, dest])
    os.rename(os.path.join(outdir, 'GRASP2fullDataset'), 
              os.path.join(outdir, 'GRASP2fullDataset.tsv'))
    os.remove(dest)

def download_nhgri_gwas_catalog(outdir):
    src = ('https://www.ebi.ac.uk/gwas/api/search/downloads/alternative')
    dest = os.path.join(outdir, 'nhgri_gwas_catalog.tsv')
    _download_file(src, dest)

def download_encode_blacklist(outdir):
    src = ('https://www.encodeproject.org/files/ENCFF001TDO/@@download/'
           'ENCFF001TDO.bed.gz')
    dest = os.path.join(outdir, 'encode_blacklist.bed.gz')
    _download_and_gunzip(src, dest)

def download_blat(outdir):
    src = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/blat/blat'
    dest = os.path.join(outdir, 'blat')
    _download_file(src, dest)
    subprocess.check_call(['chmod', '755', '{}'.format(dest)])

def download_weblogo(outdir):
    """
    Download weblogo. 

    Parameters
    ----------
    outdir : str
        Directory to save weblogo to.

    """
    src = ('http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz')
    dest = os.path.join(outdir, 'weblogo')
    _download_and_untar(src, dest, outdir)

def download_epacts(outdir):
    """
    Download and install EPACTS. bin is put in the source directory. If you have
    a specific version of R you are using with a project, you should probably
    have that R in your path. I think EPACTS may install an R package into
    whatever R it finds in your path although I'm not sure about this.

    Parameters
    ----------
    outdir : str
        Directory to save EPACTS to.

    """
    src = ('http://csg.sph.umich.edu/kang/epacts/download/EPACTS-3.2.6.tar.gz')
    dest = os.path.join(outdir, 'EPACTS-3.2.6.tar.gz')
    _download_and_untar(src, dest, outdir)
    cwd = os.getcwd()
    edir = os.path.join(outdir, 'EPACTS-3.2.6')
    os.chdir(edir)
    subprocess.check_call('./configure --prefix={}'.format(edir), shell=True)
    subprocess.check_call('make > EPACTS_make.out 2> EPACTS_make.err',
                          shell=True)
    subprocess.check_call(('make install > EPACTS_make_install.out 2> '
                           'EPACTS_make_install.err'), shell=True)
    os.chdir(cwd)

def download_rsem(outdir, lncurses=False):
    """
    Download RSEM.

    Parameters
    ----------
    outdir : str
        Directory to save RSEM.

    lncurses : bool
        Set to true to use lncurses rather than lcurses to build RSEM. See
        http://seqanswers.com/forums/showthread.php?t=6669 for more information.

    """
    src = ('http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.20.tar.gz')
    dest = os.path.join(outdir, 'rsem-1.2.20.tar.gz')
    _download_and_untar(src, dest, outdir)
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'rsem-1.2.20'))
    if lncurses:
        f = open(os.path.join('sam', 'Makefile'), 'r')
        lines = f.read().replace('lcurses', 'lncurses')
        f.close()
        f = open(os.path.join('sam', 'Makefile'), 'w')
        f.write(lines)
        f.close()

    subprocess.check_call('make')
    os.chdir(cwd)

def download_fastx_toolkit(outdir):
    """
    Download FASTX Toolkit. 

    Parameters
    ----------
    outdir : str
        Directory to save FASTX Toolkit.

    """
    src = ('https://github.com/agordon/fastx_toolkit/releases/download/'
           '0.0.14/fastx_toolkit-0.0.14.tar.bz2')
    dest = os.path.join(outdir, 'fastx_toolkit-0.0.14.tar.bz2')
    _download_and_untar(src, dest, outdir)

def download_fastqc(outdir):
    """
    Download Gencode GTF.

    Parameters
    ----------
    outdir : str
        Directory to save Gencode GTF to.

    """
    src = ('http://www.bioinformatics.babraham.ac.uk/projects/fastqc/'
           'fastqc_v0.11.2.zip')
    dest = os.path.join(outdir, 'fastqc_v0.11.2.zip')
    req = urlopen(src)
    with open(dest, 'w') as f:
        shutil.copyfileobj(req, f)
    subprocess.check_call(['unzip', '-d', outdir, dest])

def download_snpeff(outdir):
    """
    Download snpEff.

    Parameters
    ----------
    outdir : str
        Directory to save snpEff to.

    """
    url = ('http://sourceforge.net/projects/snpeff/'
           'files/snpEff_v4_1g_core.zip/download')
    dest = os.path.join(outdir, 'snpEff_v4_1g_core.zip')
    req = urlopen(url)
    with open(dest, 'w') as f:
        shutil.copyfileobj(req, f)
    subprocess.check_call(['unzip', '-d', outdir, dest])
    
def download_vcftools(outdir):
    """
    Download and compile vcftools.

    Parameters
    ----------
    outdir : str
        Directory to save vcftools to.

    """
    url = ('http://sourceforge.net/projects/vcftools/'
           'files/vcftools_0.1.12b.tar.gz/download')
    dest = os.path.join(outdir, 'vcftools_0.1.12b.tar.gz')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'vcftools_0.1.12b'))
    subprocess.check_call('make')
    os.chdir(cwd)

def download_subread(outdir):
    """
    Download Subread. Includes featureCounts.

    Parameters
    ----------
    outdir : str
        Directory to save Subread to.

    """
    url = ('http://sourceforge.net/projects/subread/files/subread-1.4.6/'
           'subread-1.4.6-Linux-x86_64.tar.gz/download')
    dest = os.path.join(outdir,
                        'subread-1.4.6-Linux-x86_64.tar.gz')
    _download_and_untar(url, dest, outdir)

def download_bcftools(outdir):
    """
    Download and compile bcftools.

    Parameters
    ----------
    outdir : str
        Directory to save bcftools to.

    """
    url = ('https://github.com/samtools/bcftools/releases/download/1.2/'
           'bcftools-1.2.tar.bz2')
    dest = os.path.join(outdir, 'bcftools-1.2.tar.bz2')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    edir = os.path.join(outdir, 'bcftools-1.2')
    os.chdir(edir)
    subprocess.check_call('make > make.out 2> make.err', shell=True)
    subprocess.check_call(('make prefix={} install > make_install.out 2> '
                           'make_install.err'.format(edir)), shell=True)
    os.chdir(cwd)

def download_htslib(outdir):
    """
    Download and compile htslib.

    Parameters
    ----------
    outdir : str
        Directory to save htslib to.

    """
    url = ('https://github.com/samtools/htslib/releases/download/1.2.1/'
           'htslib-1.2.1.tar.bz2')
    dest = os.path.join(outdir, 'htslib-1.2.1.tar.bz2')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    edir = os.path.join(outdir, 'htslib-1.2.1')
    os.chdir(edir)
    subprocess.check_call('make > make.out 2> make.err', shell=True)
    subprocess.check_call(('make prefix={} install > make_install.out 2> '
                           'make_install.err'.format(edir)), shell=True)
    os.chdir(cwd)

def download_samtools(outdir, lncurses=False):
    """
    Download and compile samtools.

    Parameters
    ----------
    outdir : str
        Directory to save samtools to.

    lncurses : bool
        Set to true to use lncurses rather than lcurses to build samtools. See
        http://seqanswers.com/forums/showthread.php?t=6669 for more information.

    """
    url = ('https://github.com/samtools/samtools/releases/download/1.2/'
           'samtools-1.2.tar.bz2')
    dest = os.path.join(outdir, 'samtools-1.2.tar.bz2')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    edir = os.path.join(outdir, 'samtools-1.2')
    os.chdir(edir)
    if lncurses:
        f = open('Makefile', 'r')
        lines = f.read().replace('lcurses', 'lncurses')
        f.close()
        f = open('Makefile', 'w')
        f.write(lines)
        f.close()

    subprocess.check_call('make > make.out 2> make.err', shell=True)
    subprocess.check_call(('make prefix={} install > make_install.out 2> '
                           'make_install.err'.format(edir)), shell=True)
    os.chdir(cwd)

def download_wgEncodeCrgMapabilityAlign100mer(outdir):
    """
    Download wgEncodeCrgMapabilityAlign100mer mappability bigwig track.

    Parameters
    ----------
    outdir : str
        Directory to save wgEncodeCrgMapabilityAlign100mer mappability bigwig
        track to.

    """
    req = urlopen(('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/'
                   'wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig'))
    dest = os.path.join(outdir, 'wgEncodeCrgMapabilityAlign100mer.bigWig')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    
def download_hg19(
    outdir,
    create_sorted=True,
    samtools_path='samtools',
    picard_path='$picard',
    twoBitToFa_path='twoBitToFa',
):
    """
    Download hg19, create index, and sequence dictionary files. 

    Parameters
    ----------
    outdir : str
        Directory to save hg19 fasta to.

    create_sorted : bool
        Create hg19_sorted.fa version that has chromosomes ordered
        lexicographically as expected by bedtools. The fasta is indexed and a
        sequence dictionary is created as well.

    samtools_path : str
        Path to Samtools executable needed to index fasta.

    """
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    req = urlopen(
        'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit')
    dest = os.path.join(outdir, 'hg19.2bit')

    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    cwd = os.getcwd()
    os.chdir(outdir)
    # req = urlopen(
    #     'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa')
    # dest = os.path.join(outdir, 'twoBitToFa')
    # with open(dest, 'w') as d:
    #     shutil.copyfileobj(req, d)
    # subprocess.check_call('chmod 755 {}'.format(os.path.join(outdir, 
    #                                                          'twoBitToFa')),
    #                       shell=True)

    # Convert from 2bit to fasta.
    subprocess.check_call('{} hg19.2bit hg19.fa'.format(twoBitToFa_path),
                          shell=True)
    # Index fasta.
    subprocess.check_call('{} faidx hg19.fa'.format(samtools_path), shell=True)
    # Create sequence dictionary.
    c = ('java -jar {} CreateSequenceDictionary R=hg19.fa '
         'O=hg19.dict'.format(picard_path))
    subprocess.check_call(c, shell=True)
    if create_sorted:
        # Sort fasta lexographically. This is the native ordering for bedtools
        # etc.
        c = ('grep \\> hg19.fa | cut -d ">" -f 2 | sort | '
             'awk \'{print $1".fa"}\' > order.txt')
        subprocess.check_call(c, shell=True)
        c = ('awk \'{if (substr($1, 1, 1) == ">") {x=substr($1, 2)}} '
             '{print $0 >> x".fa"}\' hg19.fa')
        subprocess.check_call(c, shell=True)
        c = ('xargs -d\'\\n\' cat <order.txt >> hg19_sorted.fa')
        subprocess.check_call(c, shell=True)
        subprocess.check_call('rm chr*.fa', shell=True)
        subprocess.check_call('rm order.txt', shell=True)
        # Index sorted fasta.
        subprocess.check_call('{} faidx hg19_sorted.fa'.format(samtools_path),
                              shell=True)
        # Create sequence dictionary.
        c = ('java -jar {} CreateSequenceDictionary R=hg19_sorted.fa '
             'O=hg19_sorted.dict'.format(picard_path))
        subprocess.check_call(c, shell=True)
    os.chdir(cwd)

def download_htsjdk(outdir):
    """
    Download STAR aligner.

    Parameters
    ----------
    outdir : str
        Directory to save STAR to.

    """
    url = 'https://github.com/samtools/htsjdk/tarball/master'
    dest = os.path.join(outdir, 'samtools-htsjdk-1.127-10-g18192d8.tar.gz')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'samtools-htsjdk-18192d8'))
    subprocess.check_call(['ant', 'htsjdk-jar'])
    os.chdir(cwd)

def download_star(outdir):
    """
    Download STAR aligner.

    Parameters
    ----------
    outdir : str
        Directory to save STAR to.

    """
    url = 'https://github.com/alexdobin/STAR/archive/STAR_2.4.0h.tar.gz'
    dest = os.path.join(outdir, 'STAR_2.4.0h.tar.gz')
    _download_and_untar(url, dest, outdir)

def make_star_index(
    outdir, 
    genome, 
    gtf, 
    threads=1,
    star_path='STAR'):
    """
    Make index for STAR aligner.

    Parameters
    ----------
    outdir : str
        Directory to save index to.

    genome : str
        Path to genome fasta file.

    gtf : str
        Path to GTF gene annotation.

    threads : int
        Number of threads to use to make the index.

    star_path : str
        Path to STAR executable.

    """
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    subprocess.check_call([star_path, '--runThreadN', str(threads),
                           '--runMode', 'genomeGenerate',
                           '--sjdbOverhang', '250', '--genomeDir',
                           outdir, '--genomeFastaFiles',
                           genome, '--sjdbGTFfile', gtf],
                          shell=False)
    cwd = os.getcwd()
    if os.path.realpath(cwd) != os.path.realpath(outdir):
        shutil.move('Log.out', outdir)

def download_picard(outdir):
    """
    Download Picard tools.

    Parameters
    ----------
    outdir : str
        Directory to save Picard tools to.

    """
    url = 'https://github.com/broadinstitute/picard/archive/1.131.tar.gz'
    dest = os.path.join(outdir, '1.131.tar.gz')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'picard-1.131'))
    subprocess.check_call('ant -lib lib/ant clone-htsjdk package-commands',
                          shell=True)
    os.chdir(cwd)

def download_bedtools(outdir): 
    """
    Download Bedtools.

    Parameters
    ----------
    outdir : str
        Directory to save Bedtools to.

    """
    url = ('https://github.com/arq5x/bedtools2/releases/download/v2.23.0/'
           'bedtools-2.23.0.tar.gz')
    dest = os.path.join(outdir, 'bedtools-2.23.0.tar.gz')
    _download_and_untar(url, dest, outdir)
    os.rename(os.path.join(outdir, 'bedtools2'),
              os.path.join(outdir, 'bedtools2-2.23.0'))
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'bedtools2-2.23.0'))
    subprocess.check_call('make')
    os.chdir(cwd)
    raw_input('\n\n\nYou should add\n' + 
              os.path.join(outdir, 'bedtools2-2.23.0', 'bin') + 
              '\nto your path when using this environment so\n'
              'pybedtools uses the correct bedtools installation.\n'
              'Press any key to continue.\n\n\n')

def download_r(outdir):
    """
    Download R.

    Parameters
    ----------
    outdir : str
        Directory to save R to.

    """
    rbase = 'R-3.1.1'
    url = 'http://cran.stat.ucla.edu/src/base/R-3/R-3.1.1.tar.gz'
    dest = os.path.join(outdir, '{}.tar.gz'.format(rbase))
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    os.chdir(outdir)
    shutil.move(rbase, '{}-source'.format(rbase))
    rpath = os.getcwd()
    os.chdir('{}-source'.format(rbase))
    subprocess.check_call(('./configure ' + 
                           '--enable-R-shlib ' + 
                           '--prefix={}/R-3.1.1 '.format(rpath) + 
                           '> R_configure.out ',
                           '2> R_configure.err'),
                          shell=True)
    subprocess.check_call('make > R_make.out 2> R_make.err', shell=True)
    # subprocess.check_call('make check > R_check.out 2> R_check.err', 
    #                       shell=True)
    subprocess.check_call('make install > R_install.out 2> R_install.err', 
                          shell=True)
    os.chdir(cwd)

                          
# # I execute the following python file when starting
# # the notebook server to set the library path:
# # 
# #     import os
# #     import subprocess
# # 
# #     import projectpy as ppy 
# # 
# #     spath = os.path.join(ppy.root, 'software')
# #     subprocess.check_call('export LDFLAGS="-Wl,-rpath,' + 
# #                           '{}/R-3.1.1/lib64/R/lib"'.format(spath),
# #                           shell=True)
# #     subprocess.check_call('export LD_LIBRARY_PATH="' + 
# #                           '{}/R-3.1.1/lib64/R/lib:$LD_LIBRARY_PATH"'.format(spath),
# #                           shell=True)
# 
# # In[9]:
# 
# try: 
#     import rpy2.robjects as r
# except:
#     cwd = os.getcwd()
#     os.chdir(os.path.join(ppy.root, 'software'))
#     spath = os.getcwd()
#     """
#     subprocess.check_call('export LDFLAGS="-Wl,-rpath,' + 
#                           '{}/R-3.1.1/lib64/R/lib"'.format(spath),
#                           shell=True)
#     subprocess.check_call('export LD_LIBRARY_PATH="' + 
#                           '{}/R-3.1.1/lib64/R/lib:$LD_LIBRARY_PATH"'.format(spath),
#                           shell=True)
#     """
#     req = urlopen('https://pypi.python.org/packages/'
#                   'source/r/rpy2/rpy2-2.4.2.tar.gz')
#     dest = os.path.join(ppy.root, 'software', 'rpy2-2.4.2.tar.gz')
#     with open(dest, 'w') as d:
#         shutil.copyfileobj(req, d)
#     subprocess.check_call('tar -xf {} -C {}'.format(dest, os.path.split(dest)[0]), 
#                           shell=True)
#     os.chdir('rpy2-2.4.2')
#     subprocess.check_call('python setup.py build --r-home ' + 
#                           '{}/R-3.1.1 install >& '.format(spath) + 
#                           'rpy2-2.4.2_R-3.1.1_log.txt', shell=True)
#     os.chdir(cwd)
# 
# 
# # Note that you may have to restart the notebook at this 
# # point if you are installing rpy2 for the first time so
# # that python knows the package exists.
# 
# # In[11]:
# 
# get_ipython().magic(u'load_ext rpy2.ipython')
# 
# 
# # # Bioconductor
# # 
# # I'm going to install Bioconductor and some packages
# # that I'll use.
# 
# # In[15]:
# 
# get_ipython().run_cell_magic(u'R', u'', u'\nsource("http://bioconductor.org/biocLite.R")\nbiocLite(ask=FALSE)\nbiocLite("DEXSeq", ask=FALSE)\nbiocLite("Gviz", ask=FALSE)\nbiocLite("BiocParallel", ask=FALSE)')

def download_install_rpy2(r_path, outdir):
    """
    Download and install rpy2. R must be installed and the LDFLAGS and
    LD_LIBRARY_PATH must be set. If they are not set, you can run the method to
    get instructions on how to set them.

    Parameters
    ----------
    r_path : str
        Path to R executable. The R shared library will be inferred based on
        this.

    """
    cwd = os.getcwd()
    # Path to R shared library.
    s_path = os.path.join(os.path.split(os.path.split(r_path)[0])[0], 'lib64',
                         'R', 'lib')
    lines = ('\n\n\nrpy2 has to be compiled against the version of R you\'ve \n'
             'installed here. Note that you have to set ld flags and library\n'
             'paths for installation. If you haven\'t already, paste the\n'
             'following at the command line:\n\n')
    sys.stdout.write(lines)
    command = 'export PATH="{}:$PATH"\n'.format(os.path.split(r_path)[0])
    sys.stdout.write(command)
    command = 'export LDFLAGS="-Wl,-rpath,{}"\n'.format(s_path)
    sys.stdout.write(command)
    command = 'export LD_LIBRARY_PATH="{}:$LD_LIBRARY_PATH"\n\n'.format(s_path)
    sys.stdout.write(command)
    raw_input('If these environment variables are set correctly, press any\n'
              'key to continue. Otherwise, exit and set them, then rerun.\n'
              'You need to set the PATH and LD_LIBRARY_PATH for each bash\n'
              'environment you want to use rpy2 in.\n\n')
    sys.stdout.write('Downloading and installing rpy2.\n\n')
    sys.stdout.flush()

    url = ('https://pypi.python.org/packages/source/r/rpy2/rpy2-2.7.3.tar.gz')
    dest = os.path.join(outdir, 'rpy2-2.7.3.tar.gz')
    _download_and_untar(url, dest, outdir)
    os.chdir(os.path.join(outdir, 'rpy2-2.7.3'))
    r_home = os.path.split(os.path.split(r_path)[0])[0]
    subprocess.check_call('python setup.py build --r-home ' + 
                          '{} install >& '.format(r_home) + 
                          'rpy2_install_log.txt', shell=True)
    os.chdir(cwd)

def download_gencode_gtf(outdir):
    """
    Download Gencode GTF.

    Parameters
    ----------
    outdir : str
        Directory to save Gencode GTF to.

    """
    src = ('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/'
           'gencode.v19.annotation.gtf.gz')
    dest = os.path.join(outdir, 'gencode.v19.annotation.gtf.gz')
    _download_and_gunzip(src, dest)

def download_kheradpour_motifs(outdir):
    req = urlopen('http://compbio.mit.edu/encode-motifs/motifs.txt')
    dest = os.path.join(outdir, 'motifs.txt')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)

def download_gtfToGenePred(outdir):
    req = urlopen('http://hgdownload.cse.ucsc.edu/admin/exe/'
                  'linux.x86_64/gtfToGenePred')
    dest = os.path.join(outdir, 'gtfToGenePred')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['chmod', '755', '{}'.format(dest)])

def download_liftOver(outdir):
    req = urlopen('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/'
                  'liftOver')
    dest = os.path.join(outdir, 'liftOver')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['chmod', '755', '{}'.format(dest)])

def download_bigWigAverageOverBed(outdir):
    req = urlopen('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/'
                  'bigWigAverageOverBed')
    dest = os.path.join(outdir, 'bigWigAverageOverBed')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['chmod', '755', '{}'.format(dest)])

def download_ucsc_tools(outdir):
    """Download all UCSC executables except blat."""
    import re
    url = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64'
    s = urlopen(url).read()
    
    compiled = re.compile('<a href=".*">')
    res = compiled.findall(s)
    res = [x.split('"')[-2] for x in res[2:]]
    res = [x for x in res if '/' not in x]
    to_download = ['{}/{}'.format(url, x) for x in res]
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        _download_file(src, dest)
        subprocess.check_call(['chmod', '755', '{}'.format(dest)])

def download_bedGraphToBigWig(outdir):
    req = urlopen('http://hgdownload.cse.ucsc.edu/admin/exe/'
                  'linux.x86_64/bedGraphToBigWig')
    dest = os.path.join(outdir, 'bedGraphToBigWig')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['chmod', '755', '{}'.format(dest)])

def install_bioconductor_dependencies():
    """
    Installs bioconductor and some bioconductor packages into R using rpy2.
    Packages installed are currently DESeq2 and DEXSeq.

    """
    try:
        # Have to import readline due to some weirdness with installing rpy2 in
        # a conda environment etc.
        # See https://github.com/ContinuumIO/anaconda-issues/issues/152.
        import readline
        import rpy2.robjects as robjects
    except ImportError:
        sys.stderr.write('rpy2 not installed.\n')
        sys.exit(1)
    robjects.r('source("http://bioconductor.org/biocLite.R")')
    robjects.r('biocLite(ask=FALSE)')
    robjects.r('biocLite("DESeq2", ask=FALSE)')
    robjects.r('biocLite("DEXSeq", ask=FALSE)')
    robjects.r('biocLite("Gviz", ask=FALSE)')
    robjects.r('biocLite("goseq", ask=FALSE)')
    robjects.r('biocLite("org.Hs.eg.db", ask=FALSE)')
    robjects.r('biocLite("qvalue", ask=FALSE)')

def make_dexseq_annotation(gtf, out_gff):
    try:
        # Have to import readline due to some weirdness with installing rpy2 in
        # a conda environment etc.
        # See https://github.com/ContinuumIO/anaconda-issues/issues/152.
        import readline
        import rpy2.robjects as robjects
    except ImportError:
        sys.stderr.write('rpy2 not installed.\n')
        sys.exit(1)
    robjects.r('library(DEXSeq)')
    scripts = robjects.r('system.file("python_scripts", package="DEXSeq")')
    g = scripts.items()
    scripts_path = g.next()[1]
    script = os.path.join(scripts_path, 'dexseq_prepare_annotation.py')
    command = 'python {} -r no {} {}'.format(script, gtf, out_gff)
    subprocess.check_call(command, shell=True)

def rsem_prepare_reference(
    fasta, 
    name, 
    gtf=None,
    rsem_prepare_reference_path='rsem-prepare-reference',
):
    """
    Run rsem-prepare-reference to make an RSEM reference.

    Parameters
    ----------
    fasta : str
        Path to fasta to provide to RSEM. Should be a transcriptome fasta unless
        gtf is specified. If gtf is specified, then this is a genome fasta.

    name : str
        Name of the RSEM reference.

    rsem_path : str
        Path to directory with RSEM executables.

    gtf : str
        Path to GTF file that defines genes and transcripts to provide to RSEM.

    """
    command = '{} {} {}'.format(rsem_prepare_reference_path, fasta, name)
    if gtf:
        command += ' --gtf {}'.format(gtf)
    subprocess.check_call(command, shell=True)

def download_roadmap_gwas(outdir):
    """
    Download GWAS data compiled for Roadmap project.
    """
    import re
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    url = 'http://www.broadinstitute.org/~anshul/roadmap/gwas'
    s = urlopen(url).read()
    compiled = re.compile('"EUR.*txt\.gz"')
    res = [x.strip('"') for x in compiled.findall(s)]
    to_download = ['{}/{}'.format(url, x) for x in res]
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        _download_and_gunzip(src, dest)

def download_roadmap_25_state_chromatin_model(outdir):
    """
    Download 25 state chromatin model from Roadmap Epigenomics. There is a bed
    file for each cell type as well as an annotation file with state information
    and file to convert the roadmap ID's to human readable cell types. Bed files
    are sorted so they work with bedtools -sorted option.

    Parameters
    ----------
    outdir : str
        Directory to save files to.

    """
    import re
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    pattern = 'href="E\d\d\d_25_imputed12marks_mnemonics.bed.gz"'
    compiled = re.compile(pattern)
    url = ('http://egg2.wustl.edu/roadmap/data/byFileType'
           '/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final')
    s = urlopen(url).read()
    res = compiled.findall(s)
    assert len(res) is 127
    res = [x[5:].strip('"') for x in res]
    to_download = ['{}/{}'.format(url, x) for x in res]
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        sorted_dest = '{}_sorted.bed'.format(dest.strip('.bed.gz'))
        _download_and_gunzip(src, dest)
        subprocess.check_call('sort -k 1,1 -k2,2n {} > {}'.format(
            dest.strip('.gz'), sorted_dest), shell=True)
        os.remove(dest.strip('.gz'))
    to_download = []
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/imputed12marks/'
                       'jointModel/final/EIDlegend.txt')
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/imputed12marks/'
                       'jointModel/final/annotation_25_imputed12marks.txt')
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        _download_file(src, dest)
    # Now I'll make a more useful annotation file that has better column names
    # and a few extra columns.
    import numpy as np
    import pandas as pd
    state_annot = pd.read_table(
        os.path.join(outdir, 'annotation_25_imputed12marks.txt.tsv'),
        index_col=0)
    state_annot.index.name = 'state_number'
    state_annot.columns = [x.lower().replace(' ', '_') for x in
                           state_annot.columns]
    state_annot['name'] = (state_annot.index.astype(str).values + '_' +
                           state_annot.mnemonic)
    state_annot['color_code_norm'] = state_annot['color_code'].apply(
            lambda x: np.array([int(x) for x in x.split(',')]))
    state_annot['color_code_norm'] = state_annot['color_code_norm'].apply(
        lambda x: x / 255.)
    state_annot['color_code_norm'] = state_annot.color_code_norm.apply(
        lambda x: ','.join([str(y) for y in x]))
    state_annot.to_csv(os.path.join(outdir, 'frazer_parsed_annotation.tsv'),
                       sep='\t')

def download_roadmap_18_state_chromatin_model(outdir):
    """
    Download 18 state chromatin model from Roadmap Epigenomics. There is a bed
    file for each cell type as well as an annotation file with state information
    and file to convert the roadmap ID's to human readable cell types. Bed files
    are sorted so they work with bedtools -sorted option.

    Parameters
    ----------
    outdir : str
        Directory to save files to.

    """
    import re
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    pattern = 'href="E\d\d\d_18_core_K27ac_mnemonics.bed.gz"'
    compiled = re.compile(pattern)
    url = ('http://egg2.wustl.edu/roadmap/data/byFileType/'
           'chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final')
    s = urlopen(url).read()
    res = compiled.findall(s)
    res = [x[5:].strip('"') for x in res]
    to_download = ['{}/{}'.format(url, x) for x in res]
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        sorted_dest = '{}_sorted.bed'.format(dest.strip('.bed.gz'))
        _download_and_gunzip(src, dest)
        subprocess.check_call('sort -k 1,1 -k2,2n {} > {}'.format(
            dest.strip('.gz'), sorted_dest), shell=True)
        os.remove(dest.strip('.gz'))
    to_download = []
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/core_K27ac/'
                       'jointModel/final/browserlabelmap_18_core_K27ac.tab')
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/'
                       'final/colormap_18_core_K27ac.tab')
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/imputed12marks/'
                       'jointModel/final/EIDlegend.txt')
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        _download_file(src, dest)
    # The 18 state model doesn't have a nice annotation file like the 25 state
    # model (although the table is on the website), so I'm just putting the info
    # in here and I'll write the file myself.
    columns = [u'STATE NO.', u'MNEMONIC', u'DESCRIPTION', u'COLOR NAME',
               u'COLOR CODE']
    vals = [[1, 'TssA', 'Active TSS', 'Red', '255,0,0'],
       [2, 'TssFlnk', 'Flanking TSS', 'Orange Red', '255,69,0'],
       [3, 'TssFlnkU', 'Flanking TSS Upstream', 'Orange Red', '255,69,0'],
       [4, 'TssFlnkD', 'Flanking TSS Downstream', 'Orange Red', '255,69,0'],
       [5, 'Tx', 'Strong transcription', 'Green', '0,128,0'],
       [6, 'TxWk', 'Weak transcription', 'DarkGreen', '0,100,0'],
       [7, 'EnhG1', 'Genic enhancer1', 'GreenYellow', '194,225,5'],
       [8, 'EnhG2', 'Genic enhancer2', 'GreenYellow', '194,225,5'],
       [9, 'EnhA1', 'Active Enhancer 1', 'Orange', '255,195,77'],
       [10, 'EnhA2', 'Active Enhancer 2', 'Orange', '255,195,77'],
       [11, 'EnhWk', 'Weak Enhancer', 'Yellow', '255,255,0'],
       [12, 'ZNF/Rpts', 'ZNF genes & repeats', 'Medium Aquamarine',
        '102,205,170'],
       [13, 'Het', 'Heterochromatin', 'PaleTurquoise', '138,145,208'],
       [14, 'TssBiv', 'Bivalent/Poised TSS', 'IndianRed', '205,92,92'],
       [15, 'EnhBiv', 'Bivalent Enhancer', 'DarkKhaki', '189,183,107'],
       [16, 'ReprPC', 'Repressed PolyComb', 'Silver', '128,128,128'],
       [17, 'ReprPC', 'Weak Repressed PolyComb', 'Gainsboro', '192,192,192'],
       [18, 'Quies', 'Quiescent/Low', 'White', '255,255,255']]
    import numpy as np
    import pandas as pd
    df = pd.DataFrame(vals, columns=columns)
    df.to_csv(os.path.join(outdir, 'frazer_annotation.tsv'), index=None,
              sep='\t')
    # Now I'll make a more useful annotation file that has better column names
    # and a few extra columns.
    state_annot = pd.read_table(os.path.join(outdir, 'frazer_annotation.tsv'),
                                index_col=0)
    state_annot.index.name = 'state_number'
    state_annot.columns = [x.lower().replace(' ', '_') for x in
                           state_annot.columns]
    state_annot['name'] = (state_annot.index.astype(str).values + '_' +
                           state_annot.mnemonic)
    state_annot['color_code_norm'] = state_annot['color_code'].apply(
            lambda x: np.array([int(x) for x in x.split(',')]))
    state_annot['color_code_norm'] = state_annot['color_code_norm'].apply(
        lambda x: x / 255.)
    state_annot['color_code_norm'] = state_annot.color_code_norm.apply(
        lambda x: ','.join([str(y) for y in x]))
    state_annot.to_csv(os.path.join(outdir, 'frazer_parsed_annotation.tsv'),
              sep='\t')

def download_roadmap_15_state_chromatin_model(outdir):
    """
    Download 15 state chromatin model from Roadmap Epigenomics. There is a bed
    file for each cell type as well as an annotation file with state information
    and file to convert the roadmap ID's to human readable cell types. Bed files
    are sorted so they work with bedtools -sorted option.

    Parameters
    ----------
    outdir : str
        Directory to save files to.

    """
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    src = ('http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations'
           '/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz')
    dest = os.path.join(outdir, os.path.split(src)[1])
    _download_and_untar(src, dest, outdir, remove_tarball=True)
    import glob
    fns = glob.glob(os.path.join(outdir, '*bed.gz'))
    for fn in fns:
        subprocess.check_call(['gunzip', fn])
        sorted_fn = '{}_sorted.bed'.format(fn.strip('.bed.gz'))
        subprocess.check_call('sort -k 1,1 -k2,2n {} > {}'.format(
            fn.strip('.gz'), sorted_fn), shell=True)
        os.remove(fn.strip('.gz'))
        
    to_download = []
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/imputed12marks/'
                       'jointModel/final/EIDlegend.txt')
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/coreMarks/jointModel'
                       '/final/labelmap_15_coreMarks.tab')
    to_download.append('http://egg2.wustl.edu/roadmap/data/byFileType/'
                       'chromhmmSegmentations/ChmmModels/coreMarks/jointModel'
                       '/final/colormap_15_coreMarks.tab')
    for src in to_download:
        dest = os.path.join(outdir, os.path.split(src)[1])
        _download_file(src, dest)
    # The 15 state model doesn't have a nice annotation file like the 25 state
    # model (although the table is on the website), so I'm just putting the info
    # in here and I'll write the file myself.
    columns = [u'STATE NO.', u'MNEMONIC', u'DESCRIPTION', u'COLOR NAME', 
               u'COLOR CODE']
    vals = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 
            ['TssA', 'TssAFlnk', 'TxFlnk', 'Tx', 'TxWk', 'EnhG', 'Enh',
             'ZNF/Rpts', 'Het', 'TssBiv', 'BivFlnk', 'EnhBiv', 'ReprPC',
             'ReprPCWk', 'Quies'],
            ['Active TSS', 'Flanking Active TSS', "Transcr. at gene 5' and 3'",
             'Strong transcription', 'Weak transcription', 'Genic enhancers',
             'Enhancers', 'ZNF genes & repeats', 'Heterochromatin',
             'Bivalent/Poised TSS', 'Flanking Bivalent TSS/Enh', 
             'Bivalent Enhancer', 'Repressed PolyComb', 
             'Weak Repressed PolyComb', 'Quiescent/Low'],
            ['Red', 'Orange Red', 'LimeGreen', 'Green', 'DarkGreen',
             'GreenYellow', 'Yellow', 'Medium Aquamarine', 'PaleTurquoise',
             'IndianRed', 'DarkSalmon', 'DarkKhaki', 'Silver', 'Gainsboro',
             'White'],
            ['255,0,0', '255,69,0', '50,205,50', '0,128,0', '0,100,0',
             '194,225,5', '255,255,0', '102,205,170', '138,145,208',
             '205,92,92', '233,150,122', '189,183,107', '128,128,128',
             '192,192,192', '255,255,255']
           ]
    import numpy as np
    import pandas as pd
    df = pd.DataFrame(vals, index=columns).T
    df.to_csv(os.path.join(outdir, 'frazer_annotation.tsv'), index=None,
              sep='\t')
    # Now I'll make a more useful annotation file that has better column names
    # and a few extra columns.
    state_annot = pd.read_table(os.path.join(outdir, 'frazer_annotation.tsv'),
                                index_col=0)
    state_annot.index.name = 'state_number'
    state_annot.columns = [x.lower().replace(' ', '_') for x in
                           state_annot.columns]
    state_annot['name'] = (state_annot.index.astype(str).values + '_' +
                           state_annot.mnemonic)
    state_annot['color_code_norm'] = state_annot['color_code'].apply(
            lambda x: np.array([int(x) for x in x.split(',')]))
    state_annot['color_code_norm'] = state_annot['color_code_norm'].apply(
        lambda x: x / 255.)
    state_annot['color_code_norm'] = state_annot.color_code_norm.apply(
        lambda x: ','.join([str(y) for y in x]))
    state_annot.to_csv(os.path.join(outdir, 'frazer_parsed_annotation.tsv'),
                       sep='\t')
