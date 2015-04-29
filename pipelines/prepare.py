# This notebook has commands to download some external data and software that
# are needed in various notebooks. Rather than saving these into the `output`
# directory with the output from other notebooks, I'll save them in the
# directories `software` and `external_data` to denote that they are just
# downloaded files, not files that I created. I also won't store these on
# figshare.

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

def download_rsem(outdir, lncurses=False):
    """
    Download RSEM.

    Parameters
    ----------
    outdir : str
        Directory to save RSEM.

    lncurses : bool
        Set to true to use lncurses rather than lcurses to build samtools. See
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

def download_samtools(outdir):
    """
    Download Samtools.

    Parameters
    ----------
    outdir : str
        Directory to save Samtools to.

    """
    url = ('http://sourceforge.net/projects/samtools/files/samtools/1.0/'
           'samtools-bcftools-htslib-1.0_x64-linux.tar.bz2/download')
    dest = os.path.join(outdir,
                        'samtools-bcftools-htslib-1.0_x64-linux.tar.bz2')
    _download_and_untar(url, dest, outdir)
    
def download_hg19(outdir, samtools_path):
    """
    Download hg19.

    Parameters
    ----------
    outdir : str
        Directory to save hg19 fasta to.

    samtools_path : str
        Path to Samtools executable needed to index fasta.

    """
    outdir = os.path.join(outdir, 'hg19')
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    req = urlopen(
        'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit')
    dest = os.path.join(outdir, 'hg19.2bit')

    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    req = urlopen(
        'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa')
    dest = os.path.join(outdir, 'twoBitToFa')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call('chmod 755 {}'.format(os.path.join(outdir, 
                                                             'twoBitToFa')),
                          shell=True)
    subprocess.check_call('{} {} {}'.format(
        os.path.join(outdir, 'twoBitToFa'), 
        os.path.join(outdir, 'hg19.2bit'), 
        os.path.join(outdir, 'hg19.fa')),
                          shell=True)    
    subprocess.check_call('{} faidx {}'.format(
        samtools_path,
        os.path.join(outdir, 'hg19.fa')),
                          shell=True)

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

def make_star_index(outdir, threads, genome, gtf, star_path='STARstatic'):
    """
    Make index for STAR aligner.

    Parameters
    ----------
    outdir : str
        Directory to save index to.

    threads : int
        Number of threads to use to make the index.

    genome : str
        Path to genome fasta file.

    gtf : str
        Path to GTF gene annotation.

    star_path : str
        Path to STAR executable.

    """
    dest = os.path.join(outdir, 'star_index')
    try:
        os.makedirs(dest)
    except OSError:
        pass
    subprocess.check_call([star_path, '--runThreadN', str(threads),
                           '--runMode', 'genomeGenerate',
                           '--sjdbOverhang', '250', '--genomeDir',
                           dest, '--genomeFastaFiles',
                           genome, '--sjdbGTFfile', gtf],
                          shell=False)
    shutil.move('Log.out', dest)

def download_picard(outdir):
    """
    Download Picard tools.

    Parameters
    ----------
    outdir : str
        Directory to save Picard tools to.

    """
    url = 'https://github.com/broadinstitute/picard/archive/1.128.tar.gz'
    dest = os.path.join(outdir, '1.128.tar.gz')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'picard-1.128'))
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
    url = ('https://github.com/arq5x/bedtools2/releases/'
           'download/v2.20.1/bedtools-2.20.1.tar.gz')
    dest = os.path.join(outdir, 'bedtools-2.20.1.tar.gz')
    _download_and_untar(url, dest, outdir)
    cwd = os.getcwd()
    os.chdir(os.path.join(outdir, 'bedtools2-2.20.1'))
    subprocess.check_call('make')
    os.chdir(cwd)
    raw_input('\n\n\nYou should add\n' + 
              os.path.join(outdir, 'bedtools2-2.20.1', 'bin') + 
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
    subprocess.check_call('make check > R_check.out 2> R_check.err', 
                          shell=True)
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

    url = ('https://pypi.python.org/packages/source/r/rpy2/rpy2-2.4.2.tar.gz')
    dest = os.path.join(outdir, 'rpy2-2.4.2.tar.gz')
    _download_and_untar(url, dest, outdir)
    os.chdir(os.path.join(outdir, 'rpy2-2.4.2'))
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
    dest = os.path.join(outdir, 'gencode_v19', 'gencode.v19.annotation.gtf.gz')
    _download_and_gunzip(src, dest)

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

def rsem_prepare_reference(fasta, name, rsem_path, gtf=None):
    """
    Run rsem-prepare-reference

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
    command = '{}/rsem-prepare-reference {} {}'.format(rsem_path, fasta, name)
    if gtf:
        command += ' --gtf {}'.format(gtf)
    subprocess.check_call(command, shell=True)

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
