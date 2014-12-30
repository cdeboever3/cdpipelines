# This notebook has commands to download some external data and software that
# are needed in various notebooks. Rather than saving these into the `output`
# directory with the output from other notebooks, I'll save them in the
# directories `software` and `external_data` to denote that they are just
# downloaded files, not files that I created. I also won't store these on
# figshare.

import os
import shutil
import subprocess
from urllib2 import urlopen

def _download_and_unzip(req, dest, out_dir):
    """
    Download a zipped file/directory req into out_dir and decompress it in
    out_dir.

    Parameters
    ----------
    req : str
        Zipped file.

    dest : str
        Full path to save zipped file to.

    out_dir : str
        Directory to save zipped file to and decompress to.

    """
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call('tar -xf {} -C {}'.format(dest, out_dir), shell=True)

def _download_and_untar(req, dest, out_dir):
    """
    Download a tarball req into out_dir and decompress it in out_dir.

    Parameters
    ----------
    req : str
        Tarball to download.

    dest : str
        Full path to save tarball to.

    out_dir : str
        Directory to save tarball to and decompress to.

    """
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call('tar -xf {} -C {}'.format(dest, out_dir), shell=True)

def download_samtools(out_dir):
    """
    Download Samtools.

    Parameters
    ----------
    out_dir : str
        Directory to save Samtools to.

    """
    req = urlopen('http://sourceforge.net/projects/samtools/'
                  'files/samtools/1.0/'
                  'samtools-bcftools-htslib-1.0_x64-linux.tar.bz2/download')
    dest = os.path.join(out_dir,
                        'samtools-bcftools-htslib-1.0_x64-linux.tar.bz2')
    _download_and_untar(req, dest, out_dir)
    
def download_hg19(out_dir, samtools_path):
    """
    Download hg19.

    Parameters
    ----------
    out_dir : str
        Directory to save hg19 fasta to.

    samtools_path : str
        Path to Samtools executable needed to index fasta.

    """
    out_dir = os.path.join(out_dir, 'hg19')
    try:
        os.makedirs(out_dir)
    except OSError:
        pass
    req = urlopen(
        'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit')
    dest = os.path.join(out_dir, 'hg19.2bit')

    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    req = urlopen(
        'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa')
    dest = os.path.join(out_dir, 'twoBitToFa')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call('chmod 755 {}'.format(os.path.join(out_dir, 
                                                             'twoBitToFa')),
                          shell=True)
    subprocess.check_call('{} {} {}'.format(
        os.path.join(out_dir, 'twoBitToFa'), 
        os.path.join(out_dir, 'hg19.2bit'), 
        os.path.join(out_dir, 'hg19.fa')),
                          shell=True)    
    subprocess.check_call('{} faidx {}'.format(
        samtools_path,
        os.path.join(out_dir, 'hg19.fa')),
                          shell=True)

def download_star(out_dir):
    """
    Download STAR aligner.

    Parameters
    ----------
    out_dir : str
        Directory to save STAR to.

    """
    req = urlopen('http://it-collab01.cshl.edu/shares/gingeraslab/www-data/'
                  'dobin/STAR/STARreleases/Patches/STAR_2.3.1z15.tgz')
    dest = os.path.join(out_dir, 'STAR_2.3.1z15.tgz')
    _download_and_untar(req, dest, out_dir)

def make_star_index(out_dir, threads, genome, gtf, star_path='STARstatic'):
    """
    Make index for STAR aligner.

    Parameters
    ----------
    out_dir : str
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
    dest = os.path.join(out_dir, 'star_index')
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

def download_picard(out_dir):
    """
    Download Picard tools.

    Parameters
    ----------
    out_dir : str
        Directory to save Picard tools to.

    """
    req = urlopen('http://sourceforge.net/projects/picard/files/picard-tools/'
                  '1.118/picard-tools-1.118.zip/download')
    dest = os.path.join(out_dir, 'picard-tools-1.118.zip')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['unzip', '-d', out_dir, dest])

def download_bedtools(out_dir): 
    """
    Download Bedtools.

    Parameters
    ----------
    out_dir : str
        Directory to save Bedtools to.

    """
    req = urlopen('https://github.com/arq5x/bedtools2/releases/'
                  'download/v2.20.1/bedtools-2.20.1.tar.gz')
    dest = os.path.join(out_dir, 'bedtools-2.20.1.tar.gz')
    _download_and_untar(req, dest, out_dir)
    cwd = os.getcwd()
    os.chdir(os.path.join(out_dir, 'bedtools2-2.20.1'))
    subprocess.check_call('make')
    os.chdir(cwd)

def download_r(out_dir):
    """
    Download R.

    Parameters
    ----------
    out_dir : str
        Directory to save R to.

    """
    rbase = 'R-3.1.1'
    req = urlopen('http://cran.stat.ucla.edu/src/base/R-3/R-3.1.1.tar.gz')
    dest = os.path.join(out_dir, '{}.tar.gz'.format(rbase))
    _download_and_untar(req, dest, out_dir)
    cwd = os.getcwd()
    os.chdir(out_dir)
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

# # # rpy2
# # 
# # `rpy2` has to be compiled against the version of 
# # R we've installed here. Note that you have to set 
# # ld flags and library paths for installation. I do 
# # this with the following:
# # 
# #     spath = os.getcwd()
# #     subprocess.check_call('export LDFLAGS="-Wl,-rpath,' + 
# #                           '{}/R-3.1.1/lib64/R/lib"'.format(spath),
# #                           shell=True)
# #     subprocess.check_call('export LD_LIBRARY_PATH="' + 
# #                           '{}/R-3.1.1/lib64/R/lib:$LD_LIBRARY_PATH"'.format(spath),
# #                           shell=True)
# #                           
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

def download_gencode_gtf(out_dir):
    """
    Download Gencode GTF.

    Parameters
    ----------
    out_dir : str
        Directory to save Gencode GTF to.

    """
    src = ('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/'
           'gencode.v19.annotation.gtf.gz')
    dest = os.path.join(out_dir, 'gencode_v19', 'gencode.v19.annotation.gtf.gz')
    try:
        os.makedirs(os.path.split(dest[0]))
    except OSError:
        pass
    req = urlopen(src)
    with open(dest, 'w') as f:
        shutil.copyfileobj(req, f)
    subprocess.check_call(['gunzip', dest])


def download_bedGraphToBigWig(out_dir):
    req = urlopen('http://hgdownload.cse.ucsc.edu/admin/exe/'
                  'linux.x86_64/bedGraphToBigWig')
    dest = os.path.join(out_dir, 'bedGraphToBigWig')
    with open(dest, 'w') as d:
        shutil.copyfileobj(req, d)
    subprocess.check_call(['chmod', '755', '{}'.format(dest)])
