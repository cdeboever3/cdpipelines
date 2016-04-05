import os
from setuptools import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst', format='md')
except(IOError, ImportError):
    long_description = open('README.md').read()

# I need to add a download_url before pushing to pypi:
# 
setup(
    name = 'cdpipelines',
    packages=['cdpipelines'],
    version = '0.0.1',
    author = 'Christopher DeBoever',
    author_email = 'cdeboever3@gmail.com',
    description = ('Various bioinformatics pipelines.'),
    license = 'MIT',
    keywords = ['bioinformatics'],
    url = 'https://github.com/cdeboever3/cdpipelines',
    download_url = 'https://github.com/cdeboever3/cdpipelines/tarball/0.0.1',
    long_description=long_description,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
   ]
)
