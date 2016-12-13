import os
from setuptools import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst', format='md')
except(IOError, ImportError):
    long_description = open('README.md').read()

ep = {
    'console_scripts': ['convert_bed_to_saf = cdpipelines.convert_bed_to_saf:main',
                        'make_mbased_input = cdpipelines.make_mbased_input:main',
                        'make_wasp_input = cdpipelines.make_wasp_input:main',
                        'scale_bedgraph = cdpipelines.scale_bedgraph:main']
}

setup(
    name = 'cdpipelines',
    packages=['cdpipelines'],
    entry_points=ep,
    version = '0.0.7',
    author = 'Christopher DeBoever',
    author_email = 'cdeboever3@gmail.com',
    description = ('Various bioinformatics pipelines.'),
    license = 'MIT',
    keywords = ['bioinformatics'],
    url = 'https://github.com/cdeboever3/cdpipelines',
    download_url = 'https://github.com/cdeboever3/cdpipelines/tarball/0.0.7',
    include_package_data=True,
    long_description=long_description,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
   ]
)
