import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'pipelines',
    version = '0.1.0',
    author = 'Christopher DeBoever',
    author_email = 'cdeboever3@gmail.com',
    description = ('Various computational pipelines from the Frazer Lab.'),
    packages=find_packages(),
    license = 'MIT',
    keywords = 'bioinformatics pipeline',
    url = 'https://github.com/frazer-lab/pipelines',
    long_description=read('README.md'),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
   ]
)
