from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "snpiphy",
    version = "0.1",
    author = "Daniel Bogema",
    author_email = "daniel.bogema@dpi.nsw.gov.au",
    description = ("An automated snp phylogeny pipeline"),
    license = "GPL-3.0",
    keywords = "genomics phylogenetics snp",
    url = "https://github.com/bogemad/snpiphy",
    packages=['snpiphy'],
    scripts=['bin/snpiphy'],
    install_requires=['biopython>=1.71','bcbio-gff>=0.6.4', 'pandas>=0.23.0', 'numpy>=1.14.3'],
    long_description=read('README.md'),
)
