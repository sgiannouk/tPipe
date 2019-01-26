import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='tPipe',
    version='0.2.0',
    description='Package for detecting tRNA and tsRNA molecules in small RNA-seq data.',
    long_description=read('README.md'),
    author='Stavros Giannoukakos',
    author_email='s.p.giannoukakos@hotmail.com',
    packages=['tPipe'],
    keywords=['tRNA tsRNA RNA-Seq analysis'],
    install_requires=['xlsxwriter', 'xlrd', 'natsort', 'multiqc', 'cutadapt','numpy'],
    dependency_links=['https://github.com/agordon/fastx_toolkit',
                      'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5_source.zip',
                      'https://sourceforge.net/projects/bowtie-bio/files/latest/download?source=files'],
    package_data = {
        'tPipe': ['configuration_file.txt'],
        'tPipe': ['tLibraries/reflibs_mu.fa'],
        'tPipe': ['tLibraries/sslib_mu.fa'],
        'tPipe': ['tLibraries/p.fa'],
        'tPipe': ['tLibraries/m.fa']
        },
    entry_points={
          'console_scripts': ['tpipe = tPipe.tPipe:main']
      },
)





