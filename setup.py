#!/usr/bin/python
"""
ST Pipeline is a tool to process the Spatial Transcriptomics raw data.
The data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers. 
The output contains the counts table, a stats file, a log file
and a BED file with all the transcripts.
"""

import os
import io
import glob
from setuptools import setup, find_packages
from stpipeline.version import version_number

# Get the long description from the relevant file
here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []
    
setup(
  name = 'stpipeline',
  version = version_number,
  description = __doc__.split("\n", 1)[0],
  long_description = long_description,
  keywords = 'rna-seq analysis spatial transcriptomics toolkit',
  author = 'Jose Fernandez Navarro',
  author_email = 'jose.fernandez.navarro@scilifelab.se',
  license = 'BSD',
  url = 'https://github.com/SpatialTranscriptomicsResearch/st_pipeline',
  packages = find_packages(exclude=('tests*', 'dependencies', 'utils')),
  include_package_data = False,
  package_data = {'': ['RELEASE-VERSION']},
  zip_safe = False,
  install_requires = install_requires,
  test_suite = 'tests',
  scripts = glob.glob('scripts/*.py'),
  classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: BSD :: Copyright Jose Fernandez Navarro, KTH, KI',
    'Programming Language :: Python :: 2.7',
    'Environment :: Console',
  ],
)
