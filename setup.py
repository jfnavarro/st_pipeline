#!/usr/bin/python
"""
ST Pipeline is a tool to process Spatial Transcriptomics raw data (or 
any type of single cell data whose raw data has the same configuraiton).
The data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers (UMIs). 
The output contains a counts table, a stats file, a log file
and a BED file with all the transcripts.
"""

import os
import io
import glob
import sys
from setuptools import setup, find_packages
from stpipeline.version import version_number
from distutils.core import setup, Extension

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
try:
    with io.open(os.path.join(here, 'README_SHORT'), encoding='utf-8') as f:
        long_description = f.read()
except IOError:
    long_description = ""

try:
    with open(os.path.join(here, "requirements.txt"), "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    raise SystemExit("Could not find requirements.txt file")

major, minor1, minor2, s, tmp = sys.version_info
if major != 3 or minor1 < 6:
    raise SystemExit("ST Pipeline requires Python 3.6 or bigger")

# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
# if 'setuptools.extension' in sys.modules:
#    m = sys.modules['setuptools.extension']
#    m.Extension.__dict__ = m._Extension.__dict__

setup(
    name='stpipeline',
    version=version_number,
    description="ST Pipeline: An automated pipeline for spatial mapping of unique transcripts",
    long_description=long_description,
    keywords='rna-seq analysis spatial transcriptomics toolkit',
    author='Jose Fernandez Navarro',
    author_email='jc.fernandez.navarro@gmail.com',
    license='MIT',
    url='https://github.com/SpatialTranscriptomicsResearch/st_pipeline',
    packages=find_packages(exclude=('tests*', 'utils', '*.pyx')),
    ext_modules=[
        Extension('stpipeline.common.cdistance', ['stpipeline/common/cdistance.pyx']),
        Extension('stpipeline.common.unique_events_parser', ['stpipeline/common/unique_events_parser.pyx']),
        Extension('stpipeline.common.filterInputReads', ['stpipeline/common/filterInputReads.pyx'])
    ],
    include_package_data=True,
    zip_safe=False,
    setup_requires=['setuptools_cython'],
    install_requires=install_requires,
    test_suite='tests',
    scripts=glob.glob('scripts/*.py'),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Environment :: Console',
    ],
)
