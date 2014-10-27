#!/usr/bin/python
"""
ST Pipeline is a toolkit of RNA-Seq data analysis for the Spatial Transcriptomics data.
"""

"""Based on https://github.com/pypa/sampleproject/blob/master/setup.py."""
import os
import sys
import io
import glob
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

# This is a plug-in for setuptools that will invoke py.test
# when you run python setup.py test
class PyTest(TestCommand):
    """Set up the py.test test runner."""
    
    def finalize_options(self):
        """Set options for the command line."""
        TestCommand.finalize_options(self)
        maintestfile = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + "/test/test_class.py"
        self.test_args = [maintestfile, "-q"]
        self.test_suite = True
        
    def run_tests(self):
        """Execute the test runner command."""
        # Import here, because outside the required eggs aren't loaded yet
        import pytest
        sys.exit(pytest.main(self.test_args))

# Get the long description from the relevant file
here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()


setup(
  name = 'stpipeline',
  version = '0.2.5',
  description = __doc__.split("\n", 1)[0],
  long_description = long_description,
  keywords = 'rna-seq analysis spatial transcriptomics toolkit',
  author = 'Jose Fernandez Navarro',
  author_email = 'jose.fernandez.navarro@scilifelab.se',
  license = 'Copyright Spatial Transciptomics',
  url = 'http://www.spatialtranscriptomics.com',
  packages = find_packages(exclude=('tests*', 'dependencies', 'utils')),
  include_package_data = False,
  package_data = dict(),
  zip_safe = False,
  install_requires = [
    'HTSeq',
    'setuptools',
    'pysam>=0.7.5',
    'numpy'
  ],
  tests_require = [
    'pytest',
  ],
  cmdclass=dict(
    test = PyTest,
  ),
  scripts = glob.glob('scripts/*.py'),
  classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Copyright Spatial Transciptomics',
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Environment :: Console',
  ],
)