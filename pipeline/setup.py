#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

"""ST Pipeline: processing sequencing data on Amazon EMR or a single computer trough a terminal.

ST Pipeline is a toolkit of RNA-Seq data analysis for the Spatial Transcriptomics groups
data.

"""
import glob
import os
import sys
from distutils.core import setup
from distutils.command.build import build as du_build
from distutils.command.clean import clean as du_clean

def get_arg(name):
    arg_start = "%s=" % name
    for i in xrange(len(sys.argv)):
        arg = sys.argv[i]
        if arg.startswith(arg_start):
            value = arg.replace(arg_start, "", 1)
            del sys.argv[i]
            if not value:
                raise RuntimeError("blank value specified for %s" % name)
            return value
    return None

def check_python_version():
    override = ("true" == get_arg("override_version_check"))
    if not override and sys.version_info < (2,6):
        print >>sys.stderr, "Please use a version of Python >= 2.6 (currently using vers. %s)." % ",".join( map(str, sys.version_info))
        print >>sys.stderr, "Specify setup.py override_version_check=true to override this check."
        sys.exit(1)

def get_version():
    vers = get_arg("version")
    if vers is None:
        # else, if no version specified on command line
        version_filename = os.path.join(os.path.dirname(__file__), 'VERSION')
        if os.path.exists(version_filename):
            with open(version_filename) as f:
                vers = f.read().rstrip("\n")
    else:
        from datetime import datetime
        vers = datetime.now().strftime("devel-%Y%m%d")  # _%H%M%S")
    return vers

class build(du_build):
    
    def initialize_options(self):
        du_build.initialize_options(self)
        self.scripts_path = None

    def finalize_options(self):
        du_build.finalize_options(self)
        # HACK!  Use a global variable until we find a better way to
        # pass a parameter into the build command.
        global VERSION
        self.version = VERSION
        self.override_version_check = get_arg("override_version_check") or 'false'
        
    def run(self):
        # Create (or overwrite) seal/version.py
        with open(os.path.join('main', 'version.py'), 'w') as f:
            f.write('version = "%s"' % self.version)

        # run the usual build
        du_build.run(self)

class clean(du_clean):
    
    def run(self):
        du_clean.run(self)


#############################################################################
# main
#############################################################################

# chdir to root directory (where this file is located)
os.chdir(os.path.abspath(os.path.dirname(__file__)))

check_python_version()

NAME = 'stpipeline'
DESCRIPTION = __doc__.split("\n", 1)[0]
LONG_DESCRIPTION = __doc__
URL = "http://www.spatialtranscriptomics.com"
DOWNLOAD_URL = "http://www.dowload.spatialtranscriptomics.com"
LICENSE = 'GPL'
CLASSIFIERS = [
  "Programming Language :: Python",
  "License :: OSI Approved :: GNU General Public License (GPL)",
  "Operating System :: POSIX :: Linux",
  "Topic :: Scientific/Engineering :: Bio-Informatics",
  "Intended Audience :: Science/Research",
  ]
PLATFORMS = ["Linux"]
VERSION = get_version()
AUTHOR_INFO = [
  ("Jose Fernandez Navarro", "jose.fernandez.navarro@scilifelab.se"),
  ]
MAINTAINER_INFO = [
  ("Jose Fernandez Navarro", "jose.fernandez.navarro@scilifelab.se"),
  ]
AUTHOR = ", ".join(t[0] for t in AUTHOR_INFO)
AUTHOR_EMAIL = ", ".join("<%s>" % t[1] for t in AUTHOR_INFO)
MAINTAINER = ", ".join(t[0] for t in MAINTAINER_INFO)
MAINTAINER_EMAIL = ", ".join("<%s>" % t[1] for t in MAINTAINER_INFO)

setup(name=NAME,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      classifiers=CLASSIFIERS,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      platforms=PLATFORMS,
      version=VERSION,
      packages=['main',
                'main.lib',
                'main.common',
                'main.core',
                ],
      cmdclass={
          "build": build,
          "clean": clean,
          },
      scripts=glob.glob("scripts/*"),
      )

