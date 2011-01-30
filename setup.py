#!/usr/bin/env python

from datetime import datetime

from numpy.distutils.core import setup
from numpy import get_include

DESCRIPTION=""
LONG_DESCRIPTION="""
"""

DISTNAME = 'pymaclab'
LICENSE = 'BSD'
AUTHOR = "Skipper Seabold, E.M. Scheffel"
MAINTAINER = "Skipper Seabold"
MAINTAINER_EMAIL = "jsseabold@gmail.com"
URL = ''
DOWNLOAD_URL=""
CLASSIFIERS=[]

MAJOR = 0
MINOR = 0
MICRO = 1
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
FULLVERSION = VERSION
if not ISRELEASED:
    FULLVERSION += '.dev'

def write_version_py(filename='./version.py'):
    cnt = """\
from datetime import datetime

version = '%s'
"""
    a = open(filename, 'w')
    try:
        a.write(cnt % FULLVERSION)
    finally:
        a.close()

def configuration(parent_package='', top_path=None, package_name=DISTNAME):
#    write_version_py()

    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path,
                           version=FULLVERSION)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('pymaclab')
    return config

if __name__ == '__main__':
    setup(name=DISTNAME,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          license=LICENSE,
          url=URL,
          download_url=DOWNLOAD_URL,
          long_description=LONG_DESCRIPTION,
          classifiers=CLASSIFIERS,
          platforms='any',
          configuration=configuration)
