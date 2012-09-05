#!/usr/bin/env python
# Parallel Python setup script
# For the latest version of the Parallel Python
# software visit: http://www.parallelpython.com
"""
Standard build tool for python libraries.
"""

import os.path
from distutils.core import setup

from pp import version as VERSION


LONG_DESCRIPTION = """
Parallel Python module (PP) provides an easy and efficient way to create \
parallel-enabled applications for SMP computers and clusters. PP module \
features cross-platform portability and dynamic load balancing. Thus \
application written with PP will parallelize efficiently even on \
heterogeneous and multi-platform clusters (including clusters running other \
application with variable CPU loads). Visit http://www.parallelpython.com \
for further information.
"""


setup(
        name="pp",
        url="http://www.parallelpython.com",
        version=VERSION,
        download_url="http://www.parallelpython.com/downloads/pp/pp-%s.zip" % (
            VERSION),
        author="Vitalii Vanovschi",
        author_email="support@parallelpython.com",
        py_modules=["pp", "ppauto", "ppcommon", "pptransport", "ppworker"],
        scripts=["ppserver.py"],
        description="Parallel and distributed programming for Python",
        platforms=["Windows", "Linux", "Unix"],
        long_description=LONG_DESCRIPTION,
        license="BSD-like",
        classifiers=[
        "Topic :: Software Development",
        "Topic :: System :: Distributed Computing",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Intended Audience :: Developers",
        "Development Status :: 5 - Production/Stable",
        ],
)
