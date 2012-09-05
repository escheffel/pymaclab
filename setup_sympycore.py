#!/usr/bin/env python

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
License :: OSI Approved
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

import os
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

from distutils.core import Extension, Command
from distutils.command.build_py import build_py as _build_py

expr_ext = Extension('sympycore.expr_ext',
                     sources = [os.path.join('src','expr_ext.c')],
                     )

combinatorics_ext = Extension('sympycore.arithmetic.combinatorics',
                     sources = [os.path.join('sympycore','arithmetic','combinatorics.c')],
                     )

extensions = [expr_ext,
              #combinatorics_ext
              ]

packages = ['sympycore',
            'sympycore.algebras',
            'sympycore.arithmetic',
            'sympycore.arithmetic.mpmath',
            'sympycore.arithmetic.mpmath.calculus',
            'sympycore.arithmetic.mpmath.functions',
            'sympycore.arithmetic.mpmath.libmp',
            'sympycore.arithmetic.mpmath.matrices',
            'sympycore.basealgebra',
            'sympycore.calculus',
            'sympycore.calculus.functions',
            'sympycore.functions',
            'sympycore.heads',
            'sympycore.logic',
            'sympycore.matrices',
            'sympycore.polynomials',
            'sympycore.physics',
            'sympycore.physics.sysbio',
            'sympycore.ring',
            'sympycore.sets',
            ]

packages += [p+'.tests' for p in packages \
             if os.path.exists(os.path.join(p.replace('.', os.sep), 'tests'))]

class tester(Command):
    description = "run sympycore tests"
    user_options = [('nose-args=', 'n', 'arguments to nose command'),
                    ('with-coverage', 'c', 'use nose --with-coverage flag'),
                    ('cover-package=', None, 'use nose --cover-package flag'),
                    ('detailed-errors', 'd', 'use nose --detailed-errors flag'),
                    ('nocapture', 's', 'use nose --nocapture flag'),
                    ('nose-verbose', 'v', 'use nose --verbose flag'),
                    ('match=', 'm', 'use nose --match flag'),
                    ('profile', 'p', 'use nose --profile flag'),
                    ('with-doctest', None, 'use nose --with-doctest flag'),
                    ('stop', 'x', 'use nose --stop flag')
                    ]
    def initialize_options(self):
        self.nose_args = None
        self.with_coverage = None
        self.cover_package = None
        self.detailed_errors = None
        self.nocapture = None
        self.nose_verbose = None
        self.match = None
        self.profile = None
        self.with_doctest = None
        self.stop = None
        return
    def finalize_options (self):
        if self.nose_args is None:
            self.nose_args = ''
        if self.with_coverage:
            self.nose_args += ' --with-coverage'
        if self.cover_package:
            if not self.with_coverage:
                self.nose_args += ' --with-coverage'
            self.nose_args += ' --cover-package=%s' % self.cover_package
        elif self.with_coverage:
            self.nose_args += ' --cover-package=sympycore'
        if self.detailed_errors:
            self.nose_args += ' --detailed-errors'
        if self.nocapture:
            self.nose_args += ' --nocapture'
        if self.nose_verbose:
            self.nose_args += ' --verbose'
        if self.match:
            self.nose_args += ' --match=%r' % (self.match)
        if self.profile:
            self.nose_args += ' --with-profile'
        if self.with_doctest:
            self.nose_args += ' --with-doctest'
        if self.stop:
            self.nose_args += ' --stop'
        return
    def run(self):
        import sympycore
        sympycore.test(nose_args=self.nose_args)

class build_py(_build_py):

    def find_data_files (self, package, src_dir):

        files = _build_py.find_data_files(self, package, src_dir)

        if package=='sympycore':
            revision = self._get_svn_revision(src_dir)
            if revision is not None:
                target = os.path.join(src_dir, '__svn_version__.py')
                print 'Creating ', target
                f = open(target,'w')
                f.write('version = %r\n' % (str(revision)))
                f.close()
                import atexit
                def rm_file(f=target):
                    try: os.remove(f); print 'Removed ',f
                    except OSError: pass
                    try: os.remove(f+'c'); print 'Removed ',f+'c'
                    except OSError: pass
                atexit.register(rm_file)
                files.append(target)

        if package=='sympycore.arithmetic.mpmath':
            f = os.path.join(src_dir, 'REVISION')
            if os.path.isfile(f):
                files.append(f)

        return files

    def _get_svn_revision(self, path):
        """Return path's SVN revision number.
        """
        import os, sys, re
        revision = None
        m = None
        try:
            sin, sout = os.popen4('svnversion')
            m = re.match(r'(?P<revision>\d+)', sout.read())
        except:
            pass
        if m:
            revision = int(m.group('revision'))
            return revision
        if sys.platform=='win32' and os.environ.get('SVN_ASP_DOT_NET_HACK',None):
            entries = os.path.join(path,'_svn','entries')
        else:
            entries = os.path.join(path,'.svn','entries')
        if os.path.isfile(entries):
            f = open(entries)
            fstr = f.read()
            f.close()
            if fstr[:5] == '<?xml':  # pre 1.4
                m = re.search(r'revision="(?P<revision>\d+)"',fstr)
                if m:
                    revision = int(m.group('revision'))
            else:  # non-xml entries file --- check to be sure that
                m = re.search(r'dir[\n\r]+(?P<revision>\d+)', fstr)
                if m:
                    revision = int(m.group('revision'))
        return revision

        
from distutils.core import setup
setup(name='sympycore',
      version='0.2-svn',
      author = 'Pearu Peterson, Fredrik Johansson',
      author_email = 'sympycore@googlegroups.com',
      license = 'http://sympycore.googlecode.com/svn/trunk/LICENSE',
      url = 'http://sympycore.googlecode.com',
      download_url = 'http://code.google.com/p/sympycore/downloads/',
      classifiers=filter(None, CLASSIFIERS.split('\n')),
      description = 'SympyCore: an efficient pure Python Computer Algebra System',
      long_description = '''\
SympyCore project provides a pure Python package sympycore for
representing symbolic expressions using efficient data structures as
well as methods to manipulate them. Sympycore uses a clear algebra
oriented design that can be easily extended.
''',
      platforms = ["All"],
      packages = packages,
      ext_modules = extensions,
      package_dir = {'sympycore': 'sympycore'},
      cmdclass=dict(test=tester, build_py=build_py)
      )

