"""SympyCore - An efficient pure Python Computer Algebra System.

See http://sympycore.googlecode.com/ for more information.
"""
__docformat__ = "restructuredtext"
__date__ = '2008'
__author__ = 'Pearu Peterson, Fredrik Johansson'
__license__ = 'New BSD License'

from .version import version as __version__

from .core import classes, defined_functions, DefinedFunction, Expr, Pair, IntegerList
from .core import init_module

import heads

init_module.execute()

from arithmetic import *

from ring import *

from basealgebra import *
from algebras import *
from logic import *
from sets import *
from calculus import *
from polynomials import *
from matrices import *
from physics import *
from functions import *

CollectingField = CommutativeRing

### Initialize sympycore subpackage namespaces
init_module.execute()

def profile_expr(expr):
    """ Printout the profiler information for executing ``expr``.
    """
    import sys
    import hotshot, hotshot.stats
    prof = hotshot.Profile("/tmp/sympycore_stones.prof")
    frame = sys._getframe(1)
    g = frame.f_globals
    l = frame.f_locals
    def foo():
        exec expr in g,l
    prof.runcall(foo)
    prof.close()
    stats = hotshot.stats.load("/tmp/sympycore_stones.prof")
    #stats.strip_dirs()
    stats.sort_stats('time','calls','time')
    stats.print_stats(40)
    return stats

# This function is taken from sympy (and originally sage?)
def var(s):
    """
    var('x y z') creates symbols x, y, z
    """
    import inspect
    frame = inspect.currentframe().f_back
    try:
        if not isinstance(s, list):
            s = s.split(" ")
        res = []
        for t in s:
            # skip empty strings
            if not t:
                continue
            sym = Symbol(t)
            frame.f_globals[t] = sym
            res.append(sym)
        res = tuple(res)
        if len(res) == 0:   # var('')
            res = None
        elif len(res) == 1: # var('x')
            res = res[0]
                            # otherwise var('a b ...')
        return res
    finally:
        # we should explicitly break cyclic dependencies as stated in inspect
        # doc
        del frame

FD = FDFactory(CalculusDifferentialRing)
D = DFactory(CalculusDifferentialRing)

sin = CalculusFunctionRing(heads.CALLABLE, Sin)
arcsin = CalculusFunctionRing(heads.CALLABLE, ArcSin)
cos = CalculusFunctionRing(heads.CALLABLE, Cos)
tan = CalculusFunctionRing(heads.CALLABLE, Tan)
cot = CalculusFunctionRing(heads.CALLABLE, Cot)
ln = CalculusFunctionRing(heads.CALLABLE, Ln)
log = CalculusFunctionRing(heads.CALLABLE, Log)
exp = CalculusFunctionRing(heads.CALLABLE, Exp)

class _Tester:

    def __init__(self):
        self.check_testing()
        
    def test(self, nose_args=''):
        """ Run sympycore tests using nose.
        """
        import os
        import sys
        import nose
        #mpmath = sys.modules['sympycore.arithmetic.mpmath']
        d = os.path.dirname(os.path.abspath(__file__))
        cmd = '%s %s %s %s' % (sys.executable, nose.core.__file__, d, nose_args)
        print >>sys.stderr, 'Running %r' % cmd
        s = os.system(cmd)
        if s:
            print >>sys.stderr, "TESTS FAILED"

    def check_testing(self):
        import os, sys
        try:
            import nose
        except ImportError:
            return
        if sys.platform=='win32':
            m = lambda s: s.lower()
        else:
            m = lambda s: s
        argv2 = [nose.core.__file__, os.path.dirname(os.path.abspath(__file__))]
        if map(m, sys.argv[:2]) == map(m, argv2):
            self.show_config()

    def show_config(self):
        import os, sys

        print >>sys.stderr, 'Python version: %s' % (sys.version.replace('\n',''))

        nose = None
        minimum_nose_version = (0,10,0)
        try:
            import nose
            if nose.__versioninfo__ < minimum_nose_version:
                nose = None
        except:
            pass
        if nose is None:
            msg = 'Need nose >= %d.%d.%d for tests - see ' \
              'http://somethingaboutorange.com/mrl/projects/nose' % \
              minimum_nose_version
            raise ImportError(msg)            
        else:
            print >>sys.stderr, 'nose version: %d.%d.%d' % nose.__versioninfo__
            print >>sys.stderr, 'nose is installed in %s' % (os.path.dirname(nose.__file__))

        if 'sympycore.expr_ext' in sys.modules:
            s = 'compiled'
        else:
            s = 'pure'
        print >>sys.stderr, 'sympycore version: %s (%s)' % (__version__, s)
        print >>sys.stderr, 'sympycore is installed in %s' % (os.path.dirname(__file__))

        mpmath = sys.modules['sympycore.arithmetic.mpmath']
        fn = os.path.join(os.path.dirname(mpmath.__file__), 'REVISION')
        l = []
        backend = mpmath.libmp.backend.BACKEND
        l.append ('backend=%s' % (backend))
        if os.path.isfile(fn):
            f = open(fn, 'r')
            l.append('revision=%s' % (f.read().strip()))
            f.close()
        else:
            rev = ''
        print >>sys.stderr, 'mpmath version: %s (%s)' % (mpmath.__version__, ', '.join(l))
        print >>sys.stderr, 'mpmath is installed in %s' % (os.path.dirname(mpmath.__file__))

        
        if backend=='gmpy':
            gmpy = mpmath.libmp.backend.gmpy
            print >>sys.stderr, 'gmpy version: %s' % (gmpy.version())
            print >>sys.stderr, 'gmpy is installed in %s' % (os.path.dirname(gmpy.__file__))


test = _Tester().test
del _Tester

