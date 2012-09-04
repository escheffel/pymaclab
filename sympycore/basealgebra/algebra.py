""" Provides BasicAlgebra class.
"""

__docformat__ = "restructuredtext"
__all__ = ['Algebra', 'SymbolicEquality']

from ..core import classes, Expr, defined_functions
from ..utils import LT, GT, LE, GE, NE, EQ, SYMBOL, NUMBER
from ..heads import CALLABLE, APPLY

symbolic_comparison_map = dict(
    equality = dict(__eq__=EQ, __ne__=NE),
    inequality = dict(__lt__=LT,__le__=LE, __gt__=GT,__ge__=GE)
    )

class SymbolicEquality:
    """ Contex for logical operations.

    In the ``SymbolicEquality(<Algebra class>)`` context relational
    operations return ``Logic`` instances instead of computing
    lexicographic value of relational operations.

    For example,

    >>> x = Calculus('x')
    >>> with SymbolicEquality(Calculus):
    >>>     print x==x
    ...     
    ...     
    x==x
    >>> print x==x
    True

    Add ``from __future__ import with_statement`` to the header of
    python file when using Python version 2.5.
    
    """

    def __init__(self, *classes):
        self.classes = classes

    def __enter__(self):
        for cls in self.classes:
            cls.enable_symbolic_comparison()

    def __exit__(self, type, value, tb):
        for cls in self.classes:
            cls.disable_symbolic_comparison()
        return tb is None


class Algebra(Expr):
    """ Represents an element of an algebraic structure.

    This class collects implementation specific methods of algebra
    classes.
    
    For implemented algebras, see:

      Verbatim
      CollectingField
    
    New algebras may need to redefine the following methods::
    
      __new__(cls, ...)
      convert(cls, obj, typeerror=True)
      convert_coefficient(cls, obj, typeerror=True)
      convert_exponent(cls, obj, typeerror=True)
      convert_operand(self, obj, typeerror=True)
      as_verbatim(self)
      as_algebra(self, cls, typeerror=True)

    and the following properties::

      args(self)
      func(self)

    """
    #__slots__ = ['_str_value']

    _str_value = None

    coefftypes = (int, long)
    exptypes = (int, long)

    algebra_options = dict(evaluate_addition = True)

    def __str__(self):
        h, d = self.pair
        return h.data_to_str_and_precedence(type(self), d)[0]
        s = self._str_value
        if s is None:
            head, data = self.pair
            s = head.data_to_str(type(self), data, 0.0)
            if not self.is_writable:
                self._str_value = s
        return s

    def __repr__(self):
        s = str(self)
        return '%s(%r)' % (type(self).__name__, s)

    def __nonzero__(self):
        return not not self.data

    @classmethod
    def enable_symbolic_comparison(cls, name = 'equality'):
        """ Enable returning Logic instances from relational methods.

        Argument ``name`` can be 'equality' or 'inequality'.
        """
        logic_map = symbolic_comparison_map.get(name)
        if logic_map is None:
            raise ValueError('Unknown symbolic comparison name: %r. Valid names are %s'\
                             % (name, ', '.join([repr(k) for k in symbolic_comparison_map])))
        logic_map_name = '_%s_map' % (name)
        d = getattr(cls, logic_map_name, None)
        if d is None:
            d = {}
            setattr(cls, logic_map_name, d)
            for mthname, head in logic_map.items():
                d[mthname] = getattr(cls, mthname, None)
                setattr(cls, mthname, lambda self, other, head=head: classes.Logic(head, (self, other)))

    @classmethod
    def disable_symbolic_comparison(cls, name = 'equality'):
        """ Disable returning Logic instances from relational methods.

        Argument ``name`` can be 'equality' or 'inequality'.
        """
        logic_map = symbolic_comparison_map.get(name)
        if logic_map is None:
            raise ValueError('Unknown symbolic comparison name: %r. Valid names are %s'\
                             % (name, ', '.join([repr(k) for k in symbolic_comparison_map])))
        logic_map_name = '_%s_map' % (name)
        d = getattr(cls, logic_map_name, None)
        if d is None:
            return
        for mthname, mth in d.items():
            if mth is None:
                delattr(cls, mthname)
            else:
                setattr(cls, mthname, mth)
        delattr(cls, logic_map_name)

    def as_tree(self, tab='', level=0):
        return self.as_verbatim().as_tree(tab,level)

    @classmethod
    def get_operand_algebra(cls, head, index=0):
        """ Return algebra class for index-th operand of operation with head.
        """
        return cls

    @classmethod
    def handle_get_operand_algebra_failure(cls, head, index):
        raise TypeError('Cannot convert %s-th %s operand to %s algebra' % (index, repr(head), cls.__name__))

    @classmethod
    def convert(cls, obj, typeerror=True):
        """ Convert obj to algebra element.

        With typeerror=False NotImplemented is returned when
        conversion is not possible, otherwise a TypeError is raised.
        
        Set typeerror=False when calling from operation methods like __add__,
        __mul__, etc.
        """
        tobj = type(obj)
        # check if obj is already algebra element:
        if tobj is cls:
            return obj
        
        # check if obj belongs to coefficient algebra
        r = cls.convert_number(obj, typeerror=False)
        if r is not NotImplemented:
            return cls.Number(r)

        # parse algebra expression from string:
        if isinstance(obj, (str, unicode, Verbatim)):
            return Verbatim.convert(obj).as_algebra(cls)

        # as a last resort, convert from another algebra:
        if isinstance(obj, Algebra):
            result = obj.as_algebra(cls, typeerror=typeerror)
            if result is not NotImplemented:
                return result
            
        return cls.handle_convert_failure('algebra', obj, typeerror)

    @classmethod
    def convert_exponent(cls, obj, typeerror=True):
        """ Convert obj to exponent algebra.
        """
        if isinstance(obj, cls.exptypes) or isinstance(obj, cls):
            return obj
        return cls.handle_convert_failure('exponent algebra', obj, typeerror, 'expected int|long|'+cls.__name__)

    @classmethod
    def convert_coefficient(cls, obj, typeerror=True):
        """ Convert obj to coefficient algebra.
        """
        # Use convert_number instead.
        if isinstance(obj, cls.coefftypes):
            return obj
        return cls.handle_convert_failure('coefficient algebra', obj, typeerror, 'expected int|long')

    @classmethod
    def handle_convert_failure(cls, name, obj, typeerror, hint=''):
        if typeerror:
            if hint:
                hint = ' (%s)' % (hint)
            errmsg = 'failed to convert %r to %r %s%s'\
                     % (obj.__class__.__name__, cls.__name__, name, hint)
            raise TypeError(errmsg)
        else:
            return NotImplemented

    @classmethod
    def convert_symbol(cls, obj, typeerror=True):
        """ Convert obj to algebra symbol or predefined function or
        constant.
        """
        r = cls.get_predefined_symbols(obj)
        if r is not None:
            return r
        return cls(SYMBOL, obj)

    @classmethod
    def convert_number(cls, obj, typeerror=True):
        """ Convert obj to algebra number.

        On success, the result is assumed to be used in ``Algebra.Number(<result>)``.

        On failure, the TypeError is raised if typeerror is True (default).
        Otherwise, NotImplemented is returned.
        """
        #XXX: return cls(NUMBER, obj)
        return cls.convert_coefficient(obj, typeerror)

    def convert_operand(self, obj, typeerror=True):
        """ Convert obj to operand.

        Used when operand algebra depends on expression head.
        """
        return self.convert(obj, typeerror=typeerror)

    @classmethod
    def convert_function(cls, name, args):
        """ Convert verbatim applied function call to algebra.
        """
        func = getattr(defined_functions, name.title(), None)
        arg_algebras = func.get_argument_algebras()
        return func(*[cls(a) for a,cls in zip(args, arg_algebras)])
    
    def as_verbatim(self):
        def as_verbatim(cls, head, data, target):
            return cls(head, data)
        return self.head.walk(as_verbatim, Verbatim, self.data, self)

    def as_algebra(self, cls, typeerror=True):
        """ Convert algebra to another algebra.

        This method uses default conversation via verbatim algebra that
        might not be the most efficient. For efficiency, algebras should
        redefine this method to implement direct conversation.
        """
        # todo: cache verbatim algebras
        if cls is classes.Verbatim:
            return self.as_verbatim()
        return self.as_verbatim().as_algebra(cls)

    @classmethod
    def get_predefined_symbols(cls, name):
        return getattr(defined_functions, name, None)

    @property
    def func(self):
        """ Returns a callable such that ``self.func(*self.args) == self``.
        """
        raise NotImplementedError('%s must define property func'      #pragma NO COVER
                                  % (self.__class__.__name__))                   #pragma NO COVER

    @property
    def args(self):
        """ Returns a sequence such that ``self.func(*self.args) == self``.
        """
        raise NotImplementedError('%s must define property "args"'      #pragma NO COVER
                                  % (self.__class__.__name__))        #pragma NO COVER

    @classmethod
    def Symbol(cls, obj):
        """ Construct algebra symbol directly from obj.
        """
        return cls(SYMBOL, obj)

    @classmethod
    def Number(cls, num):
        """ Construct algebra number directly from obj.
        """
        return cls(NUMBER, num)

    @classmethod
    def Apply(cls, func, *args, **kwargs):
        convert= cls.convert
        fcls = cls.get_function_algebra()
        if isinstance(func, Expr):
            pass
        elif callable(func):
            func = fcls(CALLABLE, func)
        else:
            func = fcls(SYMBOL, func)
        new_args = map(convert, args)
        for k, v in kwargs.items():
            new_args.append(cls(KWARG, (convert(k), convert(v))))
        return cls(APPLY, (func, tuple(new_args)))

    _cache_symbols = None

    @property
    def symbols(self):
        """ Return a set of atomic subexpressions in a symbolic object.
        """
        symbols = self._cache_symbols
        if symbols is None:
            def scan_for_symbols(cls, head, data, target):
                if head is SYMBOL:
                    target.add(cls(head, data)) # introduce expr argument to scan
            head, data = self.pair
            symbols = set([])
            head.scan(scan_for_symbols, type(self), data, symbols)
            self._cache_symbols = symbols
            #todo: make self read-only
        return symbols

    _cache_symbols_data = None
    @property
    def symbols_data(self):
        symbols = self._cache_symbols_data
        if symbols is None:
            symbols = set([s.data for s in self.symbols])
            self._cache_symbols_data = symbols
        return symbols

    def has_symbol(self, symbol):
        symbol = type(self)(symbol)
        return symbol in self.symbols

    def has(self, obj):
        """ Return True if self contains atomic expression obj.
        """
        convert = self.convert
        return convert(obj) in self.symbols

    def match(self, pattern, *wildcards, **options):
        """
        Pattern matching.

        Return None when expression (self) does not match with
        pattern. Otherwise return a dictionary such that

        ::

          pattern.subs_dict(self.match(pattern, *wildcards)) == self

        Don't redefine this method, redefine ``matches(..)`` method instead.

        :See:
          http://wiki.sympy.org/wiki/Generic_interface#Pattern_matching
        """
        exclude = options.get('exclude',[])
        pattern = self.convert(pattern)
        wild_expressions = []
        wild_predicates = []
        for w in wildcards:
            if type(w) in [list, tuple]:
                assert len(w)==2,`w`
                s, func = w
            else:
                if exclude:
                    s, func = w, lambda x: x not in exclude
                else:
                    s, func = w, True
            s = self.convert(s)
            wild_expressions.append(s)
            wild_predicates.append(func)
        if wild_expressions:
            wild_info = (wild_expressions, wild_predicates)
        else:
            wild_info = None
        return pattern.matches(self, {}, wild_info)

    def matches(pattern, expr, repl_dict={}, wild_info=None):
        # check if pattern has a match and return it provided that
        # the match matches with expr:
        v = repl_dict.get(pattern)
        if v is not None:
            if v==expr:
                return repl_dict
            return
        # check if pattern matches with expr:
        if pattern==expr:
            if wild_info and expr in wild_info[0]:
                return
            return repl_dict
        if wild_info:
            wild_expressions, wild_predicates = wild_info
            # check if pattern is a wild expression
            if pattern in wild_expressions:
                if expr in wild_expressions:
                    # wilds do not match other wilds
                    return
                if expr.symbols.intersection(pattern.symbols):
                    # wild does not match expressions containing wild expressions
                    return
                # wild pattern matches with expr only if predicate(expr) returns True
                predicate = wild_predicates[wild_expressions.index(pattern)]
                if (isinstance(predicate, bool) and predicate) or predicate(expr):
                    repl_dict = repl_dict.copy()
                    repl_dict[pattern] = expr
                    return repl_dict
            pattern_symbols = pattern.symbols
            for w in wild_expressions:
                if w in pattern_symbols:
                    # _matches implements implementation specific matches,
                    # pattern should contain wild symbols, otherwise
                    # _matches would always return None
                    return pattern._matches(expr, repl_dict, wild_info)
        return

    @classmethod
    def _matches_seq(cls, pattern_seq, expr_seq, repl_dict, wild_info):
        n = len(pattern_seq)
        m = len(expr_seq)
        if n!=m:
            return
        elif n==0:
            return repl_dict
        elif n==1:
            p = pattern_seq[0]
            if isinstance(p, cls):
                return p.matches(expr_seq[0], repl_dict, wild_info)
            if p==expr_seq[0]:
                return repl_dict
            return
        def matches(pattern, expr):
            if isinstance(pattern, cls):
                return pattern.matches(expr, repl_dict, wild_info)
            if pattern==expr:
                return repl_dict
            return
        def subs(expr, l):
            if isinstance(expr, cls):
                return expr.subs(l)
            return expr
        for i in xrange(n):
            p = pattern_seq[i]
            e = expr_seq[i]
            r = matches(p, e)
            if r is not None:
                new_pattern_seq = [subs(pattern_seq[j],r.items()) for j in xrange(n) if j!=i]
                new_expr_seq = [expr_seq[j] for j in xrange(n) if j!=i]
                r1 = cls._matches_seq(new_pattern_seq, new_expr_seq, r, wild_info)
                if r1 is not None:
                    return r1
        return

    def _matches(pattern, expr, repl_dict, wild_info):
        if pattern.func == expr.func:
            return pattern._matches_seq(pattern.args, expr.args, repl_dict, wild_info)
        return

    def subs(self, subexpr, newexpr=None):
        """ Substitute a sub-expression with new expression.
        
        There are three usage forms::
        
          obj.subs(subexpr, newexpr)
          obj.subs([(subexpr1, newexpr1), (subexpr2, newexpr2), ..])
          obj.subs(<dict of subexpr:newexpr>)
        """
        cls = type(self)
        if newexpr is not None:
            def subs(cls, head, data, target, subexpr=subexpr, newexpr=newexpr):
                if target==subexpr:
                    if type(newexpr) is not cls:
                        newexpr = cls(newexpr)
                        return newexpr.head.new(cls, newexpr.data)
                    return newexpr
                return target
        elif type(subexpr) is dict:
            def subs(cls, head, data, target, expr_dict=subexpr):
                r = expr_dict.get(target, target)
                if type(r) is not cls:
                    r = cls(r)
                    return r.head.new(cls, r.data)
                return r
        elif type(subexpr) is list:
            def subs(cls, head, data, target, expr_list=subexpr):
                for s, n in expr_list:
                    if s==target:
                        if type(n) is not cls:
                            n = cls(n)
                            return n.head.new(cls, n.data)
                        return n
                return target
        else:
            raise NotImplementedError(`subexpr, newexpr`)
        head, data = self.pair
        return head.walk(subs, cls, data, self)

        convert = lambda x:x
        if newexpr is None:
            if type(subexpr) is dict:
                subexpr = subexpr.iteritems()
            r = self
            for s,n in subexpr:
                r = r._subs(convert(s), convert(n))
            return r
        return self._subs(convert(subexpr), convert(newexpr))

    def _subs(self, subexpr, newexpr):
        head, data = self.pair
        t = type(data)
        cls = type(self)

        if type(subexpr) is not cls:
            if head is SYMBOL:
                if subexpr==data:
                    return self.convert_operand(newexpr)
                return self
        else:
            h, d = subexpr.pair
            if h is head and data==d:
                return self.convert_operand(newexpr)

        if head is SYMBOL:
            return self
        
        if isinstance(data, (tuple, list)):
            r = []
            for a in data:
                if hasattr(a, '_subs'):
                    r.append(a._subs(subexpr, newexpr))
                else:
                    r.append(a)
            return cls(head, t(r))

        if isinstance(data, dict):
            d = dict.__new__(t)
            for t,c in data.iteritems():
                t1 = t._subs(subexpr, newexpr) if hasattr(t,'_subs') else t
                c1 = c._subs(subexpr, newexpr) if hasattr(c,'_subs') else c
                d[t1] = c1
            return cls(head, d)

        return cls(head, data._subs(subexpr, newexpr) if hasattr(data, '_subs') else data)

from .verbatim import Verbatim
