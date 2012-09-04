#
# Created in 2007 by Fredrik Johansson
# Expression parser added by Pearu Peterson
#
""" Provides PrimitiveAlgebra class and expression parser.
"""
from __future__ import absolute_import
__docformat__ = "restructuredtext"
__all__ = ['Verbatim']

import types
import compiler
from compiler import ast

from .algebra import Algebra
from ..heads import (OR, AND, NOT, LT, LE, GT, GE, EQ, NE, BAND, BOR, BXOR,
                     INVERT, POS, NEG, ADD, SUB, MOD, MUL, DIV, FLOORDIV, POW,
                     LSHIFT, RSHIFT, IS, ISNOT, LIST, SLICE,
                     NUMBER, SYMBOL, APPLY, TUPLE, LAMBDA, TERMS, BASE_EXP_DICT,
                     IN, NOTIN, SUBSCRIPT, SPECIAL, DICT, ATTR, KWARG)
from ..heads import CALLABLE
from ..core import classes, Expr, objects

# Restrictions:
#
#    Star and double star function arguments are not implemented,
#    i.e. parsing 'f(*a)' and 'f(**b)' will fail.
#
# TODO: IfExp support: parse `a if b else c` to Verbatim(IF, (b, a, c))

EllipsisType = type(Ellipsis)
special_types = (EllipsisType, type(None), type(NotImplemented))
special_objects = set([Ellipsis, None, NotImplemented])
number_types = (int, long, float, complex)

containing_lst = set([IN, NOTIN])

atomic_lst = set([SYMBOL, NUMBER])
unary_lst = set([POS, NEG, NOT, INVERT])
binary_lst = set([AND, OR, BAND, BOR, BXOR, ADD, SUB, MUL, DIV, FLOORDIV,
                  MOD, POW, LSHIFT, RSHIFT])

convert_head_Op_map = {
    NEG : 'Neg',
    POS : 'Pos',
    ADD : 'Add',
    SUB : 'Sub',
    MUL : 'Mul',
    DIV : 'Div',
    FLOORDIV : 'FloorDiv',
    POW : 'Pow',
    MOD : 'Mod',
    AND : 'And',
    OR : 'Or',
    NOT : 'Not',
    LT : 'Lt',
    GT : 'Gt',
    LE : 'Le',
    GE : 'Ge',
    EQ : 'Eq',
    NE : 'Ne',
    }

class Verbatim(Algebra):
    """ Represents an unevaluated expression.
    """

    commutative_add = None
    commutative_mul = None
    disable_sorting = None

    _str = None

    @classmethod
    def get_value_algebra(cls):
        return cls

    def get_argument_algebra(self, index):
        return self.get_value_algebra()

    @classmethod
    def get_function_algebra(cls):
        return cls

    @classmethod
    def get_differential_algebra(cls):
        return cls
    
    @classmethod
    def convert(cls, obj):
        if isinstance(obj, (str, unicode)):
            obj = string2Verbatim(obj)
        if isinstance(obj, Verbatim):
            return obj
        if hasattr(obj, 'as_verbatim'):
            # handle low-level numbers and constants, as well as Verbatim subclasses
            try:
                return obj.as_verbatim()
            except TypeError, msg:
                if str(msg)=='unbound method as_verbatim() must be called with Verbatim instance as first argument (got nothing instead)':
                    pass
                else:
                    raise
        if isinstance(obj, slice):
            slice_args = obj.start, obj.stop, obj.step
            return cls(SLICE, tuple(map(cls.convert, slice_args)))
        elif isinstance(obj, tuple):
            return cls(TUPLE, tuple(map(cls.convert, obj)))
        elif isinstance(obj, list):
            return cls(LIST, tuple(map(cls.convert, obj)))
        elif isinstance(obj, special_types):
            return cls(SPECIAL, obj)
        elif isinstance(obj, number_types):
            return cls(NUMBER, obj)
        elif isinstance(obj, dict):
            return cls(DICT,  tuple([(cls.convert(k), cls.convert(v)) for k,v in obj.iteritems()]))
        return Verbatim(SYMBOL, obj)

    def as_verbatim(self):
        return self

    def as_algebra(self, cls, source=None):
        def as_algebra(cls, head, data, target):
            if head is NUMBER:
                return cls(data)
            if head is SYMBOL:
                return cls.convert_symbol(data)
            return head.reevaluate(cls, data)
        head, rest = self.pair
        return head.walk(as_algebra, cls, rest, self)

    def __repr__(self):
        return '%s(%r, %r)' % (type(self).__name__, self.head, self.data)

    def as_tree(self, tab='', level=0):
        if level:
            r = []
        else:
            r = [self.__class__.__name__+':']
        head, rest = self.pair
        if head in atomic_lst:
            r.append(tab + '%r[%s]' % (head, rest))
        else:
            r.append(tab + '%r[' % (head))
            if isinstance(rest, Verbatim):
                rest = rest,
            for t in rest:
                if type(t) is tuple:
                    r.append(tab + '  (')
                    for t1 in t:
                        r.append(t1.as_tree(tab=tab + '    ', level=level+1))
                    r.append(tab + '  )')
                else:
                    r.append(t.as_tree(tab=tab + '  ', level=level+1))

            r.append(tab+']')
        return '\n'.join(r)

    def __eq__(self, other):
        if type(other) is Verbatim:
            return self.pair == other.pair
        return False

    for _h in unary_lst:
        if _h.op_mth:
            exec '''\
def %s(self):
    return Verbatim(%r, self)
''' % (_h.op_mth, _h)

    for _h in binary_lst:
        if _h.op_mth:
            exec '''\
def %s(self, other):
    other = self.convert(other)
    return Verbatim(%r, (self, other))
''' % (_h.op_mth, _h)
        if _h.op_rmth:
            exec '''\
def %s(self, other):
    other = self.convert(other)
    return Verbatim(%r, (other, self))
''' % (_h.op_rmth, _h)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __pow__(self, other):
        other = self.convert(other)
        if other.head is NUMBER:
            other = other.data
        return Verbatim(POW, (self, other))

    def __rpow__(self, other):
        other = self.convert(other)
        if self.head is NUMBER:
            self = self.data
        return Verbatim(POW, (other, self))


    def __divmod__(self, other):
        other = self.convert(other)
        return Verbatim(APPLY, (Verbatim(CALLABLE, divmod), (self, other)))

    def __call__(self, *args, **kwargs):
        convert = self.convert
        args = map(convert, args)
        for k, v in kwargs.items():
            args.append(Verbatim(KWARG, (convert(k), convert(v))))
        return Verbatim(APPLY, (self, tuple(args)))

    def __getitem__(self, key):
        if type(key) is tuple:
            key = tuple(map(self.convert, key))
            return Verbatim(SUBSCRIPT, (self, key))
        key = self.convert(key)
        return Verbatim(SUBSCRIPT, (self, (key,)))

    def __getattr__(self, attr):
        # warning: with this feature hasattr() may return True when not desired.
        if not attr.startswith('_'):
            return Verbatim(ATTR, (self, self.convert(attr)))
        raise AttributeError

classes.Verbatim = Verbatim

########### string to Verbatim parser ############

node_names = []
skip_names = ['Module','Stmt','Discard']
for n, cls in ast.__dict__.items():
    if n in skip_names:
        continue
    if isinstance(cls, (type,types.ClassType)) and issubclass(cls, ast.Node):
        node_names.append(n)

node_map = dict(Add=ADD, Mul=MUL, Sub=SUB,
                Div=DIV, FloorDiv=FLOORDIV, TrueDiv=DIV,
                UnaryAdd=POS, UnarySub=NEG, Mod=MOD, Not=NOT,
                Or=OR, And=AND, Power=POW,
                Bitand=BAND,Bitor=BOR,Bitxor=BXOR,
                Tuple=TUPLE, Subscript=SUBSCRIPT,
                Invert=INVERT, LeftShift=LSHIFT, RightShift=RSHIFT,
                List=LIST, Sliceobj=SLICE
                )

callfunc_map = dict(#divmod=DIVMOD,
                    slice=SLICE)

compare_map = {'<':LT, '>':GT, '<=':LE, '>=':GE,
               '==':EQ, '!=':NE, 'in':IN, 'not in': NOTIN,
               'is':IS, 'is not':ISNOT}

class VerbatimWalker:
    """ Helper class for expression parser.
    """

    def __init__(self):
        self.stack = []

    # for composite instance:
    def start(self, head):
        self.stack.append([head, []])
    def append(self, obj):
        stack = self.stack
        if not stack:
            stack.append(obj)
        else:
            stack[-1][1].append(obj)
    def end(self, tuple_head = None):
        head, lst = self.stack.pop()
        if self.stack:
            last = self.stack[-1]
            if last[0]==head and head in [ADD, MUL]:
                # apply associativity:
                last[1].extend(lst)
                return
        if head is tuple_head:
            self.append(tuple(lst))
        elif head in unary_lst:
            assert len(lst)==1,`lst`
            self.append(Verbatim(head, lst[0]))
        elif head is POW:
            b, e = lst
            if e.head is NUMBER:
                e = e.data
            self.append(Verbatim(head, (b, e)))
        else:
            self.append(Verbatim(head, tuple(lst)))
    
    # for atomic instance:
    def add(self, *args):
        self.append(Verbatim(*args))

    for _n in node_names:
        if _n in node_map:
            continue
        exec '''\
def visit%s(self, node, *args):
    print "WARNING: using default unsupported visit%s, node=", node
    self.start(%r)
    for child in node.asList():
        if isinstance(child, ast.Node):
            self.visit(child, *args)
        else:
            self.append(child)
    self.end()
''' % (_n, _n, _n)

    for _n,_v in node_map.items():
        exec '''\
def visit%s(self, node):
    self.start(%r)
    for child in node.getChildNodes():
        self.visit(child)
    self.end()
''' % (_n, _v)

    # visitNode methods:

    def visitSubscript(self, node):
        self.start(SUBSCRIPT)
        self.visit(node.expr)
        self.start(TUPLE)
        map(self.visit, node.subs)
        self.end(tuple_head = TUPLE)
        self.end()

    def visitName(self, node):
        self.add(SYMBOL, node.name)

    def visitConst(self, node):
        if node.value in special_objects:
            self.add(SPECIAL, node.value)
        else:
            self.add(NUMBER, node.value)

    def visitCompare(self, node):
        lhs = node.expr
        op, rhs = node.ops[0]
        if len(node.ops)==1:
            self.start(compare_map[op])
            self.visit(lhs)
            self.visit(rhs)
            self.end()
            return
        n = ast.And([ast.Compare(lhs, node.ops[:1]),
                     ast.Compare(rhs, node.ops[1:])])
        self.visit(n)

    def visitLambda(self, node):
        assert not (node.kwargs or node.varargs),`node.kwargs, node.varargs` # parsing `lambda *args, **kwargs: ..` not supported
        self.start(LAMBDA)
        self.start(TUPLE)
        for n,d in zip(node.argnames, (len(node.argnames) - len(node.defaults))*[None] + list(node.defaults)):
            if d is None:
                self.visit(ast.Name(n))
            else:
                self.visit(ast.Keyword(n, d))
        self.end(tuple_head=TUPLE)
        self.visit(node.code)
        self.end()

    def visitCallFunc(self, node):
        if node.star_args is not None:
            raise NotImplementedError('parsing function star arguments')
        if node.dstar_args is not None:
            raise NotImplementedError('parsing function double star arguments')
        func = node.node
        if isinstance(func, ast.Name) and func.name in callfunc_map:
            self.start(callfunc_map[func.name])
            for child in node.args:
                self.visit(child)
            self.end()
            return
        self.start(APPLY)
        self.visit(func)
        self.start(TUPLE)
        for child in node.args:
            self.visit(child)
        self.end(tuple_head = TUPLE)
        self.end()

    def visitSlice(self, node):
        n = ast.Subscript(node.expr, compiler.consts.OP_APPLY, [ast.Sliceobj(node.asList()[2:])])
        self.visit(n)

    def visitSliceobj(self, node):
        childs = list(node.asList())
        childs.extend([None]*(3-len(childs)))
        assert len(childs)==3,`childs`
        self.start(SLICE)
        for child in childs:
            if child is None:
                self.add(SPECIAL, child)
            else:
                self.visit(child)
        self.end()

    def visitEllipsis(self, node):
        self.add(SPECIAL, Ellipsis)

    def visitDict(self, node):
        self.start(DICT)
        for k,v in node.items:
            self.start(TUPLE)
            self.visit(k)
            self.visit(v)
            self.end(tuple_head = TUPLE)
            # collect key and value to a 2-tuple:
            #self.stack[-1][1][-2] = tuple(self.stack[-1][1][-2:])
            #del self.stack[-1][1][-1]
        self.end()

    def visitGetattr(self, node):
        self.start(ATTR)
        self.visit(node.expr)
        self.visit(ast.Name(node.attrname))
        self.end()

    def visitKeyword(self, node):
        self.start(KWARG)
        self.visit(ast.Name(node.name))
        self.visit(node.expr)
        self.end()

class ast_Pair(ast.Tuple):
    pass

def string2Verbatim(expr):
    """ Parse string expr to Verbatim.
    """
    node = compiler.parse(expr)
    return compiler.walk(node, VerbatimWalker()).stack.pop()
