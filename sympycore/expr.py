# -*- coding: latin-1 -*-
"""
  This module implements extension type Expr that holds two Python
  objects, head and data, in a pair attribute.

  When adding new features to Expr class, make sure that these are
  added also to extension type Expr in src/expr_ext.c.

C Expr:

>>> from sympycore.expr_ext import *
>>> %timeit Expr(1,2)
10000000 loops, best of 3: 179 ns per loop
>>> e=Expr(1,2)
>>> %timeit h = e.head
10000000 loops, best of 3: 113 ns per loop
>>> %timeit h,d = e.pair
10000000 loops, best of 3: 141 ns per loop

Python Expr:

>>> from sympycore.expr import *
>>> %timeit Expr(1,2)
1000000 loops, best of 3: 988 ns per loop
>>> e=Expr(1,2)
>>> %timeit h = e.head
10000000 loops, best of 3: 119 ns per loop
>>> %timeit h,d = e.pair
10000000 loops, best of 3: 139 ns per loop
"""
# Author: Pearu Peterson
# Created: March 2008

__all__ = ['Expr', 'Pair']

def init_module(m):
    from .core import heads
    for n,h in heads.iterNameValue(): setattr(m, n, h)
    from .arithmetic.numbers import inttypes, numbertypes, try_power, gcd
    m.inttypes = inttypes
    m.numbertypes = numbertypes
    m.try_power = try_power
    m.gcd = gcd

class Expr(object):
    """Represents an symbolic expression in a pair form: (head, data)	
                									
The pair (head, data) is saved in an attribute ``pair``. The parts of
a pair, head and data, can be also accessed via ``head`` and ``data``
attributes, respectively. All three attributes are read-only.
  									
The head part is assumed to be an immutable object.  The data part can
be either an immutable object or Python dictionary.  In the former
case, the hash value of Expr instance is defined as::
  									
    hash((<Expr>.head, <Expr>.

Otherwise, if ``data`` contains a Python dictionary, then the hash
value is defined as::
  									
    hash((<Expr>.head, frozenset(<Expr>.data.items())

If ``data`` is a Python list, then the hash value is::

  hash((<Expr>.head, tuple(<Expr>.data)))

WARNING: the hash value of an Expr instance is computed (and cached)
when it is used as a key to Python dictionary. This means that the
instance content MUST NOT be changed after the hash is computed.  To
check if it is safe to change the ``data`` dictionary, use
``is_writable`` attribute that returns True when the hash value has
not been computed::
  									
    <Expr>.is_writable -> True or False				
  									
There are two ways to access the parts of a Expr instance from	
Python::								
  									
    a = Expr(<head>, <data>)					
    head, data = a.head, a.data     - for backward compatibility	
    head, data = a.pair             - fastest way			

When Expr constructor is called with one argument, say ``x``, then
``<Expr subclass>.convert(x)`` will be returned.

This is Python version of Expr type.
"""

    __slots__ = ['head', 'data', 'pair', '_hash']

    def __init__(self, *args):
        if len(args)==1:
            obj = self.convert(args[0])
            self.pair = obj.pair
            self._hash = obj._hash
        elif len(args)==2:
            self.pair = args
            self._hash = None
        else:
            raise TypeError("%s requires 1 or 2 arguments but got %r" % (type(self), len(args)))
        msg = self.head.is_data_ok(type(self), self.data)
        if msg:
            msg = '%s(head=%s, data=%s): %s' % (type(self).__name__, self.head, self.data, msg)
            #print msg
            raise TypeError(msg)
        
    def __repr__(self):
        return '%s%r' % (type(self).__name__, self.pair)

    def __hash__(self):
        """ Compute hash value.
        
        Different from expr_ext.Expr, an exception is raised when data
        dictionary values contain dictionaries.
        """
        h = self._hash
        if h is None:
            pair = self.pair
            obj = self.as_lowlevel()
            if obj is not pair:
                h = hash(obj)
            else:
                head, data = pair
                tdata = type(data)
                if tdata is dict:
                    h = hash((head, frozenset(data.iteritems())))
                elif tdata is list:
                    h = hash((head, tuple(data)))
                else:
                    h = hash(pair)
            self._hash = h
        return h

    @property
    def is_writable(self, _writable_types = (list, dict)):
        if self._hash is None:
            data = self.pair[1]
            tdata = type(data)
            if tdata in _writable_types:
                return True
            if tdata is Pair:
                return data.is_writable
        return False

    @property
    def head(self):
        return self.pair[0]

    @property
    def data(self):
        return self.pair[1]

    # Pickle support:
    def _sethash(self, hashvalue):
        """ Set hash value for the object.

        If hashvalue==-1, then the hash value will be reset.

        Used by pickle support in sympycore.core._reconstruct. DO NOT
        use this method directly.
        """
        if hashvalue==-1:
            self._hash = None
        else:
            self._hash = hashvalue

    def __reduce__(self):
        # see also _reconstruct function in sympycore/core.py
        version = 3
        from sympycore.core import _reconstruct
        if version==1:
            hashvalue = self._hash
            if hashvalue is None:
                hashvalue = -1
            state = (type(self), self.pair, hashvalue)
        elif version==2 or version==3:
            hashvalue = self._hash
            if hashvalue is None:
                hashvalue = -1
            cls = type(self)
            typ = type(cls)
            try:
                args = typ.__getinitargs__(cls)
            except AttributeError:
                args = None
            if args is None:
                # either metaclass does not define __getinitargs__ method
                # or cls has no metaclass
                state = (cls, self.pair, hashvalue)
            else:
                state = ((typ, args), self.pair, hashvalue)
        else:
            raise NotImplementedError('pickle state version %s' % (version))
        return  _reconstruct, (version, state)

    def __nonzero__(self):
        # Note that `not cls(MUL, [])` would return True while `cls(MUL, [])==1`.
        # So, must use as_lowlevel:
        obj = self.as_lowlevel()
        if obj is not self.pair:
            return not not obj
        return not not self.data

    def as_lowlevel(self):
        """ Return self as low-level object instance that will be used
        in comparison and in hash computation.

        By default, as_lowlevel uses heads to_lowlevel(cls, data, pair)
        method but since as_lowlevel is a most frequently called
        method then for some heads the corresponding code is copied
        here. The default return value is a pair tuple for composite
        objects and data part for atomic objects. The as_lowlevel
        method may also return an Expr instance but not self
        (otherwise infinite recurrsion will occur).

        See __hash__, __nonzero__ method for more details how the
        results of as_lowlevel method are interpreted.
        """
        head, data = pair = self.pair
        if head is NUMBER or head is SYMBOL or head is SPECIAL:
            return data
        elif head is MUL or head is DIV:
            n = len(data)
            if n==0:
                return 1
            if n==1:
                return data[0]
        elif head is ADD or head is SUB:
            n = len(data)
            if n==0:
                return 0
            if n==1:
                return data[0]
        elif head is POW:
            base, exp = data
            if exp==0 or base==1:
                return 1
            if exp==1:
                return base
        elif head is TERM_COEFF:
            term, coeff = data
            if coeff==0:
                return 0
            if coeff==1:
                return term
            if term==1:
                return coeff
        elif head is TERM_COEFF_DICT:
            n = len(data)
            if n==0:
                return 0
            if n==1:
                return type(self)(TERM_COEFF, dict_get_item(data))
        elif head is BASE_EXP_DICT:
            n = len(data)
            if n==0:
                return 1
            if n==1:
                return type(self)(POW, dict_get_item(data))
        else:
            return head.to_lowlevel(type(self), data, pair)
        return pair

    def term_coeff(self):
        head, data = self.pair
        if head is TERM_COEFF:
            return data
        if head is BASE_EXP_DICT:
            cls = type(self)
            coeff = base_exp_dict_get_coefficient(cls, data)
            if coeff is not None:
                d = data.copy()
                del d[coeff]
                r = BASE_EXP_DICT.new(cls, d)
                return r, coeff.data
                t, c = r.head.term_coeff(cls, r)
                return t, c * coeff 
            return self, 1
        if head is NUMBER:
            return self, data
        return self, 1

    def base_exp(self):
        head, data = self.pair
        if head==POW:
            return data
        return self, 1

    for _item in dict(__eq__ = '==', __ne__ = '!=',
                      __lt__ = '<', __le__ = '<=',
                      __gt__ = '>', __ge__ = '>=',
                      ).items():
        exec '''
def %s(self, other):
    if type(self) is type(other):
        other = other.as_lowlevel()
    return self.as_lowlevel() %s other
''' % _item


class Pair(Expr):
    """ Holds a pair that may contain list and dict second element.
    """
    def __init__(self, *args):
        if len(args)==2:
            self.pair = args
            self._hash = None
        else:
            raise TypeError("%s requires 2 arguments but got %r" % (type(self), len(args)))

    def __eq__(self, other):
        return self.pair == other

    def __len__(self):
        return 2
    
    def __getitem__(self, index):
        return self.pair[index]


def dict_mul_item(Algebra, d, key, value):
    c = d.get(key)
    if c is None:
        d[key] = value
    else:
        d[key] = c * value

base_exp_dict_mul_item = dict_mul_item

def term_coeff_dict_mul_item(Algebra, d, key, value):
    return dict_mul_item(Algebra, d, key, value)

def dict_add_item(Algebra, d, key, value):
    c = d.get(key)
    if c is None:
        if value:
            d[key] = value
    else:
        c = c + value
        if c:
            d[key] = c
        else:
            del d[key]

def base_exp_dict_get_coefficient(Algebra, d):
    for k, v in d.iteritems():
        if v is 1 and k.head is NUMBER:
            return k
    return

def base_exp_dict_add_item(Algebra, d, base, exp):
    """
    d is a dictonary that will be updated with base, exp pair.
    Base part must be an Algebra instance while exp part
    can either be a number instance or an expression.

    The dictionary items (base, exp) must satisfy the following
    conditions:
      1) numeric (base, exp) must be completely evaluated, that is,
      if base is numeric then exp is either 1 or rational with denominator
      part that is not a root of base.
      2) all numeric items with exp==1 must be combined to one item
      (coefficient)
      3) items with base==1 and exp==0 must not be present in the dictionary
    """
    base_head, base_data = base.pair
    value = d.get(base)
    if type(exp) is Algebra and exp.head is NUMBER:
        exp = exp.data
    if value is None:
        if base_head is NUMBER:
            assert base
            if base==1:
                return
            if exp==1:
                coeff = base_exp_dict_get_coefficient(Algebra, d)
                if coeff is None:
                    d[base] = 1
                else:
                    del d[coeff]
                    coeff = coeff * base
                    if coeff==1:
                        return
                    d[coeff] = 1
            else:
                assert not isinstance(exp, inttypes),`base, exp`
                d[base] = exp
        else:
            assert exp
            d[base] = exp
    else:
        value = value + exp
        if type(value) is Algebra and value.head is NUMBER:
            value = value.data
        if value:
            if base_head is NUMBER and isinstance(value, numbertypes):
                del d[base]
                r, l = try_power(base_data, value)
                if r!=1:
                    base_exp_dict_add_item(Algebra, d, Algebra(NUMBER, r), 1)
                elif len(l)==1 and l[0]==base_data:
                    d[base] = value
                else:
                    for b, e in l:
                        base_exp_dict_add_item(Algebra, d, Algebra(NUMBER, b), e)
            elif base_head is TERM_COEFF and isinstance(value, inttypes):
                del d[base]
                term, coeff = base_data
                base_exp_dict_add_item(Algebra, d, term, value)
                base_exp_dict_add_item(Algebra, d, Algebra(NUMBER, coeff), value)
            else:
                d[base] = value
        else:
            del d[base]

def term_coeff_dict_add_item(Algebra, d, key, value):
    khead, kdata = key.pair
    if khead is TERM_COEFF:
        term, coeff = kdata
        dict_add_item(Algebra, d, term, value * coeff)
    else:
        assert khead not in [TERM_COEFF_DICT],`Algebra, key.pair`
        dict_add_item(Algebra, d, key, value)

def term_coeff_dict_add_dict(Algebra, dict1, dict2):
    for key, value in dict2.iteritems():
        term_coeff_dict_add_item(Algebra, dict1, key, value)

def dict_get_item(d):
    return d.items()[0]

def dict_add_dict(Algebra, dict1, dict2):
    for key, value in dict2.iteritems():
        dict_add_item(Algebra, dict1, key, value)

def base_exp_dict_add_dict(Algebra, dict1, dict2):
    for key, value in dict2.iteritems():
        base_exp_dict_add_item(Algebra, dict1, key, value)

def base_exp_dict_sub_item(Algebra, dict, key, value):
    return base_exp_dict_add_item(Algebra, dict, key, -value)

def base_exp_dict_sub_dict(Algebra, dict1, dict2):
    for key, value in dict2.iteritems():
        base_exp_dict_sub_item(Algebra, dict1, key, value)

def dict_sub_dict(Algebra, dict1, dict2):
    for key, value in dict2.iteritems():
        dict_add_item(Algebra, dict1, key, -value)

def base_exp_dict_mul_dict(Algebra, d, dict1, dict2):
    for t1,c1 in dict1.iteritems():
        for t2,c2 in dict2.iteritems():
            t = t1 * t2
            t, c = t.term_coeff()
            c12 = c * c1 * c2
            base_exp_dict_add_item(Algebra, d, t, c12)

def exp_coeff_dict_mul_dict(Algebra, d, dict1, dict2):
    for t1,c1 in dict1.iteritems():
        for t2,c2 in dict2.iteritems():
            t = t1 + t2
            c12 = c1 * c2
            dict_add_item(Algebra, d, t, c12)
            
def dict_mul_dict(Algebra, d, dict1, dict2):
    for t1,c1 in dict1.iteritems():
        for t2,c2 in dict2.iteritems():
            t = t1 * t2
            #t, c = t.term_coeff()
            c12 = c1 * c2
            dict_add_item(Algebra, d, t, c12)
    
def dict_mul_value(Algebra, d, value):
    if value==1:
        return
    for t, c in d.items():
        c *= value
        if c:
            d[t] = c
        else:
            del d[t]

term_coeff_dict_mul_dict = dict_mul_dict
base_exp_dict_mul_value = dict_mul_value
term_coeff_dict_mul_value = dict_mul_value

def term_coeff_new(Algebra, data):
    term, coeff = data
    if coeff==1:
        return term
    if term==1:
        return Algebra(NUMBER, coeff)
    if coeff==0:
        return Algebra(NUMBER, 0)
    h,d = term.pair
    if h is TERM_COEFF:
        t, c = d
        return term_coeff_new(Algebra, (t, c * coeff))
    if h is NUMBER:
        return Algebra(NUMBER, d * coeff)
    return Algebra(TERM_COEFF, data)

def term_coeff_dict_new(Algebra, data):
    """
    Return canonicalized TERM_COEFF_DICT expression from data.
    """
    n = len(data)
    if n>1:
        return Algebra(TERM_COEFF_DICT, data)
    if n==0:
        return Algebra(NUMBER, 0)
    return term_coeff_new(Algebra, dict_get_item(data))

def term_coeff_dict(Algebra, expr):
    """
    Return canonicalized TERM_COEFF_DICT expression from existing
    expression:
    
      * if expr has no data, return 0
      * if expr data has one item, return TERM_COEFF expression
    """
    data = expr.data
    n = len(data)
    if n>1:
        return expr
    if n==0:
        return Algebra(NUMBER, 0)
    return term_coeff_new(Algebra, dict_get_item(data))

def pow_new(Algebra, data):
    base, exp = data
    if exp==1:
        return base
    if base==1 or exp==0:
        return Algebra(NUMBER, 1)
    #if type(exp) is Algebra and exp.head is NUMBER:
    #    data = base, exp.data
    return Algebra(POW, data)

def base_exp_dict_new(Algebra, data):
    n = len(data)
    if n==0:
        return Algebra(NUMBER, 1)
    if n==1:
        return pow_new(Algebra, dict_get_item(data))
    coeff = base_exp_dict_get_coefficient(Algebra, data)
    if coeff is None:
        return Algebra(BASE_EXP_DICT, data)
    del data[coeff]
    return term_coeff_new(Algebra, (base_exp_dict_new(Algebra, data), coeff.data))

def add_new(Algebra, data):
    n = len(data)
    if n==0:
        return Algebra(NUMBER, 0)
    if n==1:
        return data[0]
    return Algebra(ADD, data)

def add(Algebra, expr):
    data = expr.data
    n = len(data)
    if n==0:
        return Algebra(NUMBER, 0)
    if n==1:
        return data[0]
    return expr

def mul_new(Algebra, data):
    n = len(data)
    if n==0:
        return Algebra(NUMBER, 1)
    if n==1:
        return data[0]
    return Algebra(MUL, data)
