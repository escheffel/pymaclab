#
# Author: Pearu Peterson
# Created: February 2008
#
""" Provides PolynomialRing class.
"""

__docformat__ = "restructuredtext"
__all__ = "PolynomialRing"

from ..core import classes, IntegerList, Expr
from ..utils import SYMBOL, NUMBER, ADD, MUL, POW, SPARSE_POLY, SUB, DIV
from ..basealgebra.algebra import Algebra
from ..ring import CommutativeRing
from ..basealgebra.verbatim import Verbatim
from ..arithmetic.numbers import div
from ..arithmetic.number_theory import multinomial_coefficients

def cmp_symbols(x, y):
    return cmp(str(x), str(y))

class PolynomialRingFactory(type):
    """ Factory of polynomial rings with symbols and coefficient ring.
    """
    def __new__(typ, name, bases, attrdict):
        if not attrdict.has_key('ring'):
            attrdict['ring'] = classes.Calculus
        if not attrdict.has_key('variables'):
            attrdict['variables'] = ()
            attrdict['nvars'] = 0
        cls = type.__new__(typ, name, bases, attrdict)
        cls.zero = cls.Number(0)
        cls.one = cls.Number(1)
        return cls

    def __eq__(self, other):
        if isinstance(other, PolynomialRingFactory):
            return self.ring==other.ring and self.variables==other.variables
        return False

    def __ne__(self, other):
        return not self==other

    def __getitem__(self, ring_info, cache={}):
        """ Return a new polynomial ring class

        Examples::

          PolynomialRing[<seq of variables>, <coefficient ring>]
          PolynomialRing[<n>, <coefficient ring>] is
                PolynomialRing[['X%s'%i for i in range(<n>)], <coefficient ring>]
          PolynomialRing[<variable info>] is
                PolynomialRing[<variable info>, Calculus]
        """
        if isinstance(ring_info, (int, long)):
            nvars = ring_info
            assert nvars>=0,`nvars`
            variables = ['X%s'%i for i in range(nvars)]
            ring = classes.Calculus
        elif isinstance(ring_info, tuple) and isinstance(ring_info[-1], type):
            if len(ring_info)==2:
                var_info, ring = ring_info
                if isinstance(var_info,(int, long)):
                    nvars = var_info
                    assert nvars>=0,`nvars`
                    variables = ['X%s'%i for i in range(nvars)]
                elif isinstance(var_info,(tuple, list, set)):
                    variables = var_info
                elif isinstance(var_info, str):
                    variables = var_info,
                else:
                    raise TypeError(`ring_info`)
            else:
                variables = ring_info[:-1]
                ring = ring_info[-1]
        elif isinstance(ring_info, (tuple, list, set)):
            variables = ring_info
            ring = classes.Calculus
        elif isinstance(ring_info, str):
            variables = ring_info,
            ring = classes.Calculus        
        else:
            raise TypeError(`ring_info`)
        variables = tuple(sorted(variables, cmp=cmp_symbols))
        nvars = len(variables)

        r = None #cache.get((variables, ring))
        if r is None:
            name = '%s[%s, %s]' % (self.__name__, tuple(variables), ring.__name__)
            r = PolynomialRingFactory(name,
                                      (self,),
                                      dict(nvars=nvars, ring = ring,
                                           variables = variables))
            #cache[variables, ring] = r
        return r

    def is_subring(self, other):
        """ Check if self contains other as a subring, i.e. whether
        the other instances can be converted to self instances.
        """
        if not isinstance(other, PolynomialRingFactory):
            return False
        if self.ring != other.ring:
            return False
        return set(other.variables).issubset(self.variables)

AdditiveTuple = IntegerList

class PolynomialRing(CommutativeRing):
    """ Base class to polynomial rings that holds polynomial information
    using pairs ``(<exponents>: <coefficient>)`` stored in Python dictionary.

    Suitable for representing sparse multivariate polynomials.
    """

    __slots__ = ['_degree', '_ldegree']
    _degree = None
    _ldegree = None
    _str_value = None
    
    __metaclass__ = PolynomialRingFactory

    __str__ = CommutativeRing.__str__
    __repr__ = CommutativeRing.__repr__

    @classmethod
    def convert(cls, data, typeerror=True):
        if isinstance(data, list):
            if cls.nvars==1 and data:
                if type(data[0]) not in [tuple,list,set]:
                    data = [(i,c) for (i,c) in enumerate(data) if c]
            data = dict(data)
            return cls(SPARSE_POLY, data)
        elif isinstance(data, dict):
            return cls(SPARSE_POLY, data)
        r = super(CommutativeRing, cls).convert(data, typeerror=True)
        return r

    @classmethod
    def get_operand_algebra(cls, head, index=0):
        if head in [ADD, SUB, MUL]:
            return cls
        if head is POW:
            if index==0:
                return cls
            if index==1:
                return int
        return cls.handle_get_operand_algebra_failure(head, index)
    
    def __eq__(self, other):
        if type(other)==type(self):
            return self.pair == other.pair
        return CommutativeRing.__eq__(self, other)

    @classmethod
    def Symbol(cls, obj):
        """ Return symbol element of a polynomial ring. The result
        may be an instance of a super polynomial ring.

        Examples::

          r = PolynomialRing['x']
          r.Symbol('x') -> r({1:1})
          r.Symbol('y') -> PolynomialRing['x','y']({(0,1):1})
        """
        try:
            i = list(cls.variables).index(obj)
        except ValueError:
            i = None
        if i is None:
            cls = PolynomialRing[cls.variables+(obj,), cls.ring]
            i = list(cls.variables).index(obj)
        l = [0]*cls.nvars
        l[i] = 1
        return cls(SPARSE_POLY, {AdditiveTuple(l):1})

    @classmethod
    def convert_symbol(cls, obj):
        try:
            i = list(cls.variables).index(obj)
        except ValueError:
            i = None
        if i is None:
            cls = PolynomialRing[cls.variables+(obj,), cls.ring]
            i = list(cls.variables).index(obj)
        l = [0]*cls.nvars
        l[i] = 1
        return cls(SPARSE_POLY, {AdditiveTuple(l):1})

    @classmethod
    def Number(cls, obj):
        """ Return number element of a polynomial ring.

        Examples::

          r = PolynomialRing['x']
          r.Number(2) -> r({0:2})
        """
        data = {AdditiveTuple((0,)*cls.nvars): obj} if obj else {}        
        return cls(SPARSE_POLY, data)

    @classmethod
    def Add(cls, *seq):
        r = cls(SPARSE_POLY, {})
        for t in seq:
            tcls = t.__class__
            if cls==tcls:
                iadd_POLY_POLY(r, t, cls)
            elif cls.is_subring(tcls):
                t = t.as_algebra(cls)
                assert cls==t.__class__,`cls,t`
                iadd_POLY_POLY(r, t, cls)
            elif tcls.is_subring(cls):
                cls = tcls
                r = r.as_algebra(cls)
                iadd_POLY_POLY(r, t, cls)
            elif cls.ring==tcls.ring:
                cls = PolynomialRing[cls.variables+tcls.variables, cls.ring]
                r = r.as_algebra(cls)
                t = t.as_algebra(cls)
                iadd_POLY_POLY(r, t, cls)
            else:
                raise NotImplementedError(`r, t`)
        if not r.data:
            return cls.Number(0)
        return r

    @classmethod
    def Mul(cls, *seq):
        r = cls.one
        for t in seq:
            tcls = t.__class__
            if cls==tcls:
                r = mul_POLY_POLY(r, t, cls)
            elif cls.is_subring(tcls):
                t = t.as_algebra(cls)
                assert cls==t.__class__,`cls,t`
                r = mul_POLY_POLY(r, t, cls)
            elif tcls.is_subring(cls):
                cls = tcls
                r = r.as_algebra(cls)
                r = mul_POLY_POLY(r, t, cls)
            elif cls.ring==tcls.ring:
                cls = PolynomialRing[cls.variables+tcls.variables, cls.ring]
                r = r.as_algebra(cls)
                t = t.as_algebra(cls)
                r = mul_POLY_POLY(r, t, cls)
            else:
                raise NotImplementedError(`r, t`)
        if not r.data:
            return cls.Number(1)
        return r

    @classmethod
    def Pow(cls, base, exp):
        r = cls.__pow__(base, exp)
        if r is NotImplemented:
            raise NotImplementedError(`base,exp`)
        return r

    @classmethod
    def convert_coefficient(cls, obj, typeerror=True):
        if not isinstance(obj, cls.ring):
            r = cls.ring.convert_coefficient(obj, typeerror=False)
            if r is not None:
                return r
            r = cls.ring.convert(obj, typeerror=False)
            if r is not None:
                return r
        if typeerror:
            raise TypeError('%s.convert_coefficient: failed to convert %s instance'\
                            ' to coefficient algebra, expected int|long object'\
                            % (cls.__name__, obj.__class__.__name__))
        else:
            return NotImplemented

    def as_verbatim(self):
        data = self.data
        symbols = [Verbatim.convert(v) for v in self.variables]
        terms = []
        for exps in sorted(data.keys(), reverse=True):
            coeff = data[exps]
            if coeff==1:
                factors = []
            else:
                factors = [Verbatim.convert(coeff)]
            if not isinstance(exps, tuple):
                exps = exps,
            for s,e in zip(symbols, exps):
                if e:
                    if e==1:
                        factors.append(s)
                    else:
                        factors.append(s**e) #XXX: ????
            if not factors:
                terms.append(Verbatim.convert(coeff))
            elif len(factors)==1:
                terms.append(factors[0])
            else:
                terms.append(Verbatim(MUL,tuple(factors)))
        if not terms:
            return Verbatim.convert(0)
        if len(terms)==1:
            return terms[0]
        return Verbatim(ADD, tuple(terms))

    def as_algebra(self, target_cls):
        cls = self.__class__
        if cls is target_cls:
            return self
        if issubclass(target_cls, PolynomialRing):
            target_ring = target_cls.ring
            if cls.variables == target_cls.variables:
                data = dict([(e,target_ring.convert(c)) for e,c in self.data.iteritems()])
                return target_cls(data)
            if cls.ring == target_ring:
                new_cls = PolynomialRing[set(cls.variables + target_cls.variables), target_ring]
                new_variables = list(new_cls.variables)
                indices = [new_variables.index(v) for v in cls.variables]
                data = {}
                n = new_cls.nvars
                for e,c in self.data.iteritems():
                    new_e = [0] * n
                    if not isinstance(e, tuple):
                        e = e,
                    for i,j in enumerate(indices):
                        new_e[j] = e[i]
                    data[AdditiveTuple(new_e)] = c
                return new_cls(data)
        return NotImplemented

    def __neg__(self):
        head, data = self.pair
        return type(self)(head, dict([(e,-c) for e, c in data.iteritems()]))
        
    def __add__(self, other):
        cls = self.__class__
        other_cls = other.__class__
        if cls == other_cls and self.head is other.head:
            return add_POLY_POLY(self, other, cls)
        other = self.convert(other)
        if other is NotImplemented:
            return other
        other_cls = other.__class__
        if cls == other_cls and self.head is other.head:
            return add_POLY_POLY(self, other, cls)
        elif self.nvars < other.nvars:
            return other + self
        return self.head.add(cls, self, other)

    __radd__ = __add__

    def __mul__(self, other):
        cls = self.__class__
        other_cls = other.__class__
        if cls == other_cls and self.head is other.head:
            return mul_POLY_POLY(self, other, cls)
        other2 = self.convert_coefficient(other, typeerror=False)
        if other2 is not NotImplemented:
            return mul_POLY_COEFF(self, other2, cls)
        other = self.convert(other, typeerror=False)
        if other is NotImplemented:
            return other
        other_cls = other.__class__
        if cls == other_cls and self.head is other.head:
            return mul_POLY_POLY(self, other, cls)
        elif self.nvars < other.nvars:
            return other * self
        return self.head.commutative_mul(cls, self, other)

    __rmul__ = __mul__

    def __divmod__(self, other):
        cls = self.__class__
        other_cls = other.__class__
        if cls == other_cls:        
            if cls.nvars==1:
                return divmod_POLY1_POLY1_SPARSE(self, other, cls)
        return NotImplemented

    def __div__(self, other):
        cls = self.__class__
        r = self.convert_coefficient(other, typeerror=False)
        if r is not NotImplemented:
            return mul_POLY_COEFF(self, div(1, r), cls)        
        return divmod(self, other)[0]

    def __mod__(self, other):
        return divmod(self, other)[1]
    
    def __pow__(self, other):
        if isinstance(other, type(self)):
            if other.head is NUMBER:
                other = other.data
            if other.head is SPARSE_POLY:
                data = other.data
                if len(data)==1:
                    exps, coeff = data.items()[0]
                    if exps.data==[0]*len(exps):
                        other = coeff
        if isinstance(other,(int, long)):
            return pow_POLY_INT(self, other, self.__class__)
        print `self`, `other.pair`
        return NotImplemented

    def __call__(self, *args):
        # XXX: Use Horner scheme.
        r = 0
        if self.nvars==1:
            for exps, coeff in self.data.iteritems():
                r = coeff * args[0]**exps + r
        else:
            for exps, coeff in self.data.iteritems():
                r = reduce(lambda x,y: x*y, [a**e for a,e in zip(args,exps)], coeff) + r
        return r

    @property
    def degree(self):
        d = self._degree
        if d is None:
            data = self.data
            if not data:
                self._degree = d = 0
            elif self.nvars==0:
                self._degree = d = 0
            elif self.nvars==1:
                self._degree = d = max(data.keys())
            else:
                self._degree = d = map(max,zip(*data.keys()))
        return d

    @property
    def ldegree(self):
        d = self._ldegree
        if d is None:
            data = self.data
            if not data:
                self._ldegree = d = 0
            elif self.nvars==0:
                self._ldegree = d = 0
            elif self.nvars==1:
                self._ldegree = d = min(data.keys())
            else:
                self._ldegree = d = map(min,zip(*data.keys()))
        return d

    @property
    def coeff(self):
        data = self.data
        nvars = self.nvars
        if nvars==1:
            if data:
                return data[self.degree]
            return self.zero
        raise NotImplementedError(`self,nvars`)

    def diff(self, index=0): # TODO: don't use index as the order of variables is automatically sorted
        if self.nvars==0:
            return self.zero
        head, data = self.pair
        d = {}
        if self.nvars==1:
            assert index==0,`index`
            for exps, coeff in data.iteritems():
                if exps:
                    d[exps-1] = coeff * exps
        else:
            assert 0<=index<self.nvars, `index`
            for exps, coeff in data.iteritems():
                e = exps[index]
                if e:
                    d[exps.add(index,-1)] = coeff * e
        return type(self)(head, d)

    def variable_diff(self, variable, order=1):
        try:
            index = list(self.variables).index(variable)
        except ValueError:
            index = None
        if index is not None:
            return self.head.diff_index(type(self), self.data, self, index, order=order)
        raise NotImplementedError(`self.variables, variable, index`)

    def variable_integrate(self, variable, *bounds):
        """ Return integral over variable.

        poly.variable_integrate(x) -> indefinite integral
        poly.variable_integrate(x, a, b) -> definite integral
        """
        try:
            index = list(self.variables).index(variable)
        except ValueError:
            index = None
        if index is not None:
            indef_integral = self.head.integrate_indefinite_index(type(self), self.data, self, index)
            if bounds:
                low, high = bounds
                return indef_integral.variable_subs(variable, high) - indef_integral.variable_subs(variable, low)
            return indef_integral
        raise NotImplementedError(`self.variables, variable, index`)
    
    def variable_subs(self, variable, newexpr):
        """ Substitute variable with newexpr that must be instance of the same ring.
        """
        cls = type(self)
        newexpr = cls(newexpr)
        try:
            index = list(self.variables).index(variable)
        except ValueError:
            index = None
        if index is not None:
            head, data = self.pair
            result = cls.Number(0)
            variables = cls.variables
            for exps, coeff in data.iteritems():
                term = cls.Number(1)
                for i,exp in enumerate(exps):
                    if exp:
                        if i==index:
                            term *= newexpr**exp
                        else:
                            term *= cls.Symbol(variables[i])**exp
                result += term * cls.Number(coeff)
            return result
        raise NotImplementedError(`self.variables, variable, index`)

def divmod_POLY1_POLY1_SPARSE(lhs, rhs, cls):
    d2 = rhs.degree
    c2 = rhs.coeff
    if not c2:
        raise ZeroDivisionError, "polynomial division"
    other_iter = rhs.data.iteritems
    dq = {}
    dr = dict(lhs.data)
    dseq = [0] + dr.keys()
    dr_get = dr.get
    dseq_append = dseq.append
    dseq_remove = dseq.remove
    while 1:
        d1 = max(dseq)
        if d1 < d2:
            return cls(SPARSE_POLY, dq), cls(SPARSE_POLY, dr)
        c1 = dr[d1]
        e = d1 - d2
        c = div(c1, c2)
        dq[e] = c
        nc = -c
        for e0,c0 in other_iter():
            e0 = e0 + e
            c0 = nc * c0
            b = dr_get(e0)
            if b is None:
                dr[e0] = c0
                dseq_append(e0)
            else:
                b = b + c0
                if b:
                    dr[e0] = b
                else:
                    del dr[e0]
                    dseq_remove(e0)

def divmod_POLY1_POLY1_DENSE(lhs, rhs, cls):
    if not rhs.coeff:
        raise ZeroDivisionError, "polynomial division"
    n = lhs.degree
    nv = rhs.degree
    u = [0]*(n+1)
    v = [0]*(nv+1)
    for e,c in lhs.data.iteritems():
        u[e] = c
    for e,c in rhs.data.iteritems():
        v[e] = c
    r = list(u)
    q = [0] * (len(r)+1)
    for k in range(n-nv, -1, -1):
        q[k] = div(r[nv+k], v[nv])
        for j in range(nv+k-1, k-1, -1):
            r[j] -= q[k]*v[j-k]
    for j in range(nv, n+1, 1):
        r[j] = 0
    q, r = cls(q), cls(r)
    return q, r

def pow_POLY_INT(base, exp, cls):
    if exp==0:
        return cls.one
    if exp==1:
        return base
    if exp < 0:
        return NotImplemented
    data = base.data
    if len(data)==1:
        exps, coeff = data.items()[0]
        return cls(SPARSE_POLY, {exps * exp: coeff ** exp})
    nvars = cls.nvars
    d = {}
    items = data.items()
    m = len(data)
    for k,c_kn in multinomial_coefficients(m, exp).iteritems():
        new_exps = AdditiveTuple((0,)*nvars)
        new_coeff = c_kn
        for i,e in enumerate(k):
            if e:
                exps, coeff = items[i]
                new_exps += exps * e
                new_coeff *= coeff ** e
        b = d.get(new_exps)
        if b is None:
            d[new_exps] = new_coeff
        else:
            c = b + new_coeff
            if c:
                d[new_exps] = c
            else:
                del d[new_exps]
    return cls(SPARSE_POLY, d)
    

def add_POLY_POLY(lhs, rhs, cls):
    d = dict(lhs.data)
    for exps, coeff in rhs.data.iteritems():
        b = d.get(exps)
        if b is None:
            d[exps] = coeff
        else:
            c = b + coeff
            if c:
                d[exps] = c
            else:
                del d[exps]
    return cls(SPARSE_POLY, d)

def iadd_POLY_POLY(lhs, rhs, cls):
    d = lhs.data
    for exps, coeff in rhs.data.iteritems():
        b = d.get(exps)
        if b is None:
            d[exps] = coeff
        else:
            c = b + coeff
            if c:
                d[exps] = c
            else:
                del d[exps]
    

def mul_POLY_POLY(lhs, rhs, cls):
    d = {}
    for exps1, coeff1 in lhs.data.iteritems():
        for exps2, coeff2 in rhs.data.iteritems():
            exps = exps1 + exps2
            coeff = coeff1 * coeff2
            b = d.get(exps)
            if b is None:
                d[exps] = coeff
            else:
                c = b + coeff
                if c:
                    d[exps] = c
                else:
                    del d[exps]
    return cls(SPARSE_POLY, d)

def mul_POLY_COEFF(lhs, rhs, cls):
    d = {}
    for exps, coeff in lhs.data.iteritems():
        b = d.get(exps)
        if b is None:
            d[exps] = coeff * rhs
        else:
            c = b + coeff * rhs
            if c:
                d[exps] = c
            else:
                del d[exps]
    return cls(SPARSE_POLY, d)
