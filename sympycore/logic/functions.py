
from ..core import classes
from .algebra import Logic

def Eq(a, b):
    """ Return relational expression ``a==b``.
    """
    c = classes.Calculus.convert
    return Logic.Eq(c(a), c(b))

def Ne(a, b):
    """ Return relational expression ``a!=b``.
    """
    c = classes.Calculus.convert
    return Logic.Ne(c(a), c(b))

def Lt(a, b):
    """ Return relational expression ``a<b``.
    """
    c = classes.Calculus.convert
    return Logic.Lt(c(a), c(b))

def Le(a, b):
    """ Return relational expression ``a<=b``.
    """
    c = classes.Calculus.convert
    return Logic.Le(c(a), c(b))

def Gt(a, b):
    """ Return relational expression ``a>b``.
    """
    c = classes.Calculus.convert
    return Logic.Gt(c(a), c(b))

def Ge(a, b):
    """ Return relational expression ``a>=b``.
    """
    c = classes.Calculus.convert
    return Logic.Ge(c(a), c(b))

def IsElement(a, b):
    """ Return relational expression ``a in b``.
    """
    c = classes.Calculus.convert
    #XXX: b = classes.Set.convert(b)
    return Logic.IsElement(c(a), b)

def Not(a):
    """ Return boolean expression ``not a``.
    """
    return Logic.Not(Logic.convert(a))

def And(*seq):
    """ Return boolean expression ``x1 and x2 and ...``.
    """
    c = Logic.convert
    return Logic.And(*map(c, seq))

def Or(*seq):
    """ Return boolean expression ``x1 or x2 or ...``.
    """
    c = Logic.convert
    return Logic.Or(*map(c, seq))
