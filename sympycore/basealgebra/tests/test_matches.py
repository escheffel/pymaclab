
from sympycore import CollectingField as Algebra

Symbol = Algebra.Symbol
Number = Algebra.Number
Add = Algebra.Add
Mul = Algebra.Mul
Pow = Algebra.Pow
Terms = Algebra.Terms
Factors = Algebra.Factors

def test_symbol():
    p = Symbol('p')
    s = Symbol('s')
    t = Symbol('t')
    assert s.matches(s)=={}
    assert s.matches(t)==None
    assert s.matches(t,{},([s,],[True,]))=={s:t}
    assert s.matches(t,{},([s,t],[True,True]))==None

def test_number():
    s = Symbol('s')
    n = Number(2)
    assert n.matches(2)=={}
    assert n.matches(3)==None
    assert n.matches(s)==None
    assert n.matches(s+2)==None

def test_wild():
    w = Symbol('w')
    s = Symbol('s')
    wargs = [w],[True]
    assert w.matches(Number(2),{},wargs)=={w:2}
    assert w.matches(s,{},wargs)=={w:s}
    assert w.matches(w,{},wargs)==None
    assert w.matches(s+2,{},wargs)=={w:s+2}
    assert w.matches(2*s,{},wargs)=={w:2*s}
    assert w.matches(s**2,{},wargs)=={w:s**2}

def test_symbol():
    s = Symbol('s')
    assert s.matches(s)=={}
    assert s.matches(2)==None
    assert s.matches(2+s)==None
    assert s.matches(2*s)==None
    assert s.matches(s**2)==None

def test_term():
    s = Symbol('s')
    p = 2*s
    assert p.matches(2*s)=={}
    assert p.matches(3*s)==None
    assert p.matches(s)==None
    assert p.matches(Number(2))==None
    assert p.matches(s**2)==None

def _test_wild_term():
    w = Symbol('w')
    p = 2*w
    s = Symbol('s')
    t = Symbol('t')
    wargs = {},([w],[True])
    assert p.matches(Number(1),*wargs)=={w:Number(1)/2}
    assert p.matches(Number(2),*wargs)=={w:1}
    assert p.matches(2*s,*wargs)=={w:s}
    assert p.matches(3*s,*wargs)=={w:s*Number(3)/2}
    assert p.matches(t*s,*wargs)=={w:t*s/2}
    assert p.matches(s**2,*wargs)=={w:s**2/2}
    m = p.matches(2*s+2,*wargs)
    assert m is not None and m[w]==(2*(s+1))/2
    assert p.matches(2*s+4,*wargs)=={w:(s+2)*2/2}
    assert p.matches(2*s+5,*wargs)=={w:(2*s+Number(5))/2}
    assert p.matches(2*s+t,*wargs)=={w:(2*s+t)/2}
    assert p.matches(2*s-2*t,*wargs)=={w:(s-t)*2/2}

def _test_wild_symbol_term():
    w = Symbol('w')
    s = Symbol('s')
    t = Symbol('t')
    p = s+w
    wargs = {},([w],[True])
    assert p.matches(s+2,*wargs)=={w:2}
    assert p.matches(t+2,*wargs)=={w:t+2-s}

def _test_wild_wild_term():
    w1 = Symbol('w1')
    w2 = Symbol('w2')

    p = w1 + 2*w2
    s = Symbol('s')
    t = Symbol('t')
    wargs = {},([w1,w2],[True,True])
    assert p.matches(Number(2),*wargs) in [{w2:0,w1:2},{w2:1,w1:0}]
    assert p.matches(2*s+t+2,*wargs) in [{w2:1+s,w1:t},{w1:2*s+t,w2:1},{w2:s,w1:t+2},
                                         {w1:2+2*s, w2:t/2}]

def _test_wild_factor():
    w = Symbol('w')
    p = w**2
    s = Symbol('s')
    t = Symbol('t')
    wargs = {},([w],[True])
    #assert p.matches(Number(2),*wargs)=={w:Number(2)**(Number(1)/2)}
    #assert p.matches(Number(4),*wargs)=={w:2}
    #assert p.matches(Number(16),*wargs)=={w:4}
    #assert p.matches(Number(9),*wargs)=={w:3}
    #assert p.matches(Number(8),*wargs)=={w:2*Number(2)**(Number(1)/2)}
    assert p.matches(s,*wargs)==None
    assert p.matches(s**2,*wargs)=={w:s}
    assert p.matches(s**3,*wargs)==None
    #assert p.matches(s**4,*wargs)=={w:s**2}
    assert p.matches(s+2,*wargs)==None
    assert p.matches(s*2,*wargs)==None
    assert p.matches(s**2*2,*wargs)==None
    #assert p.matches(s**2*4,*wargs)=={w:2*s}
    #assert p.matches(s**2*t**2,*wargs)=={w:s*t}
    #assert p.matches(4*s**2*t**2,*wargs)=={w:2*s*t}
    #assert p.matches(s**4*t**4,*wargs)=={w:(s*t)**2}
    #assert p.matches(s**2*t**4,*wargs)=={w:s*t**2}
    assert p.matches(s**2*t**3,*wargs)==None
    #assert p.matches(s**2*t**-4,*wargs)=={w:s*t**-2}


def _test_wild_symbol_factor():
    w = Symbol('w')
    s = Symbol('s')
    t = Symbol('t')
    p = s*w
    wargs = {},([w],[True])
    assert p.matches(Number(1),*wargs)=={w:1/s}
    assert p.matches(s,*wargs)=={w:1}
    assert p.matches(2+t,*wargs)=={w:(2+t)/s}

def test_symbol2():
    x = Symbol('x')
    a,b,c,p,q = map(Symbol, 'abcpq')
    e = x
    assert e.match(x) == {}
    assert e.match(a,a) == {a: x}

    e = Number(5)
    assert e.match(c,c) == {c: 5}
    assert e.match(e) == {}
    assert e.match(e+1) == None

def _test_add():
    x,y,a,b,c = map(Symbol, 'xyabc')
    p,q,r = map(Symbol, 'pqr')

    e = a+b
    assert e.match(p+b,p) == {p: a}
    assert e.match(p+a,p) == {p: b}

    e = 1+b
    assert e.match(p+b,p) == {p: 1}

    e = a+b+c
    assert e.match(a+p+c,p) == {p: b}
    assert e.match(b+p+c,p) == {p: a}

    e = a+b+c+x
    assert e.match(a+p+x+c,p) == {p: b}
    assert e.match(b+p+c+x,p) == {p: a}
    assert e.match(b) == None
    assert e.match(b+p,p) == {p: a+c+x}
    assert e.match(a+p+c,p) == {p: b+x}
    assert e.match(b+p+c,p) == {p: a+x}

    e = 4*x+5
    assert e.match(3*x+p,p) == {p: x+5}

    assert e.match(4*x+p,(p,lambda expr: not expr.args)) == {p: 5}
    assert e.match(p*x+5,(p,lambda expr: not expr.args)) == {p: 4}
    assert e.match(p*x+q,(p,lambda expr: not expr.args),(q,lambda expr: not expr.args)) == {p: 4, q: 5}

    e = 4*x+5*y+6
    assert e.match(p*x+q*y+r,(p,lambda expr: not expr.args),
                   (q,lambda expr: not expr.args),
                   (r,lambda expr: not expr.args)) == {p: 4, q: 5, r: 6}
