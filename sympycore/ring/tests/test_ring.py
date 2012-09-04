
from sympycore import *

operation_results = '''
+(x):x
-(x):-x
(x)+(x):2*x
(x)-(x):0
(x)*(x):x**2
(x)/(x):1
(x)**(x):x**x
(x)+(y):y + x
(x)-(y):-y + x
(x)*(y):x*y
(x)/(y):x/y
(x)**(y):x**y
+(y):y
-(y):-y
(y)+(x):y + x
(y)-(x):y - x
(y)*(x):y*x
(y)/(x):y/x
(y)**(x):y**x
(y)+(y):2*y
(y)-(y):0
(y)*(y):y**2
(y)/(y):1
(y)**(y):y**y
(x)+(1):x + 1
(x)-(1):x - 1
(x)*(1):x
(x)/(1):x
(x)**(1):x
(y)+(1):y + 1
(y)-(1):y - 1
(y)*(1):y
(y)/(1):y
(y)**(1):y
(1)+(x):x + 1
(1)-(x):-x + 1
(1)*(x):x
(1)/(x):1/x
(1)**(x):1
(1)+(y):y + 1
(1)-(y):-y + 1
(1)*(y):y
(1)/(y):1/y
(1)**(y):1
(x)+(0):x
(x)-(0):x
(x)*(0):0
(x)/(0):zoo
(x)**(0):1
(y)+(0):y
(y)-(0):y
(y)*(0):0
(y)/(0):zoo
(y)**(0):1
(0)+(x):x
(0)-(x):-x
(0)*(x):0
(0)/(x):0
(0)**(x):0**x
(0)+(y):y
(0)-(y):-y
(0)*(y):0
(0)/(y):0
(0)**(y):0**y
(x)+(3):x + 3
(x)-(3):x - 3
(x)*(3):3*x
(x)/(3):1/3*x
(x)**(3):x**3
(y)+(3):y + 3
(y)-(3):y - 3
(y)*(3):3*y
(y)/(3):1/3*y
(y)**(3):y**3
(1)+(3):4
(1)-(3):-2
(1)*(3):3
(1)/(3):1/3
(1)**(3):1
(0)+(3):3
(0)-(3):-3
(0)*(3):0
(0)/(3):0
(0)**(3):0
+(3):3
-(3):-3
(3)+(x):3 + x
(3)-(x):-x + 3
(3)*(x):3*x
(3)/(x):3/x
(3)**(x):3**x
(3)+(y):y + 3
(3)-(y):-y + 3
(3)*(y):3*y
(3)/(y):3/y
(3)**(y):3**y
(3)+(1):4
(3)-(1):2
(3)*(1):3
(3)/(1):3
(3)**(1):3
(3)+(0):3
(3)-(0):3
(3)*(0):0
(3)/(0):zoo
(3)**(0):1
(3)+(3):6
(3)-(3):0
(3)*(3):9
(3)/(3):1
(3)**(3):27
(x)+(x*y):x + x*y
(x)-(x*y):x - x*y
(x)*(x*y):x**2*y
(x)/(x*y):x/y/x
(x)**(x*y):x**(x*y)
(y)+(x*y):y + x*y
(y)-(x*y):y - x*y
(y)*(x*y):y*x*y
(y)/(x*y):1/x
(y)**(x*y):y**(x*y)
(1)+(x*y):1 + x*y
(1)-(x*y):1 - x*y
(1)*(x*y):x*y
(1)/(x*y):1/y/x
(1)**(x*y):1
(0)+(x*y):x*y
(0)-(x*y):-x*y
(0)*(x*y):0
(0)/(x*y):0
(0)**(x*y):0**(x*y)
(3)+(x*y):3 + x*y
(3)-(x*y):3 - x*y
(3)*(x*y):3*x*y
(3)/(x*y):3/y/x
(3)**(x*y):3**(x*y)
+(x*y):x*y
-(x*y):-x*y
(x*y)+(x):x + x*y
(x*y)-(x):-x + x*y
(x*y)*(x):x*y*x
(x*y)/(x):x*y/x
(x*y)**(x):(x*y)**x
(x*y)+(y):y + x*y
(x*y)-(y):-y + x*y
(x*y)*(y):x*y**2
(x*y)/(y):x
(x*y)**(y):(x*y)**y
(x*y)+(1):1 + x*y
(x*y)-(1):-1 + x*y
(x*y)*(1):x*y
(x*y)/(1):x*y
(x*y)**(1):x*y
(x*y)+(0):x*y
(x*y)-(0):x*y
(x*y)*(0):0
(x*y)/(0):zoo
(x*y)**(0):1
(x*y)+(3):3 + x*y
(x*y)-(3):-3 + x*y
(x*y)*(3):3*x*y
(x*y)/(3):1/3*x*y
(x*y)**(3):(x*y)**3
(x*y)+(x*y):2*x*y
(x*y)-(x*y):0
(x*y)*(x*y):(x*y)**2
(x*y)/(x*y):1
(x*y)**(x*y):(x*y)**(x*y)
(x)+(y*x):x + y*x
(x)-(y*x):x - y*x
(x)*(y*x):x*y*x
(x)/(y*x):1/y
(x)**(y*x):x**(y*x)
(y)+(y*x):y + y*x
(y)-(y*x):y - y*x
(y)*(y*x):y**2*x
(y)/(y*x):y/x/y
(y)**(y*x):y**(y*x)
(1)+(y*x):1 + y*x
(1)-(y*x):1 - y*x
(1)*(y*x):y*x
(1)/(y*x):1/x/y
(1)**(y*x):1
(0)+(y*x):y*x
(0)-(y*x):-y*x
(0)*(y*x):0
(0)/(y*x):0
(0)**(y*x):0**(y*x)
(3)+(y*x):3 + y*x
(3)-(y*x):3 - y*x
(3)*(y*x):3*y*x
(3)/(y*x):3/x/y
(3)**(y*x):3**(y*x)
(x*y)+(y*x):y*x + x*y
(x*y)-(y*x):-y*x + x*y
(x*y)*(y*x):x*y**2*x
(x*y)/(y*x):x*y/x/y
(x*y)**(y*x):(x*y)**(y*x)
+(y*x):y*x
-(y*x):-y*x
(y*x)+(x):x + y*x
(y*x)-(x):-x + y*x
(y*x)*(x):y*x**2
(y*x)/(x):y
(y*x)**(x):(y*x)**x
(y*x)+(y):y + y*x
(y*x)-(y):-y + y*x
(y*x)*(y):y*x*y
(y*x)/(y):y*x/y
(y*x)**(y):(y*x)**y
(y*x)+(1):1 + y*x
(y*x)-(1):-1 + y*x
(y*x)*(1):y*x
(y*x)/(1):y*x
(y*x)**(1):y*x
(y*x)+(0):y*x
(y*x)-(0):y*x
(y*x)*(0):0
(y*x)/(0):zoo
(y*x)**(0):1
(y*x)+(3):3 + y*x
(y*x)-(3):-3 + y*x
(y*x)*(3):3*y*x
(y*x)/(3):1/3*y*x
(y*x)**(3):(y*x)**3
(y*x)+(x*y):y*x + x*y
(y*x)-(x*y):y*x - x*y
(y*x)*(x*y):y*x**2*y
(y*x)/(x*y):y*x/y/x
(y*x)**(x*y):(y*x)**(x*y)
(y*x)+(y*x):2*y*x
(y*x)-(y*x):0
(y*x)*(y*x):(y*x)**2
(y*x)/(y*x):1
(y*x)**(y*x):(y*x)**(y*x)

(x)+(x**2):x + x**2
(x)-(x**2):x - x**2
(x)*(x**2):x**3
(x)/(x**2):1/x
(x)**(x**2):x**x**2
(y)+(x**2):y + x**2
(y)-(x**2):y - x**2
(y)*(x**2):y*x**2
(y)/(x**2):y/x**2
(y)**(x**2):y**x**2
(1)+(x**2):1 + x**2
(1)-(x**2):1 - x**2
(1)*(x**2):x**2
(1)/(x**2):1/x**2
(1)**(x**2):1
(0)+(x**2):x**2
(0)-(x**2):-x**2
(0)*(x**2):0
(0)/(x**2):0
(0)**(x**2):0**x**2
(3)+(x**2):3 + x**2
(3)-(x**2):3 - x**2
(3)*(x**2):3*x**2
(3)/(x**2):3/x**2
(3)**(x**2):3**x**2
(x*y)+(x**2):x**2 + x*y
(x*y)-(x**2):-x**2 + x*y
(x*y)*(x**2):x*y*x**2
(x*y)/(x**2):x*y/x**2
(x*y)**(x**2):(x*y)**x**2
(y*x)+(x**2):x**2 + y*x; y*x + x**2
(y*x)-(x**2):-x**2 + y*x;y*x - x**2
(y*x)*(x**2):y*x**3
(y*x)/(x**2):y/x
(y*x)**(x**2):(y*x)**x**2
+(x**2):x**2
-(x**2):-x**2
(x**2)+(x):x + x**2
(x**2)-(x):-x + x**2
(x**2)*(x):x**3
(x**2)/(x):x
(x**2)**(x):(x**2)**x
(x**2)+(y):y + x**2
(x**2)-(y):-y + x**2
(x**2)*(y):x**2*y
(x**2)/(y):x**2/y
(x**2)**(y):(x**2)**y
(x**2)+(1):1 + x**2
(x**2)-(1):-1 + x**2
(x**2)*(1):x**2
(x**2)/(1):x**2
(x**2)**(1):x**2
(x**2)+(0):x**2
(x**2)-(0):x**2
(x**2)*(0):0
(x**2)/(0):zoo
(x**2)**(0):1
(x**2)+(3):3 + x**2
(x**2)-(3):-3 + x**2
(x**2)*(3):3*x**2
(x**2)/(3):1/3*x**2
(x**2)**(3):x**6
(x**2)+(x*y):x**2 + x*y;x*y + x**2
(x**2)-(x*y):x**2 - x*y
(x**2)*(x*y):x**3*y
(x**2)/(x*y):x**2/y/x
(x**2)**(x*y):(x**2)**(x*y)
(x**2)+(y*x):x**2 + y*x;y*x + x**2
(x**2)-(y*x):x**2 - y*x
(x**2)*(y*x):x**2*y*x
(x**2)/(y*x):x/y
(x**2)**(y*x):(x**2)**(y*x)
(x**2)+(x**2):2*x**2
(x**2)-(x**2):0
(x**2)*(x**2):x**4
(x**2)/(x**2):1
(x**2)**(x**2):(x**2)**x**2

(x)+(2*x):3*x
(x)-(2*x):-x
(x)*(2*x):2*x**2
(x)/(2*x):1/2
(x)**(2*x):x**(2*x)
(y)+(2*x):y + 2*x
(y)-(2*x):y - 2*x
(y)*(2*x):2*y*x
(y)/(2*x):1/2*y/x
(y)**(2*x):y**(2*x)
(1)+(2*x):2*x + 1
(1)-(2*x):-2*x + 1
(1)*(2*x):2*x
(1)/(2*x):1/2/x
(1)**(2*x):1
(0)+(2*x):2*x
(0)-(2*x):-2*x
(0)*(2*x):0
(0)/(2*x):0
(0)**(2*x):0**(2*x)
(3)+(2*x):2*x + 3
(3)-(2*x):-2*x + 3
(3)*(2*x):6*x
(3)/(2*x):3/2/x
(3)**(2*x):3**(2*x)
(x*y)+(2*x):2*x + x*y
(x*y)-(2*x):-2*x + x*y
(x*y)*(2*x):2*x*y*x
(x*y)/(2*x):1/2*x*y/x
(x*y)**(2*x):(x*y)**(2*x)
(y*x)+(2*x):2*x + y*x
(y*x)-(2*x):-2*x + y*x
(y*x)*(2*x):2*y*x**2
(y*x)/(2*x):1/2*y
(y*x)**(2*x):(y*x)**(2*x)
(x**2)+(2*x):2*x + x**2
(x**2)-(2*x):-2*x + x**2
(x**2)*(2*x):2*x**3
(x**2)/(2*x):1/2*x
(x**2)**(2*x):(x**2)**(2*x)

+(2*x):2*x
-(2*x):-2*x
(2*x)+(x):3*x
(2*x)-(x):x
(2*x)*(x):2*x**2
(2*x)/(x):2
(2*x)**(x):(2*x)**x
(2*x)+(y):y + 2*x
(2*x)-(y):-y + 2*x
(2*x)*(y):2*x*y
(2*x)/(y):2*x/y
(2*x)**(y):(2*x)**y
(2*x)+(1):2*x + 1
(2*x)-(1):2*x - 1
(2*x)*(1):2*x
(2*x)/(1):2*x
(2*x)**(1):2*x
(2*x)+(0):2*x
(2*x)-(0):2*x
(2*x)*(0):0
(2*x)/(0):zoo
(2*x)**(0):1
(2*x)+(3):2*x + 3
(2*x)-(3):2*x - 3
(2*x)*(3):6*x
(2*x)/(3):2/3*x
(2*x)**(3):8*x**3
(2*x)+(x*y):2*x + x*y
(2*x)-(x*y):2*x - x*y
(2*x)*(x*y):2*x**2*y
(2*x)/(x*y):2*x/y/x
(2*x)**(x*y):(2*x)**(x*y)
(2*x)+(y*x):2*x + y*x
(2*x)-(y*x):2*x - y*x
(2*x)*(y*x):2*x*y*x
(2*x)/(y*x):2/y
(2*x)**(y*x):(2*x)**(y*x)
(2*x)+(x**2):2*x + x**2
(2*x)-(x**2):2*x - x**2
(2*x)*(x**2):2*x**3
(2*x)/(x**2):2/x
(2*x)**(x**2):(2*x)**x**2
(2*x)+(2*x):4*x
(2*x)-(2*x):0
(2*x)*(2*x):4*x**2
(2*x)/(2*x):1
(2*x)**(2*x):(2*x)**(2*x)
(x)+(y + x):y + 2*x
(x)-(y + x):-y


'''

def test_operations():
    x,y,z = map(Ring, 'xyz')

    operands = [x,
                y,
                1,
                0,
                Ring(heads.NUMBER, 3),
                x*y,
                y*x,
                x*x,
                2*x,
                ]

    unary_operations = ['+', '-']
    binary_operations = ['+', '-', '*', '/', '**']

    results = {}
    for line in operation_results.split('\n'):
        line = line.strip()
        if ':' not in line: continue
        expr, result = line.split(':')
        for e in expr.split(';'):
            e = e.strip()
            results[e] = [r.strip() for r in result.split(';')]

    for op1 in operands:
        if isinstance(op1, Expr):
            for op in unary_operations:
                expr = '%s(%s)' % (op, op1) 

                try:
                    result = str(eval('%s(op1)' % (op)))
                except Exception, msg:
                    print  expr,'failed with %s' % (msg)
                    raise
                
                if expr not in results:
                    print '%s:%s' % (expr, result)
                    continue
                assert result in results[expr], `results[expr], result`
        for op2 in operands:
            if not (isinstance(op1, Expr) or isinstance(op2, Expr)):
                continue
            for op in binary_operations:
                expr = '(%s)%s(%s)' % (op1, op, op2)

                try:
                    result = str(eval('(op1)%s(op2)' % op))
                except Exception, msg:
                    print  expr,'failed with %s' % (msg)
                    raise
                
                if expr not in results:
                    print '%s:%s' % (expr, result)
                    continue
                assert result in results[expr], `results[expr], result, op1, op2`
    
