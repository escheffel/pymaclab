__all__ = ['SUB']

from .base import heads_precedence, ArithmeticHead

class SubHead(ArithmeticHead):

    """
    SubHead represents subtraction n-ary operation where operands is
    given as a n-tuple of expressions.
    """
    op_mth = '__sub__'
    op_rmth = '__rsub__'

    def new(self, cls, operands, evaluate=True):
        n = len(operands)
        if n==1:
            return operands[0]
        if n==0:
            return cls(NUMBER, 0)
        return cls(self, operands)

    def __repr__(self): return 'SUB'

    def reevaluate(self, cls, operands):
        r = operands[0] if operands else cls(NUMBER, 0)
        for op in operands[1:]:
            r -= op
        return r
    
    def data_to_str_and_precedence(self, cls, operands):
        m = len(operands)
        if m==0:
            return '0', heads_precedence.NUMBER
        if m==1:
            op = operands[0]
            return op.head.data_to_str_and_precedence(cls, op.data)
        sub_p = heads_precedence.SUB
        r = ''
        for op in operands:
            t,t_p = op.head.data_to_str_and_precedence(cls, op.data)
            if not r:
                r += '(' + t + ')' if t_p < sub_p else t
            elif t.startswith('-') and t_p > sub_p:
                r += ' + ' + t[1:]
            else:
                r += ' - (' + t + ')' if t_p <= sub_p else ' - ' + t
        return r, sub_p

    def walk(self, func, cls, data, target):
        l = []
        flag = False
        for op in data:
            o = op.head.walk(func, cls, op.data, op)
            if op is not o:
                flag = True
            l.append(o)
        if flag:
            r = SUB.new(cls, l)
            return func(cls, r.head, r.data, r)
        return func(cls, self, data, target)

SUB = SubHead()
