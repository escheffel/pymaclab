
__all__ = ['Group', 'AdditiveGroup', 'AdditiveAbelianGroup']

from ..basealgebra import Algebra
from ..core import init_module, classes
init_module.import_heads()
init_module.import_numbers()

class Group(Algebra):
    """
    Group represents algebraic group (G, <op>) where <op> is a binary
    operation that is associative, G contains identity element and
    have inverse elements.
    """
    algebra_options = dict(evaluate_addition = True,
                           evaluate_multiplication = True,
                           is_additive_group = False,
                           is_multiplicative_group = False,
                           )

    def to(self, target, *args):
        """ Convert expression to target representation.

        The following targets are recognized:

          EXP_COEFF_DICT - convert expression to exponents-coefficient
            representation, additional arguments are variables. When
            no arguments are specified, variables will be all symbols
            and non-power expressions used in the expression.

          TERM_COEFF_DICT - convert expression to term-coefficient
            representation. Note that the returned result may have
            actual head NUMBER, SYMBOL, TERM_COEFF, POW, or
            BASE_EXP_DICT instead of TERM_COEFF_DICT.
        """
        head, data = self.pair
        if target is EXP_COEFF_DICT:
            return head.to_EXP_COEFF_DICT(type(self), data, self, args or None)
        if target is TERM_COEFF_DICT:
            return head.to_TERM_COEFF_DICT(type(self), data, self)
        if target is ADD:
            return head.to_ADD(type(self), data, self)
        if target is MUL:
            return head.to_MUL(type(self), data, self)
        raise NotImplementedError('%s.to(target=%r)' % (type(self), target))

class AdditiveGroup(Group):
    """
    AdditiveGroup is a group with addition ``+`` operation, 0 is
    identity element and ``-`` denotes inverse. Multiplication ``*``
    is defined only with numbers (integers) that are considered as
    coefficients.

    Warning: AdditiveGroup allows operating with numbers although
    mathematically strictly speaking they are not defined, in general.
    E.g. the only number that belongs to additive group, is 0. Hence
    ``1 + 1`` does not make sense unless additive group represents
    specifically the set of integers.

    AdditiveGroup expressions have default heads NUMBER, SYMBOL, NEG,
    ADD, TERM_COEFF. When evaluation is disabled (see
    UnevaluatedAddition context) then MUL head will be used as well.
    """
    algebra_options = Group.algebra_options.copy()
    algebra_options.update(is_additive_group_commutative = False,
                           is_additive_group = True)

    def __pos__(self):
        return self.head.algebra_pos(type(self), self)

    def __neg__(self):
        return self.head.algebra_neg(type(self), self)

    def __add__(self, other, inplace=False):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.algebra_add_number(cls, self, other, inplace)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.algebra_add(cls, self, other, inplace)

    def __radd__(self, other):
        cls = type(self)
        other = cls.convert(other, typeerror=False)
        if other is NotImplemented: return NotImplemented
        return other.head.algebra_add(cls, other, self, False)

    def __iadd__(self, other):
        return self.__add__(other, self.is_writable)

    def __sub__(self, other, inplace=False):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.algebra_add_number(cls, self, -other, inplace)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.algebra_add(cls, self, -other, inplace)

    def __rsub__(self, other):
        return other + (-self)

    def __isub__(self, other):
        return self.__add__(-other, self.is_writable)

    def __mul__(self, other, inplace=False):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.algebra_mul_number(cls, self, other, inplace)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented:
                return NotImplemented
        if self.head is NUMBER or other.head is NUMBER:
            return self.head.algebra_mul(cls, self, other, inplace)
        raise TypeError('unsupported operand type(s) for *: %r and %r' \
                        % (self.pair, other.pair))

    def __rmul__(self, other):
        cls = type(self)
        other = cls.convert(other, typeerror=False)
        if other is NotImplemented:
            return NotImplemented
        if self.head is NUMBER or other.head is NUMBER:
            return other.head.algebra_mul(cls, other, self, False)
        raise TypeError('unsupported operand type(s) for *: %r and %r' \
                        % (type(other).__name__, cls.__name__))

    def __imul__(self, other):
        return self.__mul__(other, self.is_writable)

class AdditiveAbelianGroup(AdditiveGroup):
    """
    AdditiveAbelianGroup is an additive group with commutative
    addition.

    AdditiveAbelianGroup expressions have default heads NUMBER,
    SYMBOL, TERM_COEFF, TERM_COEFF_DICT.
    """

    algebra_options = AdditiveGroup.algebra_options.copy()
    algebra_options.update(is_additive_group_commutative = True)


class MultiplicativeGroup(Group):
    """
    MultiplicativeGroup is a Group with * operation, 1 is identity element
    and ``**(-1)`` denotes inverse. 
    """
    algebra_options = Group.algebra_options.copy()
    algebra_options.update(is_multiplicative_group_commutative = False,
                           is_multiplicative_group = True)

    def __mul__(self, other, inplace=False):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.algebra_mul_number(cls, self, other, inplace)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.algebra_mul(cls, self, other, inplace)

    def __rmul__(self, other):
        cls = type(self)
        other = cls.convert(other, typeerror=False)
        if other is NotImplemented: return NotImplemented
        return other.head.algebra_mul(cls, other, self, False)

    def __div__(self, other):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                pass
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self * other ** (-1)

    def __rdiv__(self, other):
        cls = type(self)
        other = cls.convert(other, typeerror=False)
        if other is NotImplemented: return NotImplemented
        return other * self ** (-1)

    def __pow__(self, other, inplace=False):
        cls = type(self)
        tother = type(other)
        if tother is not cls:
            if tother in numbertypes_set:
                return self.head.algebra_pow_number(cls, self, other, inplace)
            other = cls.convert(other, typeerror=False)
            if other is NotImplemented: return NotImplemented
        return self.head.algebra_pow(cls, self, other, inplace)

classes.Group = Group
classes.AdditiveGroup = AdditiveGroup
classes.AdditiveAbelianGroup = AdditiveAbelianGroup
classes.MultiplicativeGroup = MultiplicativeGroup

