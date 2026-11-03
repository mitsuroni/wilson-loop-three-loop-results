"""Finite exact Laurent jets without SymPy's general series machinery.

The output contains every nonzero coefficient through the *final* Laurent
power requested.  Poles determine the extra regular coefficients needed by
each convolution.  This is a deliberately closed rule set, not a fallback to
``series`` or a claim that an unsupported expression is regular.
"""
from __future__ import annotations

from functools import lru_cache

import sympy as sp


class JetError(ValueError):
    """The finite Laurent rule set cannot certify this expression."""


class JetLimit(JetError):
    """A finite expansion exceeded the local algebra budget."""


_MAX_ORDER = 128
_MAX_SHIFT = 256
_MAX_POWER = 256
_LEADING_SEARCH = 64


def _clean(coefficients):
    result = {}
    for power, value in coefficients.items():
        value = sp.expand(value)
        if value != 0:
            result[power] = value
    return result


def _multiply(left, right, through):
    result = {}
    for p, a in left.items():
        for q, b in right.items():
            if p + q <= through:
                result[p + q] = result.get(p + q, sp.S.Zero) + a*b
    return _clean(result)


def _exponential(argument, through):
    """exp(A) recurrence for an already certified regular argument jet."""
    if through < 0:
        return {}
    result = {0: sp.exp(argument.get(0, sp.S.Zero))}
    for n in range(1, through + 1):
        coefficient = sp.expand(sum(k*argument.get(k, sp.S.Zero)*
                                    result.get(n-k, sp.S.Zero)
                                    for k in range(1, n+1))/n)
        if coefficient != 0:
            result[n] = coefficient
    return result


def finite_laurent(expr, eps, through):
    """Return ``{epsilon_power: exact_coefficient}`` through ``through``.

    Inputs are exact SymPy expressions, never strings.  Supported epsilon
    dependence consists of sums, products, integer powers, regular exp,
    positive constant bases to regular powers, and Gamma(k+a*eps) with integer
    k.  Unknown inverse leading coefficients, essential singularities and
    unsupported functions raise :class:`JetError`.  Absent output powers up
    to ``through`` are zero; no information above it is returned.
    """
    if not isinstance(eps, sp.Symbol) or type(through) is not int:
        raise JetError('a symbolic epsilon and integer final power are required')
    if isinstance(expr, int) and not isinstance(expr, bool):
        expr = sp.Integer(expr)
    if not isinstance(expr, sp.Expr):
        raise JetError('an exact SymPy expression is required')
    if expr.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise JetError('inexact or nonfinite expression')
    if abs(through) > _MAX_ORDER:
        raise JetLimit('requested Laurent power exceeds jet budget')
    return dict(_Jets(eps).evaluate(expr, through))


class _Jets:
    def __init__(self, eps):
        self.eps = eps

    @lru_cache(maxsize=4096)
    def gamma_data(self, expression):
        argument = sp.expand(expression.args[0])
        k = argument.subs(self.eps, 0)
        a = sp.diff(argument, self.eps)
        if (k.is_Integer is not True or a.has(self.eps)
                or sp.expand(argument-k-a*self.eps) != 0):
            raise JetError('Gamma requires an affine integer-shift argument')
        k = int(k)
        if abs(k) > _MAX_SHIFT:
            raise JetLimit('Gamma integer shift exceeds jet budget')
        if k <= 0 and a.is_zero is not False:
            raise JetError('Gamma pole slope is zero or not known nonzero')
        return k, a

    @lru_cache(maxsize=8192)
    def lower(self, expression):
        """A certified lower bound; sums need not attain this bound."""
        if not expression.has(self.eps):
            return 0
        if expression == self.eps:
            return 1
        if expression.is_Add:
            return min(self.lower(a) for a in expression.args)
        if expression.is_Mul:
            return sum(self.lower(a) for a in expression.args)
        if expression.is_Pow:
            base, power = expression.args
            if power.is_Integer:
                n = int(power)
                if abs(n) > _MAX_POWER:
                    raise JetLimit('integer power exceeds jet budget')
                return n*(self.valuation(base) if n < 0 else self.lower(base))
            if not base.has(self.eps) and base.is_positive is True:
                self.regular(power*sp.log(base))
                return 0
            raise JetError('noninteger power requires a positive epsilon-independent base')
        if expression.func == sp.exp:
            self.regular(expression.args[0])
            return 0
        if expression.func == sp.gamma:
            k, _ = self.gamma_data(expression)
            return -1 if k <= 0 else 0
        raise JetError('unsupported epsilon dependence: '+str(expression.func))

    def regular(self, expression):
        # This check belongs before a lower-bound early return: exp(1/eps)
        # is unsupported even when only negative output powers are requested.
        if self.lower(expression) < 0 and self.evaluate(expression, -1):
            raise JetError('exponential has an essential Laurent singularity')

    @lru_cache(maxsize=4096)
    def valuation(self, expression):
        """Find a proved nonzero leading term; never infer generic nonzero."""
        if expression == 0:
            raise JetError('division by identically zero expression')
        lower = self.lower(expression)
        # Increase the requested jet only when cancellation has actually
        # removed all coefficients inspected so far.
        for through in range(lower, min(lower+_LEADING_SEARCH, _MAX_ORDER)+1):
            coefficients = self.evaluate(expression, through)
            for p in sorted(coefficients):
                coefficient = coefficients[p]
                if coefficient.is_zero is False:
                    return p
                if coefficient.is_zero is None:
                    raise JetError('inverse leading coefficient is not known nonzero')
        raise JetLimit('could not establish a leading coefficient within jet budget')

    def inverse(self, expression, through):
        valuation = self.valuation(expression)
        if through < -valuation:
            return {}
        order = through + valuation
        original = self.evaluate(expression, through+2*valuation)
        leading = original[valuation]
        regular = {p-valuation: c for p, c in original.items()}
        result = {0: 1/leading}
        for n in range(1, order+1):
            coefficient = sp.expand(-sum(regular.get(k, sp.S.Zero)*
                                         result.get(n-k, sp.S.Zero)
                                         for k in range(1, n+1))/leading)
            if coefficient != 0:
                result[n] = coefficient
        return {p-valuation: c for p, c in result.items()}

    def power(self, base, n, through):
        if n == 0:
            return {0: sp.S.One} if through >= 0 else {}
        if n < 0:
            lower = -self.valuation(base)
            count = -n
            factor = self.inverse(base, through-(count-1)*lower)
        else:
            lower = self.lower(base)
            count = n
            factor = self.evaluate(base, through-(count-1)*lower)
        result = {0: sp.S.One}
        for i in range(count):
            result = _multiply(result, factor, through-(count-i-1)*lower)
        return result

    def gamma(self, expression, through):
        k, a = self.gamma_data(expression)
        lower = -1 if k <= 0 else 0
        order = through-lower
        if order > _MAX_ORDER:
            raise JetLimit('required log-Gamma jet order exceeds budget')
        logarithm = {1: -sp.EulerGamma*a}
        logarithm.update({j: (-1)**j*sp.zeta(j)*a**j/sp.Integer(j)
                          for j in range(2, order+1)})
        regular = _exponential(logarithm, order)
        if k == 1:
            return regular
        if k > 1:
            shift = sp.prod(a*self.eps+i for i in range(1, k))
        else:
            shift = 1/sp.prod(a*self.eps+i for i in range(k, 1))
        return _multiply(regular, self.evaluate(shift, through), through)

    @lru_cache(maxsize=16384)
    def evaluate(self, expression, through):
        if abs(through) > _MAX_ORDER:
            raise JetLimit('required internal jet order exceeds budget')
        lower = self.lower(expression)
        if through < lower:
            return {}
        if not expression.has(self.eps):
            return {0: expression} if expression != 0 else {}
        if expression == self.eps:
            return {1: sp.S.One}
        if expression.is_Add:
            result = {}
            for term in expression.args:
                for p, coefficient in self.evaluate(term, through).items():
                    result[p] = result.get(p, sp.S.Zero)+coefficient
            return _clean(result)
        if expression.is_Mul:
            factors = expression.args
            bounds = [self.lower(a) for a in factors]
            remaining = sum(bounds)
            result = {0: sp.S.One}
            for factor, bound in zip(factors, bounds):
                coefficients = self.evaluate(factor, through-(lower-bound))
                remaining -= bound
                result = _multiply(result, coefficients, through-remaining)
            return result
        if expression.is_Pow:
            base, power = expression.args
            if power.is_Integer:
                return self.power(base, int(power), through)
            return self.evaluate(sp.exp(power*sp.log(base)), through)
        if expression.func == sp.exp:
            argument = self.evaluate(expression.args[0], max(0, through))
            if any(p < 0 for p in argument):
                raise JetError('exponential has an essential Laurent singularity')
            return _exponential(argument, through)
        if expression.func == sp.gamma:
            return self.gamma(expression, through)
        raise JetError('unsupported epsilon dependence: '+str(expression.func))
