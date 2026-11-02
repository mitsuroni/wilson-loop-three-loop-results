"""Exact ambient Gamma1 differential kernel, with no runtime CAS dependency.

Derivatives act on the independent three-point star polynomial *before*
restriction to a contour or face.  The D-dimensional metric trace is retained
as D=4-2*eps.  Pure-Python Fraction polynomials are shared with j_factor.
"""
from copy import deepcopy
from functools import lru_cache
import json
import re

from montopar.gamma1_data import OPERATORS
from montopar.geometry import dot_product_factor
from montopar.j_factor import (
    _parse_polynomial, _format_poly, _poly_add, _poly_mul, _poly_pow,
    _poly_scale,
)

ALPHAS = ('alpha1', 'alpha2', '(1-alpha1-alpha2)')
_R = {f'r{a+1}{i+1}': (a, i) for a in range(3) for i in range(3)}
_P = {'p1sq': (0, 0), 'p2sq': (1, 1), 'p1p2': (0, 1)}


def operators():
    """Return the universal N0, L0 and three M polynomials for both colours."""
    return deepcopy(OPERATORS)


@lru_cache(maxsize=4096)
def _parse_cached(expr):
    return _parse_polynomial(expr)


def parse_polynomial(expr):
    """Return an independent sparse Fraction polynomial; reject trailing input."""
    return dict(_parse_cached(str(expr)))


def _substitute(poly, replacements):
    result = {}
    powers = {}
    for monomial, coefficient in poly.items():
        term = {frozenset(): coefficient}
        for variable, power in sorted(monomial):
            key = (variable, power)
            if key not in powers:
                base = replacements.get(variable, {frozenset({(variable, 1)}): 1})
                powers[key] = _poly_pow(base, power)
            term = _poly_mul(term, powers[key])
        result = _poly_add(result, term)
    return result


def polynomial_substitute(expr, mapping):
    """Simultaneous whole-symbol polynomial substitution, then exact expansion."""
    replacements = {k: parse_polynomial(v) for k, v in mapping.items()}
    return _format_poly(_substitute(parse_polynomial(expr), replacements))


def _add(*terms):
    result = {}
    for t in terms:
        result = _poly_add(result, parse_polynomial(t))
    return _format_poly(result)


def _mul(*terms):
    result = {frozenset(): 1}
    for t in terms:
        result = _poly_mul(result, parse_polynomial(t))
    return _format_poly(result)


def _sub(a, b):
    return _add(a, _mul('-1', b))


def geometry(data):
    """Return the dimensionless star geometry at three ordered source endpoints.

    Source validation belongs to the Gamma1 parser.  The endpoint order here
    is the central-vertex leg order, never a reordering of the integration.
    """
    ends = data['endpoints']
    if len(ends) != 3 or any(p['edge'] not in (1, 2, 3, 4) for p in ends):
        raise ValueError('gamma1_geometry_requires_three_four_edge_endpoints')
    gram = [[_mul('-1/2', dot_product_factor(i, j)) for j in range(1, 5)]
            for i in range(1, 5)]
    def dot(x, y):
        return _add(*(_mul(x[i], gram[i][j], y[j]) for i in range(4) for j in range(4)))
    x = [['1' if k < p['edge'] else p['s_var'] if k == p['edge'] else '0'
          for k in range(1, 5)] for p in ends]
    n = [['1' if k == p['edge'] else '0' for k in range(1, 5)] for p in ends]
    center = [_add(*(_mul(ALPHAS[i], x[i][k]) for i in range(3))) for k in range(4)]
    r = [[_sub(x[i][k], center[k]) for k in range(4)] for i in range(3)]
    distances = [[_mul('-1', dot([_sub(x[i][k], x[j][k]) for k in range(4)],
                                  [_sub(x[i][k], x[j][k]) for k in range(4)]))
                  for j in range(3)] for i in range(3)]
    return {'q0': _add(*(_mul(ALPHAS[i], ALPHAS[j], distances[i][j])
                         for i in range(3) for j in range(i+1, 3))),
            'G': [[dot(n[i], n[j]) for j in range(3)] for i in range(3)],
            'T': [[dot(n[i], r[j]) for j in range(3)] for i in range(3)],
            'W': [[dot(r[i], r[j]) for j in range(3)] for i in range(3)]}


def _derivative_shapes(rs, ps):
    """Contractions in a symbolic Gram alphabet, before geometric substitution."""
    def a(i): return f'alpha{i+1}'
    def t(v, i): return f'T{v+1}{i+1}'
    def g(v, w): return f'G{v+1}{w+1}'
    def h(i, j): return f'({a(i)}*{int(i == j)}-{a(i)}*{a(j)})'
    def grad(r): return f'(-2*{a(r[1])}*{t(*r)})'
    def hess(r, s): return f'(-2*{h(r[1], s[1])}*{g(r[0], s[0])})'
    if not ps:
        if not rs:
            return {0: '1'}
        if len(rs) == 1:
            return {1: f'beta*{grad(rs[0])}'}
        if len(rs) == 2:
            r, s = rs
            return {2: f'beta*(beta-1)*{grad(r)}*{grad(s)}',
                    1: f'beta*{hess(r, s)}'}
        if len(rs) == 3:
            r, s, t_ = rs
            return {3: f'beta*(beta-1)*(beta-2)*{grad(r)}*{grad(s)}*{grad(t_)}',
                    2: f'beta*(beta-1)*({hess(r,s)}*{grad(t_)}+'
                       f'{hess(r,t_)}*{grad(s)}+{hess(s,t_)}*{grad(r)})'}
    elif len(ps) == 1 and len(rs) <= 1:
        i, j = ps[0]
        dotgrad = f'(4*{a(i)}*{a(j)}*W{i+1}{j+1})'
        trace = f'(-2*D*{h(i,j)})'
        if not rs:
            return {2: f'beta*(beta-1)*{dotgrad}', 1: f'beta*{trace}'}
        v, k = rs[0]
        return {3: f'beta*(beta-1)*(beta-2)*{dotgrad}*{grad(rs[0])}',
                2: f'beta*(beta-1)*({trace}*{grad(rs[0])}+'
                   f'4*{h(i,k)}*{a(j)}*{t(v,j)}+4*{h(j,k)}*{a(i)}*{t(v,i)})'}
    raise ValueError('gamma1_operator_momentum_degree_must_be_at_most_three')


@lru_cache(maxsize=512)
def _apply_cached(op, beta, geometry_key):
    geom = json.loads(geometry_key)
    # Group momentum monomials first.  This avoids repeated ambient contractions
    # for the many independent t1,t2,D coefficients in the same derivative.
    groups = {}
    allowed = set(_R) | set(_P) | {'D', 't1', 't2', 't3', 'eps', 'z1', 'z',
                                 'h12', 'h23', 'h31'}
    for mono, coefficient in parse_polynomial(op).items():
        unknown = {v for v, _ in mono} - allowed
        if unknown:
            raise ValueError('gamma1_unknown_operator_symbols: ' + ','.join(sorted(unknown)))
        rs, ps, scalar = [], [], []
        for variable, power in sorted(mono):
            if variable in _R:
                rs.extend([_R[variable]] * power)
            elif variable in _P:
                ps.extend([_P[variable]] * power)
            else:
                scalar.append((variable, power))
        key = (tuple(rs), tuple(ps))
        groups[key] = _poly_add(groups.get(key, {}), {frozenset(scalar): coefficient})
    contracted = {}
    for (rs, ps), coefficients in groups.items():
        for shift, shape in _derivative_shapes(rs, ps).items():
            contracted[shift] = _poly_add(contracted.get(shift, {}),
                                         _poly_mul(coefficients, parse_polynomial(shape)))
    replacements = {'D': '4-2*eps', 'beta': beta, 'alpha3': '1-alpha1-alpha2',
                    't3': '1-t1-t2', 'z': 'z1',
                    'h12': geom['G'][0][1], 'h23': geom['G'][1][2], 'h31': geom['G'][2][0]}
    for name in ('G', 'T', 'W'):
        for i in range(3):
            for j in range(3):
                replacements[f'{name}{i+1}{j+1}'] = geom[name][i][j]
    parsed = {k: parse_polynomial(v) for k, v in replacements.items()}
    return tuple((shift, _format_poly(result)) for shift, polynomial in sorted(contracted.items())
                 if (result := _substitute(polynomial, parsed)))


def apply_operator(op, beta, geom):
    """Return ``{shift: coefficient}`` in O q^beta = sum coeff*q^(beta-shift).

    Supports momentum degrees 0--3, including contracted second derivatives.
    Output coefficients are fully evaluated exact polynomials using ``^``.
    A zero operator returns an empty dictionary.  No derivative remains in an
    emitted coefficient; the universal operators remain separately available.
    """
    relevant = {key: geom[key] for key in ('G', 'T', 'W')}
    return dict(_apply_cached(str(op), str(beta), json.dumps(relevant, sort_keys=True)))


@lru_cache(maxsize=512)
def _apply_compact_cached(op, beta, geometry_key):
    """Evaluate every derivative, retaining products of geometry polynomials."""
    geom = json.loads(geometry_key)
    groups = {}
    allowed = set(_R) | set(_P) | {'D', 't1', 't2', 't3', 'eps', 'z1', 'z',
                                 'h12', 'h23', 'h31'}
    for mono, coefficient in parse_polynomial(op).items():
        unknown = {v for v, _ in mono} - allowed
        if unknown:
            raise ValueError('gamma1_unknown_operator_symbols: ' + ','.join(sorted(unknown)))
        rs, ps, scalar = [], [], []
        for variable, power in sorted(mono):
            if variable in _R:
                rs.extend([_R[variable]] * power)
            elif variable in _P:
                ps.extend([_P[variable]] * power)
            else:
                scalar.append((variable, power))
        key = (tuple(rs), tuple(ps))
        groups[key] = _poly_add(groups.get(key, {}), {frozenset(scalar): coefficient})
    replacements = {'D': '4-2*eps', 'beta': beta, 'alpha3': '1-alpha1-alpha2',
                    't3': '1-t1-t2', 'z': 'z1',
                    'h12': geom['G'][0][1], 'h23': geom['G'][1][2], 'h31': geom['G'][2][0]}
    for name in ('G', 'T', 'W'):
        for i in range(3):
            for j in range(3):
                replacements[f'{name}{i+1}{j+1}'] = geom[name][i][j]
    parsed = {key: parse_polynomial(value) for key, value in replacements.items()}
    identifier = re.compile(r'\b[A-Za-z_][A-Za-z0-9_]*\b')
    def inline(expression):
        # One simultaneous whole-token pass. All replacement values are already
        # independent expressions; alpha3 and D do not occur in source geometry.
        return identifier.sub(lambda match: '('+replacements[match[0]]+')'
                              if match[0] in replacements else match[0], expression)
    result = {}
    for (rs, ps), coefficient in groups.items():
        # Apply scalar h/D substitutions exactly before emitting any products;
        # this also removes zero tangent contractions in repeated-edge cases.
        scalar = _substitute(coefficient, parsed)
        shapes = _derivative_shapes(rs, ps)  # validates even a vanishing term
        if not scalar:
            continue
        coefficient_text = _format_poly(scalar)
        for shift, shape in shapes.items():
            result.setdefault(shift, []).append('('+coefficient_text+')*('+inline(shape)+')')
    return tuple((shift, '+'.join(terms)) for shift, terms in sorted(result.items()))


def apply_operator_compact(op, beta, geom):
    """Same evaluated kernel as :func:`apply_operator`, with factored coefficients.

    Every G/T/W, h, D, beta, alpha3 and t3 abbreviation is replaced inline.
    This is a self-contained ordinary polynomial, not a deferred differential
    operator. Products are retained to avoid million-character raw integrands.
    The expanded API remains useful for exact reconstruction and zero tests.
    """
    relevant = {key: geom[key] for key in ('G', 'T', 'W')}
    return dict(_apply_compact_cached(str(op), str(beta), json.dumps(relevant, sort_keys=True)))
