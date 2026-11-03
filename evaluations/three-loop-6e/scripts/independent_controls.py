"""Independent source-domain controls for the six three-loop 6e routes.

No Wilson evaluator, raw CAS primitive, historical coefficient table or result
selector is imported. Polynomial controls check the source domain, not finite
Laurent coefficients. For holomorphic all-opposite roots, fresh original-domain
QMC can directly control the finite coefficient. Numerical errors are empirical.
"""
from __future__ import annotations

import argparse
import ast
import copy
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import time

import mpmath as mp
import sympy as sp

EPS = sp.Symbol('eps')
Z = sp.Symbol('z1', positive=True)
CA, CF = sp.symbols('CA CF')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def names_for(orbit):
    names = {n: sp.Symbol(n, positive=True) for n in orbit['J']['variables']}
    names.update({f'x{i}': sp.Symbol(f'x{i}', positive=True) for i in range(1, 7)})
    names.update(eps=EPS, z1=Z, z=Z, CA=CA, CF=CF, EulerGamma=sp.EulerGamma, pi=sp.pi)
    return names


def expression(text, names):
    """A closed arithmetic parser independent of the project's source reader."""
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return names[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            return visit(node.operand) if isinstance(node.op, ast.UAdd) else -visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return a+b
            if isinstance(node.op, ast.Sub): return a-b
            if isinstance(node.op, ast.Mult): return a*b
            if isinstance(node.op, ast.Div): return a/b
            if isinstance(node.op, ast.Pow): return a**b
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
            f = {'zeta': sp.zeta, 'gamma': sp.gamma, 'exp': sp.exp, 'log': sp.log}.get(node.func.id)
            if f and len(node.args) == 1: return f(visit(node.args[0]))
        raise ValueError('Expression outside the declared arithmetic grammar')
    return visit(ast.parse(text.replace('^', '**'), mode='eval').body)


def scalar(node):
    kind = node.get('kind')
    if kind == 'rational': return sp.Rational(node['numerator'], node['denominator'])
    if kind == 'var' and node['name'] in ('z', 'z1'): return Z
    if kind != 'op': raise ValueError('Unexpected coefficient AST')
    args = [scalar(x) for x in node['args']]
    op = node['op']
    operations = {'add': lambda: args[0]+args[1], 'sub': lambda: args[0]-args[1],
                  'mul': lambda: args[0]*args[1], 'div': lambda: args[0]/args[1],
                  'neg': lambda: -args[0], 'rpow': lambda: args[0]**args[1],
                  'exp': lambda: sp.exp(args[0]), 'log': lambda: sp.log(args[0]),
                  'gamma': lambda: sp.gamma(args[0]), 'zeta': lambda: sp.zeta(args[0]),
                  'pi': lambda: sp.pi, 'euler_gamma': lambda: sp.EulerGamma}
    if op not in operations: raise ValueError('Unexpected coefficient operation')
    return operations[op]()


def source_check(orbit):
    j = orbit['J']; names = names_for(orbit)
    order = j['integration_order_outer_to_inner']
    if len(order) != 6 or set(order) != set(j['variables']): raise ValueError('Wrong variable inventory')
    if j['measure_semantics'] != 'iterated_real_interval_integral' or j['delta_constraints']:
        raise ValueError('Wrong measure interpretation')
    available = {EPS, Z}
    for name in order:
        for bound in j['bounds'][name]:
            if not expression(bound, names).free_symbols <= available:
                raise ValueError('Dependent bound uses an unavailable outer variable')
        available.add(names[name])
    factors = [expression(f, names) for f in j['propagator_polynomials'].values()]
    raw = expression(j['numerator'], names)*sp.Mul(*(f**(EPS-1) for f in factors))
    if sp.simplify(raw-expression(j['expanded'], names)) != 0:
        raise ValueError('Expanded source integrand differs from propagator product')
    preparation = orbit['color_preparation']
    original = expression(preparation['original_C']['components']['color'], names)
    basis = [expression(b, names) for b in preparation['basis']]
    coordinates = [sp.Rational(x) for x in preparation['coordinates']]
    residual = sp.cancel((original-sum(c*b for c, b in zip(coordinates, basis))).subs(CF, (CA**2-1)/(2*CA)))
    if residual != 0 or coordinates[0] != sp.Rational(preparation['a1']):
        raise ValueError('Exact full-basis color reconstruction failed')
    euler = orbit['euler_preparation']
    factor = 'Exp[-3*EulerGamma*eps]'
    if euler['factor'] != factor or euler['gamma_power'] != 3 or euler['loop_order'] != 3:
        raise ValueError('Euler preparation differs')
    if orbit['C']['components']['euler_compensation'] != factor or orbit['C']['combined'].count(factor) != 1:
        raise ValueError('Euler compensation missing or duplicated')
    multiplicity = orbit['c4_symmetry']['orbit_multiplicity']
    if type(multiplicity) is not int or multiplicity < 1: raise ValueError('Wrong orbit multiplicity')
    kappa = sp.Rational(multiplicity)*coordinates[0]/(2*(-2)**3)
    return raw, kappa, names


def integrate_polynomial_original(poly, orbit, names):
    """Integrate the declared inner variable first, retaining dependent bounds."""
    value = sp.expand(poly)
    for name in reversed(orbit['J']['integration_order_outer_to_inner']):
        variable = names[name]
        low, high = [expression(x, names) for x in orbit['J']['bounds'][name]]
        primitive = sum(c*variable**(powers[0]+1)/(powers[0]+1)
                        for powers, c in sp.Poly(value, variable).terms())
        value = sp.expand(primitive.subs(variable, high)-primitive.subs(variable, low))
    if value.free_symbols - {Z}: raise ValueError('Unintegrated source coordinate')
    return value


def cube_map(orbit, names):
    """Derive the triangular map from bounds; do not trust input_export.s_to_x."""
    mapping = {}; jacobian = sp.Integer(1)
    order = orbit['J']['integration_order_outer_to_inner']
    for index, name in enumerate(order, 1):
        low, high = [expression(x, names).subs(mapping) for x in orbit['J']['bounds'][name]]
        mapping[names[name]] = low+(high-low)*names[f'x{index}']
        jacobian *= high-low
    variables = [names[f'x{i}'] for i in range(1, 7)]
    determinant = sp.factor(sp.det(sp.Matrix([mapping[names[n]] for n in order]).jacobian(variables)))
    if sp.expand(determinant-jacobian) != 0: raise ValueError('Triangular Jacobian identity failed')
    return mapping, sp.expand(jacobian), variables


def cube_polynomial_moments(poly, mapping, jacobian, variables):
    transformed = sp.Poly(sp.expand(poly.subs(mapping, simultaneous=True)*jacobian), *variables)
    value = sum(coefficient/sp.prod(p+1 for p in powers) for powers, coefficient in transformed.terms())
    return sp.expand(value), len(transformed.terms())


def polynomial_controls(orbit):
    raw, kappa, names = source_check(orbit)
    mapping, jacobian, cube_variables = cube_map(orbit, names)
    rows = []
    for epsilon in (1, 2):
        polynomial = sp.expand(raw.subs(EPS, epsilon))
        original = integrate_polynomial_original(polynomial, orbit, names)
        independent, terms = cube_polynomial_moments(polynomial, mapping, jacobian, cube_variables)
        residual = sp.expand(original-independent)
        if residual != 0: raise ValueError('Original iterated integral and cube monomial integral differ')
        rows.append({'epsilon': epsilon, 'z': 'symbolic positive z1',
                     'exact_J': str(original), 'residual': str(residual),
                     'cube_monomials': terms, 'passed': True})
    return {'integral_id': orbit['orbit_id'], 'kappa': str(kappa),
            'method': 'Original iterated polynomial antiderivatives versus independently derived triangular cube map and product monomial moments.',
            'jacobian': str(jacobian), 'polynomial_controls': rows,
            'scope': 'At eps=1,2 only. Verifies source/domain/Jacobian manipulation; does not determine Laurent finite coefficients.'}


class IndependentHLog:
    """Log-polynomial Taylor series for 0/-1 GPLs, independent of Wilson."""
    def __init__(self, z, terms=128):
        self.z = mp.mpf(z); self.terms = terms
        if not 0 < self.z < 1: raise ValueError('Taylor control needs 0<z<1')

    @staticmethod
    def integrate_term(n, polynomial):
        if n == 0:
            return {k+1: value/(k+1) for k, value in polynomial.items()}
        out = {}
        for k, value in polynomial.items():
            for j in range(k+1):
                degree = k-j
                coefficient = value*((-1)**j)*math.factorial(k)/math.factorial(degree)/(mp.mpf(n)**(j+1))
                out[degree] = out.get(degree, mp.mpf(0))+coefficient
        return out

    @lru_cache(maxsize=None)
    def series(self, word):
        if not word: return {0: {0: mp.mpf(1)}}
        first, rest = word[0], word[1:]
        tail = self.series(rest)
        if first == 0:
            return {n: self.integrate_term(n, polynomial) for n, polynomial in tail.items()}
        if first != -1: raise ValueError('Control supports only 0/-1 letters')
        out = {}; previous = {}
        # Coefficients of tail/(1+z) satisfy h[n]=tail[n]-h[n-1].
        for power in range(self.terms):
            current = {k: -v for k, v in previous.items()}
            for k, value in tail.get(power, {}).items(): current[k] = current.get(k, mp.mpf(0))+value
            if current: out[power+1] = self.integrate_term(power+1, current)
            previous = current
        return out

    def value(self, word):
        logz = mp.log(self.z)
        return mp.fsum(self.z**n*mp.fsum(value*logz**k for k, value in polynomial.items())
                       for n, polynomial in self.series(tuple(word)).items())


def coefficient_value(payload, z):
    evaluator = IndependentHLog(str(z))
    total = mp.mpf(0)
    for term in payload['terms']:
        if term['periods']: raise ValueError('Unexpected non-scalar periods')
        gpl = term['gpl']
        if gpl['branch'] != 'positive-real' or gpl['endpoint_kind'] != 'ordinary':
            raise ValueError('Unexpected GPL branch or endpoint convention')
        if scalar(gpl['basepoint']) != 0 or scalar(gpl['endpoint']) not in (Z, 0):
            raise ValueError('Unexpected GPL path')
        letters = [scalar(a) for a in gpl['word']]
        if any(letter not in (0, -1) for letter in letters):
            raise ValueError('Control supports only exact 0/-1 letters')
        if letters and scalar(gpl['endpoint']) != Z:
            raise ValueError('Nonempty GPL must have the specified positive endpoint')
        word = [int(letter) for letter in letters]
        c = scalar(term['coefficient']).subs(Z, sp.Rational(str(z)))
        coefficient = mp.mpf(str(sp.N(c, 70)))
        total += coefficient*evaluator.value(word)
    return total


def radical_inverse(index, base):
    value = 0.0; step = 1.0/base
    while index:
        index, digit = divmod(index, base); value += digit*step; step /= base
    return value


def original_finite_qmc(orbit, coefficient, samples, shifts=8):
    """Only for holomorphic all-opposite roots: eps=0 may be taken directly."""
    raw, kappa, names = source_check(orbit)
    if any(len(expression(f, names).free_symbols-{Z}) != 2 for f in orbit['J']['propagator_polynomials'].values()):
        raise ValueError('Unexpected propagator endpoint inventory')
    # Every denominator must be of the opposite-edge polynomial type, rather
    # than an adjacent monomial whose zero-epsilon integral could diverge.
    for text in orbit['J']['propagator_polynomials'].values():
        f = expression(text, names)
        endpoints = sorted(f.free_symbols-{Z}, key=str)
        a, b = endpoints
        opposite_24 = (1-a)*(1-b)+Z*a*b
        opposite_13 = a*b+Z*(1-a)*(1-b)
        if all(sp.expand(f-model) != 0 for model in (opposite_24, opposite_13)):
            raise ValueError('Direct finite-QMC requires all-opposite holomorphic source')
    mapping, jacobian, coordinates = cube_map(orbit, names)
    integrand = raw.subs(EPS, 0).subs(mapping, simultaneous=True)*jacobian
    f = sp.lambdify([Z, *coordinates], integrand, modules='math', cse=True)
    rng = random.Random(20261103+int(orbit['orbit_id'].split('-')[1]))
    points = [tuple(radical_inverse(i, b) for b in (2, 3, 5, 7, 11, 13)) for i in range(1, samples+1)]
    estimates = []; z = 0.5
    for _ in range(shifts):
        shift = [rng.random() for _ in range(6)]
        values = []
        for point in points:
            x = tuple((a+b) % 1 for a, b in zip(point, shift))
            values.append(float(kappa)*(f(z, *x)+f(1/z, *x))/2)
        estimates.append(math.fsum(values)/samples)
    expected = float(coefficient_value(coefficient, z))
    mean = statistics.mean(estimates)
    stderr = statistics.stdev(estimates)/math.sqrt(shifts)
    residual = abs(mean-expected)
    tolerance = 6*stderr+0.005*abs(expected)
    return {'integral_id': orbit['orbit_id'], 'epsilon': 0, 'z': str(z),
            'samples_per_shift': samples, 'shifts': shifts, 'shift_estimates': estimates,
            'original_domain_estimate': mean, 'candidate_c0': expected,
            'absolute_residual': residual, 'empirical_standard_error': stderr,
            'acceptance_tolerance': tolerance, 'passed': residual <= tolerance,
            'scope': 'Fresh original six-dimensional ordered source, both reciprocal z channels, physical kappa/2; empirical randomized Halton control, no strict enclosure.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=16384)
    args = parser.parse_args()
    if args.out.exists(): raise ValueError('Choose a NEW output file')
    if not 1024 <= args.samples <= 131072: raise ValueError('Samples outside declared budget')
    started = time.monotonic(); mp.mp.dps = 75
    data = json.loads(args.input.read_text(encoding='utf-8'))
    candidates = json.loads(args.results.read_text(encoding='utf-8'))
    if 'results' in candidates and isinstance(candidates['results'], dict): candidates = candidates['results']
    orbits = {o['orbit_id']: o for o in data['orbits']}
    representatives = {'empty': '6e-01', 'single-opposite': '6e-04',
                       'same-channel-two': '6e-10', 'same-channel-three': '6e-20',
                       'double-opposite': '6e-32', 'mixed-opposite': '6e-37'}
    # Independently reconstruct every source's full color vector and kappa.
    source_rows = []
    for oid, orbit in sorted(orbits.items()):
        _, kappa, _ = source_check(orbit)
        candidate = candidates[oid]
        normalization = candidate.get('provenance', {}).get('normalization', candidate.get('normalization'))
        if normalization is not None and sp.Rational(normalization['kappa']) != kappa:
            raise ValueError('Candidate normalization differs from independently reconstructed source '+oid)
        source_rows.append({'integral_id': oid, 'kappa': str(kappa), 'passed': True})
    polynomial = []
    for oid in [*representatives.values(), '6e-38']:
        polynomial.append(polynomial_controls(orbits[oid]))
    evaluator = IndependentHLog('0.5')
    sanity = [{'word': [-1], 'residual': mp.nstr(abs(evaluator.value((-1,))-mp.log(mp.mpf('1.5'))), 8)},
              {'word': [0], 'residual': mp.nstr(abs(evaluator.value((0,))-mp.log(mp.mpf('.5'))), 8)},
              {'word': [0, -1], 'residual': mp.nstr(abs(evaluator.value((0, -1))+mp.polylog(2, mp.mpf('-.5'))), 8)}]
    if any(mp.mpf(row['residual']) > mp.mpf('1e-37') for row in sanity): raise ValueError('Independent HLOG evaluator failed known identities')
    finite = [original_finite_qmc(orbits[oid], candidates[oid]['coefficients']['0'], args.samples)
              for oid in ('6e-20', '6e-21', '6e-22', '6e-37', '6e-38')]
    negative = []
    def rejects(label, operation):
        try: operation()
        except (ValueError, KeyError, AssertionError): negative.append({'case': label, 'rejected': True})
        else: raise AssertionError('Negative control was accepted: '+label)
    bad = copy.deepcopy(orbits['6e-04']); bad['J']['integration_order_outer_to_inner'].reverse()
    rejects('unavailable dependent bound after illegal order change', lambda: source_check(bad))
    bad = copy.deepcopy(orbits['6e-01']); bad['C']['components']['euler_compensation'] = 'Exp[-6*EulerGamma*eps]'
    rejects('duplicated Euler normalization', lambda: source_check(bad))
    bad = copy.deepcopy(orbits['6e-10']); bad['color_preparation']['a1'] = '0'
    rejects('wrong projected color', lambda: source_check(bad))
    rejects('finite-domain integration incorrectly applied to pole source', lambda: original_finite_qmc(orbits['6e-04'], candidates['6e-04']['coefficients']['0'], 1024))
    bad = copy.deepcopy(orbits['6e-01']); bad['J']['bounds']['s31'] = ['0', '1']
    raw, _, names = source_check(bad)
    wrong = integrate_polynomial_original(raw.subs(EPS, 2), bad, names)
    true_raw, _, true_names = source_check(orbits['6e-01'])
    right = integrate_polynomial_original(true_raw.subs(EPS, 2), orbits['6e-01'], true_names)
    def require_same():
        if sp.expand(wrong-right) != 0: raise ValueError('Dropped dependent domain changes the integral')
    rejects('dropped ordered upper bound', require_same)
    stripped = []
    for oid in ('6e-37', '6e-38'):
        raw, kappa, names = source_check(orbits[oid])
        j2 = integrate_polynomial_original(raw.subs(EPS, 2), orbits[oid], names)
        # Cancel the common G(eps)^3 before testing an all-epsilon identity.
        # This rational factor is regular at eps=2 even though Gamma is not.
        stripped.append(sp.cancel(kappa*(j2+Z**6*j2.subs(Z, 1/Z))/(1+Z**6)))
    alias_residual = sp.cancel((stripped[0]-stripped[1]).subs(Z, sp.Rational(1, 2)))
    def wrong_full_alias():
        if alias_residual != 0: raise ValueError('Finite-only alias is not an all-epsilon identity')
    rejects('promoting 37/38 finite equality to an all-epsilon equality', wrong_full_alias)
    coverage = [{'family': family, 'polynomial_representative': oid,
                 'source_normalization_scope': 'all members individually checked',
                 'finite_original_domain_control': 'all three members' if family == 'same-channel-three' else ('both members' if family == 'mixed-opposite' else 'not run; polynomial controls do not validate finite terms')}
                for family, oid in representatives.items()]
    result = {'format': 'wilson-six-route-independent-controls/1',
              'inputs': {'ordered_6e': {'path': str(args.input), 'sha256': sha(args.input)},
                         'frozen_candidates': {'path': str(args.results), 'sha256': sha(args.results)}},
              'source_checks': source_rows, 'representative_polynomial_controls': polynomial,
              'independent_hlog_sanity': sanity, 'finite_original_domain_controls': finite,
              'negative_controls': negative, 'coverage': coverage,
              'finite_alias_scope_control': {'pair': ['6e-37', '6e-38'], 'epsilon': 2,
                                           'z': '1/2', 'I_over_G_cubed_difference': str(alias_residual),
                                           'scope': 'Nonzero exact rational difference excludes an all-epsilon alias; no Gamma pole evaluation is performed.'},
              'passed': all(r['passed'] for r in finite), 'elapsed_seconds': time.monotonic()-started,
              'formal_state': 'not_run',
              'limitations': ['No new analytic CAS integration is claimed.',
                              'Positive integer epsilon source identities do not determine Laurent finite coefficients.',
                              'Original-domain finite controls apply only to holomorphic all-opposite roots.',
                              'QMC standard errors and HLOG truncation errors are empirical, not strict intervals.',
                              'A representative structural check is not an independent finite check for every member.']}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'passed': result['passed'], 'source_checks': len(source_rows),
                      'polynomial_integrations': 2*len(polynomial), 'finite_checks': len(finite),
                      'negative_controls': len(negative), 'elapsed_seconds': result['elapsed_seconds']}))
    return 0 if result['passed'] else 1


if __name__ == '__main__': raise SystemExit(main())
