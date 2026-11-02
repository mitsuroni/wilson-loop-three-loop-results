"""C4 covariance of actual Gamma1 sources; no integral-order interchange claim.

Geometry identities are exact in all variables. Operator coefficient checks use
two exact rational points and supplement the separate unrestricted cyclic-operator
identity; they do not alone prove an identity for all z or integrated equality.
"""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path

import pytest
import sympy as sp

from montopar.gamma1_source import validate_source, c4_symmetry
from montopar.gamma1_kernel import geometry, apply_operator, operators, parse_polynomial
from montopar.gamma1_blocks import MAX_BOUND, block_specs, sectorize

SOURCE = json.loads((Path(__file__).parent / 'fixtures' / 'gamma1_orbits.json').read_text(encoding='utf-8'))
POINTS = [dict(z1=Fraction(7, 3), eps=Fraction(1, 7), alpha1=Fraction(1, 4), alpha2=Fraction(1, 3),
               t1=Fraction(2, 7), t2=Fraction(1, 5)),
          dict(z1=Fraction(11, 5), eps=Fraction(2, 9), alpha1=Fraction(2, 9), alpha2=Fraction(3, 10),
               t1=Fraction(1, 6), t2=Fraction(2, 5))]


def expression(text):
    return sp.sympify(str(text).replace('^', '**'))


def evaluate(text, values):
    total = Fraction(0)
    for mono, coefficient in parse_polynomial(str(text)).items():
        term = Fraction(coefficient)
        for variable, power in mono:
            term *= values[variable] ** power
        total += term
    return total


def rotated_source(original, rotation):
    """Independent source factory; return new-leg -> original-leg permutation."""
    old_data = validate_source(original)
    permutation = sorted(range(3), key=lambda i: ((old_data['endpoints'][i]['edge']-1+rotation) % 4+1,
                                                 old_data['endpoints'][i]['slot']))
    edges = [(old_data['endpoints'][i]['edge']-1+rotation) % 4+1 for i in permutation]
    raw = deepcopy(original); raw['orbit_id'] = original['orbit_id'] + f'-C4-{rotation}'
    counts = Counter(edges); slots = Counter(); ends = []
    for vertex, edge in zip(raw['vertices'], edges):
        vertex['line'] = edge; slots[edge] += 1; ends.append((edge, slots[edge]))
    raw['phases'] = [{'loop': f'l{i}', 'expression': '-(' + '+'.join(
        [f'v{k}' for k in range(1, edge)] + [f's{edge}{slot}*v{edge}']) + ')'}
        for i, (edge, slot) in enumerate(ends, 1)]
    raw['bounds'] = [{'var': f's{edge}{slot}', 'low': '0',
                     'high': '1' if slot == counts[edge] else f's{edge}{slot+1}'} for edge, slot in ends]
    raw['iterated_order_outer_to_inner'] = [f's{edge}{slot}' for edge in sorted(counts)
                                          for slot in range(counts[edge], 0, -1)]
    return raw, permutation


@lru_cache(maxsize=16)
def rotation_case(index, rotation):
    original = SOURCE['orbits'][index]
    raw, permutation = rotated_source(original, rotation)
    return original, raw, permutation, geometry(validate_source(original)), geometry(validate_source(raw))


def check_geometry(index, rotation):
    original, raw, perm, old_geom, new_geom = rotation_case(index, rotation)
    old_data, new_data = validate_source(original), validate_source(raw)
    a = [sp.Symbol('alpha1'), sp.Symbol('alpha2'), 1-sp.Symbol('alpha1')-sp.Symbol('alpha2')]
    z = sp.Symbol('z1')
    pullback = {sp.Symbol(p['s_var']): sp.Symbol(old_data['endpoints'][perm[i]]['s_var'])
                for i, p in enumerate(new_data['endpoints'])}
    pullback.update({sp.Symbol(f'alpha{i+1}'): a[perm[i]] for i in range(2)})
    scale = z if rotation % 2 else sp.S.One
    old_z = 1/z if rotation % 2 else z
    checks = []
    lhs = expression(new_geom['q0']).subs(pullback, simultaneous=True)
    rhs = scale * expression(old_geom['q0']).subs(z, old_z)
    assert sp.expand(lhs-rhs) == 0
    checks.append('q0')
    for name in ('G', 'T', 'W'):
        for i in range(3):
            for j in range(3):
                lhs = expression(new_geom[name][i][j]).subs(pullback, simultaneous=True)
                rhs = scale * expression(old_geom[name][perm[i]][perm[j]]).subs(z, old_z)
                assert sp.expand(lhs-rhs) == 0, (original['orbit_id'], rotation, name, i, j)
                checks.append(f'{name}{i+1}{j+1}')
    # Full alpha/t permutation preserves both simplex weights and the max base.
    t = list(sp.symbols('t1 t2 t3')); aa = list(sp.symbols('a1 a2 a3'))
    factors = [aa[i]*t[(i+1)%3]*t[(i+2)%3] for i in range(3)]
    new_factors = [aa[perm[i]]*t[perm[(i+1)%3]]*t[perm[(i+2)%3]] for i in range(3)]
    assert set(new_factors) == set(factors)
    assert sp.prod(aa[i] for i in perm) == sp.prod(aa)
    # With the max base already identified, all actual exported weights agree.
    # M's marked alpha must follow the same leg permutation as its operator.
    old_blocks = {(b['channel'], b['block_id']): b for b in block_specs(old_data)}
    weight_checks = 0
    for block in block_specs(new_data):
        old_name = block['block_id']
        if block['kind'] == 'bubble':
            old_marked = perm[block['marked_bubble_leg']-1] + 1
            old_name = f'M{old_marked}'
        expected = old_blocks[(block['channel'], old_name)]
        def weight_expr(text):
            return expression(text.replace(MAX_BOUND, 'maxbase').replace('Beta[1-eps,1-eps]', 'bubble_beta'))
        lhs = weight_expr(block['weight']).subs(pullback, simultaneous=True)
        rhs = weight_expr(expected['weight'])
        assert sp.simplify(lhs-rhs) == 0
        weight_checks += 1
    # Bound arrays retain each source's order; only same-edge chains are pulled
    # back. Rotating independent chains changes their global iterated ordering.
    s_map = {p['s_var']: old_data['endpoints'][perm[i]]['s_var'] for i, p in enumerate(new_data['endpoints'])}
    for new_var, (low, high) in new_data['bounds'].items():
        assert [low, s_map.get(high, high)] == old_data['bounds'][s_map[new_var]]
    assert new_data['order'] == raw['iterated_order_outer_to_inner']
    assert sorted([tuple(s_map[v] for v in chain) for chain in new_data['ordered_chains']]) == sorted(
        [tuple(chain) for chain in old_data['ordered_chains']])
    assert c4_symmetry(raw)['orbit_multiplicity'] == 4 and c4_symmetry(raw)['included_in_C'] is False
    assert raw['source_normalization']['euler_compensation_count'] == 0
    return {'orbit_id': original['orbit_id'], 'rotation': rotation,
            'new_to_old_leg': [i+1 for i in perm], 'geometry_identities': len(checks),
            'source_order': raw['iterated_order_outer_to_inner'],
            'pulled_back_source_order': [s_map[v] for v in new_data['order']],
            'same_edge_chains_and_bounds_preserved': True, 'simplex_and_max_weight_bases_preserved': True,
            'raw_block_weight_equalities': weight_checks}


def numeric_geometry(geom, values):
    return {key: [[str(evaluate(text, values)) for text in row] for row in geom[key]] for key in ('G', 'T', 'W')}


def check_coefficients(index, rotation, point_index):
    original, raw, perm, old_geom, new_geom = rotation_case(index, rotation)
    old_data, new_data = validate_source(original), validate_source(raw)
    point = POINTS[point_index]
    old_values = dict(point)
    for endpoint, value in zip(old_data['endpoints'], [Fraction(1, 5), Fraction(2, 5), Fraction(4, 5)]):
        old_values[endpoint['s_var']] = value
    old_values['z1'] = 1/point['z1'] if rotation % 2 else point['z1']
    alphas = [point['alpha1'], point['alpha2'], 1-point['alpha1']-point['alpha2']]
    ts = [point['t1'], point['t2'], 1-point['t1']-point['t2']]
    new_values = dict(point, alpha1=alphas[perm[0]], alpha2=alphas[perm[1]], t1=ts[perm[0]], t2=ts[perm[1]])
    for i, endpoint in enumerate(new_data['endpoints']):
        new_values[endpoint['s_var']] = old_values[old_data['endpoints'][perm[i]]['s_var']]
    old_numeric = numeric_geometry(old_geom, old_values)
    new_numeric = numeric_geometry(new_geom, new_values)
    checks = []
    for channel, ops in operators().items():
        blocks = [('N0', '3*eps', ops['N0'], ops['N0']), ('L0', '-1+3*eps', ops['L0'], ops['L0'])]
        blocks += [(f'M{i+1}', '-1+3*eps', ops['M'][i], ops['M'][perm[i]]) for i in range(3)]
        for name, beta, new_op, old_op in blocks:
            new_coefficients = apply_operator(new_op, beta, new_numeric)
            old_coefficients = apply_operator(old_op, beta, old_numeric)
            shifts = [2, 3] if name == 'N0' else [1]
            assert set(new_coefficients) <= set(shifts) and set(old_coefficients) <= set(shifts)
            for shift in shifts:
                degree = shift if name == 'N0' else 2
                lhs = evaluate(new_coefficients.get(shift, '0'), new_values)
                rhs = evaluate(old_coefficients.get(shift, '0'), old_values)
                if rotation % 2:
                    rhs *= point['z1'] ** degree
                assert lhs == rhs, (original['orbit_id'], rotation, point_index, channel, name, shift, lhs-rhs)
                checks.append(f'{channel}:{name}:shift{shift}')
    return {'orbit_id': original['orbit_id'], 'rotation': rotation, 'point': point_index,
            'exact_rational_coefficient_equalities': len(checks)}


def check_sector_maps(index, rotation):
    original, raw, perm, _, _ = rotation_case(index, rotation)
    old_data, new_data = validate_source(original), validate_source(raw)
    alpha = [sp.Symbol('alpha1'), sp.Symbol('alpha2'), 1-sp.Symbol('alpha1')-sp.Symbol('alpha2')]
    u, v = sp.symbols('u v')
    # Sector chart fields depend on the source/block domain, not its numerator.
    # Keep that check narrow: a zero polynomial avoids exporting any new C/J.
    charts = []
    for data in (old_data, new_data):
        term = deepcopy(block_specs(data)[0])
        term.update(term_id='chart-only', kernel_terms=[], q0='1')
        charts.append(sectorize(term))
    checks = []
    for new_k in range(3):
        old_k = perm[new_k]
        new_other = [j for j in range(3) if j != new_k]
        old_other = [j for j in range(3) if j != old_k]
        chart_map = {sp.Symbol(f'alpha{i+1}'): alpha[perm[i]] for i in range(2)}
        chart_map.update({u: u if perm[new_other[0]] == old_other[0] else v,
                          v: v if perm[new_other[1]] == old_other[1] else u})
        old_chart, new_chart = charts[0][old_k]['sector'], charts[1][new_k]['sector']
        assert sp.expand(expression(new_chart['Z']).subs(chart_map, simultaneous=True)-expression(old_chart['Z'])) == 0
        assert sp.cancel(expression(new_chart['jacobian']).subs(chart_map, simultaneous=True)-expression(old_chart['jacobian'])) == 0
        old_t = [expression(old_chart['t_substitution'][f't{i+1}']) for i in range(2)]
        old_t.append(1-sum(old_t))
        new_t = [expression(new_chart['t_substitution'][f't{i+1}']) for i in range(2)]
        new_t.append(1-sum(new_t))
        for i in range(3):
            assert sp.cancel(new_t[i].subs(chart_map, simultaneous=True)-old_t[perm[i]]) == 0
        assert charts[1][new_k]['bounds']['u'] == ['1', 'Infinity']
        assert charts[1][new_k]['bounds']['v'] == ['1', 'Infinity']
        checks.append({'new_sector': new_k+1, 'old_sector': old_k+1, 'uv_swapped': chart_map[u] == v})
    return {'orbit_id': original['orbit_id'], 'rotation': rotation, 'sector_charts': checks}


@pytest.mark.parametrize('index', range(4))
@pytest.mark.parametrize('rotation', range(4))
def test_all_actual_source_C4_geometries_and_original_domains(index, rotation):
    check_geometry(index, rotation)


@pytest.mark.parametrize('index', range(4))
@pytest.mark.parametrize('rotation', range(4))
@pytest.mark.parametrize('point_index', range(2))
def test_C4_raw_coefficients_at_exact_rational_points(index, rotation, point_index):
    check_coefficients(index, rotation, point_index)


@pytest.mark.parametrize('index', range(4))
@pytest.mark.parametrize('rotation', range(4))
def test_C4_sector_parameter_mapping(index, rotation):
    check_sector_maps(index, rotation)
