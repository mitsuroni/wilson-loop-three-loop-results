"""Ordered Gamma1 boundaries, finite-regulator reconstruction and rejection."""
from copy import deepcopy
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from montopar.gamma1_source import validate_source
from montopar.gamma1_kernel import geometry, apply_operator, operators
from montopar.gamma1_ibp import build_ibp, flux_data, polynomial_substitute, _inline_substitute


SOURCE = json.loads((Path(__file__).parent/'fixtures/gamma1_orbits.json').read_text(encoding='utf-8'))['orbits']


@lru_cache(maxsize=4)
def prepared(index):
    data = validate_source(SOURCE[index])
    geom = geometry(data)
    return data, geom, build_ibp(data, geom, {'g': 'Cg', 'q': 'Cq'})


@pytest.mark.parametrize('index', range(4))
def test_full_tensor_scalar_reconstruction_and_ordered_faces(index):
    data, geom, result = prepared(index)
    for channel in ('g', 'q'):
        for name in ('N0', 'L0', 'M1', 'M2', 'M3'):
            block = result['flux'][channel][name]
            assert block['coefficient_reconstruction'] == 'exact_at_each_q_power_with_finite_common_I0'
            if index != 3:
                assert block['remainder_operator'] == '0'
                inactive = 2 if index == 1 else 0
                assert block['flux_operators'][inactive] == '0'
    faces = {term['face']['face_id'] for term in result['terms'] if 'face' in term}
    assert len(faces) == (6 if index == 3 else 3)
    assert any(term['kind'] == 'contact' for term in result['terms']) == (index != 3)
    assert any(term['kind'] == 'recoil' for term in result['terms']) == (index == 3)
    for term in result['terms']:
        assert term['combined_C'] == {'g': 'Cg', 'q': 'Cq'}[term['channel']]
        assert term['variables'] == term['iterated_order_outer_to_inner']
        assert set(term['variables']) == set(term['bounds'])
        original = [v for v in data['order'] if v in term['variables']]
        assert term['variables'][:len(original)] == original
        assert len(term['variables']) == ((7 if term['integration_block_kind'] == 'triangle' else 5) - (term['kind'] != 'recoil'))
        for forbidden in ('r11', 'p1sq', 'N0', 'L0', 'delta', 'alpha3', 't3'):
            assert forbidden not in term['integrand']
        assert 'I0' in term['integrand']
        assert term['derivatives_before_restriction']


def test_334_faces_keep_remaining_third_leg_in_original_position():
    data, geom, result = prepared(1)
    for term in result['terms']:
        path_order = [v for v in term['variables'] if v.startswith('s')]
        assert path_order[-1] == 's41'
        if term['kind'] == 'contact':
            assert term['face']['normal'] == {'s31': '1', 's32': '-1', 's41': '0'}


def test_simultaneous_restriction_does_not_cascade_or_drop_parentheses():
    assert polynomial_substitute('2*x*y+x^2', {'x': 'y+1', 'y': '0'}) == '1 + 2*y + y^2'
    assert _inline_substitute('x^2+y', {'x': 'y+1', 'y': '0'}) == '(y+1)^2+(0)'


def test_flux_cache_returns_independently_owned_results():
    data, geom, _ = prepared(0)
    first = flux_data(data, geom)
    first['g']['N0']['flux_operators'][1] = '0'
    first['g']['N0']['evaluated_flux'][1].clear()
    second = flux_data(data, geom)
    assert second['g']['N0']['flux_operators'][1] != '0'
    assert second['g']['N0']['evaluated_flux'][1]


def test_wrong_parent_geometry_cannot_pass_divergence_check():
    data, geom, _ = prepared(0)
    damaged = deepcopy(geom)
    damaged['q0'] += '+s41'
    with pytest.raises(ValueError, match='coefficient_reconstruction_failed'):
        flux_data(data, damaged)


def test_unvalidated_first_third_repeat_is_explicitly_rejected():
    data, geom, _ = prepared(0)
    damaged = deepcopy(data)
    damaged['endpoints'][0]['edge'] = 4
    damaged['endpoints'][1]['edge'] = 3
    with pytest.raises(ValueError, match='noncanonical_repeated_first_third_leg'):
        flux_data(damaged, geom)


def _integrate(function, order, bounds, points):
    nodes, weights = np.polynomial.legendre.leggauss(points)
    def walk(index, values):
        if index == len(order):
            return function(values)
        variable = order[index]
        low, high = bounds[variable]
        a = values[low] if low in values else float(low)
        b = values[high] if high in values else float(high)
        return (b-a)/2*sum(weight*walk(index+1, {**values, variable: (a+b)/2+(b-a)/2*node})
                           for node, weight in zip(nodes, weights))
    return walk(0, {})


def _numeric_kernel(coefficients, beta, q0, variables):
    fixed = {'alpha1': '1/5', 'alpha2': '1/4', 't1': '1/7', 't2': '1/3',
             'eps': '1/8', 'z1': '2/3'}
    # Substitute the fixed rational parameters before SymPy compilation: the
    # fully expanded four-geometry N0 polynomial can exceed Python's AST depth
    # despite being a low-degree polynomial in the three path coordinates.
    parse = lambda text: sp.sympify(polynomial_substitute(text, fixed).replace('^', '**'))
    q = parse(q0)+1  # A finite positive constant tests the unremoved regulator.
    power = parse(beta)
    expression = sum(parse(value)*q**(power-int(shift)) for shift, value in coefficients.items())
    function = sp.lambdify(variables, expression, modules='math', cse=True)
    return lambda values: function(*(values[name] for name in variables))


@pytest.mark.parametrize('index', range(4))
def test_independent_cubature_of_raw_equals_oriented_faces_and_recoil(index):
    """Numerical contour integration detects missing contact signs/measures."""
    data, geom, result = prepared(index)
    # N0 exercises the second derivatives in the flux, including metric traces.
    for channel in ('g', 'q'):
        coefficients = apply_operator(operators()[channel]['N0'], '3*eps', geom)
        raw_function = _numeric_kernel(coefficients, '3*eps', geom['q0'], data['order'])
        original = _integrate(raw_function, data['order'], data['bounds'], 9)
        boundary = 0.0
        for term in result['terms']:
            if term['channel'] != channel or term['block_id'] != 'N0':
                continue
            order = [name for name in term['variables'] if name.startswith('s')]
            function = _numeric_kernel(term['kernel_coefficients'], term['primitive_power'], term['q0'], order)
            boundary += _integrate(function, order, term['bounds'], 9)
        assert abs(original-boundary) < 2e-11*max(1, abs(original), abs(boundary))
