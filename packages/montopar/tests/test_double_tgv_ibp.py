"""Public output/CLI contract, mutation rejection and contact-normal regression."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import pytest
from montopar.parser import parse_orbit
from montopar.double_tgv_output import build_double_tgv_result, validate_double_tgv_result
from montopar.double_tgv_ibp import ordered_faces
from montopar.cli import _compute_one

ROWS = json.loads((Path(__file__).parent / 'fixtures/double_tgv_orbits.json').read_text())


@pytest.fixture(scope='module')
def outputs():
    return {r['orbit_id']: build_double_tgv_result(parse_orbit(r), ibp_enabled=True) for r in ROWS}


@pytest.mark.parametrize('row', ROWS, ids=lambda r: r['orbit_id'])
def test_raw_unchanged_and_ibp_self_contained(row, outputs):
    out = outputs[row['orbit_id']]
    raw = build_double_tgv_result(parse_orbit(row))
    assert {k: v for k, v in out.items() if k != 'J_ibp'} == raw
    assert validate_double_tgv_result(out)
    ibp = out['J_ibp']
    assert ibp['C_raw'] == out['C']
    assert ibp['normalization']['euler_compensation_count'] == 0
    for term in ibp['terms']:
        assert len(term['ordered_variables']) in (2, 3, 4)
        assert set(term['ordered_variables']) == set(term['ordered_bounds'])
        assert len(term['alpha_measure']['independent_variables']) == 4
        assert term['regulator_group'] == out['orbit_id'] + ':all_ibp_terms'
        if ibp['physical_zero']:
            assert term['integrand'] == term['integrand_s'] == '0'
            assert term['geometric_integrand_s_before_color_screening'] != '0'


@pytest.mark.parametrize('number,expected', [(4, (9,0,0,0)), (6,(7,6,1,1)), (25,(16,8,1,1))])
def test_distinct_degeneracies(number, expected, outputs):
    assert tuple(outputs[f'4e-3g+3g-{number:02d}']['J_ibp']['counts'].values()) == expected


@pytest.mark.parametrize('kind', ['drop_normal','boundary_sign','missing_term','wrong_domain',
                                  'wrong_C','wrong_D_trace','wrong_regulator','invented_term'])
def test_validator_rejects_ibp_corruption(kind, outputs):
    out = deepcopy(outputs['4e-3g+3g-01' if kind == 'drop_normal' else '4e-3g+3g-06'])
    ibp = out['J_ibp']
    if kind == 'drop_normal':
        face = next(f for f in ibp['faces'] if f['chi_over_v14'] != '0')
        face['chi_over_v14'] = '0'
    elif kind == 'boundary_sign': ibp['terms'][0]['C_ratio_to_raw'] = '-1'
    elif kind == 'missing_term': ibp['terms'].pop()
    elif kind == 'wrong_domain': ibp['terms'][0]['ordered_bounds'] = {}
    elif kind == 'wrong_C': ibp['terms'][0]['combined_C'] = out['C']['combined']
    elif kind == 'wrong_D_trace': ibp['operator_data']['H_RR_over_A0'] = '0'
    elif kind == 'wrong_regulator': ibp['regulator']['kernel'] = 'B5'
    elif kind == 'invented_term': ibp['terms'].append(deepcopy(ibp['terms'][0]))
    with pytest.raises(ValueError, match='mismatch_J_ibp'):
        validate_double_tgv_result(out)


def test_ordered_contact_projection_is_necessary(outputs):
    # Exact smooth polynomial counterexample from r003 section 9.4.
    # TT=-z*d_s41*d_s42 on f=s41*s42, volume of the triple chain = 1/6.
    from montopar.j_factor import expand_polynomial
    ibp = outputs['4e-3g+3g-01']['J_ibp']
    normals = [f for f in ibp['faces'] if f['chi_over_v14'] != '0']
    assert len(normals) == 1
    assert normals[0]['removed'] == 's41' and normals[0]['target'] == 's42'
    assert expand_polynomial(normals[0]['chi_over_v14']) == 'z1'
    direct, tangent, normal = -Fraction(1,6), -Fraction(1,3), Fraction(1,6)
    assert direct == tangent + normal and direct != tangent


def test_cli_option_reaches_double_tgv(tmp_path, outputs):
    source = next(r for r in ROWS if r['orbit_id'].endswith('-06'))
    _compute_one(parse_orbit(source), tmp_path, tgv_ibp_enabled=True)
    file = tmp_path/'3loop'/'4e-3g+3g'/source['orbit_id']/'montopar.json'
    assert json.loads(file.read_text()) == outputs[source['orbit_id']]


def test_faces_preserve_triple_chain_contact():
    faces = ordered_faces([['a','b','c'], ['d']])
    face = next(f for f in faces if f['face_id'] == 'contact_a_b')
    assert face['normal'] == {'a':'1','b':'-1','c':'0','d':'0'}
    assert face['ordered_chains'] == [['b','c'], ['d']]
    assert face['parent_to_face'] == {'a':'b','b':'b','c':'c','d':'d'}
