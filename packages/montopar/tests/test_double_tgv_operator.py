"""Unapplied operator semantics, exact scalar conversion and additive export."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from montopar.cli import _compute_one
from montopar.double_tgv_operator import scalar_fields, to_scalar
from montopar.double_tgv_output import build_double_tgv_result, validate_double_tgv_result
from montopar.j_factor import expand_polynomial
from montopar.package_builder import build_initial_package
from montopar.parser import parse_orbit
from montopar.sign_rules import expected_source_contract

ROWS = json.loads((Path(__file__).parent / 'fixtures/double_tgv_orbits.json').read_text(encoding='utf-8'))


def row(number):
    return deepcopy(next(r for r in ROWS if r['orbit_id'].endswith(f'-{number:02}')))


def result(number, **kwargs):
    return build_double_tgv_result(parse_orbit(row(number)), operator_enabled=True, **kwargs)


@pytest.mark.parametrize('number', range(1, 27))
def test_operator_actually_recovers_all_original_scalar_fields_before_color_screening(number):
    out = result(number)
    fields = scalar_fields(out['J_operator'])
    assert fields
    for name, value in fields.items():
        assert expand_polynomial(f"({value})-({out['J_raw_s'][name]})") == '0'
    assert to_scalar(out['J_operator']) == out['J_raw_s']
    assert validate_double_tgv_result(out)
    assert out['J_operator']['combined_C'] == out['C']['combined']
    # All parent distances survive even when they vanish on the actual path.
    assert len(out['J_operator']['distance_definitions']) == 6
    assert 'd12' in out['J_operator']['kernel']['definitions'][-1]['expression']
    if out['J_raw_s']['physical_zero']:
        assert out['J_operator']['operator']['overall_coefficient'] == '-1'
        assert to_scalar(out['J_operator'])['integrand'] == '0'


@pytest.mark.parametrize('number', range(1, 27))
def test_new_export_does_not_change_any_existing_field(number):
    out = result(number)
    del out['J_operator']
    assert out == build_double_tgv_result(parse_orbit(row(number)))


@pytest.mark.parametrize('mutation', [
    'overall_sign', 'dimension', 'scope_order', 'fixed_scale', 'distance',
    'unknown_symbol', 'forward_definition', 'missing_definition', 'wrong_tensor',
    'unknown_operator', 'euler', 'C', 'outside_factor', 'domain', 'alpha_measure',
    'regulator', 'normal_response', 'left_slots', 'position', 'kernel',
])
def test_root_validator_rejects_changed_operator_contract(mutation):
    out = result(25)
    block = out['J_operator']
    if mutation == 'overall_sign': block['operator']['overall_coefficient'] = '1'
    elif mutation == 'dimension': block['operator']['dimension'] = '4'
    elif mutation == 'scope_order': block['restriction']['operation_order'].reverse()
    elif mutation == 'fixed_scale': block['held_fixed_under_position_derivatives'].remove('A0')
    elif mutation == 'distance': block['distance_definitions']['d12']['expression'] = '0'
    elif mutation == 'unknown_symbol': block['operator']['tensor_terms'][0][0]['coefficient'] = 'missing_symbol'
    elif mutation == 'forward_definition': block['kernel']['definitions'][0]['expression'] = 'dto_U'
    elif mutation == 'missing_definition': del block['kernel']['definitions'][0]
    elif mutation == 'wrong_tensor': block['operator']['tensor_terms'][0][1]['tensor'] = {'outer_legs': [0, 2]}
    elif mutation == 'unknown_operator': block['operator']['kind'] = 'arbitrary_code'
    elif mutation == 'euler': block['normalization']['euler_compensation_count'] = 1
    elif mutation == 'C': block['combined_C'] += '*2'
    elif mutation == 'outside_factor': block['normalization']['outside_kernel_factor'] = '1/A0'
    elif mutation == 'domain': block['bounds']['s11'] = ['0', '1/2']
    elif mutation == 'alpha_measure': block['alpha_measure']['jacobian'] = '2'
    elif mutation == 'regulator': block['branch_prescription']['finite_auxiliary_kernel'] = 'B0/A0+I*eta_reg'
    elif mutation == 'normal_response': block['restriction']['parent_normal_responses_retained'] = False
    elif mutation == 'left_slots': block['operator']['left_slots'] = [1, 3]
    elif mutation == 'position': block['restriction']['positions'][0]['edge_basis_components'][0] = '0'
    else: block['kernel']['power'] = '-1+2*eps'
    with pytest.raises(ValueError):
        validate_double_tgv_result(out)


@pytest.mark.parametrize('mutation', ['undefined', 'normal_order', 'outer_leg', 'D4', 'forward', 'code'])
def test_standalone_conversion_checks_its_own_payload(mutation):
    block = result(25)['J_operator']
    if mutation == 'undefined': block['operator']['tensor_terms'][0][0]['coefficient'] = 'unknown'
    elif mutation == 'normal_order': block['restriction']['operation_order'].reverse()
    elif mutation == 'outer_leg': block['operator']['tensor_terms'][0][1]['tensor']['outer_legs'] = [1, 7]
    elif mutation == 'D4': block['operator']['dimension'] = '4'
    elif mutation == 'forward': block['kernel']['definitions'][0]['expression'] = 'dto_U'
    else: block['kernel']['definitions'][0]['expression'] = "__import__('os')"
    with pytest.raises(ValueError):
        to_scalar(block)


def test_converter_uses_the_exported_tensors_and_keeps_nonzero_hessian():
    out = result(25)
    block = deepcopy(out['J_operator'])
    block['operator']['overall_coefficient'] = '1'
    wrong = scalar_fields(block)
    for name in ('pnum1', 'pnum2'):
        assert expand_polynomial(f"({wrong[name]})+({out['J_raw_s'][name]})") == '0'
        assert wrong[name] != '0'
    assert scalar_fields(result(4)['J_operator'])['pnum2'] != '0'


def test_export_does_not_depend_on_orbit_number_and_preserves_vertex_orientation():
    source = row(4)
    source['orbit_id'] = 'renamed-double-vertex'
    loops = source['internal_vertices'][0]['loops']
    source['internal_vertices'][0]['loops'] = [loops[1], loops[0], loops[2], loops[3]]
    out = build_double_tgv_result(parse_orbit(source), operator_enabled=True)
    assert to_scalar(out['J_operator']) == out['J_raw_s']
    assert out['color_contraction']['ratio_to_CA2_CF'] == '1/4'


def test_operator_forms_ibp_and_package_coexist(tmp_path):
    source = row(4)
    out = build_double_tgv_result(parse_orbit(source), operator_enabled=True,
                                 forms_enabled=True, ibp_enabled=True)
    folder = tmp_path / '3loop/4e-3g+3g' / source['orbit_id']
    folder.mkdir(parents=True)
    (folder / 'montopar.json').write_text(json.dumps(out), encoding='utf-8')
    package, _ = build_initial_package({'orbits': [source]}, tmp_path, '4e-3g+3g')
    for key in ('J_operator', 'J_ibp', 'J_compact', 'J_reduced'):
        assert package['orbits'][0][key] == out[key]


def test_compute_one_accepts_operator_option_and_rejects_other_families(tmp_path):
    source = row(25)
    source['physical_convention'] = expected_source_contract(source)
    _compute_one(parse_orbit(source), tmp_path, double_tgv_operator_enabled=True)
    out = json.loads((tmp_path / '3loop/4e-3g+3g' / source['orbit_id'] / 'montopar.json').read_text())
    assert to_scalar(out['J_operator']) == out['J_raw_s']
    with pytest.raises(ValueError, match='operator_requires_double_tgv'):
        _compute_one(SimpleNamespace(dispatch='tgv'), tmp_path, double_tgv_operator_enabled=True)
