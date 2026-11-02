"""Source rejection, complete tensor output and actual CLI Pi1 regression."""
from copy import deepcopy
import itertools
import json
from pathlib import Path
import subprocess
import sys

import pytest
from montopar.parser import parse_orbit
from montopar.self_energy import build_result, validate_result, c4_symmetry, validate_source
from montopar.package_builder import build_initial_package

FIXTURE = Path(__file__).parent / 'fixtures/self_energy_pi1_orbits.json'
RAWS = json.loads(FIXTURE.read_text(encoding='utf-8'))['orbits']


@pytest.mark.parametrize('raw', RAWS, ids=lambda o: o['orbit_id'])
def test_complete_ordered_raw(raw):
    orbit = parse_orbit(raw)
    assert orbit.dispatch == 'self_energy'
    result = build_result(orbit)
    validate_result(json.loads(json.dumps(result)))
    assert result['J_raw'] == result['J_raw_s']
    assert 'Derivative' not in result['J_raw']['integrand']
    assert not any(k.endswith('_rotated') for k in result)
    assert result['normalization']['euler_compensation_count'] == 0
    known = set()
    for var in result['J_raw']['iterated_order_outer_to_inner']:
        high = result['J_raw']['bounds'][var][1]
        assert high == '1' or high in known
        known.add(var)


@pytest.mark.parametrize('permutation', list(itertools.permutations(range(3))))
def test_propagator_array_order_does_not_change_math(permutation):
    raw = deepcopy(RAWS[4])
    original = build_result(parse_orbit(raw))
    raw['propagators'] = [raw['propagators'][i] for i in permutation]
    result = build_result(parse_orbit(raw))
    assert result['J_raw'] == original['J_raw']
    assert result['C'] == original['C']


@pytest.mark.parametrize('case', [
    'loop_order', 'family', 'operator_order', 'missing_operator', 'extra_operator',
    'wrong_se_momentum', 'extra_propagator', 'missing_propagator', 'duplicate_mu',
    'duplicate_color', 'external_color', 'internal_color', 'internal_mu',
    'duplicate_phase', 'missing_phase', 'wrong_phase', 'unknown_phase_symbol',
    'trailing_phase', 'duplicate_bound', 'missing_bound', 'wrong_bound',
    'cyclic_bound', 'wrong_endpoint', 'wrong_explicit_order', 'extra_vertex',
])
def test_bad_source_rejected(case):
    raw = deepcopy(RAWS[2])
    if case == 'loop_order': raw['nloops'] = 2
    elif case == 'family': raw['emission_type'] = '4e-Pi2'
    elif case == 'operator_order': raw['internal_vertices'][0]['se_type'] = 'vacuum-pol-2loop'
    elif case == 'missing_operator': raw['internal_vertices'] = []
    elif case == 'extra_operator': raw['internal_vertices'].append(deepcopy(raw['internal_vertices'][0]))
    elif case == 'wrong_se_momentum': raw['internal_vertices'][0]['propagator'] = 'l1'
    elif case == 'extra_propagator': raw['propagators'].append(deepcopy(raw['propagators'][0]))
    elif case == 'missing_propagator': raw['propagators'].pop()
    elif case == 'duplicate_mu': raw['vertices'][1]['lorentz'] = 'mu1'
    elif case == 'duplicate_color': raw['vertices'][1]['color'] = 'a1'
    elif case == 'external_color': raw['propagators'][0]['color'][0] = 'wrong'
    elif case == 'internal_color': raw['propagators'][1]['color'][1] = 'a1'
    elif case == 'internal_mu': raw['propagators'][2]['lorentz'][0] = 'nu1'
    elif case == 'duplicate_phase': raw['phases'].append(deepcopy(raw['phases'][0]))
    elif case == 'missing_phase': raw['phases'].pop()
    elif case == 'wrong_phase': raw['phases'][1]['expression'] = '(1-s32)*v3+s41*v4'
    elif case == 'unknown_phase_symbol': raw['phases'][0]['expression'] += '+q'
    elif case == 'trailing_phase': raw['phases'][0]['expression'] += ' 2'
    elif case == 'duplicate_bound': raw['bounds'][0] = deepcopy(raw['bounds'][1])
    elif case == 'missing_bound': raw['bounds'].pop()
    elif case == 'wrong_bound': raw['bounds'][0]['high'] = 's42'
    elif case == 'cyclic_bound': raw['bounds'][0]['high'] = 's31'
    elif case == 'wrong_endpoint': raw['propagators'][0]['lorentz'][0] = 'mu1'
    elif case == 'wrong_explicit_order': raw['iterated_order_outer_to_inner'] = ['s31', 's32', 's41', 's42']
    elif case == 'extra_vertex': raw['vertices'].append(deepcopy(raw['vertices'][0]))
    with pytest.raises(ValueError): parse_orbit(raw)


def test_unknown_internal_does_not_fall_back_to_pure():
    raw = deepcopy(RAWS[0])
    raw['emission_type'] = '2e-unknown'
    raw['internal_vertices'] = [{'type': 'unknown'}]
    assert parse_orbit(raw).dispatch == 'unsupported'


def test_two_insertions_not_a_pi2():
    raw = deepcopy(RAWS[0]); raw['internal_vertices'] *= 2
    with pytest.raises(ValueError): parse_orbit(raw)


def test_opposite_se_survives_zero_metric_limit():
    # At z=-1 the metric contraction d=0, but u*w is not zero.
    result = build_result(parse_orbit(RAWS[1]))
    assert result['J_raw']['derivative_numerator'] != '0'


def test_marked_se_bare_c4_does_not_use_pure_topology():
    assert c4_symmetry(RAWS[0])['orbit_multiplicity'] == 4
    assert c4_symmetry(RAWS[1])['orbit_multiplicity'] == 2
    assert c4_symmetry(RAWS[-1])['orbit_multiplicity'] == 2


def test_all_produced_formula_blocks_reject_tampering():
    result = json.loads(json.dumps(build_result(parse_orbit(RAWS[0]))))
    for key, field in [('C', 'combined'), ('J_raw', 'integrand'), ('J_raw_s', 'integrand'),
                       ('operator', 'C_SE'), ('normalization', 'euler_compensation_count')]:
        altered = deepcopy(result); altered[key][field] = 'wrong'
        with pytest.raises(ValueError): validate_result(altered)


@pytest.mark.parametrize('family,expected', [('2e-Pi1', 2), ('4e-Pi1', 11)])
def test_actual_cli_package_and_unrelated_unsupported_source(tmp_path, family, expected):
    raw_source = {'orbits': deepcopy(RAWS)}
    unsupported = deepcopy(RAWS[0]); unsupported['emission_type'] = '2e-Pi2'
    raw_source['orbits'].append(unsupported)
    source = tmp_path / 'source.json'
    source.write_text(json.dumps(raw_source), encoding='utf-8')
    run = subprocess.run([sys.executable, '-B', '-m', 'montopar.cli', '--source', str(source),
                          '--type', family, '--out-dir', str(tmp_path), '--package-initial'],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    package, _ = build_initial_package(raw_source, tmp_path, family)
    assert len(package['orbits']) == expected
    assert package['schema_version'] == '1.4'
    assert all(p['J'] == p['J_raw_s'] for p in package['orbits'])


def test_stale_in_memory_metadata_rejected():
    orbit = parse_orbit(RAWS[0]); orbit.n_gluons = 7
    with pytest.raises(ValueError): build_result(orbit)


def test_loop_reversal_and_delta_dummy_renaming():
    raw = deepcopy(RAWS[0])
    original = build_result(parse_orbit(raw))
    raw['phases'][0]['expression'] = '-(' + raw['phases'][0]['expression'] + ')'
    raw['propagators'][0]['lorentz'][1] = 'nu7'
    raw['propagators'][1]['lorentz'][0] = 'nu9'
    raw['propagators'][0]['color'][1] = 'b7'
    raw['propagators'][1]['color'][0] = 'b9'
    changed = build_result(parse_orbit(raw))
    assert original['C'] == changed['C']
    assert original['J_raw'] == changed['J_raw']


def test_uninterpreted_physics_metadata_is_rejected():
    raw = deepcopy(RAWS[0]); raw['overall_sign'] = 1
    with pytest.raises(ValueError): parse_orbit(raw)


@pytest.mark.parametrize('change', ['dispatch', 'integrand', 'source'])
def test_package_cannot_downgrade_or_detach_self_energy(tmp_path, change):
    raw = deepcopy(RAWS[0])
    result = json.loads(json.dumps(build_result(parse_orbit(raw))))
    if change == 'dispatch': result['dispatch'] = 'pure'
    if change == 'integrand': result['J_raw']['integrand'] = '0'
    if change == 'source': raw['phases'][0]['expression'] = '-('+raw['phases'][0]['expression']+')'
    path = tmp_path / '2loop/2e-Pi1/2e-Pi1-01/montopar.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(result), encoding='utf-8')
    with pytest.raises(ValueError): build_initial_package({'orbits': [raw]}, tmp_path, '2e-Pi1')
