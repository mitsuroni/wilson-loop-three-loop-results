"""Strict Gamma1 source admission independent of C/J formula generation."""
from collections import Counter
from copy import deepcopy
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path

import pytest

from montopar.gamma1_source import (
    OPERATOR_ID, OPERATOR_SHA256, c4_symmetry, parse_gamma1_orbit, validate_source,
)

SOURCE = json.loads((Path(__file__).parent / 'fixtures' / 'gamma1_orbits.json').read_text(encoding='utf-8'))


def from_edges(edges):
    """Construct an independent source geometry, preserving the pinned operator."""
    raw = deepcopy(SOURCE['orbits'][0])
    raw['orbit_id'] = 'independent-source-' + ''.join(map(str, edges))
    counts = Counter(edges); slots = Counter(); ends = []
    for v, edge in zip(raw['vertices'], edges):
        v['line'] = edge; slots[edge] += 1; ends.append((edge, slots[edge]))
    raw['phases'] = [{'loop': f'l{i}', 'expression': '-(' + '+'.join(
        [f'v{k}' for k in range(1, edge)] + [f's{edge}{slot}*v{edge}']) + ')'}
        for i, (edge, slot) in enumerate(ends, 1)]
    raw['bounds'] = [{'var': f's{edge}{slot}', 'low': '0',
                     'high': '1' if slot == counts[edge] else f's{edge}{slot+1}'} for edge, slot in ends]
    raw['iterated_order_outer_to_inner'] = [f's{edge}{slot}' for edge in sorted(counts)
                                          for slot in range(counts[edge], 0, -1)]
    return raw


@pytest.mark.parametrize('raw', SOURCE['orbits'], ids=lambda o: o['orbit_id'])
def test_frozen_source_preserves_all_connections_and_original_order(raw):
    before = deepcopy(raw)
    data = validate_source(raw); orbit = parse_gamma1_orbit(raw)
    assert raw == before
    assert orbit.dispatch == 'gamma1'
    assert orbit.n_gluons == 0 and orbit.pairings == []
    assert orbit.source == raw and orbit.source is not raw
    assert data['order'] == raw['iterated_order_outer_to_inner']
    assert list(data['bounds']) == data['order']
    assert [p['edge'] for p in data['endpoints']] == [v['line'] for v in raw['vertices']]
    assert [p['s_var'] for p in data['endpoints']] == [b['var'] for b in raw['bounds']]
    assert data['loops'] == ['l1', 'l2', 'l3']
    assert c4_symmetry(raw)['orbit_multiplicity'] == 4
    assert c4_symmetry(raw)['stabilizer_order'] == 1
    assert c4_symmetry(raw)['stabilizer_rotations'] == [0]
    assert c4_symmetry(raw)['included_in_C'] is False


@pytest.mark.parametrize('edges', [e for e in combinations_with_replacement(range(1, 5), 3) if len(set(e)) > 1])
def test_all_legal_sorted_geometries_are_computed_from_endpoints(edges):
    raw = from_edges(edges)
    data = validate_source(raw)
    assert [p['edge'] for p in data['endpoints']] == list(edges)
    assert parse_gamma1_orbit(raw).orbit_id == raw['orbit_id']
    assert c4_symmetry(raw)['orbit_multiplicity'] == 4


def test_coherent_momentum_index_renaming_and_storage_reordering():
    raw = deepcopy(SOURCE['orbits'][0])
    iv = raw['internal_vertices'][0]
    for i in range(3):
        ext, inner, color, inner_color, loop = f'M{i}', f'N{i}', f'A{i}', f'B{i}', f'k{i}'
        raw['vertices'][i].update(lorentz=ext, color=color)
        raw['propagators'][i].update(loop=loop, lorentz=[inner, ext], color=[inner_color, color])
        raw['phases'][i]['loop'] = loop
        iv['loops'][i] = loop; iv['incoming_momenta'][i] = '-' + loop
        iv['lorentz'][i] = inner; iv['color'][i] = inner_color
    raw['vertex_correction']['external_momentum_labels'] = deepcopy(iv['loops'])
    raw['momentum_constraints'][0].update(expression='k0+k1+k2', normalization='(2*pi)^D*delta^D(k0+k1+k2)')
    for key in ('propagators', 'phases', 'bounds'):
        raw[key].reverse()
    from montopar.sign_rules import expected_contract
    raw['physical_convention']=expected_contract(raw['vertices'])
    data = validate_source(raw)
    assert data['loops'] == ['k0', 'k1', 'k2']
    assert data['order'] == ['s31', 's42', 's41']
    assert [leg.loop for leg in parse_gamma1_orbit(raw).tgv_legs] == data['loops']


def test_equivalent_endpoint_polynomial_and_momentum_sum():
    raw = deepcopy(SOURCE['orbits'][0])
    raw['phases'][0]['expression'] = '-v1-v2-s31*v3+(v1+v2+v3+v4)'
    raw['momentum_constraints'][0]['expression'] = 'l2+l3+l1'
    validate_source(raw)


CHANGES = [
    (['emission_type'], '3e-3g'), (['sector'], 'ghost'), (['nloops'], 2), (['nloops'], True),
    (['g_power'], 4), (['orbit_id'], ''), (['sign_convention'], 'legacy'),
    (['vertices', 0, 'line'], 0), (['vertices', 0, 'line'], True),
    (['vertices', 1, 'line'], 2), (['vertices', 1, 'lorentz'], 'mu1'),
    (['vertices', 1, 'color'], 'a1'), (['vertices', 0, 'lorentz'], ['mu1']),
    (['internal_vertices', 0, 'type'], 'tgv'), (['internal_vertices', 0, 'vertex_loop_order'], 0),
    (['internal_vertices', 0, 'vertex_loop_order'], True),
    (['internal_vertices', 0, 'loops'], ['l1', 'l1', 'l3']),
    (['internal_vertices', 0, 'incoming_momenta'], ['l1', '-l2', '-l3']),
    (['internal_vertices', 0, 'lorentz'], ['nu1', 'nu1', 'nu3']),
    (['internal_vertices', 0, 'color'], ['b1', 'a1', 'b3']),
    (['internal_vertices', 0, 'operator_id'], 'Gamma1-four-dimensional'),
    (['internal_vertices', 0, 'operator_sha256'], '0' * 64),
    (['propagators', 0, 'loop'], 'l2'), (['propagators', 0, 'color'], ['a1', 'b2']),
    (['propagators', 0, 'lorentz'], ['mu1', 'nu2']),
    (['phases', 0, 'loop'], 'l2'), (['phases', 0, 'expression'], '(v1+v2+s31*v3)'),
    (['phases', 0, 'expression'], '-v1-v2-s41*v3'),
    (['phases', 0, 'expression'], '-v1-v2-s31*v3+secret'),
    (['phases', 0, 'expression'], '-v1-v2-s31*v3)'),
    (['phases', 0, 'expression'], ''),
    (['bounds', 1, 'high'], '1'), (['bounds', 0, 'low'], '1/2'),
    (['bounds', 1, 'var'], 's31'), (['bounds', 0, 'high'], 1),
    (['iterated_order_outer_to_inner'], ['s42', 's31', 's41']),
    (['iterated_order_outer_to_inner'], ['s31', 's41', 's42']),
    (['momentum_constraints', 0, 'expression'], 'l1+l2-l3'),
    (['momentum_constraints', 0, 'expression'], 'l1+l2+l3+1'),
    (['momentum_constraints', 0, 'normalization'], 'delta^D(l1+l2+l3)'),
    (['vertex_correction', 'kind'], 'connected'), (['vertex_correction', 'operator'], 'Pi1'),
    (['vertex_correction', 'subgraph_loop_order'], 2),
    (['vertex_correction', 'independent_external_momenta'], 3),
    (['vertex_correction', 'external_leg_self_energies_included'], True),
    (['vertex_correction', 'counterterm_included'], True),
    (['vertex_correction', 'internal_loop_measure_already_in_operator'], False),
    (['vertex_correction', 'source_external_sign'], '-1'),
    (['vertex_correction', 'gauge'], 'Landau'), (['vertex_correction', 'dimension'], 'D=4'),
    (['vertex_correction', 'full_color_channels'], ['CF*CA^2']),
    (['vertex_correction', 'operator_definition', 'form_factors', 'quark', 'A'], '0'),
    (['vertex_correction', 'operator_definition', 'dirac_trace_one'], 2),
    (['vertex_correction', 'operator_definition', 'dimension'], 'D=4'),
    (['source_normalization', 'normalized_wilson_trace'], 'Tr'),
    (['source_normalization', 'euler_compensation_count'], 3),
    (['source_normalization', 'euler_compensation_count'], False),
    (['source_normalization', 'orbit_multiplicity_included'], True),
    (['source_normalization', 'color_projection'], 'a1'),
    (['source_normalization', 'crossing_assembled'], True),
    (['c4_symmetry', 'orbit_multiplicity'], 1), (['c4_symmetry', 'stabilizer_order'], 4),
    (['downstream_status', 'I'], 'solved'),
]


@pytest.mark.parametrize('path,value', CHANGES)
def test_reject_modified_mathematical_source(path, value):
    raw = deepcopy(SOURCE['orbits'][0]); target = raw
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError, match='gamma1_'):
        validate_source(raw)


@pytest.mark.parametrize('key', list(SOURCE['orbits'][0]))
def test_every_source_field_is_required(key):
    raw = deepcopy(SOURCE['orbits'][0]); del raw[key]
    with pytest.raises(ValueError, match='gamma1_'):
        validate_source(raw)


@pytest.mark.parametrize('path', [[], ['vertices', 0], ['propagators', 0], ['phases', 0],
    ['bounds', 0], ['internal_vertices', 0], ['momentum_constraints', 0],
    ['vertex_correction'], ['source_normalization'], ['c4_symmetry'], ['downstream_status']])
def test_unknown_source_fields_are_not_silently_ignored(path):
    raw = deepcopy(SOURCE['orbits'][0]); target = raw
    for key in path:
        target = target[key]
    target['unexpected_factor'] = '2'
    with pytest.raises(ValueError, match='gamma1_'):
        validate_source(raw)


@pytest.mark.parametrize('key', ['vertices', 'propagators', 'phases', 'bounds', 'internal_vertices',
    'momentum_constraints', 'vertex_correction', 'source_normalization'])
def test_malformed_nested_structures_fail_cleanly(key):
    raw = deepcopy(SOURCE['orbits'][0]); raw[key] = None
    with pytest.raises(ValueError, match='gamma1_'):
        validate_source(raw)


def test_self_consistent_forged_operator_digest_is_rejected():
    raw = deepcopy(SOURCE['orbits'][0]); operator = raw['vertex_correction']['operator_definition']
    operator['form_factors']['gluon_ghost']['A'] = '0'
    digest = hashlib.sha256(json.dumps(operator, ensure_ascii=False, sort_keys=True,
                                      separators=(',', ':')).encode()).hexdigest()
    assert digest != OPERATOR_SHA256
    raw['internal_vertices'][0]['operator_sha256'] = digest
    with pytest.raises(ValueError, match='gamma1_central_operator_identity'):
        validate_source(raw)


def test_pinned_digest_and_operator_identity_are_fixture_facts():
    assert SOURCE['operator_sha256'] == OPERATOR_SHA256
    assert SOURCE['operator_id'] == OPERATOR_ID
    assert hashlib.sha256((Path(__file__).parent / 'fixtures' / 'gamma1_orbits.json').read_bytes()).hexdigest() == 'b0be12caa7aa27bdec547d2b546ee0090c1525c8aa7e1a644f914cd056c92f0b'


@pytest.mark.parametrize('edges', [(1, 1, 1), (4, 4, 4), (4, 3, 3)])
def test_reject_all_same_edge_and_nonordered_contour(edges):
    with pytest.raises(ValueError, match='gamma1_contour_edges_or_scaleless'):
        validate_source(from_edges(edges))
