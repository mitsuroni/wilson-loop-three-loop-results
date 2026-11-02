from montopar.sign_rules import validate_contract, annotate_result
"""Closed Stage 1 contract for a bare central one-loop three-gluon vertex.

The pinned DOT operator is part of the input, not a runtime CAS dependency.
Its content digest is checked against this implementation's fixed digest (not
merely against a mutable digest supplied by the input). Geometry is reconstructed
from the contour endpoints, independently of the orbit's descriptive identifier.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import json
import re

from montopar.j_factor import expand_polynomial
from montopar.models import Bound, EdgeDistribution, Orbit, Propagator, TGVLeg, Vertex
from montopar.self_energy import _poly, _position, _closure

FAMILY = '3e-Gamma1'
CONVENTION = 'wilson-central-three-gluon-1pi/1.0'
OPERATOR_ID = 'Gamma1-DOT-v3-Feynman-bare-r001'
OPERATOR_SHA256 = '8390815a76cb72b38f6b160839d6cfcf9edbd331a04c3e731c246069f89f748a'


def _fail(reason):
    raise ValueError('gamma1_' + reason)


def _record(value, keys, reason):
    if not isinstance(value, dict) or set(value) != set(keys):
        _fail(reason)


def _records(value, count, keys, reason):
    if not isinstance(value, list) or len(value) != count:
        _fail(reason)
    for item in value:
        _record(item, keys, reason)


def _label(value):
    return isinstance(value, str) and bool(re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', value))


def _same(actual, expected):
    """JSON type-sensitive equality: true must not silently stand for 1."""
    return json.dumps(actual, sort_keys=True, ensure_ascii=False) == json.dumps(expected, sort_keys=True, ensure_ascii=False)


def _symmetry(endpoints):
    def rotated(k):
        return tuple(sorted(((p['edge'] - 1 + k) % 4 + 1, p['slot']) for p in endpoints))
    stabilizer = [k for k in range(4) if rotated(k) == rotated(0)]
    return {'orbit_multiplicity': 4 // len(stabilizer),
            'stabilizer_order': len(stabilizer), 'stabilizer_rotations': stabilizer,
            'basis': 'three ordered contour endpoints; central Gamma1 operator',
            'included_in_C': False}


def validate_source(raw):
    """Validate all momentum, operator, ordering and normalization information.

    The source's exact iteration order is retained. No independent-variable
    permutation, extra loop measure, counterterm, external-leg Pi insertion,
    Euler factor or C4 multiplication is admitted by this raw-source contract.
    """
    required = {'orbit_id', 'emission_type', 'sector', 'nloops', 'g_power', 'vertices',
                'propagators', 'internal_vertices', 'phases', 'bounds',
                'iterated_order_outer_to_inner', 'momentum_constraints',
                'sign_convention', 'vertex_correction', 'source_normalization',
                'c4_symmetry', 'downstream_status'}
    _record(raw, required | {'physical_convention'}, 'unknown_or_missing_source_fields')
    if (not isinstance(raw['orbit_id'], str) or not raw['orbit_id'].strip()
            or raw['emission_type'] != FAMILY or raw['sector'] != 'dg'
            or type(raw['nloops']) is not int or raw['nloops'] != 3
            or type(raw['g_power']) is not int or raw['g_power'] != 6):
        _fail('family_or_loop_order')
    if raw['sign_convention'] != CONVENTION:
        _fail('sign_convention')

    vertices = raw['vertices']
    _records(vertices, 3, {'line', 'lorentz', 'color'}, 'external_vertices')
    edges = [v['line'] for v in vertices]
    if (any(type(e) is not int or e not in range(1, 5) for e in edges)
            or edges != sorted(edges) or len(set(edges)) < 2):
        _fail('contour_edges_or_scaleless')
    for key in ('lorentz', 'color'):
        indices = [v[key] for v in vertices]
        if not all(_label(i) for i in indices) or len(set(indices)) != 3:
            _fail('external_index_bijection')
    slots = Counter()
    try:validate_contract(raw,required=True)
    except ValueError as error:_fail(str(error))
    endpoints = []
    for vertex in vertices:
        edge = vertex['line']; slots[edge] += 1
        endpoints.append(dict(vertex, edge=edge, slot=slots[edge], s_var=f's{edge}{slots[edge]}'))
    expected_bounds = {p['s_var']: ['0', f"s{p['edge']}{p['slot'] + 1}" if p['slot'] < slots[p['edge']] else '1']
                       for p in endpoints}
    bounds = raw['bounds']
    _records(bounds, 3, {'var', 'low', 'high'}, 'ordered_bounds')
    if (any(not isinstance(b['var'], str) for b in bounds)
            or len({b['var'] for b in bounds}) != 3
            or not _same({b['var']: [b['low'], b['high']] for b in bounds}, expected_bounds)):
        _fail('ordered_bounds')
    order = [f's{edge}{slot}' for edge in sorted(slots) for slot in range(slots[edge], 0, -1)]
    if raw['iterated_order_outer_to_inner'] != order:
        _fail('original_iterated_order')

    ivs = raw['internal_vertices']
    _records(ivs, 1, {'type', 'id', 'loops', 'lorentz', 'color', 'incoming_momenta',
                      'vertex_loop_order', 'operator_id', 'operator_sha256'}, 'central_vertex')
    vertex = ivs[0]
    if (vertex['type'] != 'gamma1' or not _label(vertex['id'])
            or type(vertex['vertex_loop_order']) is not int or vertex['vertex_loop_order'] != 1
            or vertex['operator_id'] != OPERATOR_ID or vertex['operator_sha256'] != OPERATOR_SHA256):
        _fail('central_operator_identity')
    loops = vertex['loops']
    if (not isinstance(loops, list) or len(loops) != 3 or not all(_label(l) for l in loops)
            or len(set(loops)) != 3 or set(loops) & (set(order) | {'v1', 'v2', 'v3', 'v4'})
            or vertex['incoming_momenta'] != ['-' + l for l in loops]):
        _fail('three_external_momenta')
    for key in ('lorentz', 'color'):
        inner = vertex[key]
        if not isinstance(inner, list) or len(inner) != 3 or not all(_label(i) for i in inner):
            _fail('central_indices')
        if len(set(inner + [p[key] for p in endpoints])) != 6:
            _fail('index_bijection')

    props = raw['propagators']
    _records(props, 3, {'loop', 'lorentz', 'color'}, 'three_propagators')
    expected_segments = []
    for i, loop in enumerate(loops):
        expected_segments.append((loop, tuple(sorted([(endpoints[i]['lorentz'], endpoints[i]['color']),
                                                     (vertex['lorentz'][i], vertex['color'][i])]))))
    actual_segments = []
    for prop in props:
        if (not _label(prop['loop']) or any(not isinstance(prop[k], list) or len(prop[k]) != 2
                or not all(_label(i) for i in prop[k]) for k in ('lorentz', 'color'))):
            _fail('propagator_indices')
        actual_segments.append((prop['loop'], tuple(sorted(zip(prop['lorentz'], prop['color'])))))
    if Counter(actual_segments) != Counter(expected_segments):
        _fail('segment_connectivity')

    phases = raw['phases']
    _records(phases, 3, {'loop', 'expression'}, 'phases')
    if any(not _label(p['loop']) for p in phases) or {p['loop'] for p in phases} != set(loops):
        _fail('phase_bijection')
    by_phase = {p['loop']: p['expression'] for p in phases}
    for i, loop in enumerate(loops):
        if not isinstance(by_phase[loop], str):
            _fail('phase_expression')
        try:
            phase = _poly(by_phase[loop], [*order, 'v1', 'v2', 'v3', 'v4'])
            residual = _closure(f'({phase})+({_position(endpoints[i])})')
        except (ValueError, IndexError, ZeroDivisionError) as error:
            raise ValueError('gamma1_phase_polynomial') from error
        if residual != '0':
            _fail('phase_endpoint_mismatch')

    constraints = raw['momentum_constraints']
    _records(constraints, 1, {'type', 'expression', 'normalization'}, 'central_delta')
    total = '+'.join(loops); constraint = constraints[0]
    if (constraint['type'] != 'delta' or not isinstance(constraint['expression'], str)
            or constraint['normalization'] != f'(2*pi)^D*delta^D({total})'):
        _fail('central_delta')
    try:
        delta_poly = _poly(constraint['expression'], loops)
    except (ValueError, IndexError, ZeroDivisionError) as error:
        raise ValueError('gamma1_delta_polynomial') from error
    if delta_poly != expand_polynomial(total):
        _fail('central_delta')

    correction = raw['vertex_correction']
    fixed = {'kind': '1pi', 'operator': 'Gamma1', 'subgraph_loop_order': 1,
             'skeleton_loop_order': 2, 'total_loop_order': 3, 'independent_external_momenta': 2,
             'external_momentum_labels': loops, 'external_leg_self_energies_included': False,
             'counterterm_included': False, 'internal_loop_measure_already_in_operator': True,
             'source_external_sign': '+1', 'gauge': 'Feynman', 'dimension': 'D=4-2*eps',
             'renormalization': 'bare', 'full_color_channels': ['CF*CA^2', 'CF*CA*TF*nf']}
    _record(correction, set(fixed) | {'operator_definition'}, 'operator_contract_fields')
    if not _same({k: correction[k] for k in fixed}, fixed):
        _fail('operator_convention')
    operator = correction['operator_definition']
    if not isinstance(operator, dict) or operator.get('id') != OPERATOR_ID:
        _fail('embedded_operator_identity')
    try:
        digest = hashlib.sha256(json.dumps(operator, ensure_ascii=False, sort_keys=True,
                                          separators=(',', ':'), allow_nan=False).encode('utf-8')).hexdigest()
    except (TypeError, ValueError) as error:
        raise ValueError('gamma1_invalid_embedded_operator') from error
    if digest != OPERATOR_SHA256:
        _fail('embedded_operator_digest')

    normalization = {'normalized_wilson_trace': 'Tr/N', 'euler_compensation_count': 0,
                     'color_projection': 'none', 'orbit_multiplicity_included': False,
                     'crossing_assembled': False}
    _record(raw['source_normalization'], normalization, 'source_normalization_fields')
    if not _same(raw['source_normalization'], normalization):
        _fail('source_normalization')
    symmetry = _symmetry(endpoints)
    expected_c4 = {k: symmetry[k] for k in ('orbit_multiplicity', 'stabilizer_order')}
    if not _same(raw['c4_symmetry'], expected_c4):
        _fail('c4_source_multiplicity')
    status = {'momentum_source': 'complete_in_stated_operator_convention',
              'montopar': 'not_implemented', 'C_J': 'not_available', 'I': 'not_solved'}
    if not _same(raw['downstream_status'], status):
        _fail('immutable_source_status')
    return {'endpoints': endpoints, 'loops': deepcopy(loops), 'order': order,
            'bounds': {var: expected_bounds[var] for var in order},
            'ordered_chains': [[p['s_var'] for p in endpoints if p['edge'] == edge] for edge in sorted(slots)]}


def parse_gamma1_orbit(raw):
    data = validate_source(raw)
    by_edge = {p['edge']: [q['s_var'] for q in data['endpoints'] if q['edge'] == p['edge']]
               for p in data['endpoints']}
    return Orbit(raw['orbit_id'], FAMILY, 0,
                 EdgeDistribution({e: len(v) for e, v in by_edge.items()}, by_edge), [],
                 [Bound(b['var'], b['low'], b['high']) for b in raw['bounds']],
                 [Vertex(**v) for v in raw['vertices']],
                 [Propagator(p['loop'], tuple(p['lorentz']), tuple(p['color'])) for p in raw['propagators']],
                 source=deepcopy(raw), dispatch='gamma1',
                 tgv_legs=[TGVLeg(p['edge'], p['s_var'], loop) for p, loop in zip(data['endpoints'], data['loops'])])


def c4_symmetry(raw):
    """Return the source multiplicity without multiplying it into raw C."""
    return _symmetry(validate_source(raw)['endpoints'])
