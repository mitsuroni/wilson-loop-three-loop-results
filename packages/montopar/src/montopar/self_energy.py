"""Bare one-loop 1PI gluon self energy, with the complete transverse tensor.

No CAS runtime dependency. Formula/phase contract: docs/self-energy-pi1.md.
Only one insertion on 2 or 4 Wilson endpoints is supported in this release.
"""
from montopar.sign_rules import real_phase, annotate_result, validate_contract
from montopar.sign_rules import real_phase, annotate_result, validate_contract
from collections import Counter, defaultdict
from copy import deepcopy
import re

from montopar.geometry import dot_product_factor
from montopar.j_factor import (
    _Parser, _tokenize, expand_polynomial, substitute_s_to_x,
)
from montopar.models import Bound, EdgeDistribution, Orbit, Propagator, Vertex

SCHEMA = '1.4'
FAMILIES = {'2e-Pi1': (2, 2), '4e-Pi1': (4, 3)}
SE_COLOR = '((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)'


def _poly(expression, allowed):
    """Use the existing exact polynomial parser, requiring full consumption."""
    tokens = _tokenize(expression)
    names = {t for t in tokens if re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', t)}
    if not names <= set(allowed):
        raise ValueError('self_energy_unknown_polynomial_symbol')
    parser = _Parser(tokens)
    parser.parse_expr()
    if parser.i != len(tokens):
        raise ValueError('self_energy_trailing_polynomial_tokens')
    return expand_polynomial(expression)


def _position(endpoint):
    return '+'.join([f'v{k}' for k in range(1, endpoint['edge'])]
                    + [f"{endpoint['s_var']}*v{endpoint['edge']}"])


def _closure(expression):
    return expand_polynomial(substitute_s_to_x(expression, {'v4': '(-v1-v2-v3)'}))


def validate_source(raw):
    """Bind every endpoint, delta bridge, momentum, phase and ordered bound."""
    validate_contract(raw)
    required = {'orbit_id', 'emission_type', 'sector', 'nloops', 'vertices', 'propagators',
                'internal_vertices', 'phases', 'bounds'}
    if not required <= set(raw) or not set(raw) <= required | {'g_power', 'iterated_order_outer_to_inner', 'physical_convention'}:
        raise ValueError('self_energy_unknown_or_missing_source_fields')
    family = raw.get('emission_type')
    if family not in FAMILIES:
        raise ValueError('unsupported_self_energy_family: only 2e-Pi1 and 4e-Pi1')
    size, loops = FAMILIES[family]
    if raw.get('nloops') != loops or raw.get('sector') != 'dg':
        raise ValueError('self_energy_loop_order_or_sector_mismatch')
    if 'g_power' in raw and raw['g_power'] != 2*loops:
        raise ValueError('self_energy_coupling_order_mismatch')
    ivs = raw.get('internal_vertices', [])
    if len(ivs) != 1 or set(ivs[0]) != {'type', 'propagator', 'se_type'}:
        raise ValueError('self_energy_requires_one_unambiguous_internal_operator')
    iv = ivs[0]
    if iv['type'] != 'se' or iv['se_type'] != 'vacuum-pol':
        raise ValueError('self_energy_requires_one_loop_1pi_operator')
    vertices = raw.get('vertices', [])
    if len(vertices) != size or any(type(v.get('line')) is not int or v['line'] not in range(1, 5) for v in vertices):
        raise ValueError('self_energy_endpoint_count_or_edge_mismatch')
    if [v['line'] for v in vertices] != sorted(v['line'] for v in vertices):
        raise ValueError('self_energy_vertices_require_contour_order')
    if any(set(v) != {'line', 'lorentz', 'color'} for v in vertices):
        raise ValueError('self_energy_ambiguous_external_vertex')
    for field in ('lorentz', 'color'):
        if len({v[field] for v in vertices}) != size:
            raise ValueError('self_energy_duplicate_external_index')
    slots = Counter()
    endpoints = {}
    for rank, v in enumerate(vertices):
        slots[v['line']] += 1
        endpoints[v['lorentz']] = dict(v, edge=v['line'], slot=slots[v['line']],
                                     s_var=f"s{v['line']}{slots[v['line']]}", rank=rank)
    bounds = raw.get('bounds', [])
    expected_bounds = {p['s_var']: ['0', f"s{p['edge']}{p['slot']+1}" if p['slot'] < slots[p['edge']] else '1']
                       for p in endpoints.values()}
    if (len(bounds) != size or len({b['var'] for b in bounds}) != size
            or {b['var']: [str(b['low']), str(b['high'])] for b in bounds} != expected_bounds):
        raise ValueError('self_energy_endpoint_ordered_domain_mismatch')
    order = [f's{edge}{slot}' for edge in sorted(slots) for slot in range(slots[edge], 0, -1)]
    if 'iterated_order_outer_to_inner' in raw and raw['iterated_order_outer_to_inner'] != order:
        raise ValueError('self_energy_explicit_order_mismatch')
    props = raw.get('propagators', [])
    if len(props) != size//2+1:
        raise ValueError('self_energy_propagator_count_mismatch')
    grouped = defaultdict(list)
    visited, internal_mu, internal_color = [], [], []
    external_colors = {v['color'] for v in vertices}
    for prop in props:
        if (set(prop) != {'loop', 'lorentz', 'color'} or len(prop['lorentz']) != 2
                or len(prop['color']) != 2 or not isinstance(prop['loop'], str)):
            raise ValueError('self_energy_malformed_propagator')
        ext = []
        for mu, color in zip(prop['lorentz'], prop['color']):
            if mu in endpoints:
                endpoint = endpoints[mu]
                if color != endpoint['color']:
                    raise ValueError('self_energy_external_color_mismatch')
                visited.append(mu); ext.append(endpoint)
            else:
                internal_mu.append(mu); internal_color.append(color)
        if (prop['loop'] == iv['propagator'] and len(ext) != 1
                or prop['loop'] != iv['propagator'] and len(ext) != 2):
            raise ValueError('self_energy_internal_connectivity_mismatch')
        grouped[prop['loop']].append((prop, ext))
    if (Counter(visited) != Counter(endpoints.keys()) or len(internal_mu) != 2
            or len(set(internal_mu)) != 2 or len(set(internal_color)) != 2
            or set(internal_color) & external_colors):
        raise ValueError('self_energy_endpoint_or_internal_index_bijection_failed')
    se_loop = iv['propagator']
    if len(grouped[se_loop]) != 2 or any(len(v) != 1 for k, v in grouped.items() if k != se_loop):
        raise ValueError('self_energy_momentum_bridge_mismatch')
    phases = raw.get('phases', [])
    if (len(phases) != len(grouped) or len({p['loop'] for p in phases}) != len(phases)
            or {p['loop'] for p in phases} != set(grouped)):
        raise ValueError('self_energy_phase_momentum_bijection_failed')
    phase_by = {p['loop']: p['expression'] for p in phases}
    lines = []
    for loop, pieces in grouped.items():
        ends = sorted([p for _, ext in pieces for p in ext], key=lambda p: p['rank'])
        if len(ends) != 2 or ends[0]['edge'] == ends[1]['edge']:
            raise ValueError('self_energy_requires_distinct_edge_pair')
        phase = _poly(phase_by[loop], [*expected_bounds, 'v1', 'v2', 'v3', 'v4'])
        expected = f'({_position(ends[1])})-({_position(ends[0])})'
        if _closure(f'({phase})-({expected})') == '0':
            orientation = 1
        elif _closure(f'({phase})+({expected})') == '0':
            orientation = -1
        else:
            raise ValueError('self_energy_phase_endpoint_geometry_mismatch')
        lines.append({'kind': 'se' if loop == se_loop else 'bare', 'loop': loop,
                      'endpoints': ends, 'propagators': [deepcopy(p) for p, _ in pieces],
                      'phase_Phi': phase_by[loop], 'phase_orientation': orientation})
    lines.sort(key=lambda p: (p['kind'] != 'se', p['loop']))
    return {'lines': lines, 'bounds': expected_bounds, 'order': order,
            'internal_delta': {'lorentz': sorted(internal_mu), 'color': sorted(internal_color)},
            'endpoints': list(endpoints.values())}


def parse_self_energy_orbit(raw):
    validate_source(raw)
    bounds = [Bound(b['var'], str(b['low']), str(b['high'])) for b in raw['bounds']]
    by_edge = defaultdict(list)
    for b in bounds:
        by_edge[int(b.var[1])].append(b.var)
    return Orbit(orbit_id=raw['orbit_id'], emission_type=raw['emission_type'],
                 n_gluons=len(raw['vertices'])//2-1, pairings=[], bounds=bounds,
                 edge_distribution=EdgeDistribution({k: len(v) for k, v in by_edge.items()}, dict(by_edge)),
                 vertices=[Vertex(**v) for v in raw['vertices']],
                 propagators=[Propagator(p['loop'], tuple(p['lorentz']), tuple(p['color'])) for p in raw['propagators']],
                 source=deepcopy(raw), dispatch='self_energy')


def color_data(data):
    pairs = sorted(tuple(p['rank'] for p in line['endpoints']) for line in data['lines'])
    crossed = pairs == [(0, 2), (1, 3)]
    trace = 'CF' if len(pairs) == 1 else 'CF*(CF-CA/2)' if crossed else 'CF^2'
    return {'wilson_trace': trace, 'self_energy_color': SE_COLOR,
            'full_color': f'({trace})*({SE_COLOR})', 'endpoint_pairs_in_contour_order': pairs,
            'method': 'normalized fundamental trace; Ta*Tb*Ta=(CF-CA/2)*Tb',
            'TF': '1/2', 'nf': 'number of massless quark flavours',
            'projection': 'none'}


def c4_symmetry(raw, data=None):
    data = validate_source(raw) if data is None else data
    def fingerprint(rotation):
        return tuple(sorted((line['kind'], tuple(sorted(((p['edge']-1+rotation)%4+1, p['slot'])
                                                     for p in line['endpoints']))) for line in data['lines']))
    stab = [r for r in range(4) if fingerprint(r) == fingerprint(0)]
    return {'orbit_multiplicity': 4//len(stab), 'stabilizer_order': len(stab),
            'stabilizer_rotations': stab,
            'definition': 'C4 acting on labelled SE/bare endpoint pairs and ordered slots; no mirror quotient'}


def geometry_for_line(line):
    a, b = line['endpoints']
    # X = x_a-x_b. The tensor contains X twice, so reversing loop routing is harmless.
    vector = [expand_polynomial(f"({a['s_var'] if k == a['edge'] else '1' if k < a['edge'] else '0'})"
                                f"-({b['s_var'] if k == b['edge'] else '1' if k < b['edge'] else '0'})")
              for k in range(1, 5)]
    def dot(edge):
        return expand_polynomial('+'.join(f'({v})*({dot_product_factor(edge,k)})' for k, v in enumerate(vector, 1)))
    f = expand_polynomial('+'.join(f'({vector[i-1]})*({vector[j-1]})*({dot_product_factor(i,j)})'
                                   for i in range(1, 5) for j in range(i+1, 5)))
    return {'loop': line['loop'], 'kind': line['kind'], 'edges': [a['edge'], b['edge']],
            's_vars': [a['s_var'], b['s_var']], 'X_coefficients': vector,
            'f': f, 'd': dot_product_factor(a['edge'], b['edge']),
            'u': dot(a['edge']), 'w': dot(b['edge'])}


def build_result(orbit):
    raw = orbit.source
    data = validate_source(raw)
    parsed = parse_self_energy_orbit(raw)
    if orbit != parsed:
        raise ValueError('self_energy_orbit_metadata_disagrees_with_source')
    loop_order = raw['nloops']
    color = color_data(data)
    geometries = [geometry_for_line(line) for line in data['lines']]
    se = geometries[0]
    f, d, u, w = (se[k] for k in ('f', 'd', 'u', 'w'))
    metric_num = expand_polynomial(f'2*(1+eps)*({d})*({f})')
    derivative_num = expand_polynomial(f'-({d})*({f})+(1-2*eps)*({u})*({w})')
    full_num = expand_polynomial(f'({metric_num})+({derivative_num})')
    se_factor = f'({f})^(-2+2*eps)/(2*(1+eps))'
    bare_factor = '*'.join(f"({p['d']})*({p['f']})^(-1+eps)" for p in geometries[1:]) or '1'
    component = lambda numerator: f'({numerator})*({se_factor})*({bare_factor})'
    j = {'representation': 'ordered_s', 'variables': data['order'], 'bounds': data['bounds'],
         'iterated_order_outer_to_inner': data['order'], 'source_bounds_array': deepcopy(raw['bounds']),
         'order_basis': 'ascending contour edges; descending slots; checked against source integrals.md at publication',
         'jacobian': '1', 'integrand': component(full_num),
         'tensor_decomposition': {'metric': component(metric_num), 'derivative': component(derivative_num),
                                  'identity': 'integrand = metric + derivative; common C'},
         'se_numerator': full_num, 'metric_numerator': metric_num, 'derivative_numerator': derivative_num,
         'geometry_polynomials': geometries,
         'power_prescription': 'each f power is the boundary value (f+I0)^power in A_ref>0; continue jointly with C'}
    sign = real_phase(len(raw["vertices"]), 3*2, 1, 2, 1, 2)
    components = {'wilson_color': color['wilson_trace'], 'self_energy_color': SE_COLOR,
                  'coupling': f'(alpha_s/Pi)^{loop_order}',
                  'gamma_fourier': f'(Pi^eps*Gamma[1-eps])^{loop_order}',
                  'rational': f'{sign}/({2**(loop_order+2)}*eps*(1-2*eps))',
                  'scale': f'(-2*v14+I0)^({loop_order}*eps)'}
    c = {'components': components, 'factor_order': list(components), 'component_semantics': 'multiplicative',
         'combined': '*'.join(f'({v})' for v in components.values()),
         'phase_ledger': {'source_external_sign': '+1', 'wilson_vertices': f'I^{len(raw["vertices"])}',
                          'se_two_propagators': '(-I)^2', 'Pi1_scalar_phase': 'I',
                          'se_metric_k_squared_conversion': '-1', 'se_Fourier_phase': 'I',
                          'additional_bare_coordinate_propagator': '-1' if loop_order == 3 else '1',
                          'coordinate_product_before_extracting_vref': '1',
                          'vref_over_Aref_each_exchange': '-1/2',
                          'overall_dimensionless_sign': str(sign)}}
    result={'schema_version': SCHEMA, 'orbit_id': raw['orbit_id'], 'emission_type': raw['emission_type'],
            'dispatch': 'self_energy', 'loop_order': loop_order,
            'geometry': {'n_lines': 4, 'z_definition': 'z1=v12/v14', 'reference_scale': 'v14',
                         'vij_definition': 'v_i dot v_j', 'A_ref': '-2*v14', 'lightlike_edges': True,
                         'line_invariants': 'f=X^2/(2*v14); d=vi.vj/v14; u=vi.X/v14; w=vj.X/v14'},
            'representative': {'scheme': 'single_source_representative', 'rotation': 0,
                               'source_topology': deepcopy(raw), 'canonical_topology': deepcopy(raw),
                               'effective_lines': data['lines'], 'internal_delta': data['internal_delta']},
            'operator': {'name': 'Pi1', 'subgraph_loop_order': 1, 'insertion_count': 1, 'kind': '1pi',
                         'renormalization': 'bare', 'dimension': '4-2*eps', 'gauge': 'Feynman',
                         'source': 'hep-ph/9605348 Appendix C (C.1),(C.3); xi_D=0',
                         'scalar': 'I*g^2/(2*(4*Pi)^(2-eps))*Gamma[eps]*Gamma[1-eps]^2/Gamma[2-2*eps]*C_SE*(-l^2-I0)^(-eps)',
                         'C_SE': SE_COLOR, 'tensor': 'delta_ab*(l^2*g_mu_nu-l_mu*l_nu)', 'epsilon_precision': 'exact'},
            'color_contraction': color, 'c4_symmetry': c4_symmetry(raw, data),
            'normalization': {'scheme': 'montopar-raw-CJ', 'loop_order': loop_order,
                              'color_projection': 'none', 'euler_compensation_count': 0,
                              'orbit_multiplicity_included': False, 'crossing_assembled': False,
                              'normalized_wilson_trace': 'Tr/N', 'extra_factorial': 'none'},
            'C': c, 'J_raw': j, 'J_raw_s': deepcopy(j),
            'branch_prescription': {'dimensionful_distance': '-X^2+I0', 'momentum': '-l^2-I0',
                                    'starting_region': {'z1': '>0', 'v14': '<0'},
                                    'power_convention': 'common Feynman boundary value; joint analytic continuation',
                                    'complex_continuation_verified': False},
            'translation_status': 'translated', 'evaluation_status': 'not_evaluated', 'formal': 'not_run'}
    return annotate_result(result, raw)


def validate_result(result):
    if result.get('schema_version') != SCHEMA or result.get('dispatch') != 'self_energy':
        raise ValueError('regenerate_self_energy_result')
    expected = build_result(parse_self_energy_orbit(result['representative']['source_topology']))
    # Normalize tuple/list JSON round trips without permitting unvalidated extra fields.
    import json
    if json.loads(json.dumps(expected)) != result:
        raise ValueError('self_energy_result_does_not_reconstruct_from_source')
