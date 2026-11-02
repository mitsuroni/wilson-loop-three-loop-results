"""Unapplied two-vertex operator on an independent-position scalar kernel.

The descriptor is finite: identity/outer products, four derivative slots and
six distances. Conversion differentiates that kernel before its path pullback.
No CAS dependency, executable expression evaluator or general tensor IR.
"""
from copy import deepcopy
import hashlib
import json
import re

from montopar.double_tgv import (
    ALPHAS, _add, _mul, _dot, _scale, _vadd, _point, _tangent,
)
from montopar.j_factor import expand_polynomial, differentiate_polynomial

SCHEMA = 'double-tgv-operator-1.0'
PREFIX = 'dto_'
NAMES = re.compile(r'\b[A-Za-z_][A-Za-z0-9_]*\b')
PAIRS = [(i, j) for i in range(1, 5) for j in range(i + 1, 5)]
DISTANCES = [f'd{i}{j}' for i, j in PAIRS]
GRAM = [f'g{i}{j}' for i in range(1, 5) for j in range(i, 5)]
ORDER = ['differentiate_independent_positions', 'restrict_to_paths']
FIXED = ['u1', 'u2', 'u3', 'u4', 'A0', *ALPHAS, 'eps', 'eta_reg']
DOMAIN_KEYS = ('variables', 'bounds', 'iterated_order_outer_to_inner',
               'ordered_chains', 'alpha_vars', 'alpha_measure',
               'delta_constraints', 'integration_order', 'jacobian')


def _tensor_terms():
    def term(coefficient, outer=None):
        return {'coefficient': coefficient,
                'tensor': 'identity' if outer is None else {'outer_legs': list(outer)}}
    return [
        [term('-g12'), term('2', (1, 2)), term('-1', (2, 1))],
        [term('g12'), term('1', (1, 2)), term('-2', (2, 1))],
        [term('-g34'), term('2', (3, 4)), term('-1', (4, 3))],
        [term('g34'), term('1', (3, 4)), term('-2', (4, 3))],
    ]


def _definitions():
    return [
        {'symbol': PREFIX + 'SigmaL', 'expression': 'alpha1+alpha2+alpha5'},
        {'symbol': PREFIX + 'SigmaR', 'expression': 'alpha3+alpha4+alpha5'},
        {'symbol': PREFIX + 'U', 'expression': 'dto_SigmaL*dto_SigmaR-alpha5^2'},
        {'symbol': PREFIX + 'B0', 'expression':
         'alpha1*alpha2*dto_SigmaR*d12+alpha3*alpha4*dto_SigmaL*d34+'
         'alpha5*(alpha1*alpha3*d13+alpha1*alpha4*d14+'
         'alpha2*alpha3*d23+alpha2*alpha4*d24)'},
    ]


def build_operator(orbit, raw):
    """Export the parent object; do not evaluate its derivatives here."""
    legs = orbit.double_tgv_legs
    if any(l.s_var.startswith(PREFIX) or l.s_var in set(DISTANCES + GRAM + ALPHAS +
               ['D', 'eps', 'z1', 'A0', 'v14', 'eta_reg']) for l in legs):
        raise ValueError('double_tgv_operator_reserved_variable')
    source = json.dumps(orbit.source, sort_keys=True, separators=(',', ':'),
                        ensure_ascii=False).encode('utf-8')
    block = {
        'schema_version': SCHEMA,
        'representation': 'ordered_s_parent_operator',
        'source': {'orbit_id': orbit.orbit_id,
                   'topology_sha256': hashlib.sha256(source).hexdigest(),
                   'legs': deepcopy(raw['representative']['legs'])},
        'independent_positions': [f'X{i}' for i in range(1, 5)],
        'held_fixed_under_position_derivatives': FIXED.copy(),
        'position_metric': 'D-dimensional Lorentz metric; signature (+,-,...)',
        'distance_definitions': {f'd{i}{j}': {'positions': [i, j],
                                          'expression': f'-(X{i}-X{j})^2'}
                                 for i, j in PAIRS},
        'coefficient_definitions': {f'g{i}{j}': {'tangents': [i, j],
                                              'expression': f'u{i} dot u{j}'}
                                   for i in range(1, 5) for j in range(i, 5)},
        'kernel': {
            'definitions': _definitions(), 'B0': 'dto_B0', 'U': 'dto_U',
            'argument': 'dto_B0/A0+I*eta_reg*dto_U',
            'power': '-1+3*eps',
            'alpha_weight': '(alpha1*alpha2*alpha3*alpha4*alpha5)^(-eps)*dto_U^(-1-2*eps)',
        },
        'operator': {
            'kind': 'contract_two_first_order_vertices', 'overall_coefficient': '-1',
            'left_slots': [1, 2], 'right_slots': [3, 4],
            'tensor_terms': _tensor_terms(), 'dimension': '4-2*eps',
            'derivative_convention': 'D-dimensional gradient in independent X_slot',
            'outer_product_convention': '(u tensor v)w = u*(v dot w)',
            'trace_convention': 'identity:identity=D; do not use the edge Gram rank',
        },
        'restriction': {
            'operation_order': ORDER.copy(),
            'positions': [{'slot': i + 1, 'edge': leg.edge, 's_var': leg.s_var,
                           'edge_basis_components': _point(leg),
                           'tangent_edge_basis_components': _tangent(leg)}
                          for i, leg in enumerate(legs)],
            'parent_normal_responses_retained': True,
        },
        'normalization': {
            'relative_to': 'C * J_raw_s', 'relative_C': '1',
            'outside_kernel_factor': '1/((1-3*eps)*A0)',
            'reference_scale': 'A0=-2*v14', 'color_projection': 'none',
            'euler_compensation_count': 0, 'orbit_multiplicity_included': False,
            'crossing_assembled': False,
        },
        'C_raw': deepcopy(raw['C']), 'combined_C': raw['C']['combined'],
        'color_contraction': deepcopy(raw['color_contraction']),
        'c4_symmetry': deepcopy(raw['c4_symmetry']),
        'physical_zero': raw['J_raw_s']['physical_zero'],
        'zero_policy': 'retain the parent geometric operator before color screening; scalar effective J is zero',
        'branch_prescription': {
            'eta_reg': 'positive dimensionless coordinate auxiliary regulator; kappa=A0*eta_reg',
            'finite_auxiliary_kernel': 'B0/A0+I*eta_reg*U',
            'physical_kernel': 'common boundary value as eta_reg tends to zero from above',
            'starting_region': {'A0': '>0', 'z1': '>0'},
            'independent_term_limits_or_derivative_integral_exchange_proved': False,
        },
        'equivalence': {
            'to': 'J_raw_s', 'operation': 'apply parent operator, then restrict and use the second chain rule',
            'finite_auxiliary_regulator_identity': True,
            'physical_complex_boundary_limit_verified': False,
        },
        'scalar_conversion': 'montopar.double_tgv_operator.to_scalar',
        'formal': 'not_run',
    }
    for key in DOMAIN_KEYS:
        block[key] = deepcopy(raw['J_raw_s'][key])
    return block


def _polynomial(expression, symbols, values=None):
    """Use the existing exact polynomial grammar, with an explicit scope."""
    if not isinstance(expression, str) or set(NAMES.findall(expression)) - set(symbols):
        raise ValueError('double_tgv_operator_undefined_polynomial_symbol')
    values = values or {}
    expression = NAMES.sub(lambda m: '(' + values[m[0]] + ')' if m[0] in values else m[0], expression)
    try:
        return expand_polynomial(expression)
    except (ValueError, ZeroDivisionError, IndexError) as error:
        raise ValueError('double_tgv_operator_invalid_polynomial') from error


def _parent_fields(block):
    values = {}
    for definition in block['kernel']['definitions']:
        symbol = definition['symbol']
        if not re.fullmatch(r'dto_[A-Za-z0-9_]+', symbol) or symbol in values:
            raise ValueError('double_tgv_operator_invalid_definition')
        values[symbol] = _polynomial(definition['expression'],
                                     set(ALPHAS + DISTANCES) | values.keys(), values)
    return {key: _polynomial(block['kernel'][key], values.keys(), values)
            for key in ('B0', 'U')}


def _validate_block(block):
    """Check the supported finite operator contract before any conversion."""
    if not isinstance(block, dict) or block.get('schema_version') != SCHEMA:
        raise ValueError('unsupported_double_tgv_operator_schema')
    try:
        if (block['representation'] != 'ordered_s_parent_operator'
                or block['independent_positions'] != [f'X{i}' for i in range(1, 5)]
                or block['held_fixed_under_position_derivatives'] != FIXED
                or block['restriction']['operation_order'] != ORDER
                or block['restriction']['parent_normal_responses_retained'] is not True):
            raise ValueError('double_tgv_operator_parent_scope_required')
        op = block['operator']
        if (op['kind'] != 'contract_two_first_order_vertices'
                or op['left_slots'] != [1, 2] or op['right_slots'] != [3, 4]
                or op['dimension'] != '4-2*eps'
                or op['derivative_convention'] != 'D-dimensional gradient in independent X_slot'
                or op['outer_product_convention'] != '(u tensor v)w = u*(v dot w)'
                or len(op['tensor_terms']) != 4):
            raise ValueError('unsupported_double_tgv_operator')
        _polynomial(op['overall_coefficient'], [])
        for terms in op['tensor_terms']:
            if not isinstance(terms, list) or not terms:
                raise ValueError('invalid_double_tgv_operator_tensor')
            for term in terms:
                _polynomial(term['coefficient'], GRAM)
                tensor = term['tensor']
                if tensor != 'identity' and (not isinstance(tensor, dict)
                        or set(tensor) != {'outer_legs'} or len(tensor['outer_legs']) != 2
                        or any(type(i) is not int or i not in range(1, 5) for i in tensor['outer_legs'])):
                    raise ValueError('unsupported_double_tgv_operator_tensor')
        norm = block['normalization']
        if (norm['relative_to'] != 'C * J_raw_s' or norm['relative_C'] != '1'
                or norm['outside_kernel_factor'] != '1/((1-3*eps)*A0)'
                or norm['reference_scale'] != 'A0=-2*v14'
                or norm['euler_compensation_count'] != 0 or norm['color_projection'] != 'none'
                or norm['orbit_multiplicity_included'] is not False
                or norm['crossing_assembled'] is not False
                or block['combined_C'] != block['C_raw']['combined']):
            raise ValueError('double_tgv_operator_normalization_mismatch')
        kernel = block['kernel']
        if (kernel['argument'] != 'dto_B0/A0+I*eta_reg*dto_U'
                or kernel['power'] != '-1+3*eps'
                or kernel['alpha_weight'] != '(alpha1*alpha2*alpha3*alpha4*alpha5)^(-eps)*dto_U^(-1-2*eps)'
                or block['branch_prescription']['finite_auxiliary_kernel'] != 'B0/A0+I*eta_reg*U'
                or block['equivalence']['finite_auxiliary_regulator_identity'] is not True):
            raise ValueError('double_tgv_operator_kernel_or_prescription_mismatch')
        if block['distance_definitions'] != {
                f'd{i}{j}': {'positions': [i, j], 'expression': f'-(X{i}-X{j})^2'} for i, j in PAIRS}:
            raise ValueError('double_tgv_operator_distance_scope_mismatch')
        if block['coefficient_definitions'] != {
                f'g{i}{j}': {'tangents': [i, j], 'expression': f'u{i} dot u{j}'}
                for i in range(1, 5) for j in range(i, 5)}:
            raise ValueError('double_tgv_operator_tangent_scope_mismatch')
        positions = block['restriction']['positions']
        if len(positions) != 4 or [p['slot'] for p in positions] != [1, 2, 3, 4]:
            raise ValueError('double_tgv_operator_position_slots_mismatch')
        names = [p['s_var'] for p in positions]
        if len(set(names)) != 4 or set(names) != set(block['variables']) or set(names) != set(block['bounds']):
            raise ValueError('double_tgv_operator_domain_mismatch')
        if (sorted(block['iterated_order_outer_to_inner']) != sorted(names)
                or sorted(v for chain in block['ordered_chains'] for v in chain) != sorted(names)):
            raise ValueError('double_tgv_operator_order_mismatch')
        for p in positions:
            edge = p['edge']
            if type(edge) is not int or edge not in (1, 2, 3, 4):
                raise ValueError('double_tgv_operator_unknown_edge')
            point = ['1' if k < edge - 1 else p['s_var'] if k == edge - 1 else '0' for k in range(4)]
            tangent = ['1' if k == edge - 1 else '0' for k in range(4)]
            if p['edge_basis_components'] != point or p['tangent_edge_basis_components'] != tangent:
                raise ValueError('double_tgv_operator_pullback_mismatch')
        if block['alpha_vars'] != ALPHAS or block['jacobian'] != '1':
            raise ValueError('double_tgv_operator_requires_original_simplex_and_paths')
        from montopar.double_tgv_output import alpha_simplex5, CONSTRAINT
        if (block['alpha_measure'] != alpha_simplex5()
                or block['delta_constraints'] != {CONSTRAINT: '1'}
                or block['integration_order'] != 'path_parameters_outer_alpha_simplex_inner'):
            raise ValueError('double_tgv_operator_simplex_mismatch')
        return _parent_fields(block)
    except (KeyError, TypeError, AttributeError) as error:
        raise ValueError('incomplete_double_tgv_operator') from error


def scalar_fields(block):
    """Apply the exported operator, then pull back; retain both chain-rule terms.

    All vector components below represent external contour vectors only.
    Their Gram contractions do not supply the separate symbolic D trace.
    """
    parent = _validate_block(block)
    positions = block['restriction']['positions']
    points = [p['edge_basis_components'] for p in positions]
    tangents = [p['tangent_edge_basis_components'] for p in positions]
    gram = {f'g{i}{j}': _dot(tangents[i - 1], tangents[j - 1])
            for i in range(1, 5) for j in range(i, 5)}

    def substitute(expression, values):
        return _polynomial(expression, set(values) | {'eps'}, values)

    # Differentiate the unrestricted invariant polynomial before using points(s).
    coefficients = {}
    remainder = parent['B0']
    for i, j in PAIRS:
        distance = f'd{i}{j}'
        coefficient = differentiate_polynomial(parent['B0'], distance)
        if set(NAMES.findall(coefficient)) - set(ALPHAS):
            raise ValueError('double_tgv_operator_requires_distance_linear_parent_kernel')
        coefficients[i - 1, j - 1] = coefficient
        remainder = _add(remainder, _mul('-1', coefficient, distance))
    if remainder != '0' or set(NAMES.findall(parent['U'])) - set(ALPHAS):
        raise ValueError('double_tgv_operator_invalid_parent_kernel')

    def coefficient(i, j):
        return coefficients[tuple(sorted((i, j)))]

    gradients = [_vadd(*(_scale(_mul('-2', coefficient(i, j)),
                          _vadd(points[i], _scale('-1', points[j])))
                         for j in range(4) if j != i)) for i in range(4)]

    def apply(terms, vector):
        values = []
        for term in terms:
            tensor = term['tensor']
            value = vector if tensor == 'identity' else _scale(
                _dot(tangents[tensor['outer_legs'][1] - 1], vector),
                tangents[tensor['outer_legs'][0] - 1])
            values.append(_scale(substitute(term['coefficient'], gram), value))
        return _vadd(*values)

    tensors = block['operator']['tensor_terms']
    left = _vadd(apply(tensors[0], gradients[0]), apply(tensors[1], gradients[1]))
    right = _vadd(apply(tensors[2], gradients[2]), apply(tensors[3], gradients[3]))

    def trace_pair(terms1, terms2):
        terms = []
        for a in terms1:
            for b in terms2:
                A, B = a['tensor'], b['tensor']
                if A == 'identity' and B == 'identity':
                    value = block['operator']['dimension']
                elif A == 'identity' or B == 'identity':
                    pair = (B if A == 'identity' else A)['outer_legs']
                    value = _dot(tangents[pair[0] - 1], tangents[pair[1] - 1])
                else:
                    p, q = A['outer_legs']
                    r, t = B['outer_legs']
                    value = _mul(_dot(tangents[p - 1], tangents[r - 1]),
                                 _dot(tangents[q - 1], tangents[t - 1]))
                terms.append(_mul(substitute(a['coefficient'], gram),
                                  substitute(b['coefficient'], gram), value))
        return _add(*terms)

    overall = block['operator']['overall_coefficient']
    t1 = _mul(overall, _dot(left, right))
    t2 = _mul(overall, _add(*(_mul('2', coefficient(i, j), trace_pair(tensors[i], tensors[j]))
                              for i in (0, 1) for j in (2, 3))))
    # B0/vref=-sum c_ij dot(Xi-Xj,Xi-Xj); A0/vref=-2.
    bhat = _add(*(_mul('1/2', c,
                       _dot(_vadd(points[i], _scale('-1', points[j])),
                            _vadd(points[i], _scale('-1', points[j]))))
                 for (i, j), c in coefficients.items()))
    return {'B5': bhat, 'U': parent['U'], 'T1_hat': t1, 'T2_hat': t2,
            'pnum1': _mul('-1/8', '2-3*eps', t1), 'pnum2': _mul('-1/4', t2),
            'K_plus': trace_pair(tensors[0], tensors[2]),
            'K_minus': trace_pair(tensors[0], tensors[3])}


def to_scalar(block):
    """Convert J_operator to the existing ordered scalar block, with the same C.

    This applies no IBP and exchanges no divergent integrals or limits.
    It computes the geometric kernel even when color later sets effective J=0.
    """
    fields = scalar_fields(block)
    from montopar.double_tgv_output import _block
    out = _block(fields, block['variables'], deepcopy(block['bounds']),
                 block['iterated_order_outer_to_inner'], '1', block['physical_zero'])
    out.update({'representation': 'ordered_s',
                'ordered_chains': deepcopy(block['ordered_chains'])})
    return out
