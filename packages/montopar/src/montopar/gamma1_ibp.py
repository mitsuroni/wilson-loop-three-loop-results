"""Contour IBP for the full one-loop three-gluon vertex.

The differential-operator decomposition is polynomial and takes place before
restriction to the contour.  Each exported face is an ordinary scalar kernel
on the inherited ordered domain.  All terms retain one common boundary value.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import json
import re

from montopar.double_tgv_ibp import ordered_faces
from montopar.j_factor import (
    _parse_polynomial, _format_poly, _poly_add, _poly_mul, _poly_pow,
    differentiate_polynomial,
)
from montopar.gamma1_kernel import operators, apply_operator, apply_operator_compact


def _add(*expressions):
    result = {}
    for expression in expressions:
        result = _poly_add(result, _parse_polynomial(expression))
    return _format_poly(result)


def _mul(*expressions):
    result = {frozenset(): 1}
    for expression in expressions:
        result = _poly_mul(result, _parse_polynomial(expression))
    return _format_poly(result)


def polynomial_substitute(expression, substitutions):
    """Simultaneous exact substitution (including compound replacements)."""
    parsed = {name: _parse_polynomial(value) for name, value in substitutions.items()}
    result = {}
    for monomial, coefficient in _parse_polynomial(expression).items():
        term = {frozenset(): coefficient}
        for name, exponent in sorted(monomial):
            replacement = parsed.get(name, {frozenset({(name, 1)}): 1})
            term = _poly_mul(term, _poly_pow(replacement, exponent))
        result = _poly_add(result, term)
    return _format_poly(result)


def _tangent_basis(data, geom):
    edges = [endpoint['edge'] for endpoint in data['endpoints']]
    if len(edges) != 3 or len(set(edges)) == 1:
        raise ValueError('gamma1_ibp_requires_three_nonscaleless_legs')
    mapping = {'r11': 'd1', 'r12': 'e1', 'r21': 'e2', 'r22': 'd2',
               'r31': 'e3', 'r32': '-d3-e3',
               'h12': geom['G'][0][1], 'h23': geom['G'][1][2],
               'h31': geom['G'][2][0]}
    if edges[1] == edges[2]:
        mapping.update(r21='-d2-d3', r22='d2', r31='-d2-d3', r32='d2')
        directions = [1, 2]
    elif edges[0] == edges[1]:
        mapping.update(r11='d1', r12='d2', r21='d1', r22='d2')
        directions = [0, 1]
    elif edges[0] == edges[2]:
        # The validated canonical source orders all external legs by edge.
        # Do not silently use a different cyclic convention for an unvalidated
        # source; its color sign and path slots would need a separate bridge.
        raise ValueError('gamma1_ibp_noncanonical_repeated_first_third_leg')
    else:
        directions = [0, 1, 2]
    return mapping, directions


def decompose_operator(operator, substitutions, directions):
    """Split O = sum(delta_i O_i) + R without dividing by geometry.

    Each monomial is assigned to its first nonzero delta in the specified
    order.  This is exactly the successive zero-substitution construction in
    section 8 of the parameterization theory.
    """
    transformed = _parse_polynomial(polynomial_substitute(operator, substitutions))
    flux = [{}, {}, {}]
    remainder = {}
    for monomial, coefficient in transformed.items():
        powers = dict(monomial)
        selected = next((i for i in directions if powers.get(f'd{i+1}', 0)), None)
        if selected is None:
            remainder[monomial] = coefficient
            continue
        variable = f'd{selected+1}'
        powers[variable] -= 1
        if not powers[variable]:
            del powers[variable]
        reduced = frozenset(powers.items())
        flux[selected][reduced] = flux[selected].get(reduced, 0) + coefficient
    reconstructed = dict(remainder)
    for i, quotient in enumerate(flux):
        reconstructed = _poly_add(reconstructed, _poly_mul(
            {frozenset({(f'd{i+1}', 1)}): 1}, quotient))
    if _poly_add(reconstructed, transformed, -1):
        raise ValueError('gamma1_ibp_operator_decomposition_failed')
    basis_flux = [_format_poly(value) for value in flux]
    basis_remainder = _format_poly(remainder)
    inverse = {'d1': 'r11', 'd2': 'r22', 'd3': '-r31-r32',
               'e1': 'r12', 'e2': 'r21', 'e3': 'r31'}
    return {'flux_operators': [polynomial_substitute(value, inverse) for value in basis_flux],
            'remainder_operator': polynomial_substitute(basis_remainder, inverse),
            'tangent_basis_flux_operators': basis_flux,
            'tangent_basis_remainder_operator': basis_remainder,
            'operator_reconstruction': 'exact_polynomial_identity'}


def _clean_map(mapping):
    return {int(shift): value for shift, value in mapping.items() if value != '0'}


def reconstruct_coefficients(raw, evaluated_flux, evaluated_remainder, beta, data, geom):
    """Verify the divergence identity at each distinct power of q + I0.

    I0 stays constant and is never set to zero.  The checks do not replace q
    with its contour polynomial to cancel neighboring powers, so they also
    prove the identity for a finite common I0.
    """
    reconstructed = _clean_map(evaluated_remainder)
    for endpoint, coefficients in zip(data['endpoints'], evaluated_flux):
        variable = endpoint['s_var']
        dq = differentiate_polynomial(geom['q0'], variable)
        for shift, coefficient in coefficients.items():
            shift = int(shift)
            derivative = differentiate_polynomial(coefficient, variable)
            reconstructed[shift] = _add(reconstructed.get(shift, '0'), derivative)
            raised = _mul(coefficient, _add(beta, str(-shift)), dq)
            reconstructed[shift+1] = _add(reconstructed.get(shift+1, '0'), raised)
    for shift in set(raw) | set(reconstructed):
        difference = _add(raw.get(shift, '0'), _mul('-1', reconstructed.get(shift, '0')))
        if difference != '0':
            raise ValueError(f'gamma1_ibp_coefficient_reconstruction_failed:shift={shift}')
    return 'exact_at_each_q_power_with_finite_common_I0'


def _compute_flux_data(data, geom):
    substitutions, directions = _tangent_basis(data, geom)
    result = {}
    for channel, channel_operators in operators().items():
        result[channel] = {}
        for label, operator, beta in [
            ('N0', channel_operators['N0'], '3*eps'),
            ('L0', channel_operators['L0'], '-1+3*eps'),
            *[(f'M{i+1}', value, '-1+3*eps') for i, value in enumerate(channel_operators['M'])],
        ]:
            entry = decompose_operator(operator, substitutions, directions)
            entry['primitive_power'] = beta
            entry['evaluated_flux'] = [
                _clean_map(apply_operator(value, beta, geom))
                for value in entry['flux_operators']]
            entry['evaluated_remainder'] = _clean_map(apply_operator(entry['remainder_operator'], beta, geom))
            if len(directions) == 2 and entry['remainder_operator'] != '0':
                raise ValueError('gamma1_ibp_repeated_tangent_has_remainder')
            raw = _clean_map(apply_operator(operator, beta, geom))
            entry['coefficient_reconstruction'] = reconstruct_coefficients(
                raw, entry['evaluated_flux'], entry['evaluated_remainder'], beta, data, geom)
            result[channel][label] = entry
    result['decomposition'] = {'directions': [f'd{i+1}' for i in directions],
                               'tangent_substitutions': substitutions,
                               'derivatives_before_restriction': True}
    return result


@lru_cache(maxsize=16)
def _cached_flux_data(data_key, geometry_key):
    return _compute_flux_data(json.loads(data_key), json.loads(geometry_key))


def flux_data(data, geom):
    """Return independently owned scalar fluxes, keyed by complete inputs.

    Result validation may replay the same immutable source many times. Cache
    only completed exact checks, including the complete geometry in the key;
    a changed q or contraction therefore cannot inherit an earlier success.
    """
    data_key = json.dumps(data, sort_keys=True, separators=(',', ':'))
    geometry_key = json.dumps(geom, sort_keys=True, separators=(',', ':'))
    return deepcopy(_cached_flux_data(data_key, geometry_key))


def scalar_integrand(weight, coefficients, beta, q0):
    """A complete evaluable scalar expression, with no operator placeholders."""
    terms = [f'({value})*(({q0})+I0)^({beta}-{int(shift)})'
             for shift, value in sorted(coefficients.items()) if value != '0']
    return f'({weight})*(' + '+'.join(terms) + ')' if terms else '0'


def _inline_substitute(expression, substitutions):
    """Whole-token, parenthesized, simultaneous substitution without expansion."""
    changes = {key: value for key, value in substitutions.items() if key != value}
    return re.sub(r'\b[A-Za-z_][A-Za-z0-9_]*\b',
                  lambda match: '('+changes[match[0]]+')' if match[0] in changes else match[0],
                  expression)


def _compact_flux(flux, geom):
    """Retain compact scalar polynomials after exact expansion has checked zeros."""
    output = deepcopy(flux)
    for channel in ('g', 'q'):
        for label, entry in flux[channel].items():
            target = output[channel][label]
            beta = entry['primitive_power']
            compact = [apply_operator_compact(operator, beta, geom)
                       for operator in entry['flux_operators']]
            target['evaluated_flux'] = [
                {shift: compact[i][shift] for shift in expanded}
                for i, expanded in enumerate(entry['evaluated_flux'])]
            remainder = apply_operator_compact(entry['remainder_operator'], beta, geom)
            target['evaluated_remainder'] = {shift: remainder[shift]
                                             for shift in entry['evaluated_remainder']}
            target['scalar_polynomial_representation'] = 'fully_inline_factored'
    return output


def build_ibp(data, geom, channels_C):
    """Construct every nonzero weighted face and recoil term separately."""
    from montopar.gamma1_blocks import block_specs, kernel_terms

    flux = flux_data(data, geom)
    compact_flux = _compact_flux(flux, geom)
    terms = []
    variables = [endpoint['s_var'] for endpoint in data['endpoints']]
    faces = ordered_faces(data['ordered_chains'])
    for block in block_specs(data):
        channel, label = block['channel'], block['block_id']
        entry = flux[channel][label]
        compact_entry = compact_flux[channel][label]
        beta = block['primitive_power']
        if _add(beta, _mul('-1', entry['primitive_power'])) != '0':
            raise ValueError('gamma1_ibp_primitive_power_mismatch')
        parameters = [name for name in block['iterated_order_outer_to_inner'] if name not in variables]
        for face in faces:
            coefficients = {}
            for i, variable in enumerate(variables):
                for shift, coefficient in entry['evaluated_flux'][i].items():
                    coefficients[shift] = _add(coefficients.get(shift, '0'), _mul(face['normal'][variable], coefficient))
            coefficients = {shift: polynomial_substitute(value, face['parent_to_face'])
                            for shift, value in coefficients.items()}
            coefficients = _clean_map(coefficients)
            if not coefficients:
                continue
            # The fully expanded parent flux determines whether a coefficient
            # vanishes. Emit the equivalent factored expression only for those
            # surviving powers; textual nonzero tests are never used here.
            compact_coefficients = {}
            for shift in coefficients:
                components = [f'({face["normal"][variable]})*({compact_entry["evaluated_flux"][i][shift]})'
                              for i, variable in enumerate(variables)
                              if face['normal'][variable] != '0'
                              and shift in compact_entry['evaluated_flux'][i]]
                compact_coefficients[shift] = _inline_substitute('+'.join(components), face['parent_to_face'])
            coefficients = compact_coefficients
            remaining = [name for name in data['order'] if name != face['removed']]
            bounds = {name: ['0', chain[i+1] if i+1 < len(chain) else '1']
                      for chain in face['ordered_chains'] for i, name in enumerate(chain)}
            order = remaining + parameters
            q0 = polynomial_substitute(geom['q0'], face['parent_to_face'])
            terms.append({
                'term_id': f'ibp_{len(terms)+1:03d}', 'kind': face['kind'],
                'channel': channel, 'block_id': label, 'origin': 'full_vertex_operator_flux',
                'integration_block_kind': block['kind'],
                'combined_C': channels_C[channel], 'face': deepcopy(face),
                'variables': order, 'bounds': {**{name: bounds[name] for name in remaining},
                                             **{name: block['bounds'][name] for name in parameters}},
                'iterated_order_outer_to_inner': order, 'jacobian': '1',
                'q0': q0, 'primitive_power': beta, 'kernel_coefficients': coefficients,
                'kernel_terms': kernel_terms(coefficients, beta),
                'weight': block['weight'],
                'integrand': scalar_integrand(block['weight'], coefficients, beta, q0),
                'regulator_group': 'one_parent_all_channels_and_blocks_common_boundary_value',
                'derivatives_before_restriction': True,
            })
        coefficients = compact_entry['evaluated_remainder']
        if coefficients:
            order = list(block['iterated_order_outer_to_inner'])
            terms.append({
                'term_id': f'ibp_{len(terms)+1:03d}', 'kind': 'recoil',
                'channel': channel, 'block_id': label, 'origin': 'full_vertex_operator_remainder',
                'integration_block_kind': block['kind'],
                'combined_C': channels_C[channel], 'variables': order,
                'bounds': deepcopy(block['bounds']), 'iterated_order_outer_to_inner': order,
                'jacobian': '1', 'q0': geom['q0'], 'primitive_power': beta,
                'kernel_coefficients': deepcopy(coefficients), 'weight': block['weight'],
                'kernel_terms': kernel_terms(coefficients, beta),
                'integrand': scalar_integrand(block['weight'], coefficients, beta, geom['q0']),
                'regulator_group': 'one_parent_all_channels_and_blocks_common_boundary_value',
                'derivatives_before_restriction': True,
            })
    return {
        'schema_version': 'gamma1-ibp-1.0',
        'representation': 'ordered_face_sum_with_recoil_by_color_and_block',
        'flux': compact_flux, 'terms': terms, 'term_count': len(terms),
        'scalar_polynomial_representation': 'fully_inline_factored_after_exact_expanded_checks',
        'common_C_by_channel': deepcopy(channels_C),
        'equality': 'sum_c C_c*J_c = sum_t combined_C_t*integral_t',
        'conditions': {
            'first': 'finite common I0=I*eta and alpha/t interior cutoffs, then contour IBP',
            'boundary_value': 'combine all terms before joint regulator removal and meromorphic continuation',
            'individual_term_convergence_proved': False,
            'original_order_preserved': True,
        },
        'evaluation_status': 'not_evaluated', 'formal': 'not_run',
    }
