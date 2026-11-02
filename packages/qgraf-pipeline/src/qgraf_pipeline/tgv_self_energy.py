from .sign_rules import attach_contract
"""A bare three-gluon vertex with one one-loop 1PI insertion on a leg.

Topology labels in Diagram are graph edge IDs. Exported momenta impose the
two-point momentum constraint explicitly: both segments of the marked leg
carry the same l_r. This branch is not the one-loop 1PI vertex Gamma1.
"""
from collections import Counter

FAMILY = '3e-3g+Pi1'
CONVENTION = 'wilson-self-energy-1pi-insertion/2.0'


def is_tgv_pi1(orbit):
    return (orbit.emission_type == FAMILY and orbit.internal == {'3g': 1, 'Pi1': 1}
            and orbit.eikonal_count == 3 and orbit.g_power == 6)


def insertion_position(orbit, info):
    from .topology import vert_kind
    if not is_tgv_pi1(orbit):
        raise ValueError('unsupported_tgv_self_energy')
    d = orbit.representative
    kinds = Counter(vert_kind(v) for v in d.verts)
    if (kinds != {'eikonal': 3, 'three_gluon': 1, 'two_gluon': 1}
            or len(d.props) != 4 or info.n_wl_wl_gluons != 0
            or sorted(map(len, info.internal_connections)) != [1, 2]
            or info.internal_gluon_links != [(0, 1)]
            or sorted(p for c in info.internal_connections for p in c) != [1, 2, 3]):
        raise ValueError('invalid_tgv_self_energy_topology')
    # The two-point vertex has one direct contour attachment; the TGV has two.
    return next(c[0] for c in info.internal_connections if len(c) == 1)


def build_definition(orbit, info):
    from .definition_exporter import _coord_to_origin_expr, _negate_phase
    r = insertion_position(orbit, info)
    eps = [(edge, slot) for edge in range(1, 5)
           for slot in range(1, info.edge_counts[edge-1]+1)]
    props = []
    for j in range(1, 4):
        if j == r:
            props += [dict(loop=f'l{j}', lorentz=[f'mu{j}', 'rho'], color=[f'a{j}', 'c']),
                      dict(loop=f'l{j}', lorentz=['sigma', f'nu{j}'], color=['d', f'b{j}'])]
        else:
            props.append(dict(loop=f'l{j}', lorentz=[f'mu{j}', f'nu{j}'], color=[f'a{j}', f'b{j}']))
    bounds = [dict(var=f's{edge}{slot}', low='0',
                   high='1' if slot == info.edge_counts[edge-1] else f's{edge}{slot+1}')
              for edge, slot in eps]
    return dict(
        orbit_id=orbit.orbit_id, emission_type=FAMILY, sector='dg', nloops=3, g_power=6,
        vertices=[dict(line=edge, lorentz=f'mu{j}', color=f'a{j}')
                  for j, (edge, slot) in enumerate(eps, 1)],
        propagators=props,
        internal_vertices=[
            dict(type='tgv', id='v1', loops=['l1', 'l2', 'l3'],
                 lorentz=['nu1', 'nu2', 'nu3'], color=['b1', 'b2', 'b3'],
                 incoming_momenta=['-l1', '-l2', '-l3'], vertex_loop_order=0),
            dict(type='se', id='se1', propagator=f'l{r}', se_type='vacuum-pol',
                 lorentz=['rho', 'sigma'], color=['c', 'd'])],
        phases=[dict(loop=f'l{j}', expression=_negate_phase(_coord_to_origin_expr(*ep)))
                for j, ep in enumerate(eps, 1)],
        bounds=bounds,
        iterated_order_outer_to_inner=[f's{edge}{slot}' for edge in range(1, 5)
                                          for slot in range(info.edge_counts[edge-1], 0, -1)],
        momentum_constraints=[dict(type='delta', expression='l1+l2+l3',
                                   normalization='(2*pi)^D*delta^D(l1+l2+l3)')],
        sign_convention=CONVENTION,
        self_energy_insertion=dict(kind='1pi', operator='Pi1', insertion_count=1,
            subgraph_loop_order=1, skeleton_loop_order=2, total_loop_order=3,
            marked_contour_position=r, marked_endpoint=list(eps[r-1]),
            shared_momentum=f'l{r}', source_external_sign='+1', dyson_coefficient='1',
            tensor='Pi1=-i*(k^2*g-k*k)*P1', not_vertex_correction=True,
            scalar_contains_color=True, renormalization='bare'),
        downstream_status=dict(montopar='not_implemented', C_J='not_available', I='not_solved'))


def emit_integral(orbit, info):
    from .emitter import _emit_measure, _emit_eikonal_chain, _emit_coord_to_origin
    r = insertion_position(orbit, info)
    eps = [(edge, slot) for edge in range(1, 5)
           for slot in range(1, info.edge_counts[edge-1]+1)]
    chain, _ = _emit_eikonal_chain(info.edge_counts)
    out = [f'I_{{({orbit.orbit_id})}} = '+_emit_measure(info.edge_counts)+r'\;'+chain,
           r'\times\left[\prod_{j=1}^{3}\int\frac{d^D\ell_j}{(2\pi)^D}\right]'
           r'(2\pi)^D\delta^{(D)}(\ell_1+\ell_2+\ell_3)',
           r'\times V^{b_1b_2b_3}_{\nu_1\nu_2\nu_3}(-\ell_1,-\ell_2,-\ell_3)']
    for j, ep in enumerate(eps, 1):
        if j == r:
            out.append(rf'\times\frac{{-ig_{{\mu_{j}\rho}}\delta^{{a_{j}c}}}}{{\ell_{j}^2+i0}}'
                       rf'\Pi^{{(1)\,\rho\sigma}}_{{cd}}(\ell_{j})'
                       rf'\frac{{-ig_{{\sigma\nu_{j}}}\delta^{{db_{j}}}}}{{\ell_{j}^2+i0}}')
        else:
            out.append(rf'\times\frac{{-ig_{{\mu_{j}\nu_{j}}}\delta^{{a_{j}b_{j}}}}}{{\ell_{j}^2+i0}}')
        out.append(rf'e^{{-i\ell_{j}\cdot\bigl('+_emit_coord_to_origin(*ep)+r'\bigr)}')
    out.append(r'\quad\text{(bare TGV; one Pi1 leg insertion; positive Dyson connection)}')
    return '\n'.join(out)


def emit_drawing(orbit, info, positions):
    from .tikz_emitter import _push_inward
    r = insertion_position(orbit, info)
    eps = sorted(positions)
    xy = [positions[p] for p in eps]
    center = _push_inward((sum(p[0] for p in xy)/3, sum(p[1] for p in xy)/3), .30)
    out = []
    for j, (x, y) in enumerate(xy, 1):
        out.append(rf'  \draw ({x:.2f},{y:.2f}) -- ({center[0]:.2f},{center[1]:.2f});')
        if j == r:
            bx, by = ((x+center[0])/2, (y+center[1])/2)
            out.append(rf'  \draw[fill=white] ({bx:.2f},{by:.2f}) circle (0.12cm);')
    out.append(rf'  \filldraw ({center[0]:.2f},{center[1]:.2f}) circle (1.5pt);')
    return out
