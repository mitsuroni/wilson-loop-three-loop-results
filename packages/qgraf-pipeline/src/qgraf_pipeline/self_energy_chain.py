from .sign_rules import attach_contract
"""Two one-loop 1PI insertions in series, distinct from a two-loop 1PI blob."""
FAMILY = '2e-Pi1+Pi1'


def is_pi1_chain(orbit):
    return (orbit.emission_type == FAMILY and orbit.internal == {'Pi1': 2}
            and orbit.eikonal_count == 2 and orbit.g_power == 6)


def generate_pi1_chain(et, start_idx):
    from .generator import _canonical_edge_distributions, _build_two_internal_diagram
    if et.name != FAMILY or et.internal != {'Pi1': 2} or et.eikonal_count != 2 or et.wl_wl_pairs != 0 or et.g_power != 6:
        raise ValueError('unsupported_self_energy_chain')
    # Identical two-point blobs have one ordered Dyson word, not two label permutations.
    # The canonical edge generator excludes same-edge placements; a separate
    # negative test checks that an explicit same-edge chain is also filtered.
    return [_build_two_internal_diagram(start_idx+i, ec, [0], [1], 'Pi1', 'Pi1', 1,
                et_name=FAMILY) for i,ec in enumerate(_canonical_edge_distributions(2))]


def build_chain_definition(orbit, info):
    from .definition_exporter import _coord_diff_expr, _negate_phase
    if not is_pi1_chain(orbit):
        raise ValueError('unsupported_self_energy_chain')
    eps=[(edge,slot) for edge in range(1,5) for slot in range(1,info.edge_counts[edge-1]+1)]
    if len(eps)!=2 or info.n_wl_wl_gluons or len(info.internal_connections)!=2 or len(info.internal_gluon_links)!=1:
        raise ValueError('invalid_self_energy_chain_topology')
    a,b=eps
    return {
        'orbit_id':orbit.orbit_id,'emission_type':FAMILY,'sector':'dg','nloops':3,'g_power':6,
        'vertices':[{'line':edge,'lorentz':f'mu{i}','color':f'a{i}'} for i,(edge,slot) in enumerate(eps,1)],
        'propagators':[
            {'loop':'l1','lorentz':['mu1','nu1'],'color':['a1','b1']},
            {'loop':'l1','lorentz':['nu2','nu3'],'color':['b2','b3']},
            {'loop':'l1','lorentz':['nu4','mu2'],'color':['b4','a2']}],
        'internal_vertices':[
            {'type':'se','id':'se1','propagator':'l1','se_type':'vacuum-pol','lorentz':['nu1','nu2'],'color':['b1','b2']},
            {'type':'se','id':'se2','propagator':'l1','se_type':'vacuum-pol','lorentz':['nu3','nu4'],'color':['b3','b4']}],
        'phases':[{'loop':'l1','expression':_negate_phase(_coord_diff_expr(*a,*b))}],
        'bounds':[{'var':f's{edge}{slot}','low':'0','high':str(1) if slot==info.edge_counts[edge-1] else f's{edge}{slot+1}'} for edge,slot in eps],
        'self_energy_chain':{'kind':'one_particle_reducible','insertion_count':2,'ordered_insertions':['se1','se2'],
            'subgraph_loop_orders':[1,1],'dyson_coefficient':'1','source_external_sign':'+1',
            'sign_convention':'wilson-self-energy-1pi-insertion/2.0',
            'not_two_loop_1pi':True}}


def emit_chain_integral(orbit, info):
    from .emitter import _emit_measure, _emit_eikonal_chain, _emit_coord_diff
    definition=build_chain_definition(orbit,info)
    eps=[(edge,slot) for edge in range(1,5) for slot in range(1,info.edge_counts[edge-1]+1)]
    chain,_=_emit_eikonal_chain(info.edge_counts)
    # Each Pi symbol is the same one-loop tensor as the existing single-insertion source.
    # (-Pi1)^2 supplies a positive external sign; no new diagram factorial is added.
    return '\n'.join([
        f"I_{{({orbit.orbit_id})}} = {_emit_measure(info.edge_counts)}\\;{chain}",
        r'\times \int\frac{d^D\ell}{(2\pi)^D}\,'
        r'\frac{-ig_{\mu_1\nu_1}\delta^{a_1b_1}}{\ell^2+i0}\,'
        r'\Pi^{(1)\,\nu_1\nu_2}_{b_1b_2}(\ell)\,'
        r'\frac{-ig_{\nu_2\nu_3}\delta^{b_2b_3}}{\ell^2+i0}\,'
        r'\Pi^{(1)\,\nu_3\nu_4}_{b_3b_4}(\ell)\,'
        r'\frac{-ig_{\nu_4\mu_2}\delta^{b_4a_2}}{\ell^2+i0}\,'
        r'e^{-i\ell\cdot\bigl('+_emit_coord_diff(*eps[0],*eps[1])+r'\bigr)}',
        r'\quad\text{(two one-loop 1PI insertions in series; Dyson coefficient 1)}'])
