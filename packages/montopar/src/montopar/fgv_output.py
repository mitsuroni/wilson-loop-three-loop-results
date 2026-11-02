"""Self-contained single-representative FGV outputs and compatibility exports."""
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path
import json
from montopar.mixed_ibp import ordered_chains
from montopar.fgv import compute_Q_prime, compute_lorentzZ
from montopar.j_factor import expand_polynomial

ALPHAS = [f'alpha{i}' for i in range(1,5)]
MONOMIAL = '(alpha1*alpha2*alpha3*alpha4)^(-eps)'
POWER = '2 - 3*eps'

def alpha_simplex4():
    return {'kind':'simplex','variables':list(ALPHAS),
        'constraint':'alpha1+alpha2+alpha3+alpha4=1',
        'nonnegative_variables':list(ALPHAS),
        'independent_variables':ALPHAS[:3],
        'substitution':{'alpha4':'1-alpha1-alpha2-alpha3'},
        'bounds':{'alpha1':['0','1'],'alpha2':['0','1-alpha1'],
                  'alpha3':['0','1-alpha1-alpha2']},
        'iterated_order_outer_to_inner':ALPHAS[:3], 'jacobian':'1',
        'weight_location':'integrand'}

def branch_prescription():
    return {'propagator':'-I/(p^2+I0)', 'coordinate_distance':'-(Xi-Xj)^2+I0',
        'star_denominator':'-Sum[i<j,alpha_i*alpha_j*(Xi-Xj)^2]+I0',
        'reference_power':'(-2*v14+I0)^(3*eps)',
        'parameter_starting_region':{'z1':'>0','v14':'<0'},
        'power_convention':'principal branch, continued from the stated starting region',
        'complex_continuation_verified':False}

def attach_measure(block):
    block['alpha_vars']=list(ALPHAS)
    block['alpha_measure']=alpha_simplex4()
    block['delta_constraints']={'alpha1+alpha2+alpha3+alpha4':'1'}
    block['integration_order']='path_parameters_outer_alpha_simplex_inner'
    return block

def _validate_orbit(orbit):
    if orbit.dispatch!='fourv' or orbit.emission_type!='4e-4g' or orbit.source.get('nloops')!=3:
        raise ValueError('expected_three_loop_4e4g_source')
    if len(orbit.fgv_legs)!=4 or len(orbit.propagators)!=4 or len(orbit.vertices)!=4:
        raise ValueError('expected_four_distinct_fgv_legs')
    if len({t.s_var for t in orbit.fgv_legs})!=4:
        raise ValueError('repeated_fgv_parameter')
    if [p.color[0] for p in orbit.propagators]!=[v.color for v in orbit.vertices]:
        raise ValueError('fgv_color_word_not_in_contour_order')
    if [p.lorentz[0] for p in orbit.propagators]!=[v.lorentz for v in orbit.vertices]:
        raise ValueError('fgv_lorentz_word_not_in_contour_order')
    edge={t.s_var:t.edge for t in orbit.fgv_legs}
    chains=ordered_chains(orbit.bounds)
    if set(edge)!={b.var for b in orbit.bounds}:
        raise ValueError('fgv_endpoint_domain_mismatch')
    if any(len({edge[v] for v in chain})!=1 for chain in chains):
        raise ValueError('fgv_order_crosses_edges')
    if len(chains)!=len(set(edge.values())):
        raise ValueError('fgv_same_edge_split_into_independent_chains')
    # The color formula is for this contour word, including same-edge slot order.
    for chain in chains:
        actual=[t.s_var for t in orbit.fgv_legs if t.edge==edge[chain[0]]]
        if actual!=chain:
            raise ValueError('fgv_contour_word_disagrees_with_bounds')
    return chains

def finalize_fgv(result, orbit, hc):
    chains=_validate_orbit(orbit)
    q=expand_polynomial(compute_Q_prime(orbit.fgv_legs))
    lz=expand_polynomial(compute_lorentzZ(orbit.fgv_legs))
    iszero=lz=='0'
    for key in list(result):
        if key.endswith('_rotated'): del result[key]
    result['schema_version']='1.2'
    result['emission_type']='4e-4g'
    result['representative']={
        'scheme':'single_source_contour_representative', 'rotation':0,
        'leg_permutation':[1,2,3,4],
        'legs':[asdict(t) for t in orbit.fgv_legs],
        'canonical_topology':deepcopy(orbit.source),
        'source_topology':deepcopy(orbit.source),
        's_to_x':dict(hc.s_to_x),'jacobian':hc.jacobian}
    result['normalization']={
        'scheme':'montopar-raw-CJ','loop_order':3,
        'euler_compensation_count':0,'color_projection':'none',
        'color_trace':'normalized fundamental SU(N) trace',
        'orbit_multiplicity_included':False,'crossing_assembled':False,
        'note':'Euler compensation and physical I assembly belong to downstream preparation'}
    result['branch_prescription']=branch_prescription()
    result['evaluation_status']='algebraic_zero' if iszero else 'not_evaluated'
    result['formal']='not_run'
    result['zero_reason']='all_three_Lorentz_channels_vanish' if iszero else None
    result['J_raw_s']=attach_measure({
        'representation':'ordered_s','integrand':'0' if iszero else
            f'({lz}) * ({MONOMIAL}) * ({q})^(-({POWER}))',
        'variables':[b.var for b in orbit.bounds],
        'bounds':{b.var:[b.low,b.high] for b in orbit.bounds},
        'ordered_chains':chains,
        'iterated_order_outer_to_inner':[v for c in chains for v in reversed(c)],
        'lorentzZ':lz,'Q_prime':q,'alpha_monomial':MONOMIAL,
        'denom_power':POWER,'jacobian':'1','physical_zero':iszero})
    cube=result['J_raw']
    cube.update({'representation':'unit_cube','variables':list(hc.x_vars),
        'bounds':{x:[0,1] for x in hc.x_vars},
        'iterated_order_outer_to_inner':list(hc.x_vars),
        's_to_x':dict(hc.s_to_x),'jacobian':hc.jacobian,'physical_zero':iszero})
    attach_measure(cube)
    validate_fgv_result(result)
    return result

def validate_fgv_result(result):
    if result.get('schema_version')!='1.2' or result.get('dispatch')!='fourv':
        raise ValueError('regenerate_self_contained_fgv_source')
    if any(k.endswith('_rotated') for k in result):
        raise ValueError('fgv_requires_single_representative')
    for key in ('representative','normalization','branch_prescription','J_raw_s','C','J_raw'):
        if key not in result: raise ValueError('missing_fgv_'+key)
    if result['normalization'].get('euler_compensation_count')!=0:
        raise ValueError('fgv_raw_must_not_contain_euler_compensation')
    for name in ('J_raw_s','J_raw'):
        j=result[name]
        if j.get('alpha_measure')!=alpha_simplex4() or j.get('alpha_vars')!=ALPHAS:
            raise ValueError('fgv_alpha_simplex_incomplete')
        if j.get('delta_constraints')!={'alpha1+alpha2+alpha3+alpha4':'1'}:
            raise ValueError('fgv_delta_constraint_missing')
        if set(j.get('iterated_order_outer_to_inner',[]))!=set(j.get('variables',[])):
            raise ValueError('fgv_iterated_order_incomplete')
    cs=result['C']; comp=cs['components']
    expected={'color','coupling','ft_prefactor','v14_scale','twoscale','i_factor'}
    if set(comp)!=expected or cs.get('component_semantics')!='multiplicative':
        raise ValueError('fgv_C_components_incomplete')
    if len(cs.get('factor_order',[]))!=len(expected) or set(cs['factor_order'])!=expected:
        raise ValueError('fgv_C_factor_order_invalid')
    if result['branch_prescription']!=branch_prescription():
        raise ValueError('fgv_branch_prescription_incomplete')
    if cs['combined']!=' * '.join(f'({comp[k]})' for k in cs['factor_order']):
        raise ValueError('fgv_C_component_product_mismatch')
    return True

def write_compatibility_exports(result, directory):
    """Derive legacy filenames from one validated current result, never old state."""
    validate_fgv_result(result)
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    j=deepcopy(result['J_raw']); c=deepcopy(result['C']); oid=result['orbit_id']
    state={'schema_version':'1.2','orbit_id':oid,'integral_id':oid,
        'C':c['combined'],'J':j['integrand'],'J1_integrand':j['integrand'],
        'J_raw_s':deepcopy(result['J_raw_s']),
        **{k:deepcopy(j[k]) for k in ('variables','bounds','alpha_vars','alpha_measure','delta_constraints','iterated_order_outer_to_inner','integration_order')},
        'x_vars':j['variables'],'geometry':deepcopy(result['geometry']),
        'reference_scale':'v14','normalization':deepcopy(result['normalization']),
        'representative':deepcopy(result['representative']),
        'branch_prescription':deepcopy(result['branch_prescription']),
        'evaluation_status':result['evaluation_status'],'formal':'not_run'}
    final={'version':'1.2','geometry':deepcopy(result['geometry']),
        'momentum_integral':{'source':deepcopy(result['representative']['source_topology']),'expression':{}},
        'results':{'C':c,'J':j,'J_raw_s':deepcopy(result['J_raw_s']),
            'color':c['components']['color'],'lorentzZ':j['lorentzZ'],
            'Q_prime':j['Q_prime'],'alpha_monomial':j['alpha_monomial'],
            'denom_power':j['denom_power'],'jacobian':j['jacobian'],
            'alpha_vars':list(ALPHAS),'alpha_measure':alpha_simplex4(),
            'physical_zero':j['physical_zero'],'representation':'parametric',
            'normalization':deepcopy(result['normalization']),
            'branch_prescription':deepcopy(result['branch_prescription'])},
        'provenance':{'generator':'montopar.fgv_output.write_compatibility_exports',
                      'source':'validated single-representative montopar.json'}}
    for name,payload in (('final.json',final),('state.json',state),('export.dimreduce.json',state)):
        target=directory/name;temp=target.with_suffix('.json.tmp')
        temp.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf8');temp.replace(target)