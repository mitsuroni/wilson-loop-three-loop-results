from montopar.sign_rules import validate_contract, annotate_result
"""Bare Pi1-D0-Pi1 chain, with the complete D-dimensional tensor and source path.

Theory: docs/self-energy-pi1-chain.md. This is not a two-loop 1PI operator.
No CAS runtime dependency; raw output has neither color projection nor Euler factor.
"""
from collections import Counter, defaultdict
from copy import deepcopy
import json
import re

from montopar.models import Bound, EdgeDistribution, Orbit, Propagator, Vertex
from montopar.self_energy import (_poly, _closure, _position, geometry_for_line,
                                 c4_symmetry as _c4, SE_COLOR)
from montopar.j_factor import expand_polynomial

FAMILY = '2e-Pi1+Pi1'
SCHEMA = '1.5'
GAMMA_RATIO = 'Gamma[1+eps]^2*Gamma[1-eps]*Gamma[1-3*eps]/(Gamma[1-2*eps]^2*Gamma[1+2*eps])'
SIGN_CONVENTION = 'wilson-self-energy-1pi-insertion/2.0'


def _identifier(x):
    return isinstance(x, str) and re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', x) is not None


def _validate_source(raw):
    validate_contract(raw)
    required = {'orbit_id','emission_type','sector','nloops','g_power','vertices','propagators',
                'internal_vertices','phases','bounds','self_energy_chain'}
    if not required <= set(raw) or not set(raw) <= required | {'iterated_order_outer_to_inner','physical_convention'}:
        raise ValueError('chain_unknown_or_missing_source_fields')
    if raw['emission_type'] != FAMILY or not isinstance(raw['orbit_id'],str):
        raise ValueError('unsupported_self_energy_chain_family')
    if type(raw['nloops']) is not int or type(raw['g_power']) is not int or (raw['nloops'],raw['g_power'],raw['sector']) != (3,6,'dg'):
        raise ValueError('chain_loop_coupling_or_sector_mismatch')
    vertices = raw['vertices']
    if len(vertices)!=2 or any(set(v)!={'line','lorentz','color'} for v in vertices):
        raise ValueError('chain_requires_two_endpoints')
    if any(type(v['line']) is not int or v['line'] not in range(1,5) for v in vertices) or not vertices[0]['line']<vertices[1]['line']:
        raise ValueError('chain_requires_two_distinct_ordered_edges')
    if any(not _identifier(v[k]) for v in vertices for k in ['lorentz','color']):
        raise ValueError('chain_bad_external_index')
    if any(len({v[k] for v in vertices})!=2 for k in ['lorentz','color']):
        raise ValueError('chain_duplicate_external_index')
    endpoints=[dict(v,edge=v['line'],slot=1,s_var=f"s{v['line']}1",rank=i) for i,v in enumerate(vertices)]
    bounds={v['s_var']:['0','1'] for v in endpoints};order=list(bounds)
    if len(raw['bounds'])!=2 or any(set(b)!={'var','low','high'} for b in raw['bounds']) or {b['var']:[str(b['low']),str(b['high'])] for b in raw['bounds']}!=bounds:
        raise ValueError('chain_original_domain_mismatch')
    if 'iterated_order_outer_to_inner' in raw and raw['iterated_order_outer_to_inner']!=order:
        raise ValueError('chain_original_order_mismatch')
    ivs=raw['internal_vertices'];props=raw['propagators']
    if len(ivs)!=2 or len(props)!=3:
        raise ValueError('chain_requires_two_insertions_three_propagators')
    if any(set(v)!={'type','id','propagator','se_type','lorentz','color'} for v in ivs):
        raise ValueError('chain_insertion_fields_mismatch')
    if any(v['type']!='se' or v['se_type']!='vacuum-pol' or not _identifier(v['id']) for v in ivs) or len({v['id'] for v in ivs})!=2:
        raise ValueError('chain_requires_two_distinct_one_loop_operators')
    if any(set(p)!={'loop','lorentz','color'} for p in props):
        raise ValueError('chain_propagator_fields_mismatch')
    momenta={p['loop'] for p in props}|{v['propagator'] for v in ivs}
    if len(momenta)!=1 or not all(_identifier(m) for m in momenta):
        raise ValueError('chain_requires_one_common_momentum')
    loop=next(iter(momenta))
    # Treat each (Lorentz,color) pair as a port. The path must alternate P,S,P,S,P.
    edges=[];mu_to_color={};color_to_mu={};ports=[]
    for kind,objects in [('propagator',props),('insertion',ivs)]:
        for idx,obj in enumerate(objects):
            if any(not isinstance(obj[k],list) or len(obj[k])!=2 or not all(_identifier(t) for t in obj[k]) for k in ['lorentz','color']):
                raise ValueError('chain_malformed_port_pair')
            pair=list(zip(obj['lorentz'],obj['color']))
            if pair[0]==pair[1]:raise ValueError('chain_self_connected_component')
            for mu,color in pair:
                if (mu in mu_to_color and mu_to_color[mu]!=color) or (color in color_to_mu and color_to_mu[color]!=mu):
                    raise ValueError('chain_lorentz_color_connectivity_disagrees')
                mu_to_color[mu]=color;color_to_mu[color]=mu
            ports.extend(pair);edges.append((kind,obj.get('id',idx),pair))
    ext=[(v['lorentz'],v['color']) for v in vertices];counts=Counter(ports)
    if set(ext)-set(counts) or len(counts)!=6 or any(n!=(1 if p in ext else 2) for p,n in counts.items()):
        raise ValueError('chain_internal_index_multiplicity_mismatch')
    neighbors=defaultdict(list)
    for n,(_,_,pair) in enumerate(edges):
        for port in pair:neighbors[port].append(n)
    current=ext[0];used=[];path=[]
    while current!=ext[1]:
        candidates=[i for i in neighbors[current] if i not in used]
        if len(candidates)!=1:raise ValueError('chain_disconnected_or_branched_path')
        idx=candidates[0];kind,label,pair=edges[idx]
        path.append({'kind':kind,'id':label,'ports':[list(current),list(pair[1] if current==pair[0] else pair[0])]})
        current=pair[1] if current==pair[0] else pair[0];used.append(idx)
        if len(used)>5:raise ValueError('chain_cycle')
    if len(used)!=5 or [p['kind'] for p in path]!=['propagator','insertion','propagator','insertion','propagator']:
        raise ValueError('chain_not_one_alternating_path')
    ordered_ids=[p['id'] for p in path if p['kind']=='insertion']
    metadata=raw['self_energy_chain']
    expected={'kind':'one_particle_reducible','insertion_count':2,'ordered_insertions':ordered_ids,
              'subgraph_loop_orders':[1,1],'dyson_coefficient':'1','source_external_sign':'+1',
              'sign_convention':SIGN_CONVENTION,'not_two_loop_1pi':True}
    if metadata!=expected or type(metadata['insertion_count']) is not int or any(type(n) is not int for n in metadata['subgraph_loop_orders']):
        raise ValueError('chain_operator_order_sign_or_Dyson_metadata_mismatch')
    phases=raw['phases']
    if len(phases)!=1 or set(phases[0])!={'loop','expression'} or phases[0]['loop']!=loop:
        raise ValueError('chain_requires_one_endpoint_phase')
    phase=_poly(phases[0]['expression'],[*bounds,'v1','v2','v3','v4'])
    expected_phi=f'({_position(endpoints[1])})-({_position(endpoints[0])})'
    if _closure(f'({phase})-({expected_phi})')=='0':orientation=1
    elif _closure(f'({phase})+({expected_phi})')=='0':orientation=-1
    else:raise ValueError('chain_phase_endpoint_geometry_mismatch')
    line={'kind':'self_energy_chain','loop':loop,'endpoints':endpoints,'propagators':deepcopy(props),
          'phase_Phi':phases[0]['expression'],'phase_orientation':orientation,'ordered_insertions':ordered_ids}
    return {'lines':[line],'endpoints':endpoints,'bounds':bounds,'order':order,'contraction_path':path}


def validate_source(raw):
    try:return _validate_source(raw)
    except (KeyError,TypeError,IndexError) as exc:
        raise ValueError('chain_malformed_source') from exc


def parse_chain_orbit(raw):
    data=validate_source(raw)
    return Orbit(orbit_id=raw['orbit_id'],emission_type=FAMILY,n_gluons=0,pairings=[],
        bounds=[Bound(b['var'],str(b['low']),str(b['high'])) for b in raw['bounds']],
        edge_distribution=EdgeDistribution({v['edge']:1 for v in data['endpoints']},{v['edge']:[v['s_var']] for v in data['endpoints']}),
        vertices=[Vertex(**v) for v in raw['vertices']],
        propagators=[Propagator(p['loop'],tuple(p['lorentz']),tuple(p['color'])) for p in raw['propagators']],
        source=deepcopy(raw),dispatch='self_energy_chain')


def c4_symmetry(raw):
    result=_c4(raw,validate_source(raw))
    result['definition']='C4 on the two chain endpoints; identical Pi1 labels do not create extra multiplicity'
    return result


def build_result(orbit):
    raw=orbit.source;data=validate_source(raw)
    if orbit!=parse_chain_orbit(raw):raise ValueError('chain_orbit_metadata_disagrees_with_source')
    se=geometry_for_line(data['lines'][0]);f,d,u,w=(se[k] for k in ['f','d','u','w'])
    numerator=expand_polynomial(f'(1+4*eps)*({d})*({f})+(1-3*eps)*({u})*({w})')
    j={'representation':'ordered_s','variables':data['order'],'bounds':data['bounds'],
       'iterated_order_outer_to_inner':data['order'],'source_bounds_array':deepcopy(raw['bounds']),
       'order_basis':'ascending contour edges; frozen chain source order','jacobian':'1',
       'integrand':f'({numerator})*({f})^(-2+3*eps)/(2*(1+2*eps))',
       'se_numerator':numerator,'geometry_polynomials':[se],
       'power_prescription':'each f power is the common Feynman boundary value (f+I0)^power; continue jointly with C'}
    color={'wilson_trace':'CF','self_energy_color':f'({SE_COLOR})^2','full_color':f'CF*({SE_COLOR})^2',
           'flavor_channels':[
               {'monomial':'CA^2*CF','coefficient':'(10-6*eps)^2/(3-2*eps)^2'},
               {'monomial':'CA*CF*TF*nf','coefficient':'-16*(10-6*eps)*(1-eps)/(3-2*eps)^2'},
               {'monomial':'CF*(TF*nf)^2','coefficient':'64*(1-eps)^2/(3-2*eps)^2'}],
           'TF':'1/2','nf':'number of massless quark flavours','projection':'none',
           'method':'normalized Tr/N; path contracts all internal deltas to t^a t^a = CF'}
    components={'wilson_color':'CF','self_energy_color':color['self_energy_color'],
                'coupling':'(alpha_s/Pi)^3','gamma_fourier':'(Pi^eps*Gamma[1-eps])^3',
                'gamma_chain_ratio':GAMMA_RATIO,'rational':'-1/(128*eps^2*(1-2*eps)^2)',
                'scale':'(-2*v14+I0)^(3*eps)'}
    c={'components':components,'factor_order':list(components),'component_semantics':'multiplicative',
       'combined':'*'.join(f'({v})' for v in components.values()),
       'phase_ledger':{'source_external_sign':'+1','wilson_vertices':'I^2','three_propagators':'(-I)^3',
           'Pi1_scalar_phases':'I^2','momentum_prefactor_phase':'I','q_conversion':'-1','Fourier_phase':'I',
           'coordinate_sign_before_scale':'1','vref_over_Aref':'-1/2','overall_dimensionless_sign':'-1','dyson_coefficient':'1'}}
    result={'schema_version':SCHEMA,'orbit_id':raw['orbit_id'],'emission_type':FAMILY,'dispatch':'self_energy_chain','loop_order':3,
      'geometry':{'n_lines':4,'z_definition':'z1=v12/v14','reference_scale':'v14','vij_definition':'v_i dot v_j',
                  'A_ref':'-2*v14','lightlike_edges':True,'line_invariants':'f=X^2/(2*v14); d=vi.vj/v14; u=vi.X/v14; w=vj.X/v14'},
      'representative':{'scheme':'single_source_representative','rotation':0,'source_topology':deepcopy(raw),
                        'canonical_topology':deepcopy(raw),'effective_lines':data['lines'],'contraction_path':data['contraction_path']},
      'operator':{'name':'Pi1+Pi1','kind':'one_particle_reducible','insertion_count':2,'subgraph_loop_orders':[1,1],
          'insertions':deepcopy(raw['internal_vertices']),'dyson_coefficient':'1','renormalization':'bare','dimension':'4-2*eps','gauge':'Feynman',
          'Pi1_source':'hep-ph/9605348 Appendix C (C.1),(C.3); xi_D=0','theory':'docs/self-energy-pi1-chain.md, SC1-SC24',
          'Pi1_scalar':'I*g^2/(2*(4*Pi)^(2-eps))*Gamma[eps]*Gamma[1-eps]^2/Gamma[2-2*eps]*C_SE*(-l^2-I0)^(-eps)',
          'C_SE':SE_COLOR,'tensor_identity':'T*g*T=l^2*T; T=l^2*g-l*l','epsilon_precision':'exact'},
      'color_contraction':color,'c4_symmetry':c4_symmetry(raw),
      'normalization':{'scheme':'montopar-raw-CJ','loop_order':3,'color_projection':'none','euler_compensation_count':0,
          'orbit_multiplicity_included':False,'crossing_assembled':False,'normalized_wilson_trace':'Tr/N','extra_factorial':'none',
          'gamma_residual':GAMMA_RATIO,'gamma_residual_Euler_linear_coefficient':'0'},
      'C':c,'J_raw':j,'J_raw_s':deepcopy(j),
      'branch_prescription':{'dimensionful_distance':'-X^2+I0','momentum':'-l^2-I0',
          'starting_region':{'z1':'>0','v14':'<0','eps':'0<Re(eps)<1/3'},
          'power_convention':'common Feynman boundary value; joint analytic continuation','complex_continuation_verified':False},
      'translation_status':'translated','evaluation_status':'not_evaluated','formal':'not_run'}
    return annotate_result(result, raw)


def validate_result(result):
    if result.get('schema_version')!=SCHEMA or result.get('dispatch')!='self_energy_chain':
        raise ValueError('regenerate_self_energy_chain_result')
    expected=build_result(parse_chain_orbit(result['representative']['source_topology']))
    if json.loads(json.dumps(expected))!=result:raise ValueError('chain_result_does_not_reconstruct_from_source')
