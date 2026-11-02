from montopar.sign_rules import three_slot_trace_phase, annotate_result, validate_contract
"""Tree TGV with one bare Pi1 leg, exact transverse tensor and ordered C/J.

Formula contract: docs/tgv-pi1.md. No CAS runtime, Euler compensation or
orbit multiplicity in raw C. Optional contour IBP retains raw and all faces.
"""
from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
import json
import re

from montopar.geometry import dot_product_factor
from montopar.j_factor import expand_polynomial, substitute_s_to_x
from montopar.models import Bound, EdgeDistribution, Orbit, Propagator, TGVLeg, Vertex
from montopar.self_energy import SE_COLOR, _poly, _position, _closure

SCHEMA = '1.7'
SIGN = 'wilson-self-energy-1pi-insertion/2.0'
ALPHAS = ['alpha1', 'alpha2', '(1-alpha1-alpha2)']


def add(*terms):
    return expand_polynomial('+'.join(f'({t})' for t in terms) or '0')


def mul(*terms):
    return expand_polynomial('*'.join(f'({t})' for t in terms) or '1')


def sub(a, b):
    return add(a, mul('-1', b))


def validate_source(raw):
    """Reconstruct the marked leg, indices, momentum delta, phases and order."""
    validate_contract(raw)
    required = {'orbit_id', 'emission_type', 'sector', 'nloops', 'g_power', 'vertices',
                'propagators', 'internal_vertices', 'phases', 'bounds',
                'iterated_order_outer_to_inner', 'momentum_constraints',
                'sign_convention', 'self_energy_insertion'}
    if not required <= set(raw):
        raise ValueError('tgv_pi1_missing_source_fields')
    if not set(raw) <= required | {'downstream_status','physical_convention'}:
        raise ValueError('tgv_pi1_unknown_source_fields')
    if raw['emission_type'] != '3e-3g+Pi1' or raw['nloops'] != 3 or raw['g_power'] != 6 or raw['sector'] != 'dg':
        raise ValueError('tgv_pi1_family_or_loop_order')
    if raw['sign_convention'] != SIGN:
        raise ValueError('tgv_pi1_requires_sign_convention_2_0')
    vertices = raw['vertices']
    if len(vertices) != 3 or any(set(v) != {'line', 'lorentz', 'color'} for v in vertices):
        raise ValueError('tgv_pi1_external_vertices')
    edges = [v['line'] for v in vertices]
    if any(type(e) is not int or e not in range(1,5) for e in edges) or edges != sorted(edges) or len(set(edges)) == 1:
        raise ValueError('tgv_pi1_contour_edges_or_scaleless')
    for key in ('lorentz', 'color'):
        if len({v[key] for v in vertices}) != 3:
            raise ValueError('tgv_pi1_duplicate_external_index')
    slots = Counter()
    ends = []
    for v in vertices:
        slots[v['line']] += 1
        ends.append(dict(v, edge=v['line'], slot=slots[v['line']], s_var=f"s{v['line']}{slots[v['line']]}"))
    expected = {p['s_var']: ['0', f"s{p['edge']}{p['slot']+1}" if p['slot'] < slots[p['edge']] else '1'] for p in ends}
    bounds = raw['bounds']
    if (len(bounds) != 3 or any(set(b) != {'var','low','high'} for b in bounds)
            or len({b['var'] for b in bounds}) != 3
            or {b['var']: [str(b['low']),str(b['high'])] for b in bounds} != expected):
        raise ValueError('tgv_pi1_ordered_bounds')
    order = raw['iterated_order_outer_to_inner']
    if len(order) != 3 or set(order) != set(expected):
        raise ValueError('tgv_pi1_iterated_order')
    if any(high != '1' and order.index(high) > order.index(var) for var, (_,high) in expected.items()):
        raise ValueError('tgv_pi1_parent_must_be_outer')
    ivs = raw['internal_vertices']
    if len(ivs) != 2 or sorted(v.get('type','') for v in ivs) != ['se','tgv']:
        raise ValueError('tgv_pi1_requires_one_tree_tgv_and_one_pi1')
    tgv = next(v for v in ivs if v['type'] == 'tgv')
    se = next(v for v in ivs if v['type'] == 'se')
    if (set(tgv) != {'type','id','loops','lorentz','color','incoming_momenta','vertex_loop_order'}
            or set(se) != {'type','id','propagator','se_type','lorentz','color'}):
        raise ValueError('tgv_pi1_ambiguous_internal_operator')
    loops = tgv.get('loops', [])
    if (len(loops) != 3 or len(set(loops)) != 3 or tgv.get('vertex_loop_order') != 0
            or tgv.get('incoming_momenta') != ['-'+l for l in loops]
            or len(tgv.get('lorentz',[])) != 3 or len(tgv.get('color',[])) != 3):
        raise ValueError('tgv_pi1_tgv_momenta_or_indices')
    if se.get('se_type') != 'vacuum-pol' or len(se.get('lorentz',[])) != 2 or len(se.get('color',[])) != 2:
        raise ValueError('tgv_pi1_operator_indices')
    ins = raw['self_energy_insertion']
    fixed = {'kind':'1pi','operator':'Pi1','insertion_count':1,'subgraph_loop_order':1,
             'skeleton_loop_order':2,'total_loop_order':3,'source_external_sign':'+1',
             'dyson_coefficient':'1','tensor':'Pi1=-i*(k^2*g-k*k)*P1',
             'not_vertex_correction':True,'scalar_contains_color':True,'renormalization':'bare'}
    if any(ins.get(k) != v for k,v in fixed.items()):
        raise ValueError('tgv_pi1_operator_or_dyson_contract')
    marked = ins.get('marked_contour_position')
    if type(marked) is not int or marked not in (1,2,3):
        raise ValueError('tgv_pi1_marked_position')
    r = marked-1
    if ins.get('marked_endpoint') != [ends[r]['edge'],ends[r]['slot']] or ins.get('shared_momentum') != loops[r] or se.get('propagator') != loops[r]:
        raise ValueError('tgv_pi1_marked_momentum_or_endpoint')
    # Each segment connects the exact external/SE/TGV index pair at equal momentum.
    # Reversing a segment's two ends is harmless; changing the pairing is not.
    pairs = []
    for i,l in enumerate(loops):
        ext = (ends[i]['lorentz'], ends[i]['color'])
        inner = (tgv['lorentz'][i], tgv['color'][i])
        if i == r:
            left = (se['lorentz'][0], se['color'][0])
            right = (se['lorentz'][1], se['color'][1])
            pairs += [(l, tuple(sorted([ext,left]))), (l,tuple(sorted([right,inner])))]
        else:
            pairs.append((l,tuple(sorted([ext,inner]))))
    props = raw['propagators']
    if len(props) != 4 or any(set(p) != {'loop','lorentz','color'} or len(p['lorentz']) != 2 or len(p['color']) != 2 for p in props):
        raise ValueError('tgv_pi1_four_segments')
    actual = [(p['loop'],tuple(sorted(zip(p['lorentz'],p['color'])))) for p in props]
    if Counter(actual) != Counter(pairs):
        raise ValueError('tgv_pi1_segment_connectivity')
    for key in ('lorentz','color'):
        all_indices = [v[key] for v in vertices] + tgv[key] + se[key]
        if len(set(all_indices)) != 8 or any(not isinstance(x,str) or not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*',x) for x in all_indices):
            raise ValueError('tgv_pi1_index_bijection')
    phases = raw['phases']
    if len(phases) != 3 or {p['loop'] for p in phases} != set(loops):
        raise ValueError('tgv_pi1_phase_bijection')
    by_phase = {p['loop']:p['expression'] for p in phases}
    for i,l in enumerate(loops):
        phase = _poly(by_phase[l], [*expected,'v1','v2','v3','v4'])
        if _closure(add(phase,_position(ends[i]))) != '0':
            raise ValueError('tgv_pi1_phase_endpoint_mismatch')
    constraints = raw['momentum_constraints']
    total = '+'.join(loops)
    if (len(constraints) != 1 or constraints[0].get('type') != 'delta'
            or _poly(constraints[0].get('expression',''),loops) != expand_polynomial(total)
            or constraints[0].get('normalization') != f'(2*pi)^D*delta^D({total})'):
        raise ValueError('tgv_pi1_central_delta')
    return {'endpoints':ends, 'loops':loops, 'marked_leg':r+1, 'order':deepcopy(order),
            'bounds':{var:expected[var] for var in order},
            'ordered_chains':[[p['s_var'] for p in ends if p['edge']==e] for e in sorted(slots)]}


def parse_tgv_pi1_orbit(raw):
    data = validate_source(raw)
    by_edge = {p['edge']:[v['s_var'] for v in data['endpoints'] if v['edge']==p['edge']] for p in data['endpoints']}
    return Orbit(raw['orbit_id'],raw['emission_type'],0,
                 EdgeDistribution({e:len(v) for e,v in by_edge.items()},by_edge), [],
                 [Bound(b['var'],str(b['low']),str(b['high'])) for b in raw['bounds']],
                 [Vertex(**v) for v in raw['vertices']],
                 [Propagator(p['loop'],tuple(p['lorentz']),tuple(p['color'])) for p in raw['propagators']],
                 source=deepcopy(raw),dispatch='tgv_pi1',
                 tgv_legs=[TGVLeg(p['edge'],p['s_var'],l) for p,l in zip(data['endpoints'],data['loops'])])


def geometry(data):
    """Reuse the four-edge Gram convention; all polynomials use A0=-2*v14."""
    ends = data['endpoints']
    gram = [[mul('-1/2',dot_product_factor(i,j)) for j in range(1,5)] for i in range(1,5)]
    dot = lambda x,y: add(*(mul(x[i],gram[i][j],y[j]) for i in range(4) for j in range(4)))
    x = [['1' if k < p['edge'] else p['s_var'] if k == p['edge'] else '0' for k in range(1,5)] for p in ends]
    u = [['1' if k==p['edge'] else '0' for k in range(1,5)] for p in ends]
    center = [add(*(mul(ALPHAS[i],x[i][k]) for i in range(3))) for k in range(4)]
    rr = [[sub(x[i][k],center[k]) for k in range(4)] for i in range(3)]
    qq = [[mul('-1',dot([sub(x[i][k],x[j][k]) for k in range(4)],
                         [sub(x[i][k],x[j][k]) for k in range(4)])) for j in range(3)] for i in range(3)]
    g = [[dot(u[i],u[j]) for j in range(3)] for i in range(3)]
    t = [[dot(u[i],rr[j]) for j in range(3)] for i in range(3)]
    w = [[dot(rr[i],rr[j]) for j in range(3)] for i in range(3)]
    q0 = add(*(mul(ALPHAS[i],ALPHAS[j],qq[i][j]) for i in range(3) for j in range(i+1,3)))
    l = add(*(mul(g[i][j],sub(mul(ALPHAS[i],t[k][i]),mul(ALPHAS[j],t[k][j]))) for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]))
    r = data['marked_leg']-1; s=(r+1)%3; v=(r+2)%3
    e = add(mul(add('1',mul('-1',ALPHAS[r]),ALPHAS[s]),g[r][v],t[s][r]),
            mul('-1',add('1',mul('-1',ALPHAS[r]),ALPHAS[v]),g[r][s],t[v][r]),
            mul('4-2*eps',sub(ALPHAS[v],ALPHAS[s]),g[s][v],t[r][r]))
    h = add(mul('-1',ALPHAS[s],t[s][r],t[v][s]),mul(ALPHAS[v],t[v][r],t[s][v]),
            mul(g[s][v],sub(mul(ALPHAS[s],w[r][s]),mul(ALPHAS[v],w[r][v]))))
    f = mul(t[r][r],h)
    return {'q0':q0,'L':l,'E':e,'F':f,'H':h,'G':g,'T':t,'W':w,
            'K1':sub(mul(add('1+2*eps',ALPHAS[r]),l),mul(ALPHAS[r],e)),
            'K2':mul('-2*(2-3*eps)',ALPHAS[r],f)}


def c4_symmetry(raw):
    data = validate_source(raw)
    def fingerprint(k):
        return tuple(sorted(((p['edge']-1+k)%4+1,p['slot'],i==data['marked_leg']-1) for i,p in enumerate(data['endpoints'])))
    stab = [k for k in range(4) if fingerprint(k)==fingerprint(0)]
    return {'orbit_multiplicity':4//len(stab),'stabilizer_order':len(stab),'stabilizer_rotations':stab,
            'basis':'three ordered contour endpoints with the marked leg'}


def alpha_measure():
    return {'original_simplex':{'variables':['alpha1','alpha2','alpha3'], 'constraint':'alpha1+alpha2+alpha3=1',
                                'nonnegative':True},
            'eliminated_variable':'alpha3','substitution':{'alpha3':'1-alpha1-alpha2'},
            'variables':['alpha1','alpha2'],'bounds':{'alpha1':['0','1'],'alpha2':['0','1-alpha1']},'jacobian':'1'}


def build_result(orbit, ibp_enabled=False):
    raw=orbit.source; data=validate_source(raw)
    if orbit != parse_tgv_pi1_orbit(raw):
        raise ValueError('tgv_pi1_metadata_disagrees_with_source')
    geom=geometry(data); r=data['marked_leg']-1
    weight='*'.join(f'({al})^({"-2*eps" if i==r else "-eps"})' for i,al in enumerate(ALPHAS))
    measure=alpha_measure()
    order=data['order']+measure['variables']
    j={'representation':'ordered_s_alpha_simplex','variables':order,
       'bounds':{**data['bounds'],**measure['bounds']},'iterated_order_outer_to_inner':order,
       'source_bounds_array':deepcopy(raw['bounds']),'source_order_outer_to_inner':data['order'],
       'jacobian':'1','alpha_measure':measure,'alpha_weight':weight,
       'integrand':f'({weight})*(({geom["K1"]})/({geom["q0"]}+I0)^(2-3*eps)+({geom["K2"]})/({geom["q0"]}+I0)^(3-3*eps))',
       'kernel_terms':[{'numerator':geom['K1'],'denominator_power':'2-3*eps'},
                       {'numerator':geom['K2'],'denominator_power':'3-3*eps'}],
       'geometry_polynomials':geom,'dimension':'4-2*eps',
       'power_prescription':'common (q0+I0) boundary value; continue jointly with C'}
    components={'wilson_color':'CA*CF','self_energy_color':SE_COLOR,'coupling':'(alpha_s/Pi)^3',
                'gamma_fourier':'Gamma[1-eps]^2*Gamma[2-3*eps]/Gamma[2-2*eps]',
                'rational':f'{three_slot_trace_phase(raw)}/(64*eps*(1+eps))','scale':'(Pi*mu^2*(-2*v14+I0))^(3*eps)'}
    c={'components':components,'factor_order':list(components),'component_semantics':'multiplicative',
       'combined':'*'.join(f'({v})' for v in components.values()),
       'phase_ledger':{'source_external_sign':'+1','dyson_coefficient':'1','wilson_and_trace':'I^3*(-I/2)=-1/2',
                       'coordinate_legs':'-1','tree_vertex_derivative':'-I','central_wick':'-I','net':str(three_slot_trace_phase(raw))}}
    result={'schema_version':SCHEMA,'orbit_id':raw['orbit_id'],'emission_type':raw['emission_type'],
            'dispatch':'tgv_pi1','loop_order':3,
            'source_sha256':hashlib.sha256(json.dumps(raw,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest(),
            'representative':{'scheme':'single_source_representative','rotation':0,
                              'source_topology':deepcopy(raw),'canonical_topology':deepcopy(raw),'marked_leg':r+1},
            'geometry':{'n_lines':4,'z_definition':'z1=v12/v14','reference_scale':'v14',
                        'A_ref':'-2*v14','vij_definition':'v_i dot v_j','lightlike_edges':True,
                        'dimensionless_gram':'u_i.u_j/A0; A0=-2*v14'},
            'operator':{'name':'Pi1','subgraph_loop_order':1,'skeleton_loop_order':2,'insertion_count':1,
                        'kind':'1pi','renormalization':'bare','dimension':'4-2*eps','gauge':'Feynman',
                        'sign_convention':SIGN,'scalar_contains_color':True,'tensor':'k^2*g-k*k','epsilon_precision':'exact'},
            'color_contraction':{'wilson_trace':'-I*CA*CF/2','self_energy_color':SE_COLOR,
                                 'full_color':f'CA*CF*({SE_COLOR})','channels':['CF*CA^2','CF*CA*TF*nf'],'projection':'none'},
            'c4_symmetry':c4_symmetry(raw),
            'normalization':{'scheme':'montopar-raw-CJ','loop_order':3,'color_projection':'none',
                             'euler_compensation_count':0,'orbit_multiplicity_included':False,'crossing_assembled':False,
                             'extra_factorial':'none','normalized_wilson_trace':'Tr/N'},
            'C':c,'J_raw':j,'J_raw_s':deepcopy(j),
            'branch_prescription':{'dimensionful_distance':'-X^2+I0','momentum':'-l^2-I0',
                                   'starting_region':{'z1':'>0','v14':'<0'},
                                   'power_convention':'common Feynman boundary value; joint analytic continuation',
                                   'complex_continuation_verified':False},
            'translation_status':'translated','evaluation_status':'not_evaluated','formal':'not_run'}
    if ibp_enabled:
        from montopar.tgv_pi1_ibp import build_ibp
        result['J_ibp']=build_ibp(data,geom,weight,c['combined'])
    return annotate_result(result, raw)


def validate_result(result):
    if result.get('schema_version') != SCHEMA or result.get('dispatch') != 'tgv_pi1':
        raise ValueError('regenerate_tgv_pi1_result')
    expected=build_result(parse_tgv_pi1_orbit(result['representative']['source_topology']),ibp_enabled='J_ibp' in result)
    if json.loads(json.dumps(expected)) != result:
        raise ValueError('tgv_pi1_result_does_not_reconstruct_from_source')
