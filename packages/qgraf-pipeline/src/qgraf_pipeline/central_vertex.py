from .sign_rules import attach_contract
"""Complete effective one-loop central three-gluon source, with explicit momentum conservation."""
from collections import Counter
from copy import deepcopy
import hashlib
import json
from .gamma1_operator import OPERATOR

FAMILY = '3e-Gamma1'
CONVENTION = 'wilson-central-three-gluon-1pi/1.0'

def operator_hash():
    return hashlib.sha256(json.dumps(OPERATOR,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def is_gamma1(orbit):
    return (orbit.emission_type==FAMILY and orbit.internal=={'Gamma1':1}
            and orbit.eikonal_count==3 and orbit.g_power==6)

def _from_edges(orbit_id, edges):
    from .definition_exporter import _coord_to_origin_expr, _negate_phase
    if len(edges)!=3 or any(type(e) is not int or e not in range(1,5) for e in edges) or edges!=sorted(edges) or len(set(edges))<2:
        raise ValueError('gamma1_invalid_contour_edges')
    slots=Counter();ends=[]
    for edge in edges:
        slots[edge]+=1;ends.append((edge,slots[edge]))
    rotations={tuple(sorted((((e-1+r)%4+1,k) for e,k in ends))) for r in range(4)}
    if len(rotations)!=4:raise ValueError('gamma1_unexpected_stabilizer')
    return attach_contract(dict(
        orbit_id=orbit_id,emission_type=FAMILY,sector='dg',nloops=3,g_power=6,
        vertices=[dict(line=e,lorentz=f'mu{j}',color=f'a{j}') for j,(e,k) in enumerate(ends,1)],
        propagators=[dict(loop=f'l{j}',lorentz=[f'mu{j}',f'nu{j}'],color=[f'a{j}',f'b{j}']) for j in range(1,4)],
        internal_vertices=[dict(type='gamma1',id='v1',loops=['l1','l2','l3'],
          lorentz=['nu1','nu2','nu3'],color=['b1','b2','b3'],
          incoming_momenta=['-l1','-l2','-l3'],vertex_loop_order=1,
          operator_id=OPERATOR['id'],operator_sha256=operator_hash())],
        phases=[dict(loop=f'l{j}',expression=_negate_phase(_coord_to_origin_expr(e,k))) for j,(e,k) in enumerate(ends,1)],
        bounds=[dict(var=f's{e}{k}',low='0',high='1' if k==slots[e] else f's{e}{k+1}') for e,k in ends],
        iterated_order_outer_to_inner=[f's{e}{k}' for e in range(1,5) for k in range(slots[e],0,-1)],
        momentum_constraints=[dict(type='delta',expression='l1+l2+l3',normalization='(2*pi)^D*delta^D(l1+l2+l3)')],
        sign_convention=CONVENTION,
        vertex_correction=dict(kind='1pi',operator='Gamma1',subgraph_loop_order=1,skeleton_loop_order=2,
          total_loop_order=3,independent_external_momenta=2,external_momentum_labels=['l1','l2','l3'],
          external_leg_self_energies_included=False,counterterm_included=False,
          internal_loop_measure_already_in_operator=True,source_external_sign='+1',
          gauge='Feynman',dimension='D=4-2*eps',renormalization='bare',
          full_color_channels=['CF*CA^2','CF*CA*TF*nf'],operator_definition=deepcopy(OPERATOR)),
        source_normalization=dict(normalized_wilson_trace='Tr/N',euler_compensation_count=0,
          color_projection='none',orbit_multiplicity_included=False,crossing_assembled=False),
        c4_symmetry=dict(orbit_multiplicity=4,stabilizer_order=1),
        downstream_status=dict(momentum_source='complete_in_stated_operator_convention',
          montopar='not_implemented',C_J='not_available',I='not_solved')))

def validate_definition(raw):
    if not isinstance(raw,dict) or raw.get('emission_type')!=FAMILY:
        raise ValueError('gamma1_wrong_family')
    try:
        edges=[v['line'] for v in raw['vertices']]
        expected=_from_edges(raw['orbit_id'],edges)
    except (KeyError,TypeError) as error:
        raise ValueError('gamma1_missing_source_fields') from error
    if raw!=expected:
        raise ValueError('gamma1_source_contract_mismatch')
    return raw

def build_definition(orbit, info):
    from .topology import vert_kind
    if not is_gamma1(orbit):raise ValueError('unsupported_central_vertex')
    if (orbit.representative is None or
        Counter(vert_kind(v) for v in orbit.representative.verts)!={'eikonal':3,'three_gluon':1} or
        len(orbit.representative.props)!=3 or info.n_wl_wl_gluons!=0 or
        info.internal_connections!=[(1,2,3)] or info.internal_gluon_links):
        raise ValueError('gamma1_invalid_source_topology')
    edges=[e for e in range(1,5) for _ in range(info.edge_counts[e-1])]
    result=_from_edges(orbit.orbit_id,edges)
    validate_definition(result)
    return result

def emit_integral(orbit, info):
    from .emitter import _emit_measure,_emit_eikonal_chain,_emit_coord_to_origin
    source=build_definition(orbit,info)
    chain,_=_emit_eikonal_chain(info.edge_counts)
    ends=[];slots=Counter()
    for v in source['vertices']:
        edge=v['line'];slots[edge]+=1;ends.append((edge,slots[edge]))
    lines=[rf'\mathcal A_{{({orbit.orbit_id})}} = '+_emit_measure(info.edge_counts)+
       r'\;\frac{1}{N_c}\operatorname{Tr}\bigl['+chain+r'\bigr]',
       r'\times\left[\prod_{j=1}^{3}\int\frac{d^D\ell_j}{(2\pi)^D}\right]'
       r'(2\pi)^D\delta^{(D)}(\ell_1+\ell_2+\ell_3)',
       r'\times \mathcal V^{(1)\,b_1b_2b_3}_{\nu_1\nu_2\nu_3}(-\ell_1,-\ell_2,-\ell_3)']
    for j,ep in enumerate(ends,1):
        lines.append(rf'\times\frac{{-ig_{{\mu_{j}\nu_{j}}}\delta^{{a_{j}b_{j}}}}}{{\ell_{j}^2+i0}}'
                     rf'e^{{-i\ell_{j}\cdot\bigl('+_emit_coord_to_origin(*ep)+r'\bigr)}')
    lines.append(r'\quad\text{(bare 1PI; full D; }\mathcal V_1=i\Gamma_{\rm DOT}^{(1)}'
                 r'\text{; DOT hep-ph/9605348v3, (2.4), (3.4)--(3.9), (3.17)--(3.22); 1999 erratum)}')
    return '\n'.join(lines)
