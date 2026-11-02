"""One bare two-loop 1PI gluon self energy on two Wilson endpoints.

Full D-dimensional Feynman-gauge scalar, no CAS runtime dependency.
The incidence adapter reuses Pi1 source geometry validation only; no Pi1
operator or prefactor is used to compute Pi2. Raw outputs retain all colors.
"""
from montopar.sign_rules import real_phase, annotate_result, validate_contract
from montopar.sign_rules import real_phase, annotate_result, validate_contract
from copy import deepcopy
import json
import re

from montopar.models import Bound, EdgeDistribution, Orbit, Propagator, Vertex
from montopar.self_energy import validate_source as _validate_incidence
from montopar.self_energy import geometry_for_line, c4_symmetry as _c4
from montopar.self_energy_chain import GAMMA_RATIO
from montopar.j_factor import expand_polynomial

FAMILY = '2e-Pi2'
SCHEMA = '1.6'
R_RATIO = 'Gamma[1+2*eps]*Gamma[1-2*eps]^2/(Gamma[1+eps]^2*Gamma[1-eps]*Gamma[1-3*eps])'
COEFFICIENTS = {
    'AA': {
        'operator_monomial': 'CA^2', 'full_monomial': 'CA^2*CF',
        'A': '(2*eps^2-4*eps-1)/(2*eps^3*(2*eps-3))',
        'B': '(2*eps-1)*(6*eps^5-eps^4-44*eps^3+79*eps^2-30*eps-2)/(2*eps^3*(eps+1)*(2*eps-3)*(3*eps-2)*(3*eps-1))'},
    'Af': {
        'operator_monomial': 'CA*TF*nf', 'full_monomial': 'CA*CF*TF*nf',
        'A': '2*(eps-1)*(2*eps^2-3*eps+2)/(eps^3*(2*eps-3))',
        'B': '-2*(eps-1)*(2*eps-1)*(2*eps^4-11*eps^3-10*eps^2+17*eps-4)/(eps^3*(eps+1)*(2*eps-3)*(3*eps-2)*(3*eps-1))'},
    'Ff': {
        'operator_monomial': 'CF*TF*nf', 'full_monomial': 'CF^2*TF*nf',
        'A': '-4*(eps-1)*(2*eps^2-eps+2)/(eps^3*(2*eps-3))',
        'B': '8*(eps-1)*(2*eps-1)^2*(eps^2-2*eps+2)/(eps^3*(2*eps-3)*(3*eps-2)*(3*eps-1))'},
}


def coefficient_table():
    return {name: dict(item, p=f"({item['A']})+({item['B']})*({R_RATIO})",
                       U_times_p=f"({item['A']})*({GAMMA_RATIO})+({item['B']})")
            for name, item in COEFFICIENTS.items()}


def scalar_polynomial():
    return '+'.join(f"({v['operator_monomial']})*({v['p']})" for v in coefficient_table().values())


def _identifier(value):
    return isinstance(value, str) and re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', value) is not None


def _validate_source(raw):
    validate_contract(raw)
    required = {'orbit_id','emission_type','sector','nloops','vertices','propagators','internal_vertices','phases','bounds'}
    if not isinstance(raw, dict) or not required <= set(raw) or not set(raw) <= required | {'g_power','iterated_order_outer_to_inner','physical_convention'}:
        raise ValueError('pi2_unknown_or_missing_source_fields')
    if raw['emission_type'] != FAMILY or not isinstance(raw['orbit_id'], str) or not re.fullmatch(r'2e-Pi2-[0-9]+', raw['orbit_id']):
        raise ValueError('pi2_family_or_id_mismatch')
    if type(raw['nloops']) is not int or raw['nloops'] != 3 or raw['sector'] != 'dg':
        raise ValueError('pi2_requires_total_three_loops')
    if 'g_power' in raw and (type(raw['g_power']) is not int or raw['g_power'] != 6):
        raise ValueError('pi2_requires_g6')
    for key, count in [('vertices',2),('propagators',2),('internal_vertices',1),('phases',1),('bounds',2)]:
        if not isinstance(raw[key], list) or len(raw[key]) != count or any(not isinstance(v,dict) for v in raw[key]):
            raise ValueError('pi2_source_cardinality_or_shape')
    iv=raw['internal_vertices'][0]
    if set(iv) != {'type','propagator','se_type'} or iv['type'] != 'se' or iv['se_type'] != 'vacuum-pol-2loop' or not _identifier(iv['propagator']):
        raise ValueError('pi2_requires_one_bare_two_loop_1pi_operator')
    for v in raw['vertices']:
        if set(v) != {'line','lorentz','color'} or not all(_identifier(v[k]) for k in ['lorentz','color']):
            raise ValueError('pi2_malformed_endpoint')
    for prop in raw['propagators']:
        if set(prop) != {'loop','lorentz','color'} or not _identifier(prop['loop']):
            raise ValueError('pi2_malformed_propagator')
        if any(not isinstance(prop[k],list) or len(prop[k])!=2 or not all(_identifier(v) for v in prop[k]) for k in ['lorentz','color']):
            raise ValueError('pi2_malformed_propagator_indices')
    if any(set(b) != {'var','low','high'} or not all(isinstance(b[k],str) for k in ['var','low','high']) for b in raw['bounds']):
        raise ValueError('pi2_malformed_bounds')
    phase=raw['phases'][0]
    if set(phase) != {'loop','expression'} or not _identifier(phase['loop']) or not isinstance(phase['expression'], str):
        raise ValueError('pi2_malformed_phase')
    # Geometry-only adapter: one two-port insertion has the same endpoint,
    # delta-bridge, momentum and domain incidence at either subgraph loop order.
    geometry_input=deepcopy(raw)
    geometry_input.update(emission_type='2e-Pi1',nloops=2)
    if 'g_power' in geometry_input:geometry_input['g_power']=4
    geometry_input['internal_vertices'][0]['se_type']='vacuum-pol'
    data=_validate_incidence(geometry_input)
    for line in data['lines']:line['kind']='self_energy_pi2'
    return data


def validate_source(raw):
    try:return _validate_source(raw)
    except (KeyError,TypeError,IndexError,AttributeError) as exc:
        raise ValueError('pi2_malformed_source') from exc


def parse_pi2_orbit(raw):
    data=validate_source(raw)
    return Orbit(orbit_id=raw['orbit_id'],emission_type=FAMILY,n_gluons=0,pairings=[],
        bounds=[Bound(b['var'],b['low'],b['high']) for b in raw['bounds']],
        edge_distribution=EdgeDistribution({v['edge']:1 for v in data['endpoints']},{v['edge']:[v['s_var']] for v in data['endpoints']}),
        vertices=[Vertex(**v) for v in raw['vertices']],
        propagators=[Propagator(p['loop'],tuple(p['lorentz']),tuple(p['color'])) for p in raw['propagators']],
        source=deepcopy(raw),dispatch='self_energy_pi2')


def c4_symmetry(raw):
    result=_c4(raw,validate_source(raw))
    result['definition']='C4 on the two endpoints of one two-loop 1PI operator; no extra factorial'
    return result


def build_result(orbit):
    raw=orbit.source;data=validate_source(raw)
    if orbit != parse_pi2_orbit(raw):raise ValueError('pi2_orbit_metadata_disagrees_with_source')
    se=geometry_for_line(data['lines'][0]);f,d,u,w=(se[k] for k in ['f','d','u','w'])
    numerator=expand_polynomial(f'(1+4*eps)*({d})*({f})+(1-3*eps)*({u})*({w})')
    j={'representation':'ordered_s','variables':data['order'],'bounds':data['bounds'],
       'iterated_order_outer_to_inner':data['order'],'source_bounds_array':deepcopy(raw['bounds']),
       'order_basis':'ascending contour edges; original ordered Pi2 source','jacobian':'1',
       'integrand':f'({numerator})*({f})^(-2+3*eps)/(2*(1+2*eps))',
       'se_numerator':numerator,'geometry_polynomials':[se],
       'power_prescription':'each f power is the common Feynman boundary value (f+I0)^power; continue jointly with C'}
    table=coefficient_table();p2=scalar_polynomial()
    color={'wilson_trace':'CF','self_energy_scalar':p2,'full_color':f'CF*({p2})',
           'flavor_channels':[{'channel':k,'monomial':v['full_monomial'],'coefficient':v['p']} for k,v in table.items()],
           'ordered_basis_coordinates':[
               {'flavor_power':0,'basis':['CA^2*CF','CA*CF^2','CF^3'],'coordinates':[table['AA']['p'],'0','0']},
               {'flavor_power':1,'flavor_factor':'TF*nf','basis':['CA*CF','CF^2'],'coordinates':[table['Af']['p'],table['Ff']['p']]},
               {'flavor_power':2,'flavor_factor':'(TF*nf)^2','basis':['CF'],'coordinates':['0']}],
           'TF':'1/2','nf':'number of massless quark flavours','projection':'none',
           'method':'normalized Tr/N gives CF; retain all three 1PI channels including Ff second coordinate'}
    components={'wilson_color':'CF','self_energy_scalar':p2,'coupling':'(alpha_s/Pi)^3',
        'gamma_fourier':'(Pi^eps*Gamma[1-eps])^3','gamma_residual':GAMMA_RATIO,
        'rational':f'{real_phase(2, 3*2, 3, 2, 1, 2)}/(32*(1-2*eps)^2)','scale':'(-2*v14+I0)^(3*eps)'}
    c={'components':components,'factor_order':list(components),'component_semantics':'multiplicative',
       'combined':'*'.join(f'({v})' for v in components.values()),
       'phase_ledger':{'source_external_sign':'+1','wilson_vertices':'I^2','two_propagators':'(-I)^2',
           'Pi2_scalar_phase_relative_to_P2':'-I','momentum_prefactor_phase_relative_to_P2':'-I',
           'q_conversion':'-1','Fourier_phase':'I','coordinate_sign_before_scale_relative_to_P2':'-1',
           'vref_over_Aref':'-1/2','overall_dimensionless_sign_relative_to_P2':'1','extra_factorial':'none'}}
    result={'schema_version':SCHEMA,'orbit_id':raw['orbit_id'],'emission_type':FAMILY,'dispatch':'self_energy_pi2','loop_order':3,
       'geometry':{'n_lines':4,'z_definition':'z1=v12/v14','reference_scale':'v14','vij_definition':'v_i dot v_j',
          'A_ref':'-2*v14','lightlike_edges':True,'line_invariants':'f=X^2/(2*v14); d=vi.vj/v14; u=vi.X/v14; w=vj.X/v14'},
       'representative':{'scheme':'single_source_representative','rotation':0,'source_topology':deepcopy(raw),
          'canonical_topology':deepcopy(raw),'effective_lines':data['lines'],'internal_delta':data['internal_delta']},
       'operator':{'name':'Pi2','kind':'1pi','insertion_count':1,'subgraph_loop_order':2,
          'renormalization':'bare','dimension':'4-2*eps','gauge':'Feynman','xi_RUVV':1,'xi_DOT':0,
          'source':'arXiv:1703.08532v2 ancillary gluonB.h a^2; G10 from FORCER reduceG',
          'source_url':'https://arxiv.org/src/1703.08532v2/anc/gluonB.h',
          'source_sha256':'7d33d4f63e9a776a736c1d75a22376a0cb8f97b68fe4cb3f3ce9f8edea056ee3',
          'theory':'docs/self-energy-pi2.md, P2.1-P2.20',
          'scalar':f'-I*g^4/(4*Pi)^(4-2*eps)*(Gamma[1+eps]*Gamma[1-eps]^2/((1-2*eps)*Gamma[1-2*eps]))^2*({p2})*(-l^2-I0)^(-2*eps)',
          'tensor':'delta_ab*(l^2*g_mu_nu-l_mu*l_nu)','P2':p2,'R_G10':R_RATIO,'U':GAMMA_RATIO,
          'channels':table,'epsilon_precision':'exact','not_Pi1_squared':True},
       'color_contraction':color,'c4_symmetry':c4_symmetry(raw),
       'normalization':{'scheme':'montopar-raw-CJ','loop_order':3,'color_projection':'none','euler_compensation_count':0,
          'orbit_multiplicity_included':False,'crossing_assembled':False,'normalized_wilson_trace':'Tr/N','extra_factorial':'none',
          'gamma_residual':GAMMA_RATIO,'gamma_residual_Euler_linear_coefficient':'0',
          'shared_I_convention':{'status':'definition_only_not_prepared','operator_factor_B2':f'4*eps^2*({p2})',
             'external_factor':'(alpha_s/Pi)^3*CF*B2','common_R':'-G^3*U/(128*eps^2*(1-2*eps)^2)',
             'G':'Exp[-EulerGamma*eps]*Gamma[1-eps]','I':'(m/2)*common_R*(J(z1)+z1^(3*eps)*J(1/z1))/(1+z1^(3*eps))',
             'a1_projection':'separate from shared I; Ff remains the second flavor-linear coordinate'}},
       'C':c,'J_raw':j,'J_raw_s':deepcopy(j),
       'branch_prescription':{'dimensionful_distance':'-X^2+I0','momentum':'-l^2-I0',
          'starting_region':{'z1':'>0','v14':'<0','eps':'0<Re(eps)<1/3'},
          'power_convention':'common Feynman boundary value; joint analytic continuation','complex_continuation_verified':False},
       'translation_status':'translated','evaluation_status':'not_evaluated','formal':'not_run'}
    return annotate_result(result, raw)


def validate_result(result):
    try:
        if result.get('schema_version') != SCHEMA or result.get('dispatch') != 'self_energy_pi2':
            raise ValueError('regenerate_self_energy_pi2_result')
        expected=build_result(parse_pi2_orbit(result['representative']['source_topology']))
        if json.loads(json.dumps(expected)) != result:
            raise ValueError('pi2_result_does_not_reconstruct_from_source')
    except (KeyError,TypeError,AttributeError) as exc:
        raise ValueError('pi2_malformed_result') from exc
