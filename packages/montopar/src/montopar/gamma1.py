from montopar.sign_rules import three_slot_trace_phase, annotate_result
"""Central one-loop three-gluon vertex: complete full-D two-colour C/J."""
from copy import deepcopy
import hashlib,json
from montopar.gamma1_source import validate_source,parse_gamma1_orbit,c4_symmetry
from montopar.gamma1_kernel import geometry,operators
from montopar.gamma1_blocks import block_specs,evaluate_block,sectorize,MAX_BOUND

SCHEMA='1.8'
THEORY_ID='Gamma1-gram-free-CJ-r002'
THEORY_DATA_SHA256='032dfb99302ddc56ce83181b770a435ca66126c8005bd24adc29784eb95869d2'

def build_result(orbit,ibp_enabled=False):
    raw=orbit.source;data=validate_source(raw)
    if orbit!=parse_gamma1_orbit(raw):raise ValueError('gamma1_metadata_disagrees_with_source')
    geom=geometry(data)
    components={'wilson_color':'CA*CF','coupling':'(alpha_s/Pi)^3',
                'gamma_fourier':'Gamma[1-3*eps]','rational':f'{three_slot_trace_phase(raw)}/128',
                'scale':'(Pi*mu^2*(-2*v14+I0))^(3*eps)'}
    common='*'.join(f'({v})' for v in components.values())
    channel_colors={'g':'CA','q':'TF*nf'}
    channels={ch:{'color':color,'full_color':'CF*CA^2' if ch=='g' else 'CF*CA*TF*nf',
                  'combined':f'({common})*({color})',
                  'bare_coupling_form':f'{three_slot_trace_phase(raw)}*g0^6*({"CF*CA^2" if ch=="g" else "CF*CA*TF*nf"})*Gamma[1-3*eps]*(-2*v14+I0)^(3*eps)/(8192*Pi^(6-3*eps))'}
              for ch,color in channel_colors.items()}
    c={'components':components,'factor_order':list(components),'component_semantics':'multiplicative',
       'combined':common,'channels':channels,'combined_semantics':'common factor; each term includes its internal colour',
       'coupling_dictionary':'g0=mu^eps*g_b; alpha_s=g_b^2/(4*Pi); scale mu is counted once',
       'phase_ledger':{'source':'V1=i*Gamma_DOT','wilson_and_propagators':'I^3*(-I)^3=1',
                       'normalized_color_trace':'-I*CF*CA/2','canonical_internal_loop':'1/I',
                       'triangle_wick_fourier':'I^3*(-1)^6*(-I)^3=1',
                       'bubble_wick_fourier':'I^3*(-1)^5*(-I)=1','net_rational_color_trace':'-1/2',
                       'DOT_eta_extra_factor':False}}
    terms=[]
    for block in block_specs(data):
        term=evaluate_block(block,geom)
        term.update({'term_id':f'raw_{block["channel"]}_{block["block_id"]}',
                     'combined_C':channels[block['channel']]['combined'],
                     'color_coefficient':channel_colors[block['channel']]})
        terms.append(term)
    j={'representation':'ordered_channel_block_sum','terms':terms,'term_count':len(terms),
       'source_bounds_array':deepcopy(raw['bounds']),'source_order_outer_to_inner':deepcopy(data['order']),
       'geometry_polynomials':geom,'max_bound':MAX_BOUND,'operators':operators(),
       'dimension':'4-2*eps','equality':'amplitude=sum_terms combined_C * integral(integrand)',
       'common_C_form':'amplitude=C.combined*sum_terms color_coefficient*integral(integrand)',
       'power_prescription':'common (q0+I0) boundary value; complete sum continued jointly',
       'stage3_consumption':'requires explicit channel/block handling; no single-integrand shortcut'}
    sector_terms=[t for term in terms for t in sectorize(term)]
    result={'schema_version':SCHEMA,'orbit_id':raw['orbit_id'],'emission_type':raw['emission_type'],
            'dispatch':'gamma1','loop_order':3,
            'source_sha256':hashlib.sha256(json.dumps(raw,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest(),
            'representative':{'scheme':'single_source_representative','rotation':0,
                              'source_topology':deepcopy(raw),'canonical_topology':deepcopy(raw)},
            'geometry':{'n_lines':4,'z_definition':'z1=v12/v14','reference_scale':'v14','A_ref':'-2*v14',
                        'vij_definition':'v_i dot v_j','lightlike_edges':True},
            'operator':{'name':'Gamma1','subgraph_loop_order':1,'skeleton_loop_order':2,'kind':'1pi',
                        'renormalization':'bare','dimension':'4-2*eps','gauge':'Feynman',
                        'input_operator_id':raw['internal_vertices'][0]['operator_id'],
                        'input_operator_sha256':raw['internal_vertices'][0]['operator_sha256'],
                        'parametric_theory':THEORY_ID,'theory_data_sha256':THEORY_DATA_SHA256,
                        'epsilon_precision':'exact','external_leg_Pi1_included':False,'counterterm_included':False,
                        'basis_connection':'all four generic-D scalar-master coefficients match DOT in each colour channel'},
            'normalization':{'scheme':'montopar-raw-CJ','loop_order':3,'color_projection':'none',
                             'euler_compensation_count':0,'orbit_multiplicity_included':False,'crossing_assembled':False,
                             'extra_factorial':'none','normalized_wilson_trace':'Tr/N'},
            'color_contraction':{'wilson_trace':'-I*CA*CF/2','channels':['CF*CA^2','CF*CA*TF*nf'],
                                 'projection':'none','distinct_channel_integrals':True},
            'c4_symmetry':c4_symmetry(raw),'C':c,'J_raw':j,'J_raw_s':deepcopy(j),
            'J_sector':{'representation':'ordered_channel_rectangular_sector_sum','terms':sector_terms,
                        'term_count':len(sector_terms),'equality':'same complete amplitude as J_raw',
                        'sectors_are_domain_partitions_not_graph_multiplicities':True},
            'branch_prescription':{'dimensionful_distance':'-X^2+I0','momentum':'p^2+I0',
                                   'starting_region':{'z1':'>0','v14':'<0'},
                                   'power_convention':'common Feynman boundary value; joint meromorphic continuation',
                                   'complex_continuation_verified':False},
            'translation_status':'translated','evaluation_status':'not_evaluated','formal':'not_run'}
    if ibp_enabled:
        from montopar.gamma1_ibp import build_ibp
        result['J_ibp']=build_ibp(data,geom,{ch:info['combined'] for ch,info in channels.items()})
        result['J_ibp']['sector_terms']=[t for term in result['J_ibp']['terms'] for t in sectorize(term)]
        result['J_ibp']['sector_term_count']=len(result['J_ibp']['sector_terms'])
        result['J_ibp']['sector_equality']='same complete IBP sum, with each triangle t simplex partitioned once'
    return annotate_result(result, raw)

def validate_result(result):
    if result.get('schema_version')!=SCHEMA or result.get('dispatch')!='gamma1':
        raise ValueError('regenerate_gamma1_result')
    try:
        expected=build_result(parse_gamma1_orbit(result['representative']['source_topology']),ibp_enabled='J_ibp' in result)
    except (KeyError,TypeError) as error:raise ValueError('gamma1_result_missing_source') from error
    if json.loads(json.dumps(expected))!=result:
        raise ValueError('gamma1_result_does_not_reconstruct_from_source')
