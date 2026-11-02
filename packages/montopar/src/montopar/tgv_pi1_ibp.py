"""Single-layer contour IBP of the complete TGV-Pi1 tensor.

Flux_i = b0_i/(1-3 eps) Q^(-1+3 eps) + b1_i Q^(-2+3 eps).
The longitudinal self-energy contribution is an exact marked-path derivative.
For a repeated edge the metric tree vertex is entirely tangential. For three
distinct edges a translation identity leaves an explicit recoil integral.
"""
from copy import deepcopy
from montopar.tgv_pi1 import ALPHAS, add, mul, sub, alpha_measure
from montopar.j_factor import differentiate_polynomial, substitute_s_to_x, expand_polynomial
from montopar.double_tgv_ibp import ordered_faces


def flux_data(data,geom):
    g,t=geom['G'],geom['T']; r=data['marked_leg']-1; s=(r+1)%3; v=(r+2)%3
    edges=[p['edge'] for p in data['endpoints']]
    if edges[1]==edges[2]:
        c=['0',mul('-1',g[0][1]),g[0][2]]; recoil='0'
    elif edges[0]==edges[1]:
        c=[mul('-1',g[0][2]),g[1][2],'0']; recoil='0'
    elif edges[0]==edges[2]:
        c=[g[0][1],'0',mul('-1',g[1][2])]; recoil='0'
    else:
        c=[g[1][2],g[2][0],g[0][1]]
        recoil=mul('4*(1+eps)',add(*(mul(g[i][j],ALPHAS[i],t[k][i]) for i,j,k in [(0,1,2),(1,2,0),(2,0,1)])))
    b0=[mul('1+eps',ci) for ci in c]; b1=['0']*3
    b0[r]=add(b0[r],mul('3/2-eps',g[s][v],sub(ALPHAS[s],ALPHAS[v])))
    b1[r]=mul('-1',geom['H'])
    # Coefficient-by-coefficient reconstruction, valid even for finite I0.
    div1=[];div2=[]
    for i,p in enumerate(data['endpoints']):
        var=p['s_var'];dq=differentiate_polynomial(geom['q0'],var)
        if differentiate_polynomial(b0[i],var) != '0':
            raise ValueError('tgv_pi1_flux0_must_be_path_constant')
        div1.append(sub(differentiate_polynomial(b1[i],var),mul(b0[i],dq)))
        div2.append(mul('-(2-3*eps)',b1[i],dq))
    if sub(geom['K1'],add(*div1,recoil)) != '0' or sub(geom['K2'],add(*div2)) != '0':
        raise ValueError('tgv_pi1_ibp_does_not_reconstruct_raw')
    return {'b0_numerators':b0,'b0_denominator':'1-3*eps','b1_numerators':b1,
            'recoil_numerator':recoil,'scalar_path_coefficients':c,
            'identity':'raw_kernel = sum_i d/ds_i flux_i + recoil',
            'coefficient_reconstruction':'exact polynomial identities at both denominator powers'}


def build_ibp(data,geom,weight,combined_c):
    flux=flux_data(data,geom);measure=alpha_measure();variables=[p['s_var'] for p in data['endpoints']]
    terms=[]
    restrict=lambda value,mapping:expand_polynomial(substitute_s_to_x(value,mapping))
    for face in ordered_faces(data['ordered_chains']):
        mapping=face['parent_to_face']
        n0=restrict(add(*(mul(face['normal'][var],flux['b0_numerators'][i]) for i,var in enumerate(variables))),mapping)
        n1=restrict(add(*(mul(face['normal'][var],flux['b1_numerators'][i]) for i,var in enumerate(variables))),mapping)
        if n0=='0' and n1=='0':
            continue
        q=restrict(geom['q0'],mapping)
        remaining=[v for v in data['order'] if v!=face['removed']]
        bounds={v:['0',chain[i+1] if i+1<len(chain) else '1'] for chain in face['ordered_chains'] for i,v in enumerate(chain)}
        order=remaining+measure['variables']
        terms.append({'term_id':f'ibp_{len(terms)+1:03d}','kind':face['kind'],
                      'origin':'complete_tensor_flux','face':face,'combined_C':combined_c,
                      'variables':order,'bounds':{**{v:bounds[v] for v in remaining},**measure['bounds']},
                      'iterated_order_outer_to_inner':order,'jacobian':'1','alpha_measure':deepcopy(measure),
                      'q0':q,'flux_numerators':[n0,n1],
                      'integrand':f'({weight})*(({n0})/(1-3*eps)/({q}+I0)^(1-3*eps)+({n1})/({q}+I0)^(2-3*eps))',
                      'regulator_group':'one_parent_common_boundary_value'})
    if flux['recoil_numerator']!='0':
        order=data['order']+measure['variables']
        terms.append({'term_id':f'ibp_{len(terms)+1:03d}','kind':'recoil','origin':'metric_tree_translation_remainder',
                      'combined_C':combined_c,'variables':order,'bounds':{**data['bounds'],**measure['bounds']},
                      'iterated_order_outer_to_inner':order,'jacobian':'1','alpha_measure':deepcopy(measure),
                      'numerator':flux['recoil_numerator'],
                      'integrand':f'({weight})*({flux["recoil_numerator"]})/({geom["q0"]}+I0)^(2-3*eps)',
                      'regulator_group':'one_parent_common_boundary_value'})
    return {'schema_version':'tgv-pi1-ibp-1.0','representation':'ordered_face_sum_with_recoil',
            'flux':flux,'terms':terms,'term_count':len(terms),'common_C':combined_c,
            'equality':'raw C*J = sum_t combined_C_t * integral_t',
            'conditions':{'first':'finite common I0=I*eta and alpha interior cutoff, then ordinary contour IBP',
                          'boundary_value':'combine terms before joint regulator removal and meromorphic continuation',
                          'individual_term_convergence_proved':False},
            'evaluation_status':'not_evaluated'}
