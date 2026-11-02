"""Complete scalar block measures for Gamma1; shared by raw and contour IBP.

The finite-simplex form retains its max bound. Explicit rectangular sectors
are exported independently without changing the original Wilson s order.
"""
from copy import deepcopy
from montopar.gamma1_kernel import operators, apply_operator_compact
from montopar.tgv_pi1 import alpha_measure
from montopar.j_factor import substitute_s_to_x

ALPHAS=['alpha1','alpha2','(1-alpha1-alpha2)']
P='(alpha1*alpha2*(1-alpha1-alpha2))'
MAX_BOUND='Max[alpha1*t2*(1-t1-t2),alpha2*t1*(1-t1-t2),(1-alpha1-alpha2)*t1*t2]'

def block_specs(data):
    """Every block is one integral, with its own actual dimension and bounds."""
    alpha=alpha_measure()
    base_order=list(data['order'])+alpha['variables']
    base_bounds={**deepcopy(data['bounds']),**alpha['bounds']}
    out=[]
    for channel,ops in operators().items():
        for name,beta,weight in [
            ('N0','3*eps',f'-({P})^(-eps)*({MAX_BOUND})^(-1-eps)/(3*eps*(1+eps))'),
            ('L0','-1+3*eps',f'({P})^(-eps)*({MAX_BOUND})^(-eps)/eps')]:
            order=base_order+['t1','t2']
            out.append({'channel':channel,'block_id':name,'kind':'triangle','operator':ops[name],
                        'primitive_power':beta,'weight':weight,'variables':order,
                        'bounds':{**base_bounds,'t1':['0','1'],'t2':['0','1-t1']},
                        'iterated_order_outer_to_inner':order,'alpha_measure':deepcopy(alpha),
                        'internal_measure':{'variables':['t1','t2'],'eliminated':'t3',
                                            'substitution':'t3=1-t1-t2','jacobian':'1'}})
        for r,op in enumerate(ops['M']):
            if op=='0':continue
            out.append({'channel':channel,'block_id':f'M{r+1}','kind':'bubble','marked_bubble_leg':r+1,
                        'operator':op,'primitive_power':'-1+3*eps',
                        'weight':f'4*Beta[1-eps,1-eps]/eps*({P})^(-eps)*({ALPHAS[r]})^(-eps)',
                        'variables':list(base_order),'bounds':deepcopy(base_bounds),
                        'iterated_order_outer_to_inner':list(base_order),'alpha_measure':deepcopy(alpha),
                        'internal_measure':{'integrated':'Beta[1-eps,1-eps]',
                                            'origin':'integral_0^1 [x*(1-x)]^(-eps) dx'}})
    return out

def kernel_terms(coefficients,beta):
    if beta not in ('3*eps','-1+3*eps'):
        raise ValueError('gamma1_unknown_primitive_power')
    offset=0 if beta=='3*eps' else 1
    return [{'numerator':value,'denominator_power':f'{shift+offset}-3*eps','derivative_shift':shift}
            for shift,value in sorted(coefficients.items()) if value!='0']

def scalar_integrand(weight,terms,q0):
    body='+'.join(f'({t["numerator"]})/(({q0})+I0)^({t["denominator_power"]})' for t in terms) or '0'
    return f'({weight})*({body})'

def evaluate_block(block,geom):
    result=deepcopy(block)
    terms=kernel_terms(apply_operator_compact(block['operator'],block['primitive_power'],geom),block['primitive_power'])
    result.update({'kernel_terms':terms,'q0':geom['q0'],
                   'integrand':scalar_integrand(block['weight'],terms,geom['q0']),
                   'jacobian':'1','regulator_group':'gamma1_complete_parent',
                   'ambient_derivatives_evaluated':True,
                   'numerator_representation':'self_contained_factored_polynomial'})
    return result

def sectorize(term):
    """Split one already evaluated triangle term into its exact three sectors."""
    if term.get('integration_block_kind',term['kind'])!='triangle':return [deepcopy(term)]
    out=[]
    for k in range(3):
        others=[j for j in range(3) if j!=k]
        cs=['1']*3;cs[others[0]]='u';cs[others[1]]='v'
        Z='+'.join(f'({ALPHAS[i]})*({cs[i]})' for i in range(3))
        substitutions={f't{i+1}':f'(({ALPHAS[i]})*({cs[i]})/({Z}))' for i in range(2)}
        if term['block_id']=='N0':
            weight=f'-({P})^(-2*eps)*({Z})^(-1+2*eps)/(u*v)^(1+eps)/(3*eps*(1+eps))'
        elif term['block_id']=='L0':
            weight=f'({P})^(1-2*eps)*({Z})^(-3+2*eps)/(u*v)^eps/eps'
        else:raise ValueError('gamma1_unrecognized_triangle_block')
        result=deepcopy(term)
        if 'operator' in result:
            result['operator_before_sector_substitution']=result.pop('operator')
        terms=[dict(t,numerator=substitute_s_to_x(t['numerator'],substitutions)) for t in term['kernel_terms']]
        order=[v for v in term['variables'] if v not in ('t1','t2')]+['u','v']
        bounds={v:deepcopy(term['bounds'][v]) for v in order if v not in ('u','v')}
        bounds.update({'u':['1','Infinity'],'v':['1','Infinity']})
        result.update({'term_id':term['term_id']+f'_sector{k+1}',
                       'sector':{'index':k+1,'Z':Z,'t_substitution':substitutions,
                                 'criterion':'alpha_k/t_k is maximal','jacobian':f'{P}/({Z})^3'},
                       'variables':order,'bounds':bounds,'iterated_order_outer_to_inner':order,
                       'weight':weight,'kernel_terms':terms,
                       'integrand':scalar_integrand(weight,terms,term['q0']),
                       'internal_measure':{'variables':['u','v'],'sector_index':k+1,
                                           'jacobian_included_in_weight':True}})
        if 'kernel_coefficients' in result:
            result['kernel_coefficients']={t['derivative_shift']:t['numerator'] for t in terms}
        out.append(result)
    return out
