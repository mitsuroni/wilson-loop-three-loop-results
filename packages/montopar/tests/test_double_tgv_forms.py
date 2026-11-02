"""Exact source connections and consumer contracts for optional initial forms."""
from copy import deepcopy
import json
from pathlib import Path
import re

import pytest

from montopar.double_tgv_forms import polynomial_fields
from montopar.double_tgv_output import build_double_tgv_result,validate_double_tgv_result
from montopar.j_factor import expand_polynomial as expand
from montopar.package_builder import build_initial_package
from montopar.parser import parse_orbit

ROWS=json.loads((Path(__file__).parent/'fixtures/double_tgv_orbits.json').read_text(encoding='utf-8'))


def row(n):
    return deepcopy(next(r for r in ROWS if r['orbit_id'].endswith(f'-{n:02}')))


def result(n,**kw):
    return build_double_tgv_result(parse_orbit(row(n)),forms_enabled=True,**kw)


def sub(expression,mapping):
    return re.sub(r'\b[A-Za-z_][A-Za-z0-9_]*\b',lambda m:'('+mapping[m[0]]+')' if m[0] in mapping else m[0],expression)


@pytest.mark.parametrize('n',range(1,27))
def test_compact_fields_reconstruct_every_raw_numerator_and_domain(n):
    out=result(n)
    block=out['J_compact']
    for key,value in polynomial_fields(block).items():
        assert expand('('+value+')-('+out['J_raw_s'][key]+')')=='0'
    for key in ('variables','bounds','iterated_order_outer_to_inner','alpha_measure','ordered_chains'):
        assert block[key]==out['J_raw_s'][key]
    assert block['combined_C']==out['C']['combined']
    if out['J_raw_s']['physical_zero']:
        assert block['integrand']=='0'
        assert block['geometric_integrand_before_color_screening']!='0'
    # Optional mode is additive, not a change of the old definition.
    without={k:v for k,v in out.items() if k not in ('J_compact','J_reduced')}
    assert without==build_double_tgv_result(parse_orbit(row(n)))


@pytest.mark.parametrize('n',[n for n in range(1,27) if n not in (4,19)])
def test_central_coefficients_and_projective_B_reconstruct_original(n):
    out=result(n)
    reduced=out['J_reduced']
    assert reduced['method']=='central_parameter_gauss'
    assert reduced['dimension_after']==7
    expressions={'P1':'0','P2':'0','B5':'dtg_h*dtg_Q0+dtg_c*dtg_Q1'}
    for term in reduced['terms']:
        k=term['source']['denominator_shift'];m=term['source']['central_power']
        key='P1' if k==2 else 'P2'
        expressions[key]+=f"+dtg_c^{m}*({term['coefficient']})"
        assert term['combined_C']==f"({out['C']['combined']}) * ({term['relative_C']})"
        assert term['iterated_order_outer_to_inner'][:4]==out['J_raw_s']['iterated_order_outer_to_inner']
        assert term['hypergeometric']['normalization']=='ordinary Gauss 2F1'
    variables=reduced['terms'][0]['variables']+['dtg_c']
    fields=polynomial_fields({'definitions':reduced['definitions'],'fields':expressions,'variables':variables})
    mapping={'alpha1':'dtg_a*dtg_r','alpha2':'dtg_a*(1-dtg_r)',
             'alpha3':'(1-dtg_a)*dtg_t','alpha4':'(1-dtg_a)*(1-dtg_t)','alpha5':'dtg_c'}
    for key,oldkey in [('P1','pnum1'),('P2','pnum2'),('B5','B5')]:
        expected=sub(out['J_raw_s'][oldkey],mapping)
        assert expand(f"({fields[key]})-({expected})")=='0'
    # The retained shared definitions cannot reference the integrated parameter.
    assert not any(re.search(r'\bdtg_c\b',d['expression']) for d in reduced['definitions'])
    assert reduced['equivalence']['finite_auxiliary_regulator_identity'] is False


@pytest.mark.parametrize('n',[4,19])
def test_same_edge_terms_keep_all_nine_boundaries_and_color(n):
    out=result(n);block=out['J_reduced']
    assert block['method']=='paired_same_edge_beta'
    assert block['dimension_after']==4
    assert len(block['terms'])==9
    assert [t['source']['boundary_coefficient'] for t in block['terms']]==[4,-2,-2,-2,1,1,-2,1,1]
    for term in block['terms']:
        assert term['combined_C'].startswith('('+out['C']['combined']+')')
        assert set(term['variables'])==set(term['bounds'])==set(term['iterated_order_outer_to_inner'])
    assert 'Gamma' not in str(block['definitions'])
    assert block['normalization']['euler_compensation_count']==0


def test_same_edge_criterion_is_geometric_not_orbit_id_or_edge_occupancy():
    source=row(4);source['orbit_id']='renamed-adjacent-pair'
    assert build_double_tgv_result(parse_orbit(source),forms_enabled=True)['J_reduced']['method']=='paired_same_edge_beta'
    assert result(6)['J_reduced']['method']=='central_parameter_gauss'
    source=row(4)
    loops=source['internal_vertices'][0]['loops']
    source['internal_vertices'][0]['loops']=[loops[1],loops[0],loops[2],loops[3]]
    changed=build_double_tgv_result(parse_orbit(source),forms_enabled=True)
    assert all(t['source']['orientation']==-1 for t in changed['J_reduced']['terms'])
    assert changed['C']['components']['color']=='(1/4)*CA^2*CF'


@pytest.mark.parametrize('mutation',['hessian','scope','domain','C','missing_compact','missing_reduced','regulator','contact','gauss'])
def test_corrupted_optional_forms_are_rejected(mutation):
    out=result(4 if mutation=='contact' else 25)
    if mutation=='hessian':out['J_compact']['fields']['T2_hat']='0'
    elif mutation=='scope':out['J_compact']['definitions'][0]['symbol']='alpha1'
    elif mutation=='domain':out['J_reduced']['terms'][0]['bounds']['s11']=['0','1/2']
    elif mutation=='C':out['J_reduced']['terms'][0]['relative_C']='1'
    elif mutation=='missing_compact':del out['J_compact']
    elif mutation=='missing_reduced':del out['J_reduced']
    elif mutation=='regulator':out['J_reduced']['equivalence']['finite_auxiliary_regulator_identity']=True
    elif mutation=='contact':out['J_reduced']['terms'].pop(0)
    else:out['J_reduced']['terms'][0]['hypergeometric']['c']='3-eps'
    with pytest.raises(ValueError):validate_double_tgv_result(out)


def test_polynomial_definition_consumer_rejects_forward_reference():
    block={'definitions':[{'symbol':'dtg_a','expression':'dtg_b+1'}],'fields':{'B5':'dtg_a'}}
    with pytest.raises(ValueError):polynomial_fields(block)


def test_forms_and_ibp_coexist_and_survive_initial_package(tmp_path):
    source=row(4);out=result(4,ibp_enabled=True)
    folder=tmp_path/'3loop'/'4e-3g+3g'/source['orbit_id'];folder.mkdir(parents=True)
    (folder/'montopar.json').write_text(json.dumps(out),encoding='utf-8')
    package,_=build_initial_package({'orbits':[source]},tmp_path,'4e-3g+3g')
    for key in ('J_compact','J_reduced','J_ibp'):
        assert package['orbits'][0][key]==out[key]
