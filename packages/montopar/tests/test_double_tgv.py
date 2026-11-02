"""Regression tests for the corrected double-TGV branch and source contract."""
from copy import deepcopy
from pathlib import Path
import json
import itertools
import re
import pytest
from montopar.parser import parse_orbit
from montopar.double_tgv import compute_numerators, compute_B5
from montopar.double_tgv_output import (build_double_tgv_result, validate_double_tgv_result,
    color_data, double_tgv_c4_symmetry, alpha_simplex5)
from montopar.j_factor import expand_polynomial as expand, substitute_s_to_x
from montopar.package_builder import build_initial_package

ROWS = json.loads((Path(__file__).parent/'fixtures/double_tgv_orbits.json').read_text(encoding='utf8'))
def raw(n):
    return deepcopy(next(r for r in ROWS if r['orbit_id'].endswith(f'-{n:02d}')))
def zero(expr):
    return expand(expr) == '0'

def test_04_has_exact_chain_and_scale_coefficients():
    orbit = parse_orbit(raw(4));num = compute_numerators(orbit.double_tgv_legs)
    b = compute_B5(orbit.double_tgv_legs)
    h = '-2*alpha5*z1^2*(alpha1-alpha2)*(alpha3-alpha4)'
    assert zero(f"({num['T2_hat']})-({h})")
    assert zero(f"({num['T1_hat']})+2*({b})*({h})")
    # Polynomial grammar supports rational factors, not expression division.
    target = f'-1/2*(1-3*eps)*alpha5*z1^2*(alpha1-alpha2)*(alpha3-alpha4)*({b})'
    assert zero(f"({num['pnum1']})+({b})*({num['pnum2']})-({target})")

def test_25_retains_dimensionally_regulated_trace():
    orbit = parse_orbit(raw(25));num = compute_numerators(orbit.double_tgv_legs)
    assert zero(f"({num['K_plus']})-((7-2*eps)*z1^2+10*z1+1)")
    assert zero(f"({num['K_minus']})-((2+2*eps)*z1^2+8*z1-1)")
    assert 'eps' in num['T2_hat']

@pytest.mark.parametrize('n,expected', [(4,'-1/4'),(5,'0'),(6,'1/4')])
def test_same_edge_distribution_does_not_determine_color(n, expected):
    assert color_data(parse_orbit(raw(n)))['ratio_to_CA2_CF'] == expected

def test_all_propagator_permutations_preserve_oriented_pairing():
    source=raw(25);ref=parse_orbit(source)
    for order in itertools.permutations(range(4)):
        candidate=deepcopy(source);candidate['propagators']=[source['propagators'][i] for i in order]
        actual=parse_orbit(candidate)
        assert actual.double_tgv_legs == ref.double_tgv_legs
        assert color_data(actual)==color_data(ref)

@pytest.mark.parametrize('order,alpha_map,sign', [([1,0,2,3],{'alpha1':'alpha2','alpha2':'alpha1'},-1),
    ([0,1,3,2],{'alpha3':'alpha4','alpha4':'alpha3'},-1),
    ([2,3,0,1],{'alpha1':'alpha3','alpha2':'alpha4','alpha3':'alpha1','alpha4':'alpha2'},1)])
def test_vertex_orientation_moves_color_and_numerator_together(order, alpha_map, sign):
    source=raw(25);ref=parse_orbit(source);base=compute_numerators(ref.double_tgv_legs)
    loops=source['internal_vertices'][0]['loops'];source['internal_vertices'][0]['loops']=[loops[i] for i in order]
    changed=parse_orbit(source);nums=compute_numerators(changed.double_tgv_legs)
    for key in ('T1_hat','T2_hat','pnum1','pnum2'):
        restored=re.sub(r'\balpha[1-5]\b', lambda m: alpha_map.get(m[0],m[0]), nums[key])
        assert zero(f"({restored})-({sign})*({base[key]})")
    assert color_data(changed)['ratio_to_CA2_CF'] == ('1/4' if sign==-1 else '-1/4')

@pytest.mark.parametrize('n,expected',[(4,4),(19,2),(20,2),(21,2),(25,2),(26,1)])
def test_c4_stabilizers_include_internal_vertex_exchange(n, expected):
    assert double_tgv_c4_symmetry(parse_orbit(raw(n)))['orbit_multiplicity']==expected

@pytest.mark.parametrize('n', [4,5,25])
def test_export_contains_raw_ordered_domain_and_multiplicative_C(n):
    out=build_double_tgv_result(parse_orbit(raw(n)))
    assert out['schema_version']=='1.3' and not any(k.endswith('_rotated') for k in out)
    c=out['C'];assert c['combined']==' * '.join(f"({c['components'][k]})" for k in c['factor_order'])
    assert 'I^6' not in str(c['components'])
    for key in ('J_raw_s','J_raw'):
        j=out[key];assert j['alpha_measure']==alpha_simplex5()
        assert set(j['variables'])==set(j['iterated_order_outer_to_inner'])
    if n==5:
        assert out['J_raw']['integrand']=='0' and out['J_raw']['B5']!='0'
        assert out['C']['components']['color']=='0'
    assert validate_double_tgv_result(out)

@pytest.mark.parametrize('mutation', ['old_version','euler','trace','measure','domain','ledger','rotated','color','order'])
def test_corrupted_exports_are_rejected(mutation):
    out=build_double_tgv_result(parse_orbit(raw(4)))
    if mutation=='old_version':out['schema_version']='1.0'
    elif mutation=='euler':out['normalization']['euler_compensation_count']=1
    elif mutation=='trace':out['J_raw_s']['pnum2']='0'
    elif mutation=='measure':out['J_raw_s']['delta_constraints']={}
    elif mutation=='domain':out['J_raw_s']['bounds']['s31']=['0','1']
    elif mutation=='ledger':out['C']['components']['v14_scale']='v14^2'
    elif mutation=='rotated':out['J_raw_rotated']={}
    elif mutation=='color':out['C']['components']['color']='CF^3'
    else:out['representative']['legs'].reverse()
    with pytest.raises(ValueError):validate_double_tgv_result(out)

@pytest.mark.parametrize('mutation', ['loops','phase','split','color'])
def test_ambiguous_sources_are_rejected(mutation):
    source=raw(4)
    if mutation=='loops':source['internal_vertices'][0]['loops'][1]='l1'
    elif mutation=='phase':source['phases'][0]['expression']='-(v1+v2+2*s31*v3)'
    elif mutation=='split':source['bounds'][0]['high']='1'
    else:source['propagators'][0]['color'][0]='unbound_color'
    with pytest.raises(ValueError):build_double_tgv_result(parse_orbit(source))

def test_explicit_two_vertices_preserve_cyclic_orientation():
    source=raw(25);ref=parse_orbit(source)
    source['internal_vertices']=[{'type':'tgv','id':'left','loops':['q','l1','l2']},
                                 {'type':'tgv','id':'right','loops':['l4','q','l3']}]
    assert parse_orbit(source).double_tgv_legs==ref.double_tgv_legs

def test_package_rejects_legacy_double_tgv_files(tmp_path):
    source=raw(4);folder=tmp_path/'3loop'/'4e-3g+3g'/source['orbit_id'];folder.mkdir(parents=True)
    (folder/'montopar.json').write_text(json.dumps({'dispatch':'double_tgv','C':{},'J_raw':{}}))
    with pytest.raises(ValueError):build_initial_package({'orbits':[source]},tmp_path,'4e-3g+3g')

def test_corrected_package_keeps_double_vertex_core(tmp_path):
    source=raw(4);out=build_double_tgv_result(parse_orbit(source))
    folder=tmp_path/'3loop'/'4e-3g+3g'/source['orbit_id'];folder.mkdir(parents=True)
    (folder/'montopar.json').write_text(json.dumps(out))
    package,_=build_initial_package({'orbits':[source]},tmp_path,'4e-3g+3g')
    assert package['schema_version']=='1.3'
    assert package['orbits'][0]['integral_structure']['core']['kind']=='double_tgv'
    assert package['orbits'][0]['J_raw_s']==out['J_raw_s']
