from copy import deepcopy
from pathlib import Path
import json, subprocess, sys
import pytest
from montopar.parser import parse_orbit
from montopar.self_energy_pi2 import build_result, validate_result, validate_source, c4_symmetry
from montopar.package_builder import build_initial_package

RAWS=json.loads((Path(__file__).parent/'fixtures/self_energy_pi2_orbits.json').read_text(encoding='utf-8'))['orbits']

@pytest.mark.parametrize('raw',RAWS,ids=lambda r:r['orbit_id'])
def test_full_pi2(raw):
    result=json.loads(json.dumps(build_result(parse_orbit(raw))));validate_result(result)
    assert result['schema_version']=='1.6' and result['dispatch']=='self_energy_pi2'
    assert result['operator']['kind']=='1pi' and result['operator']['subgraph_loop_order']==2
    assert result['operator']['insertion_count']==1 and result['J_raw']==result['J_raw_s']
    assert result['normalization']['euler_compensation_count']==0
    assert result['color_contraction']['flavor_channels'][2]['monomial']=='CF^2*TF*nf'
    assert result['color_contraction']['ordered_basis_coordinates'][1]['coordinates'][1]!='0'
    assert result['normalization']['shared_I_convention']['status']=='definition_only_not_prepared'

@pytest.mark.parametrize('mode',['reverse_props','reverse_ports','dummy_indices','reverse_momentum','optional_g6'])
def test_equivalent_source_encodings(mode):
    raw=deepcopy(RAWS[0]);before=build_result(parse_orbit(raw))
    if mode=='reverse_props':raw['propagators'].reverse()
    if mode=='reverse_ports':
        for p in raw['propagators']:p['lorentz'].reverse();p['color'].reverse()
    if mode=='dummy_indices':
        for p in raw['propagators']:
            for f,prefix in [('lorentz','nu'),('color','b')]:p[f]=['renamed_'+v if v.startswith(prefix) else v for v in p[f]]
    if mode=='reverse_momentum':raw['phases'][0]['expression']='-('+raw['phases'][0]['expression']+')'
    if mode=='optional_g6':raw['g_power']=6
    after=build_result(parse_orbit(raw))
    assert after['C']==before['C'] and after['J_raw']==before['J_raw']

@pytest.mark.parametrize('case',[
    'nloops','bool_nloops','g_power','float_g','sector','id','family','Pi1','extra_insertion','chain',
    'missing_prop','extra_prop','momentum','duplicate_mu','duplicate_color','color_crosswire','internal_repeat',
    'endpoint_disconnected','same_edge','unsorted','phase','phase_momentum','phase_trailing','phase_extra','phase_nonstring',
    'duplicate_bound','wrong_bound','bound_extra','order','vertex_extra','operator_extra','prop_extra','bad_port','extra_field','missing_field'])
def test_invalid_sources(case):
    r=deepcopy(RAWS[0])
    if case=='nloops':r['nloops']=2
    elif case=='bool_nloops':r['nloops']=True
    elif case=='g_power':r['g_power']=4
    elif case=='float_g':r['g_power']=6.0
    elif case=='sector':r['sector']='tree'
    elif case=='id':r['orbit_id']='2e-Pi1-01'
    elif case=='family':r['emission_type']='4e-Pi2'
    elif case=='Pi1':r['internal_vertices'][0]['se_type']='vacuum-pol'
    elif case=='extra_insertion':r['internal_vertices'].append(deepcopy(r['internal_vertices'][0]))
    elif case=='chain':r['self_energy_chain']={'insertion_count':2}
    elif case=='missing_prop':r['propagators'].pop()
    elif case=='extra_prop':r['propagators'].append(deepcopy(r['propagators'][0]))
    elif case=='momentum':r['propagators'][0]['loop']='l2'
    elif case=='duplicate_mu':r['vertices'][1]['lorentz']='mu1'
    elif case=='duplicate_color':r['vertices'][1]['color']='a1'
    elif case=='color_crosswire':r['propagators'][0]['color'].reverse()
    elif case=='internal_repeat':r['propagators'][1]['lorentz'][0]='nu1'
    elif case=='endpoint_disconnected':r['propagators'][1]['lorentz'][1]='hidden'
    elif case=='same_edge':r['vertices'][1]['line']=3
    elif case=='unsorted':r['vertices'].reverse()
    elif case=='phase':r['phases'][0]['expression']+=' + v1'
    elif case=='phase_momentum':r['phases'][0]['loop']='l2'
    elif case=='phase_trailing':r['phases'][0]['expression']+=' 2'
    elif case=='phase_extra':r['phases'][0]['phase_factor']=2
    elif case=='phase_nonstring':r['phases'][0]['expression']=['wrong']
    elif case=='duplicate_bound':r['bounds'][1]=deepcopy(r['bounds'][0])
    elif case=='wrong_bound':r['bounds'][0]['high']='s41'
    elif case=='bound_extra':r['bounds'][0]['jacobian']='2'
    elif case=='order':r['iterated_order_outer_to_inner']=['s41','s31']
    elif case=='vertex_extra':r['vertices'][0]['weight']=2
    elif case=='operator_extra':r['internal_vertices'][0]['renormalization']='MSbar'
    elif case=='prop_extra':r['propagators'][0]['power']=2
    elif case=='bad_port':r['propagators'][0]['lorentz']='mu1nu1'
    elif case=='extra_field':r['source_external_sign']=1
    elif case=='missing_field':del r['phases']
    with pytest.raises(ValueError):parse_orbit(r)

@pytest.mark.parametrize('key,field',[('C','combined'),('C','components'),('C','phase_ledger'),('J_raw','integrand'),('J_raw_s','bounds'),
    ('operator','channels'),('operator','kind'),('operator','scalar'),('color_contraction','full_color'),('color_contraction','ordered_basis_coordinates'),
    ('normalization','euler_compensation_count'),('normalization','shared_I_convention'),('c4_symmetry','orbit_multiplicity')])
def test_result_tampering_rejected(key,field):
    r=json.loads(json.dumps(build_result(parse_orbit(RAWS[0]))));r[key][field]='wrong'
    with pytest.raises(ValueError):validate_result(r)

def test_orbit_metadata_must_match():
    orbit=parse_orbit(RAWS[0]);orbit.n_gluons=1
    with pytest.raises(ValueError):build_result(orbit)

def test_c4():
    assert [c4_symmetry(r)['orbit_multiplicity'] for r in RAWS]==[4,2]

def test_actual_cli_and_package(tmp_path):
    source=tmp_path/'source.json';source.write_text(json.dumps({'orbits':RAWS}),encoding='utf-8')
    result=subprocess.run([sys.executable,'-B','-m','montopar.cli','--source',str(source),'--type','2e-Pi2','--out-dir',str(tmp_path),'--package-initial'],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    package,_=build_initial_package({'orbits':RAWS},tmp_path,'2e-Pi2')
    assert package['schema_version']=='1.6' and len(package['orbits'])==2
    for o in package['orbits']:
        assert o['integral_structure']['core']['subgraph_loop_order']==2
        assert o['J']==o['J_raw_s'] and o['topology']==o['representative']['source_topology']

@pytest.mark.parametrize('change',['dispatch','schema','source','P2','Ff','R','J'])
def test_package_rejects_stale_or_detached_result(tmp_path,change):
    raw=deepcopy(RAWS[0]);r=json.loads(json.dumps(build_result(parse_orbit(raw))))
    if change=='dispatch':r['dispatch']='self_energy'
    elif change=='schema':r['schema_version']='1.4'
    elif change=='source':raw['phases'][0]['expression']='-('+raw['phases'][0]['expression']+')'
    elif change=='P2':r['operator']['P2']='1'
    elif change=='Ff':r['color_contraction']['flavor_channels'].pop()
    elif change=='R':r['operator']['R_G10']='1'
    elif change=='J':r['J_raw']['se_numerator']='0'
    p=tmp_path/'3loop/2e-Pi2'/raw['orbit_id']/'montopar.json';p.parent.mkdir(parents=True);p.write_text(json.dumps(r),encoding='utf-8')
    with pytest.raises(ValueError):build_initial_package({'orbits':[raw]},tmp_path,'2e-Pi2')
