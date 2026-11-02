from copy import deepcopy
from pathlib import Path
import itertools,json,subprocess,sys
import pytest
from montopar.parser import parse_orbit
from montopar.self_energy_chain import build_result,validate_result,validate_source,c4_symmetry
from montopar.package_builder import build_initial_package

RAWS=json.loads((Path(__file__).parent/'fixtures/self_energy_chain_orbits.json').read_text(encoding='utf-8'))['orbits']


@pytest.mark.parametrize('raw',RAWS,ids=lambda o:o['orbit_id'])
def test_full_chain_and_source_reconstruction(raw):
    orbit=parse_orbit(raw);assert orbit.dispatch=='self_energy_chain'
    result=json.loads(json.dumps(build_result(orbit)));validate_result(result)
    assert result['schema_version']=='1.5' and result['J_raw']==result['J_raw_s']
    assert result['operator']['insertion_count']==2 and result['operator']['kind']=='one_particle_reducible'
    assert len(result['representative']['contraction_path'])==5
    assert result['normalization']['euler_compensation_count']==0
    assert 'Gamma[1-3*eps]' in result['C']['components']['gamma_chain_ratio']
    assert len(result['color_contraction']['flavor_channels'])==3


@pytest.mark.parametrize('perm',list(itertools.permutations(range(3))))
def test_propagator_permutation_and_port_reversal(perm):
    raw=deepcopy(RAWS[0]);original=build_result(parse_orbit(raw))
    raw['propagators']=[raw['propagators'][i] for i in perm]
    for p in raw['propagators']:
        p['lorentz'].reverse();p['color'].reverse()
    raw['internal_vertices'].reverse()
    changed=build_result(parse_orbit(raw))
    assert changed['J_raw']==original['J_raw'] and changed['C']==original['C']


def test_dummy_renaming_and_routing_reversal():
    raw=deepcopy(RAWS[1]);before=build_result(parse_orbit(raw))
    for p in raw['propagators']+raw['internal_vertices']:
        for field,prefix in [('lorentz','nu'),('color','b')]:
            p[field]=['renamed_'+x if x.startswith(prefix) else x for x in p[field]]
    raw['phases'][0]['expression']='-('+raw['phases'][0]['expression']+')'
    after=build_result(parse_orbit(raw));assert before['C']==after['C'] and before['J_raw']==after['J_raw']


@pytest.mark.parametrize('case',[
    'missing_middle','extra_prop','mixed_momenta','Pi2','missing_operator','duplicate_operator_id',
    'wrong_g','wrong_n','wrong_sector','bad_sign','factorial','wrong_order','wrong_counts','unknown_metadata',
    'duplicate_external_mu','duplicate_external_color','internal_external_collision','color_crosswire',
    'phase_momentum','phase_significant_change','phase_unknown','phase_trailing','duplicate_phase',
    'wrong_bound','duplicate_bound','same_edge','bad_explicit_order','closed_components','internal_self_loop',
    'extra_field','missing_chain','renormalized','bad_family','stale_subgraph_order'])
def test_bad_chain_sources_rejected(case):
    raw=deepcopy(RAWS[0])
    if case=='missing_middle':raw['propagators'].pop(1)
    elif case=='extra_prop':raw['propagators'].append(deepcopy(raw['propagators'][1]))
    elif case=='mixed_momenta':raw['propagators'][1]['loop']='l2'
    elif case=='Pi2':raw['internal_vertices'][0]['se_type']='vacuum-pol-2loop'
    elif case=='missing_operator':raw['internal_vertices'].pop()
    elif case=='duplicate_operator_id':raw['internal_vertices'][1]['id']='se1'
    elif case=='wrong_g':raw['g_power']=4
    elif case=='wrong_n':raw['nloops']=2
    elif case=='wrong_sector':raw['sector']='tree'
    elif case=='bad_sign':raw['self_energy_chain']['source_external_sign']='-1'
    elif case=='factorial':raw['self_energy_chain']['dyson_coefficient']='1/2'
    elif case=='wrong_order':raw['self_energy_chain']['ordered_insertions'].reverse()
    elif case=='wrong_counts':raw['self_energy_chain']['insertion_count']=1
    elif case=='unknown_metadata':raw['self_energy_chain']['ignored_factor']='2'
    elif case=='duplicate_external_mu':raw['vertices'][1]['lorentz']='mu1'
    elif case=='duplicate_external_color':raw['vertices'][1]['color']='a1'
    elif case=='internal_external_collision':raw['propagators'][1]['lorentz'][0]='mu1'
    elif case=='color_crosswire':raw['propagators'][1]['color'].reverse()
    elif case=='phase_momentum':raw['phases'][0]['loop']='l2'
    elif case=='phase_significant_change':raw['phases'][0]['expression']+=' + v1'
    elif case=='phase_unknown':raw['phases'][0]['expression']+=' + hidden'
    elif case=='phase_trailing':raw['phases'][0]['expression']+=' 2'
    elif case=='duplicate_phase':raw['phases'].append(deepcopy(raw['phases'][0]))
    elif case=='wrong_bound':raw['bounds'][0]['high']='s41'
    elif case=='duplicate_bound':raw['bounds'][1]=deepcopy(raw['bounds'][0])
    elif case=='same_edge':raw['vertices'][1]['line']=3
    elif case=='bad_explicit_order':raw['iterated_order_outer_to_inner']=['s41','s31']
    elif case=='closed_components':
        raw['propagators']=[{'loop':'l1','lorentz':['mu1','mu2'],'color':['a1','a2']},
          {'loop':'l1','lorentz':['nu1','nu2'],'color':['b1','b2']},
          {'loop':'l1','lorentz':['nu3','nu4'],'color':['b3','b4']}]
    elif case=='internal_self_loop':raw['internal_vertices'][0]['lorentz']=['nu1','nu1']
    elif case=='extra_field':raw['hidden_sign']=-1
    elif case=='missing_chain':del raw['self_energy_chain']
    elif case=='renormalized':raw['internal_vertices'][0]['renormalization']='MSbar'
    elif case=='bad_family':raw['emission_type']='2e-Pi2'
    elif case=='stale_subgraph_order':raw['self_energy_chain']['subgraph_loop_orders']=[2,0]
    with pytest.raises(ValueError):parse_orbit(raw)


def test_c4_identical_insertions_do_not_double_count():
    assert [c4_symmetry(r)['orbit_multiplicity'] for r in RAWS]==[4,2]
    raw=deepcopy(RAWS[1]);raw['internal_vertices'][0]['id']='first';raw['internal_vertices'][1]['id']='last'
    raw['self_energy_chain']['ordered_insertions']=['first','last'];assert c4_symmetry(raw)['orbit_multiplicity']==2


@pytest.mark.parametrize('key,field',[('C','combined'),('C','components'),('J_raw','integrand'),('J_raw_s','integrand'),
    ('operator','insertion_count'),('color_contraction','full_color'),('normalization','euler_compensation_count'),('c4_symmetry','orbit_multiplicity')])
def test_output_reconstruction_rejects_altered_physics(key,field):
    result=json.loads(json.dumps(build_result(parse_orbit(RAWS[0]))));result[key][field]='wrong'
    with pytest.raises(ValueError):validate_result(result)


def test_stale_orbit_rejected():
    orbit=parse_orbit(RAWS[0]);orbit.n_gluons=1
    with pytest.raises(ValueError):build_result(orbit)


def test_actual_cli_and_initial_package(tmp_path):
    source=tmp_path/'source.json';raw={'orbits':deepcopy(RAWS)}
    unsupported=deepcopy(RAWS[0]);unsupported['emission_type']='2e-Pi2';raw['orbits'].append(unsupported)
    source.write_text(json.dumps(raw),encoding='utf-8')
    run=subprocess.run([sys.executable,'-B','-m','montopar.cli','--source',str(source),'--type','2e-Pi1+Pi1','--out-dir',str(tmp_path),'--package-initial'],capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
    package,_=build_initial_package(raw,tmp_path,'2e-Pi1+Pi1')
    assert package['schema_version']=='1.5' and len(package['orbits'])==2
    for o in package['orbits']:
        assert o['integral_structure']['core']['insertion_count']==2
        assert o['topology']['self_energy_chain']['dyson_coefficient']=='1'
        assert o['J']==o['J_raw_s']


@pytest.mark.parametrize('change',['dispatch','source','gamma_ratio','numerator'])
def test_package_rejects_downgrade_and_detachment(tmp_path,change):
    raw=deepcopy(RAWS[0]);result=json.loads(json.dumps(build_result(parse_orbit(raw))))
    if change=='dispatch':result['dispatch']='self_energy'
    elif change=='source':raw['phases'][0]['expression']='-('+raw['phases'][0]['expression']+')'
    elif change=='gamma_ratio':result['C']['components']['gamma_chain_ratio']='1'
    elif change=='numerator':result['J_raw']['se_numerator']='0'
    p=tmp_path/f"3loop/2e-Pi1+Pi1/{raw['orbit_id']}/montopar.json";p.parent.mkdir(parents=True);p.write_text(json.dumps(result),encoding='utf-8')
    with pytest.raises(ValueError):build_initial_package({'orbits':[raw]},tmp_path,'2e-Pi1+Pi1')
