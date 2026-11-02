"""Bind production output to independent theory fixtures and real CLI use."""
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import sympy as sp

from montopar.parser import parse_orbit
from montopar.tgv_pi1 import build_result, validate_source, validate_result, geometry
from montopar.package_builder import build_initial_package

FIXTURES=Path(__file__).parent/'fixtures'
SOURCE=json.loads((FIXTURES/'tgv_pi1_orbits.json').read_text(encoding='utf-8'))
GEOMETRY=json.loads((FIXTURES/'tgv_pi1_geometry.json').read_text(encoding='utf-8'))


def symbolic(expr):
    return sp.sympify(expr.replace('^','**'))


@pytest.mark.parametrize('raw',SOURCE['orbits'],ids=lambda o:o['orbit_id'])
def test_direct_kernel_against_independent_theory(raw):
    fixture=next(g for g in GEOMETRY if g['orbit_id']==raw['orbit_id'])
    result=build_result(parse_orbit(raw),ibp_enabled=True)
    geom=result['J_raw']['geometry_polynomials']
    names={'a':'alpha1','b':'alpha2','e':'eps','z':'z1',**{f's{i+1}':v for i,v in enumerate(fixture['s_aliases'])}}
    for key,expected in [('q0',fixture['q']),('L',fixture['L']),('E',fixture['E']),('F',fixture['F'])]:
        comparison=symbolic(expected).subs({sp.Symbol(k):sp.Symbol(v) for k,v in names.items()},simultaneous=True)
        assert sp.expand(symbolic(geom[key])-comparison)==0
    assert len(result['J_raw']['variables'])==5
    assert result['J_raw']['source_order_outer_to_inner']==raw['iterated_order_outer_to_inner']
    assert result['normalization']['euler_compensation_count']==0
    assert result['normalization']['color_projection']=='none'
    assert result['c4_symmetry']['orbit_multiplicity']==4
    assert result['evaluation_status']=='not_evaluated'
    assert result['J_ibp']['flux']['coefficient_reconstruction'].startswith('exact')
    assert any(t['kind']=='contact' for t in result['J_ibp']['terms']) == (raw['vertices'][0]['line']==raw['vertices'][1]['line'] or raw['vertices'][1]['line']==raw['vertices'][2]['line'])
    if len(set(v['line'] for v in raw['vertices']))==3:
        assert any(t['kind']=='recoil' for t in result['J_ibp']['terms'])
    else:
        assert not any(t['kind']=='recoil' for t in result['J_ibp']['terms'])
    validate_result(json.loads(json.dumps(result)))


@pytest.mark.parametrize('path,value',[
    (['nloops'],2),(['g_power'],4),(['sign_convention'],'legacy'),
    (['self_energy_insertion','operator'],'Pi2'),(['self_energy_insertion','marked_contour_position'],1),
    (['self_energy_insertion','source_external_sign'],'-1'),(['self_energy_insertion','dyson_coefficient'],'1/2'),
    (['self_energy_insertion','scalar_contains_color'],False),(['self_energy_insertion','renormalization'],'renormalized'),
    (['internal_vertices',0,'incoming_momenta'],['l1','-l2','-l3']),
    (['internal_vertices',0,'vertex_loop_order'],1),(['internal_vertices',1,'propagator'],'l1'),
    (['propagators',2,'loop'],'l1'),(['propagators',3,'lorentz'],['rho','nu3']),
    (['propagators',0,'color'],['a1','b2']),(['phases',0,'expression'],'s31*v1'),
    (['bounds',1,'high'],'1'),(['bounds',0,'low'],'1/2'),
    (['iterated_order_outer_to_inner'],['s31','s41','s42']),
    (['momentum_constraints',0,'expression'],'l1+l2-l3'),
    (['momentum_constraints',0,'normalization'],'delta^D(l1+l2+l3)'),
])
def test_reject_changed_mathematical_source(path,value):
    raw=deepcopy(SOURCE['orbits'][0]);target=raw
    for key in path[:-1]:target=target[key]
    target[path[-1]]=value
    with pytest.raises(ValueError):parse_orbit(raw)


def test_propagator_array_order_and_index_renaming_are_semantic():
    raw=deepcopy(SOURCE['orbits'][0])
    raw['propagators'].reverse()
    # Renaming both ends of one internal contraction remains a legitimate source.
    for v in raw['internal_vertices']:
        v['lorentz']=['renamed' if x=='rho' else x for x in v['lorentz']]
    for p in raw['propagators']:
        p['lorentz']=['renamed' if x=='rho' else x for x in p['lorentz']]
    validate_source(raw)
    assert geometry(validate_source(raw))==geometry(validate_source(SOURCE['orbits'][0]))


def test_result_validation_rejects_dropped_tensor_and_repeated_euler():
    result=build_result(parse_orbit(SOURCE['orbits'][0]),ibp_enabled=True)
    for path,value in [(['J_raw','kernel_terms',1,'numerator'],'0'),
                       (['normalization','euler_compensation_count'],1),
                       (['J_ibp','terms',0,'combined_C'],'1')]:
        changed=deepcopy(result);target=changed
        for key in path[:-1]:target=target[key]
        target[path[-1]]=value
        with pytest.raises(ValueError):validate_result(changed)


@pytest.mark.parametrize('number',range(12))
@pytest.mark.parametrize('rotation',range(4))
def test_c4_scale_and_marked_leg_are_preserved(number,rotation):
    raw=deepcopy(SOURCE['orbits'][number]); old=deepcopy(raw)
    old_data=validate_source(old); permutation=sorted(range(3),key=lambda i:(old['vertices'][i]['line']-1+rotation)%4+1)
    raw['vertices']=[dict(old['vertices'][i],line=(old['vertices'][i]['line']-1+rotation)%4+1) for i in permutation]
    tgv=raw['internal_vertices'][0]
    for key in ('loops','lorentz','color'):
        tgv[key]=[old['internal_vertices'][0][key][i] for i in permutation]
    tgv['incoming_momenta']=['-'+l for l in tgv['loops']]
    slots={};newvars=[];bounds=[];phase=[]
    for v,l in zip(raw['vertices'],tgv['loops']):
        edge=v['line'];slots[edge]=slots.get(edge,0)+1
        var=f's{edge}{slots[edge]}';newvars.append(var)
        phase.append({'loop':l,'expression':'-('+'+'.join([*(f'v{k}' for k in range(1,edge)),f'{var}*v{edge}'])+')'})
    for v,var in zip(raw['vertices'],newvars):
        edge=v['line'];slot=int(var[2:])
        bounds.append({'var':var,'low':'0','high':f's{edge}{slot+1}' if slot<slots[edge] else '1'})
    raw['bounds']=bounds;raw['phases']=phase
    raw['iterated_order_outer_to_inner']=[f's{edge}{slot}' for edge in sorted(slots) for slot in range(slots[edge],0,-1)]
    ins=raw['self_energy_insertion'];ins['marked_contour_position']=permutation.index(old_data['marked_leg']-1)+1
    ri=ins['marked_contour_position']-1
    ins['marked_endpoint']=[raw['vertices'][ri]['line'],int(newvars[ri][2:])]
    total='+'.join(tgv['loops'])
    raw['momentum_constraints']=[{'type':'delta','expression':total,'normalization':f'(2*pi)^D*delta^D({total})'}]
    from montopar.sign_rules import expected_contract
    raw['physical_convention']=expected_contract(raw['vertices'])
    expected=geometry(old_data);actual=geometry(validate_source(raw))
    a1,a2,z=sp.symbols('alpha1 alpha2 z1');oldalphas=[a1,a2,1-a1-a2]
    mapping={sp.Symbol(newvars[i]):sp.Symbol(old_data['endpoints'][permutation[i]]['s_var']) for i in range(3)}
    mapping.update({a1:oldalphas[permutation[0]],a2:oldalphas[permutation[1]]})
    for key,degree in [('q0',1),('L',2),('E',2),('F',3)]:
        lhs=symbolic(actual[key]).subs(mapping,simultaneous=True)
        rhs=symbolic(expected[key]) if rotation%2==0 else z**degree*symbolic(expected[key]).subs(z,1/z)
        assert sp.expand(lhs-rhs)==0


def test_real_cli_and_initial_package(tmp_path):
    env={**os.environ,'PYTHONPATH':str(Path(__file__).parents[1]/'src')}
    cmd=[sys.executable,'-m','montopar.cli','--source',str(FIXTURES/'tgv_pi1_orbits.json'),
         '--type','3e-3g+Pi1','--tgv-ibp','--package-initial','--out-dir',str(tmp_path)]
    completed=subprocess.run(cmd,env=env,capture_output=True,text=True,timeout=120)
    assert completed.returncode==0,completed.stdout+completed.stderr
    package=json.loads((tmp_path/'3loop/3e-3g+Pi1/montopar-initial.json').read_text(encoding='utf-8'))
    assert package['schema_version']=='1.7' and len(package['orbits'])==12
    assert all(o['dispatch']=='tgv_pi1' and 'J_ibp' in o for o in package['orbits'])
    path=tmp_path/'3loop/3e-3g+Pi1/3e-3g+Pi1-01/montopar.json'
    result=json.loads(path.read_text(encoding='utf-8'));result['C']['combined']='1'
    path.write_text(json.dumps(result),encoding='utf-8')
    with pytest.raises(ValueError):build_initial_package(SOURCE,tmp_path,'3e-3g+Pi1')
