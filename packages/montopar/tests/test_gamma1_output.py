"""Production routes, exact scalar domains, sector bridges, and result integrity."""
from copy import deepcopy
from functools import lru_cache
import json,os,re,subprocess,sys
from pathlib import Path
import pytest
import sympy as s

from montopar.parser import parse_orbit
from montopar.gamma1 import build_result,validate_result
from montopar.gamma1_blocks import sectorize
from montopar.package_builder import build_initial_package

FIXTURES=Path(__file__).parent/'fixtures'
SOURCE=json.loads((FIXTURES/'gamma1_orbits.json').read_text(encoding='utf-8'))

@lru_cache(None)
def result(i):return build_result(parse_orbit(SOURCE['orbits'][i]),ibp_enabled=True)

def expression(text):
    return s.sympify(text.replace('^','**').replace('[','(').replace(']',')'),
                    locals={'Gamma':s.gamma,'Beta':s.beta,'Max':s.Max,'Pi':s.pi,'Infinity':s.oo})

@pytest.mark.parametrize('i',range(4))
def test_full_output_scope_and_order(i):
    out=result(i);raw=SOURCE['orbits'][i]
    assert out['schema_version']=='1.8' and out['dispatch']=='gamma1'
    assert out['normalization']['euler_compensation_count']==0
    assert not out['normalization']['orbit_multiplicity_included']
    assert out['C']['channels']['g']['full_color']=='CF*CA^2'
    assert out['C']['channels']['q']['full_color']=='CF*CA*TF*nf'
    assert out['J_raw']['term_count']==7
    assert out['J_sector']['term_count']==15
    for term in out['J_raw']['terms']:
        assert term['variables']==term['iterated_order_outer_to_inner']
        assert term['variables'][:3]==raw['iterated_order_outer_to_inner']
        assert len(term['variables'])==(7 if term['kind']=='triangle' else 5)
        assert set(term['variables'])==set(term['bounds'])
        assert term['combined_C']==out['C']['channels'][term['channel']]['combined']
        assert 'partial' not in term['integrand']
        assert all(not re.search(r'\b(r[123][123]|p1sq|p2sq|p1p2|D|h12|h23|h31|t3|alpha3)\b',k['numerator']) for k in term['kernel_terms'])
    for term in out['J_sector']['terms']+out['J_ibp']['sector_terms']:
        assert 'Max[' not in term['integrand']
        if term.get('integration_block_kind',term['kind'])=='triangle':
            assert not re.search(r'\bt[123]\b',term['integrand'])
            assert term['variables'][-2:]==['u','v']
            if 'kernel_coefficients' in term:
                assert term['kernel_coefficients']=={k['derivative_shift']:k['numerator'] for k in term['kernel_terms']}
    assert out['evaluation_status']=='not_evaluated' and out['formal']=='not_run'
    validate_result(json.loads(json.dumps(out)))

@pytest.mark.parametrize('channel',['g','q'])
def test_c_from_bare_formula(channel):
    out=result(0)
    gb=s.Symbol('gb',positive=True);e=s.Rational(1,7)
    replacements={s.Symbol('alpha_s'):gb**2/(4*s.pi),s.Symbol('mu'):1,
                  s.Symbol('v14'):-s.Rational(1,2),s.Symbol('I0'):0,s.Symbol('eps'):e,s.Symbol('g0'):gb}
    actual=expression(out['C']['channels'][channel]['combined']).subs(replacements)
    expected=expression(out['C']['channels'][channel]['bare_coupling_form']).subs(replacements)
    assert s.simplify(actual/expected)==1

def test_sector_structured_coefficients_use_the_same_parenthesized_substitution():
    parent={'term_id':'probe','kind':'contact','integration_block_kind':'triangle','block_id':'N0',
            'primitive_power':'3*eps','weight':'1','q0':'1+alpha1',
            'variables':['s31','alpha1','alpha2','t1','t2'],
            'bounds':{'s31':['0','1'],'alpha1':['0','1'],'alpha2':['0','1-alpha1'],
                      't1':['0','1'],'t2':['0','1-t1']},
            'kernel_coefficients':{2:'t1^2+2*t1*t2'},
            'kernel_terms':[{'derivative_shift':2,'numerator':'t1^2+2*t1*t2','denominator_power':'2-3*eps'}]}
    for term in sectorize(parent):
        numerator=term['kernel_terms'][0]['numerator']
        assert term['kernel_coefficients'][2]==numerator
        mapping={s.Symbol(k):expression(v) for k,v in term['sector']['t_substitution'].items()}
        expected=expression(parent['kernel_terms'][0]['numerator']).subs(mapping,simultaneous=True)
        assert s.cancel(expression(numerator)-expected)==0

@pytest.mark.parametrize('which',['raw','ibp'])
@pytest.mark.parametrize('block',['N0','L0'])
def test_explicit_sector_jacobian_identity(which,block):
    out=result(0)
    source=out['J_raw']['terms'] if which=='raw' else out['J_ibp']['terms']
    parent=next(t for t in source if t['channel']=='g' and t['block_id']==block)
    # This verifies the actual exported integrands with all numerators and
    # a finite regulator, not only a shared symbolic weight formula.
    point={'alpha1':s.Rational(1,5),'alpha2':s.Rational(1,3),'u':s.Rational(7,5),'v':s.Rational(9,5),
           'eps':s.Rational(1,8),'z1':s.Rational(7,4),'I0':s.Rational(1,17)}
    path_bounds={b['var']:b['high'] for b in SOURCE['orbits'][0]['bounds']}
    for variable in SOURCE['orbits'][0]['iterated_order_outer_to_inner']:
        high=point[path_bounds[variable]] if path_bounds[variable] in point else s.S.One
        point[variable]=s.Rational(2,5)*high
    sub={s.Symbol(k):v for k,v in point.items()}
    for sector in sectorize(parent):
        ts={s.Symbol(k):expression(v).subs(sub) for k,v in sector['sector']['t_substitution'].items()}
        jac=expression(sector['sector']['jacobian']).subs(sub)
        def value(term,replacements):
            q=expression(term['q0']).subs(replacements)+replacements[s.Symbol('I0')]
            return expression(term['weight']).subs(replacements)*sum(
                expression(t['numerator']).subs(replacements)*q**(-expression(t['denominator_power']).subs(replacements))
                for t in term['kernel_terms'])
        expected=value(parent,{**sub,**ts})*jac
        actual=value(sector,sub)
        assert abs(complex(s.N(actual-expected,35)))<1e-25

@pytest.mark.parametrize('path,value',[
    (['C','channels','q','combined'],'1'),(['normalization','euler_compensation_count'],1),
    (['J_raw','terms',0,'integrand'],'0'),(['J_sector','terms',0,'weight'],'1'),
    (['J_ibp','terms',0,'combined_C'],'1'),(['J_ibp','sector_terms',0,'integrand'],'0'),
    (['representative','source_topology','g_power'],4),(['c4_symmetry','orbit_multiplicity'],8),
])
def test_reject_changed_result(path,value):
    out=deepcopy(result(0));target=out
    for key in path[:-1]:target=target[key]
    target[path[-1]]=value
    with pytest.raises(ValueError):validate_result(out)

def test_runtime_without_site_packages():
    src=Path(__file__).parents[1]/'src'
    code=("import sys,json;from pathlib import Path;sys.path.insert(0,"+repr(str(src))+");"
          "from montopar.parser import parse_orbit;from montopar.gamma1 import build_result;"
          "r=json.loads(Path("+repr(str(FIXTURES/'gamma1_orbits.json'))+").read_text(encoding='utf-8'))['orbits'][0];"
          "a=build_result(parse_orbit(r),ibp_enabled=True);"
          "assert len(a['J_raw']['terms'])==7;assert 'sympy' not in sys.modules;print('NO_CAS_RUNTIME_OK')")
    proc=subprocess.run([sys.executable,'-S','-c',code],capture_output=True,text=True,timeout=120)
    assert proc.returncode==0,proc.stdout+proc.stderr
    assert 'NO_CAS_RUNTIME_OK' in proc.stdout

def test_real_cli_and_package_reject_wrong_source(tmp_path):
    env={**os.environ,'PYTHONPATH':str(Path(__file__).parents[1]/'src')}
    command=[sys.executable,'-m','montopar.cli','--source',str(FIXTURES/'gamma1_orbits.json'),
             '--type','3e-Gamma1','--gamma1-ibp','--package-initial','--out-dir',str(tmp_path)]
    proc=subprocess.run(command,env=env,capture_output=True,text=True,timeout=600)
    assert proc.returncode==0,proc.stdout+proc.stderr
    package=json.loads((tmp_path/'3loop/3e-Gamma1/montopar-initial.json').read_text(encoding='utf-8'))
    assert package['schema_version']=='1.8' and len(package['orbits'])==4
    assert all(o['dispatch']=='gamma1' and 'J_ibp' in o and 'J_sector' in o for o in package['orbits'])
    path=tmp_path/'3loop/3e-Gamma1/3e-Gamma1-01/montopar.json'
    changed=json.loads(path.read_text(encoding='utf-8'));changed['dispatch']='pure'
    path.write_text(json.dumps(changed),encoding='utf-8')
    with pytest.raises(ValueError,match='gamma1_result_dispatch_mismatch'):
        build_initial_package(SOURCE,tmp_path,'3e-Gamma1')
