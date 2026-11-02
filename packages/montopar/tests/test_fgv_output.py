"""Regression/negative tests for the FGV source contract and regenerated exports."""
import json
from copy import deepcopy
from pathlib import Path
import pytest
from montopar.cli import _compute_one
from montopar.parser import parse_orbit
from montopar.fgv_output import validate_fgv_result, write_compatibility_exports
from montopar.package_builder import build_initial_package, c4_symmetry
from montopar.j_factor import expand_polynomial

SOURCE=json.loads((Path(__file__).parent/'fixtures/4e-4g_orbits.json').read_text(encoding='utf8'))

def generate(tmp_path, index=8):
    raw=deepcopy(SOURCE['orbits'][index]);_compute_one(parse_orbit(raw),tmp_path)
    p=tmp_path/'3loop/4e-4g'/raw['orbit_id']/'montopar.json'
    return raw,json.loads(p.read_text(encoding='utf8')),p

@pytest.mark.parametrize('index',range(9))
def test_current_canonical_self_contained_and_compatibility(tmp_path,index):
    raw,d,p=generate(tmp_path,index)
    assert validate_fgv_result(d)
    assert not any(k.endswith('_rotated') for k in d)
    assert d['C']['components']['color']=='CA^2*CF'
    assert d['normalization']['euler_compensation_count']==0
    assert d['J_raw']['physical_zero']==(index in (0,2,3))
    assert len(d['J_raw_s']['variables'])==4
    assert len(d['J_raw_s']['alpha_measure']['independent_variables'])==3
    target=tmp_path/'compat';write_compatibility_exports(d,target)
    for filename in ('state.json','export.dimreduce.json'):
        exp=json.loads((target/filename).read_text())
        assert exp['J']==d['J_raw']['integrand']
        assert exp['C']==d['C']['combined']
        assert exp['alpha_measure']==d['J_raw']['alpha_measure']
        assert exp['J_raw_s']==d['J_raw_s']
    final=json.loads((target/'final.json').read_text())
    assert final['results']['J']==d['J_raw']
    package,_=build_initial_package({'orbits':[raw]},tmp_path,'4e-4g')
    assert package['schema_version']=='1.2'
    assert package['orbits'][0]['J_raw_s']==d['J_raw_s']

@pytest.mark.parametrize('mutation',('simplex','delta','order','factor','phase','euler','rotated','factor_order','branch'))
def test_rejects_incomplete_or_inconsistent_semantics(tmp_path,mutation):
    _,d,_=generate(tmp_path)
    if mutation=='factor_order': d['C']['factor_order'].append('color')
    if mutation=='branch': d['branch_prescription']={}
    if mutation=='simplex': d['J_raw_s']['alpha_measure']={}
    if mutation=='delta': d['J_raw']['delta_constraints']={}
    if mutation=='order': d['J_raw_s']['iterated_order_outer_to_inner']=[]
    if mutation=='factor': d['C']['components']['twoscale']='1'
    if mutation=='phase': d['C']['combined']='-('+d['C']['combined']+')'
    if mutation=='euler': d['normalization']['euler_compensation_count']=1
    if mutation=='rotated': d['J_raw_rotated']=deepcopy(d['J_raw'])
    with pytest.raises(ValueError):validate_fgv_result(d)

@pytest.mark.parametrize('mutation',('old','split','loop','word'))
def test_rejects_bad_sources_and_stale_packages(tmp_path,mutation):
    raw=deepcopy(SOURCE['orbits'][1])
    if mutation=='old':
        folder=tmp_path/'3loop/4e-4g'/raw['orbit_id'];folder.mkdir(parents=True)
        (folder/'montopar.json').write_text(json.dumps({'dispatch':'fourv','C':{'combined':'old'},'J_raw':{'integrand':'old'}}))
        with pytest.raises(ValueError,match='regenerate'):build_initial_package({'orbits':[raw]},tmp_path,'4e-4g')
        return
    if mutation=='split':raw['bounds'][0]['high']='1'
    if mutation=='loop':raw['nloops']=2
    if mutation=='word':raw['propagators'][0],raw['propagators'][1]=raw['propagators'][1],raw['propagators'][0]
    with pytest.raises(ValueError):_compute_one(parse_orbit(raw),tmp_path)

def test_09_wrap_around_and_c4_multiplicity(tmp_path):
    _,d,_=generate(tmp_path)
    # Coefficient extraction is performed independently in the external Gram audit;
    # here keep a minimal regression for the historically corrupted monomials.
    q=d['J_raw']['Q_prime']
    assert 'x1*alpha1*alpha4' in q
    assert ' + x4*alpha1*alpha4' not in q
    assert d['J_raw']['lorentzZ']=='-1/4 - z1 - 1/4*z1^2'
    assert [c4_symmetry(r)['orbit_multiplicity'] for r in SOURCE['orbits']]==[4,4,4,4,4,4,2,4,1]