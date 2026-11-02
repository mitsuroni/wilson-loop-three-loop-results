"""Reject the re-audit counterexamples before constructing any integral."""
from copy import deepcopy
import itertools
import json
from pathlib import Path
import pytest
from montopar.parser import parse_orbit
from montopar.double_tgv_output import build_double_tgv_result, validate_double_tgv_result
ROWS=json.loads((Path(__file__).parent/'fixtures/double_tgv_orbits.json').read_text(encoding='utf8'))
def raw():return deepcopy(next(r for r in ROWS if r['orbit_id'].endswith('-04')))
def mutate(r,kind):
    if kind=='duplicate_endpoint':
        r['propagators'][1]['lorentz'][0]=r['propagators'][0]['lorentz'][0]
        r['propagators'][1]['color'][0]=r['propagators'][0]['color'][0]
    elif kind=='extra_propagator':r['propagators'].append({'loop':'extra','lorentz':['nu5','nu6'],'color':['b5','b6']})
    elif kind=='extra_vertex':r['internal_vertices'].append({'type':'fgv','loops':['l1','l2','l3','l4']})
    elif kind=='duplicate_internal_mu':r['propagators'][1]['lorentz'][1]=r['propagators'][0]['lorentz'][1]
    elif kind=='duplicate_internal_color':r['propagators'][1]['color'][1]=r['propagators'][0]['color'][1]
    elif kind=='external_internal_color_collision':r['propagators'][1]['color'][1]=r['vertices'][0]['color']
    elif kind=='malformed_connector':
        r['internal_vertices']=[{'type':'tgv','loops':['l1','l2','q']},{'type':'tgv','loops':['l3','l4','q']}]
        r['propagators'].append({'loop':'q','lorentz':['nu1','nu1'],'color':['b1','b1']})
    elif kind=='extra_phase':r['phases'].append({'loop':'unbound','expression':'v1'})
    elif kind=='vertex_annotations':r['internal_vertices'][0]['color']=['b2','b1','b3','b4']
    elif kind=='propagator_arity':r['propagators'][0]['lorentz'].append('nu_extra')
    elif kind=='propagator_annotation':r['propagators'][0]['power']=2
    elif kind=='phase_duplicate':r['phases'][1]=deepcopy(r['phases'][0])
    elif kind=='wrong_vertex_loop':r['internal_vertices'][0]['loops'][0]='missing'
    elif kind=='extra_tgv':r['internal_vertices'].append({'type':'tgv','loops':['l1','l2','q']})
    return r
KINDS=['duplicate_endpoint','extra_propagator','extra_vertex','duplicate_internal_mu','duplicate_internal_color','external_internal_color_collision','malformed_connector','extra_phase','vertex_annotations','propagator_arity','propagator_annotation','phase_duplicate','wrong_vertex_loop','extra_tgv']
@pytest.mark.parametrize('kind',KINDS)
def test_bad_sources_fail_before_build(kind):
    with pytest.raises(ValueError,match='double_tgv_'):parse_orbit(mutate(raw(),kind))
@pytest.mark.parametrize('kind',KINDS)
def test_result_validator_rejects_mutated_source(kind):
    result=build_double_tgv_result(parse_orbit(raw()))
    mutate(result['representative']['canonical_topology'],kind)
    with pytest.raises(ValueError,match='double_tgv_'):validate_double_tgv_result(result)
@pytest.mark.parametrize('index',range(26))
def test_all_current_source_rows_parse(index):
    assert len(parse_orbit(deepcopy(ROWS[index])).double_tgv_legs)==4

def test_explicit_loop_only_vertices_build_same_mathematics():
    source=raw();ref=build_double_tgv_result(parse_orbit(source))
    source['internal_vertices']=[{'type':'tgv','id':'L','loops':['q','l1','l2']},{'type':'tgv','id':'R','loops':['l4','q','l3']}]
    out=build_double_tgv_result(parse_orbit(source))
    for key in ('C','J_raw_s','J_raw','color_contraction','c4_symmetry'):assert out[key]==ref[key]
    assert validate_double_tgv_result(out)

def test_endpoint_direction_and_label_renaming():
    source=raw();ref=build_double_tgv_result(parse_orbit(source))
    for p in source['propagators']:
        p['lorentz'].reverse();p['color'].reverse()
    out=build_double_tgv_result(parse_orbit(source))
    for key in ('C','J_raw_s','J_raw'):assert out[key]==ref[key]

def test_vertex_phase_and_bound_array_permutations():
    source=raw();ref=build_double_tgv_result(parse_orbit(source))
    source['vertices'].reverse();source['phases'].reverse();source['propagators'].reverse()
    out=build_double_tgv_result(parse_orbit(source))
    for key in ('C','J_raw_s','J_raw'):assert out[key]==ref[key]

def test_changed_orbit_identity_rejected():
    orbit=parse_orbit(raw());orbit.orbit_id='wrong'
    with pytest.raises(ValueError,match='metadata_disagrees'):build_double_tgv_result(orbit)
