from copy import deepcopy
from collections import Counter
import pytest
from qgraf_pipeline.internal_types import enumerate_emission_types
from qgraf_pipeline.generator import generate_diagrams, _build_two_internal_diagram
from qgraf_pipeline.classifier import classify_orbits
from qgraf_pipeline.topology import classify_topology
from qgraf_pipeline.definition_exporter import _build_orbit_definition
from qgraf_pipeline.emitter import emit_integral
from qgraf_pipeline.tikz_emitter import _emit_one_diagram
from qgraf_pipeline.self_energy_chain import FAMILY, generate_pi1_chain, build_chain_definition


def chain_orbits():
    types=[t for t in enumerate_emission_types(3) if t.name==FAMILY]
    return classify_orbits(generate_diagrams(3,types),types)


def test_distinct_from_pi2_and_not_duplicated():
    types=enumerate_emission_types(3)
    assert len([t for t in types if t.name==FAMILY])==1
    t=next(t for t in types if t.name==FAMILY)
    assert (t.eikonal_count,t.g_power,t.internal,t.wl_wl_pairs)==(2,6,{'Pi1':2},0)
    assert next(t for t in types if t.name=='2e-Pi2').internal=={'Pi2':1}
    assert all(t.name!=FAMILY for t in enumerate_emission_types(2))
    assert all(t.name!=FAMILY for t in enumerate_emission_types(4))


def test_two_geometries_one_word_each():
    orbits=chain_orbits()
    assert [o.edge_counts for o in orbits]==[(0,0,1,1),(0,1,0,1)]
    assert [o.n_diagrams for o in orbits]==[1,1]
    for o in orbits:
        d=o.representative;info=classify_topology(d)
        assert len(d.props)==3 and len(d.verts)==4
        assert len(info.internal_connections)==2 and len(info.internal_gluon_links)==1
        assert o.is_2pi  # contour web cut; not 1PI of the corrected propagator


def test_same_edge_chain_scaleless():
    d=_build_two_internal_diagram(1,(0,0,0,2),[0],[1],'Pi1','Pi1',1,et_name=FAMILY)
    info=classify_topology(d)
    assert info.is_scaleless and info.scaleless_reason=='B'  # earlier all-endpoints-on-one-edge filter


@pytest.mark.parametrize('index',[0,1])
def test_chain_index_contractions_and_shared_momentum(index):
    d=_build_orbit_definition(chain_orbits()[index])
    assert d['nloops']==3 and d['g_power']==6
    assert len(d['propagators'])==3 and len(d['internal_vertices'])==2
    assert {p['loop'] for p in d['propagators']}=={'l1'} and len(d['phases'])==1
    for key in ['lorentz','color']:
        counts=Counter(i for p in d['propagators']+d['internal_vertices'] for i in p[key])
        external={v[key] for v in d['vertices']}
        assert all(counts[i]==1 for i in external)
        assert all(n==2 for i,n in counts.items() if i not in external)
    assert d['self_energy_chain']['ordered_insertions']==['se1','se2']
    assert d['self_energy_chain']['dyson_coefficient']=='1'
    assert d['self_energy_chain']['source_external_sign']=='+1'
    assert all(v['se_type']=='vacuum-pol' for v in d['internal_vertices'])
    assert len(d['bounds'])==2 and all(b['low']=='0' and b['high']=='1' for b in d['bounds'])


@pytest.mark.parametrize('index',[0,1])
def test_formula_and_drawing_have_two_insertions(index):
    orbit=chain_orbits()[index];formula=emit_integral(orbit)
    assert formula.count(r'\Pi^{(1)')==2
    assert formula.count(r'\ell^2+i0')==3
    assert formula.count(r'\int\frac{d^D\ell}')==1
    assert '\\Pi^{(2)' not in formula and ' = -' not in formula
    assert _emit_one_diagram(orbit).count('draw[fill=white]')==2


def test_reject_disconnected_chain():
    orbit=chain_orbits()[0];info=deepcopy(classify_topology(orbit.representative));info.internal_gluon_links=[]
    with pytest.raises(ValueError,match='invalid_self_energy_chain_topology'):build_chain_definition(orbit,info)


def test_reject_unsupported_mixed_loop_chain():
    t=next(t for t in enumerate_emission_types(3) if t.name==FAMILY);t=deepcopy(t);t.internal={'Pi1':1,'Pi2':1}
    with pytest.raises(ValueError,match='unsupported_self_energy_chain'):generate_pi1_chain(t,1)
