from collections import Counter
from copy import deepcopy
from itertools import combinations_with_replacement
import pytest
from qgraf_pipeline.internal_types import enumerate_emission_types
from qgraf_pipeline.generator import generate_diagrams
from qgraf_pipeline.classifier import classify_orbits
from qgraf_pipeline.topology import classify_topology
from qgraf_pipeline.definition_exporter import _build_orbit_definition
from qgraf_pipeline.emitter import emit_integral
from qgraf_pipeline.central_vertex import build_definition,validate_definition,OPERATOR,operator_hash

@pytest.fixture(scope='module')
def orbits():
    ts=[t for t in enumerate_emission_types(3) if t.name=='3e-Gamma1']
    return classify_orbits(generate_diagrams(3,ts),ts)

def test_independent_c4_enumeration(orbits):
    keys=set()
    for ends in combinations_with_replacement(range(1,5),3):
        if len(set(ends))==1:continue
        keys.add(min(tuple(sorted((e-1+r)%4+1 for e in ends)) for r in range(4)))
    actual={tuple(e for e,n in enumerate(o.edge_counts,1) for _ in range(n)) for o in orbits}
    assert len(keys)==len(actual)==4
    assert {min(tuple(sorted((e-1+r)%4+1 for e in ends)) for r in range(4)) for ends in actual}==keys

@pytest.mark.parametrize('index',range(4))
def test_complete_momentum_indices_and_order(orbits,index):
    o=orbits[index];d=_build_orbit_definition(o);validate_definition(d)
    assert d['g_power']==6 and d['nloops']==3
    assert [p['loop'] for p in d['propagators']]==['l1','l2','l3']
    assert d['momentum_constraints']==[dict(type='delta',expression='l1+l2+l3',normalization='(2*pi)^D*delta^D(l1+l2+l3)')]
    assert d['internal_vertices'][0]['incoming_momenta']==['-l1','-l2','-l3']
    for field in ['lorentz','color']:
        c=Counter(v[field] for v in d['vertices'])
        c.update(x for p in d['propagators']+d['internal_vertices'] for x in p[field])
        assert set(c.values())=={2}
    order=d['iterated_order_outer_to_inner'];bounds={b['var']:b for b in d['bounds']}
    for i,v in enumerate(order):assert bounds[v]['high']=='1' or bounds[v]['high'] in order[:i]
    assert len(d['phases'])==3
    assert d['vertex_correction']['internal_loop_measure_already_in_operator']
    assert not d['vertex_correction']['external_leg_self_energies_included']
    assert not d['vertex_correction']['counterterm_included']
    assert d['internal_vertices'][0]['operator_sha256']==operator_hash()
    text=emit_integral(o)
    assert text.count('^2+i0')==3 and r'\ell^2+i0' not in text
    assert r'\prod_{j=1}^{3}\int' in text and r'\delta^{(D)}(\ell_1+\ell_2+\ell_3)' in text
    assert text.count(r'\cdot\bigl(')==3 and r'\operatorname{Tr}' in text
    assert 'hep-ph/9605348v3' in text

@pytest.mark.parametrize('field', ['propagators','internal_vertices','phases','bounds','momentum_constraints','vertex_correction','sign_convention'])
def test_missing_source_fields_rejected(orbits,field):
    d=_build_orbit_definition(orbits[0]);del d[field]
    with pytest.raises(ValueError):validate_definition(d)

@pytest.mark.parametrize('mutation', ['wrong_loop','wrong_colour','wrong_phase','wrong_bound','wrong_order','wrong_delta','wrong_operator','wrong_shape','wrong_gauge','counterterm','extra_field'])
def test_corrupt_source_rejected(orbits,mutation):
    d=_build_orbit_definition(orbits[0])
    if mutation=='wrong_loop':d['propagators'][1]['loop']='l1'
    elif mutation=='wrong_colour':d['propagators'][0]['color'][1]='b2'
    elif mutation=='wrong_phase':d['phases'][0]['expression']='s31*v3'
    elif mutation=='wrong_bound':d['bounds'][1]['high']='1'
    elif mutation=='wrong_order':d['iterated_order_outer_to_inner'].reverse()
    elif mutation=='wrong_delta':d['momentum_constraints'][0]['expression']='l1+l2-l3'
    elif mutation=='wrong_operator':d['internal_vertices'][0]['type']='tgv'
    elif mutation=='wrong_shape':d['vertex_correction']['operator_definition']['form_factors']['quark']['A']='0'
    elif mutation=='wrong_gauge':d['vertex_correction']['gauge']='Landau'
    elif mutation=='counterterm':d['vertex_correction']['counterterm_included']=True
    else:d['unexpected']=True
    with pytest.raises(ValueError):validate_definition(d)

def test_wrong_topology_rejected(orbits):
    o=deepcopy(orbits[0]);info=classify_topology(o.representative)
    info.internal_connections=[(1,2)]
    with pytest.raises(ValueError):build_definition(o,info)

def test_formula_origin_is_specific():
    assert OPERATOR['reference']['arxiv']=='hep-ph/9605348v3'
    assert OPERATOR['reference']['erratum_doi']=='10.1103/PhysRevD.59.109901'
    assert OPERATOR['form_factors']['gluon_ghost']['S']==OPERATOR['form_factors']['quark']['S']=='0'
    assert set(OPERATOR['tensor_basis'])==set('ABCSFH')
