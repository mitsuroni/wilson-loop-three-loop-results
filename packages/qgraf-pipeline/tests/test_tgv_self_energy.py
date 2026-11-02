from collections import Counter
from copy import deepcopy
from itertools import product
import pytest

from qgraf_pipeline.internal_types import enumerate_emission_types
from qgraf_pipeline.generator import generate_diagrams, _build_two_internal_diagram
from qgraf_pipeline.classifier import classify_orbits, _build_type_map
from qgraf_pipeline.topology import classify_topology
from qgraf_pipeline.definition_exporter import _build_orbit_definition
from qgraf_pipeline.emitter import emit_integral
from qgraf_pipeline.tikz_emitter import _emit_one_diagram
from qgraf_pipeline.tgv_self_energy import FAMILY, CONVENTION, build_definition


def orbits():
    ts=[t for t in enumerate_emission_types(3) if t.name==FAMILY]
    return classify_orbits(generate_diagrams(3,ts),ts)


def test_distinct_operator_and_loop_count():
    ts=enumerate_emission_types(3)
    t=next(t for t in ts if t.name==FAMILY)
    assert (t.eikonal_count,t.g_power,t.wl_wl_pairs)==(3,6,0)
    assert t.internal=={'3g':1,'Pi1':1} and t.category=='loop'
    assert next(t for t in ts if t.name=='3e-Gamma1').internal=={'Gamma1':1}
    assert _build_type_map(ts)[(3,True,frozenset())].name=='3e-Gamma1'
    assert all(t.name!=FAMILY for n in (1,2,4) for t in enumerate_emission_types(n))


def test_independent_C4_count_and_marked_positions():
    # Enumerate edge multisets and one distinguished ordered slot independently.
    keys=set()
    for edges in product(range(1,5),repeat=3):
        if tuple(sorted(edges))!=edges or len(set(edges))==1:
            continue
        ends=[]; counts=Counter()
        for edge in edges:
            counts[edge]+=1;ends.append((edge,counts[edge]))
        for marked in ends:
            rotations=[]
            for r in range(4):
                rotated=tuple(sorted((((a-1+r)%4+1,b) for a,b in ends)))
                mark=((marked[0]-1+r)%4+1,marked[1])
                rotations.append((rotated,mark))
            keys.add(min(rotations))
    actual=set()
    for o in orbits():
        d=_build_orbit_definition(o)
        ends=tuple((a,b) for a in range(1,5) for b in range(1,o.edge_counts[a-1]+1))
        rotations=[]
        for r in range(4):
            mark=d['self_energy_insertion']['marked_endpoint']
            rotations.append((tuple(sorted((((a-1+r)%4+1,b) for a,b in ends))),
                              ((mark[0]-1+r)%4+1,mark[1])))
        actual.add(min(rotations))
    assert actual==keys and len(keys)==12
    assert Counter(o.edge_counts for o in orbits())=={
        (0,0,1,2):3,(0,0,2,1):3,(0,1,0,2):3,(0,1,1,1):3}


@pytest.mark.parametrize('index',range(12))
def test_closed_indices_shared_momentum_order_and_formula(index):
    o=orbits()[index];d=_build_orbit_definition(o)
    assert d['sign_convention']==CONVENTION
    ins=d['self_energy_insertion'];r=ins['marked_contour_position']
    assert ins['source_external_sign']=='+1' and ins['dyson_coefficient']=='1'
    assert ins['not_vertex_correction'] and ins['subgraph_loop_order']==1
    assert Counter(p['loop'] for p in d['propagators'])==Counter(['l1','l2','l3',f'l{r}'])
    for field in ['lorentz','color']:
        counts=Counter(v[field] for v in d['vertices'])
        counts.update(x for p in d['propagators']+d['internal_vertices'] for x in p[field])
        assert set(counts.values())=={2}
    assert d['momentum_constraints'][0]['expression']=='l1+l2+l3'
    order=d['iterated_order_outer_to_inner'];bounds={b['var']:b for b in d['bounds']}
    assert len(order)==len(set(order))==len(bounds)==3
    for i,var in enumerate(order):
        hi=bounds[var]['high']
        assert hi=='1' or hi in order[:i]
    assert len(d['phases'])==3 and {p['loop'] for p in d['phases']}=={'l1','l2','l3'}
    formula=emit_integral(o)
    assert formula.count(r'\Pi^{(1)')==1 and formula.count('^2+i0')==4
    assert r'\delta^{(D)}(\ell_1+\ell_2+\ell_3)' in formula
    assert ' = -' not in formula and r'\Gamma^{(1)' not in formula
    assert _emit_one_diagram(o).count('draw[fill=white]')==1
    assert d['downstream_status']['montopar']=='not_implemented'
    assert o.is_2pi and o.n_diagrams==1


def test_three_insertions_remain_distinct_for_repeated_edge():
    same=[_build_orbit_definition(o) for o in orbits() if o.edge_counts==(0,0,1,2)]
    assert {tuple(x['self_energy_insertion']['marked_endpoint']) for x in same}=={(3,1),(4,1),(4,2)}


def test_mirror_is_metadata_not_identification():
    os=orbits();byid={o.orbit_id:o for o in os}
    assert len([o for o in os if o.mirror_of is None])==2
    for o in os:
        if o.mirror_of:assert byid[o.mirror_of].mirror_of==o.orbit_id


def test_single_edge_scaleless_and_untagged_topology():
    d=_build_two_internal_diagram(1,(0,0,0,3),[0,1],[2],'3g','Pi1',1,et_name=FAMILY)
    assert classify_topology(d).is_scaleless
    ts=enumerate_emission_types(3)
    ds=generate_diagrams(3,[t for t in ts if t.name==FAMILY])
    for d in ds:d.emission_type=''
    assert len(classify_orbits(ds,ts))==12
    assert {o.emission_type for o in classify_orbits(ds,ts)}=={FAMILY}


@pytest.mark.parametrize('mutation',['missing_link','duplicate_endpoint','extra_prop','wrong_operator'])
def test_reject_bad_source_topology(mutation):
    o=deepcopy(orbits()[0]);info=classify_topology(o.representative)
    if mutation=='missing_link':info.internal_gluon_links=[]
    if mutation=='duplicate_endpoint':info.internal_connections=[(1,2),(2,)]
    if mutation=='extra_prop':o.representative.props.append(deepcopy(o.representative.props[0]))
    if mutation=='wrong_operator':o.internal={'Gamma1':1}
    with pytest.raises(ValueError):build_definition(o,info)


def test_pdf_summary_counts_mirrors_across_edge_groups():
    from qgraf_pipeline.pdf_emitter import _emit_type_document
    text = _emit_type_document(orbits(), FAMILY, 3)
    assert text.count('3 cross-group partners') == 2
    assert text.count('1 self-mirror, 1 mirror pairs') == 2
