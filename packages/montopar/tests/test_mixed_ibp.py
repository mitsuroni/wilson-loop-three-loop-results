"""Regression cases for endpoint interleaving and self-contained mixed exports."""
import pytest
from montopar.models import Orbit, Bound, TGVLeg, Pairing, EdgeDistribution
from montopar.mixed_ibp import ordered_chains, oriented_faces, compute_mixed_ibp
from montopar.emitter import _build_c_mixed_block, _build_jraw_s_mixed
from montopar.c_factor import compute_c_mixed

def orbit(interleaved=False):
    # Pure endpoint follows (01) or separates (02) the repeated TGV legs.
    variables=["s31","s32","s41","s42","s43"]
    bounds=[Bound("s31","0","s32"),Bound("s32","0","1"),
            Bound("s41","0","s42"),Bound("s42","0","s43"),Bound("s43","0","1")]
    return Orbit("test","5e-3g",1,EdgeDistribution({3:2,4:3},{3:variables[:2],4:variables[2:]}),
        [Pairing(3,4,"pure","s32","s42" if interleaved else "s43")],bounds,[],[],
        dispatch="mixed",tgv_legs=[TGVLeg(3,"s31","a"),TGVLeg(4,"s41","b"),
                                  TGVLeg(4,"s43" if interleaved else "s42","c")])

def constraints(result):
    return {tuple({**t["fixed_parameters"],**t["identified_parameters"]}.items())[0]:t["coefficient"] for t in result["terms"]}

def test_following_endpoint_replaces_spurious_upper_cusp():
    r=compute_mixed_ibp(orbit(),"C0")
    assert constraints(r)=={("s41","0"):"z1",("s41","s42"):"-2*z1",("s42","s43"):"z1"}
    assert len(r["terms"])==3
    for t in r["terms"]:
        assert t["ordered_bounds"]["s31"]==["0","s32"]
        assert len(t["variables"])==4
        assert t["jacobian"]!="1"
        assert "f1" not in t["integrand"]

def test_interleaving_creates_two_adjacent_contacts():
    r=compute_mixed_ibp(orbit(True),"C0")
    assert constraints(r)=={("s41","0"):"z1",("s43","1"):"z1",
                            ("s41","s42"):"-z1",("s42","s43"):"-z1"}
    assert all(t["C_ref"]=="C_star" for t in r["terms"])
    assert r["operator"]["pnum_recoil"]=="0"

@pytest.mark.parametrize("bounds",[
    [Bound("a","0","b"),Bound("b","0","a")],
    [Bound("a","0","c"),Bound("b","0","c"),Bound("c","0","1")],
    [Bound("a","x","1")],
    [Bound("a","0","missing")],
    [Bound("a","0","1"),Bound("a","0","1")],
])
def test_malformed_domains_rejected(bounds):
    with pytest.raises(ValueError):ordered_chains(bounds)

def test_cannot_split_ordered_edge_into_cube():
    o=orbit()
    o.bounds=[Bound(b.var,"0","1") for b in o.bounds]
    with pytest.raises(ValueError,match="split"):compute_mixed_ibp(o,"C0")

def test_endpoint_overlap_rejected():
    o=orbit();o.pairings[0].s_a="s31"
    with pytest.raises(ValueError,match="distinct"):compute_mixed_ibp(o,"C0")

@pytest.mark.parametrize("n",range(5))
def test_eikonal_plus_fourier_phase_is_independent_of_n(n):
    c=_build_c_mixed_block(compute_c_mixed(n),n)
    assert c["phase_provenance"]["coordinate_propagators"]==f"(-1)^{n+3}"
    assert c["components"]["i_factor"]=="I"
    assert c["component_semantics"]=="multiplicative"
    assert (1j)**(2*n+3)*(-1)**n==(1j)**3

def test_raw_s_requires_own_f_polynomials():
    with pytest.raises(ValueError,match="f_polynomials"):
        _build_jraw_s_mixed("z1","f1^(-1+eps)","P","S",[],[])
