"""Tests for montopar.models."""
from montopar.models import (
    Vertex, Propagator, Bound, Pairing,
    EdgeDistribution, Orbit, HypercubeResult, CFactors, JFactors,
    TGVLeg, FGVLeg, TGVJFactors, FGVJFactors,
)


def test_vertex_creation():
    v = Vertex(line=3, lorentz="mu1", color="a1")
    assert v.line == 3
    assert v.lorentz == "mu1"
    assert v.color == "a1"


def test_propagator_creation():
    p = Propagator(loop="l1", lorentz=("mu1", "mu3"), color=("a1", "a3"))
    assert p.loop == "l1"
    assert p.lorentz == ("mu1", "mu3")


def test_pairing_creation():
    p = Pairing(edge_a=3, edge_b=4, loop="l1", s_a="s31", s_b="s41")
    assert p.edge_a == 3
    assert p.s_a == "s31"


def test_edge_distribution():
    ed = EdgeDistribution(counts={3: 2, 4: 2}, s_vars={3: ["s31", "s32"], 4: ["s41", "s42"]})
    assert ed.counts[3] == 2
    assert ed.s_vars[4] == ["s41", "s42"]


def test_hypercube_result():
    hc = HypercubeResult(
        s_to_x={"s31": "x1*x3", "s32": "x1"},
        jacobian="x1*x2",
        x_vars=["x1", "x2", "x3", "x4"],
    )
    assert hc.s_to_x["s32"] == "x1"
    assert len(hc.x_vars) == 4


def test_c_factors():
    c = CFactors(
        color="-1/2*((CA - 2*CF)*CF)",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-eps]^2/(4*Pi^(2-eps))^2",
        v14_scale="(-2*v14+I0)^(2*eps)",
        twoscale="(-2)^(-2)",
    )
    assert c.color.startswith("-1/2")


def test_j_factors():
    j = JFactors(
        lorentz_z="z1^2",
        jacobian="x1*x2",
        f_product="(x2*x4*z1 - x1*x2*x3*x4*z1)^(-1+eps) * (x2*z1 - x1*x2*z1)^(-1+eps)",
        f_polynomials={"f1": "x2*x4*z1 - x1*x2*x3*x4*z1", "f2": "x2*z1 - x1*x2*z1"},
        lorentz_v14_power="2",
    )
    assert "f1" in j.f_polynomials


def test_tgv_leg():
    leg = TGVLeg(edge=3, s_var="s31", loop="l1")
    assert leg.edge == 3
    assert leg.s_var == "s31"
    assert leg.loop == "l1"


def test_fgv_leg():
    leg = FGVLeg(edge=1, s_var="s11", loop="l1")
    assert leg.edge == 1
    assert leg.s_var == "s11"


def test_orbit_dispatch_defaults():
    """Orbit backward compat: dispatch/tgv_legs/fgv_legs default correctly."""
    orbit = Orbit(
        orbit_id="4e-01",
        emission_type="4e",
        n_gluons=2,
        edge_distribution=EdgeDistribution(counts={3: 2, 4: 2}, s_vars={3: ["s31", "s32"], 4: ["s41", "s42"]}),
        pairings=[],
        bounds=[],
        vertices=[],
        propagators=[],
    )
    assert orbit.dispatch == "pure"
    assert orbit.tgv_legs == []
    assert orbit.fgv_legs == []


def test_orbit_with_tgv_legs():
    leg1 = TGVLeg(edge=1, s_var="s11", loop="l1")
    leg2 = TGVLeg(edge=3, s_var="s31", loop="l1")
    orbit = Orbit(
        orbit_id="3e-3g-01",
        emission_type="3e-3g",
        n_gluons=3,
        edge_distribution=EdgeDistribution(counts={}, s_vars={}),
        pairings=[],
        bounds=[],
        vertices=[],
        propagators=[],
        dispatch="tgv",
        tgv_legs=[leg1, leg2],
    )
    assert orbit.dispatch == "tgv"
    assert len(orbit.tgv_legs) == 2
    assert orbit.tgv_legs[0].edge == 1


def test_c_factors_derivative():
    """CFactors backward compat: derivative defaults to None."""
    c = CFactors(
        color="CF",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-eps]^2",
        v14_scale="1",
        twoscale="1",
    )
    assert c.derivative is None

    c_tgv = CFactors(
        color="CF",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-eps]^2",
        v14_scale="1",
        twoscale="1",
        derivative="D[v14]",
    )
    assert c_tgv.derivative == "D[v14]"


def test_tgv_j_factors():
    j = TGVJFactors(
        jacobian="x2",
        pnum="alpha1*x1 + alpha2*x2",
        alpha_monomial="(alpha1*alpha2*alpha3)^(-eps)",
        S_polynomial="alpha1*x1 + alpha2*x2 + z1",
        denom_power="2 - 2*eps",
        x_vars=["x1", "x2", "x3"],
    )
    assert j.jacobian == "x2"
    assert len(j.x_vars) == 3


def test_fgv_j_factors():
    j = FGVJFactors(
        jacobian="x1",
        lorentzZ="-(1+4*z1+z1^2)/4",
        alpha_monomial="(alpha1*alpha2*alpha3*alpha4)^(-eps)",
        Q_prime_polynomial="alpha1*x1 + z1",
        denom_power="2 - 3*eps",
        x_vars=["x1", "x2", "x3", "x4"],
    )
    assert j.physical_zero is False
    assert j.denom_power == "2 - 3*eps"

    j_zero = FGVJFactors(
        jacobian="x1",
        lorentzZ="0",
        alpha_monomial="(alpha1*alpha2*alpha3*alpha4)^(-eps)",
        Q_prime_polynomial="0",
        denom_power="2 - 3*eps",
        x_vars=["x1", "x2", "x3", "x4"],
        physical_zero=True,
    )
    assert j_zero.physical_zero is True
