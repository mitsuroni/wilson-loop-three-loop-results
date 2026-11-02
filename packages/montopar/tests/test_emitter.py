"""Tests for montopar.emitter."""
from montopar.models import (
    Orbit, EdgeDistribution, Pairing, Bound, Vertex, Propagator,
    CFactors, JFactors, HypercubeResult,
    TGVJFactors, FGVJFactors, MixedJFactors, TGVIBPResult, TGVIBPTerm,
)
from montopar.c_factor import compute_c_star
from montopar.emitter import (
    emit_montopar_json,
    build_kinematics_map,
    build_c_components,
    build_j_components,
    combine_c,
)


def test_build_kinematics_map():
    """N=4 kinematics map."""
    km = build_kinematics_map()
    assert km["SP[v1,v4]"] == "v14"
    assert km["SP[v3,v4]"] == "z1*v14"
    assert km["SP[v1,v2]"] == "z1*v14"
    assert km["SP[v1,v3]"] == "(-1 - z1)*v14"
    assert km["SP[v2,v3]"] == "v14"
    assert km["SP[v2,v4]"] == "(-1 - z1)*v14"


def test_build_c_components():
    c = CFactors(
        color="-1/2*((CA - 2*CF)*CF)",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-eps]^2/(4*Pi^(2-eps))^2",
        v14_scale="(-2*v14+I0)^(2*eps)",
        twoscale="(-2)^(-2)",
    )
    components = build_c_components(c, n=2)
    assert components["color"] == c.color
    assert components["coupling"] == c.coupling
    assert "ft_prefactor" in components
    assert "ft_v14_scale" in components


def test_combine_c():
    c = CFactors(
        color="-1/2*((CA - 2*CF)*CF)",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-eps]^2/(4*Pi^(2-eps))^2",
        v14_scale="(-2*v14+I0)^(2*eps)",
        twoscale="(-2)^(-2)",
    )
    combined = combine_c(c)
    assert "alphaS" in combined
    assert "Gamma" in combined


def test_build_j_components():
    j = JFactors(
        lorentz_z="z1^2",
        jacobian="x1*x2",
        f_product="f1^(-1+eps) * f2^(-1+eps)",
        f_polynomials={"f1": "x2*x4*z1 - x1*x2*x3*x4*z1", "f2": "x2*z1 - x1*x2*z1"},
        lorentz_v14_power="2",
    )
    components = build_j_components(j)
    assert components["lorentz_z"] == "z1^2"
    assert components["jacobian"] == "x1*x2"
    assert "f1^(-1+eps)" in components["f_product"]


def test_emit_montopar_json_pure():
    """montopar.json for pure dispatch should have orbit_id, dispatch, C, J_raw."""
    orbit = Orbit(
        orbit_id="4e-01", emission_type="4e", n_gluons=2,
        edge_distribution=EdgeDistribution(counts={3: 2, 4: 2}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[],
        source={"orbit_id": "4e-01"},
    )
    c = CFactors(
        color="-1/2*((CA - 2*CF)*CF)",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-eps]^2/(4*Pi^(2-eps))^2",
        v14_scale="(-2*v14+I0)^(2*eps)",
        twoscale="(-2)^(-2)",
    )
    j = JFactors(
        lorentz_z="z1^2", jacobian="x1*x2",
        f_product="f1^(-1+eps) * f2^(-1+eps)",
        f_polynomials={"f1": "expr1", "f2": "expr2"},
        lorentz_v14_power="2",
    )
    hc = HypercubeResult(
        s_to_x={}, jacobian="x1*x2", x_vars=["x1", "x2", "x3", "x4"],
    )

    result = emit_montopar_json(orbit, c, j, hc)

    # Top-level keys
    assert result["orbit_id"] == "4e-01"
    assert result["dispatch"] == "pure"
    assert result["n_gluons"] == 2
    # Geometry
    assert result["geometry"]["n_lines"] == 4
    assert result["geometry"]["reference_scale"] == "v14"
    # C block
    assert "components" in result["C"]
    assert "combined" in result["C"]
    assert result["C"]["components"]["color"] == c.color
    # J_raw block
    assert "integrand" in result["J_raw"]
    assert "expanded" in result["J_raw"]
    assert "variables" in result["J_raw"]
    assert result["J_raw"]["variables"] == ["x1", "x2", "x3", "x4"]
    # J_raw.expanded should not have f symbols
    assert "f1" not in result["J_raw"]["expanded"]
    assert "f2" not in result["J_raw"]["expanded"]


def test_emit_montopar_json_tgv():
    """montopar.json for TGV dispatch should have derivative in C, alpha_vars in J_raw."""
    orbit = Orbit(
        orbit_id="3e-3g-01", emission_type="3e-3g", n_gluons=3,
        dispatch="tgv",
        edge_distribution=EdgeDistribution(counts={3: 1, 4: 2}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[], source={},
    )
    c = CFactors(
        color="I/2*CA*CF", coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-2*eps]/(64*Pi^(4-2*eps))",
        v14_scale="(-2*v14+I0)^(2*eps)", twoscale="1",
        derivative="2*(1-2*eps)",
    )
    j = TGVJFactors(
        jacobian="x2", pnum="P",
        alpha_monomial="(alpha1*alpha2*alpha3)^(-eps)",
        S_polynomial="S", denom_power="2 - 2*eps",
        x_vars=["x1", "x2", "x3"],
    )
    hc = HypercubeResult(s_to_x={}, jacobian="x2", x_vars=["x1", "x2", "x3"])

    result = emit_montopar_json(orbit, c, j, hc)

    assert result["dispatch"] == "tgv"
    assert result["C"]["components"]["derivative"] == "2*(1-2*eps)"
    assert "2*(1-2*eps)" in result["C"]["combined"]
    assert result["J_raw"]["alpha_vars"] == ["alpha1", "alpha2", "alpha3"]
    assert "P" in result["J_raw"]["integrand"]


def test_emit_montopar_json_tgv_with_ibp():
    orbit = Orbit(
        orbit_id="3e-3g-test", emission_type="3e-3g", n_gluons=3,
        dispatch="tgv",
        edge_distribution=EdgeDistribution(counts={3: 1, 4: 2}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[], source={},
    )
    c = CFactors(
        color="I/2*CA*CF", coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-2*eps]/(64*Pi^(4-2*eps))",
        v14_scale="(-2*v14+I0)^(2*eps)", twoscale="1",
        derivative="2*(1-2*eps)",
    )
    j = TGVJFactors(
        jacobian="1", pnum="P", alpha_monomial="A",
        S_polynomial="S", denom_power="2 - 2*eps", x_vars=["x1", "x2", "x3"],
    )
    hc = HypercubeResult(s_to_x={}, jacobian="1", x_vars=["x1", "x2", "x3"])
    term = TGVIBPTerm(
        term_id="boundary_s1_0", kind="cusp_boundary", coefficient="-z1",
        c_ref="C_star", combined_C="(-z1)*(Cstar)", integrand="K",
        variables=["x1", "x2"], alpha_vars=["alpha1", "alpha2", "alpha3"],
        bounds={"x1": [0, 1], "x2": [0, 1]}, delta_constraints={},
        denom_power="1 - 2*eps", fixed_parameters={"s1": "0"},
    )
    ibp = TGVIBPResult(
        method="contour_ibp", status="complete", case="nondegenerate",
        C_star=compute_c_star(), terms=[term], sum_expression="term(boundary_s1_0)",
        certificates={"pnum_identity": True}, derivative_vector=["z1", "-1-z1", "1"],
    )

    result = emit_montopar_json(orbit, c, j, hc, ibp=ibp)

    assert result["J_ibp"]["status"] == "complete"
    assert result["J_ibp"]["terms"][0]["C_ref"] == "C_star"
    assert result["J_ibp"]["terms"][0]["fixed_parameters"] == {"s1": "0"}
    assert "256*Pi" in result["J_ibp"]["prefactors"]["C_star"]


def test_emit_montopar_json_fgv_physical_zero():
    """FGV physical zero should produce C=0, J_raw with physical_zero=True."""
    orbit = Orbit(
        orbit_id="4e-4g-01", emission_type="4e-4g", n_gluons=4,
        dispatch="fourv",
        edge_distribution=EdgeDistribution(counts={3: 1, 4: 3}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[], source={},
    )
    c = CFactors(color="0", coupling="0", ft_base="0", v14_scale="0", twoscale="0")
    j = FGVJFactors(
        jacobian="0", lorentzZ="0", alpha_monomial="0",
        Q_prime_polynomial="0", denom_power="0",
        x_vars=[], physical_zero=True,
    )
    hc = HypercubeResult(s_to_x={}, jacobian="0", x_vars=[])

    result = emit_montopar_json(orbit, c, j, hc)

    assert result["dispatch"] == "fourv"
    assert result["C"]["combined"] == "0"
    assert result["J_raw"]["physical_zero"] is True
    assert result["J_raw"]["integrand"] == "0"


def test_emit_montopar_json_fgv_nonzero():
    """FGV non-zero should have lorentzZ in J_raw integrand."""
    orbit = Orbit(
        orbit_id="4e-4g-02", emission_type="4e-4g", n_gluons=4,
        dispatch="fourv",
        edge_distribution=EdgeDistribution(counts={3: 1, 4: 3}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[], source={},
    )
    c = CFactors(
        color="CA^2*CF", coupling="(4*Pi*alphaS)^3",
        ft_base="Gamma[1-3*eps]/(64*Pi^(6-3*eps))",
        v14_scale="(-2*v14+I0)^(3*eps)", twoscale="1",
    )
    j = FGVJFactors(
        jacobian="x1*x2", lorentzZ="-(1+4*z1+z1^2)/4",
        alpha_monomial="(alpha1*alpha2*alpha3*alpha4)^(-eps)",
        Q_prime_polynomial="Qp", denom_power="2 - 3*eps",
        x_vars=["x1", "x2", "x3", "x4"],
    )
    hc = HypercubeResult(s_to_x={}, jacobian="x1*x2", x_vars=["x1", "x2", "x3", "x4"])

    result = emit_montopar_json(orbit, c, j, hc)

    assert result["dispatch"] == "fourv"
    assert result["J_raw"]["variables"] == ["x1", "x2", "x3", "x4"]
    assert "Qp" in result["J_raw"]["integrand"]


def test_emit_mixed():
    """emit_montopar_json with mixed dispatch."""
    orbit = Orbit(
        orbit_id="test-mixed", emission_type="5e-3g", n_gluons=4,
        edge_distribution=EdgeDistribution(counts={3: 2, 4: 3}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[],
        dispatch="mixed", tgv_legs=[],
    )
    c = CFactors(
        color="CF", coupling="(4*Pi*alphaS)^3",
        ft_base="Gamma[1-eps]^1*Gamma[1-2*eps]/((4*Pi^(2-eps))^1*64*Pi^(4-2*eps))",
        v14_scale="(-2*v14+I0)^(3*eps)", twoscale="(-2)^(-1)",
        derivative="2*(1-2*eps)",
    )
    j = MixedJFactors(
        lorentz_z="z1", jacobian="x1*x2^2",
        f_product="f1^(-1+eps)", f_polynomials={"f1": "x2*z1 - x1*x2*z1"},
        pnum="N", alpha_monomial="(alpha1*alpha2*alpha3)^(-eps)",
        S_polynomial="S", denom_power="2 - 2*eps",
        x_vars=["x1", "x2", "x3", "x4", "x5"], n_pure=1,
    )
    hc = HypercubeResult(
        s_to_x={}, jacobian="x1*x2^2",
        x_vars=["x1", "x2", "x3", "x4", "x5"],
    )
    result = emit_montopar_json(orbit, c, j, hc)
    assert result["dispatch"] == "mixed"
    assert result["n_pure"] == 1
    assert result["C"]["components"]["derivative"] == "2*(1-2*eps)"
    assert result["J_raw"]["pure"]["lorentz_z"] == "z1"
    assert result["J_raw"]["tgv"]["pnum"] == "N"
    assert "alpha_vars" in result["J_raw"]
