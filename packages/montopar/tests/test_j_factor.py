"""Tests for montopar.j_factor."""
from montopar.models import (
    Orbit, EdgeDistribution, Pairing, Bound, Vertex, Propagator,
    HypercubeResult, JFactors, TGVLeg, FGVLeg, TGVJFactors, FGVJFactors,
    MixedJFactors,
)
from montopar.j_factor import (
    compute_j,
    compute_f_p,
    compute_lorentz_factor,
    dot_product_factor,
    substitute_s_to_x,
    compute_j_tgv,
    compute_j_fgv,
    compute_j_mixed,
    differentiate_polynomial,
    expand_polynomial,
    linear_combination_polynomials,
)


def test_polynomial_power_is_accepted_on_reparse():
    assert expand_polynomial("z1^2 + (1+z1)^2") == "1 + 2*z1 + 2*z1^2"


def test_differentiate_polynomial_exactly():
    result = differentiate_polynomial("1/4*x1^2*z1 - x1*x2", "x1")
    assert result == "-x2 + 1/2*x1*z1"


def test_linear_combination_polynomials_exactly():
    result = linear_combination_polynomials([
        ("1/4*(1+z1)", "x1-x2"),
        ("-1/4", "x1+x2"),
    ])
    assert result == "-1/2*x2 + 1/4*x1*z1 - 1/4*x2*z1"


def test_dot_product_adjacent():
    """Adjacent edges: (3,4) → z, (4,1) → 1."""
    assert dot_product_factor(3, 4) == "z1"
    assert dot_product_factor(4, 1) == "1"
    assert dot_product_factor(1, 2) == "z1"
    assert dot_product_factor(2, 3) == "1"


def test_dot_product_opposite():
    """Opposite edges: (1,3) → -(1+z1), (2,4) → -(1+z1)."""
    assert dot_product_factor(1, 3) == "-(1+z1)"
    assert dot_product_factor(2, 4) == "-(1+z1)"


def test_compute_f_p_adjacent():
    """f_p for adjacent edges: (1-s_a)*s_b * (v_a·v_b/v_ref)."""
    result = compute_f_p(edge_a=3, edge_b=4, s_a="s31", s_b="s41")
    assert result == "(1-s31)*s41*z1"


def test_compute_f_p_adjacent_2_3():
    """f_p for adjacent edges {2,3}: (1-s2)*s3 * 1."""
    result = compute_f_p(edge_a=2, edge_b=3, s_a="s21", s_b="s31")
    assert result == "(1-s21)*s31*1"


def test_compute_f_p_opposite_1_3():
    """f_p for opposite edges {1,3}: z*(1-s1)*(1-s3) + s1*s3."""
    result = compute_f_p(edge_a=1, edge_b=3, s_a="s11", s_b="s31")
    assert result == "z1*(1-s11)*(1-s31) + s11*s31"


def test_compute_f_p_opposite_2_4():
    """f_p for opposite edges {2,4}: (1-s2)*(1-s4) + z*s2*s4."""
    result = compute_f_p(edge_a=2, edge_b=4, s_a="s21", s_b="s41")
    assert result == "(1-s21)*(1-s41) + z1*s21*s41"


def test_substitute_s_to_x():
    """Replace s variables with their x expressions."""
    s_to_x = {"s31": "x1*x3", "s32": "x1", "s41": "x2*x4", "s42": "x2"}
    expr = "(1-s31)*s41*z1"
    result = substitute_s_to_x(expr, s_to_x)
    # (1-x1*x3)*x2*x4*z1
    assert "x1*x3" in result
    assert "x2*x4" in result
    assert "s31" not in result
    assert "s41" not in result


def test_compute_lorentz_factor():
    """ℓ_n = product of (v_a·v_b/v_ref) for each pairing."""
    pairings = [
        Pairing(edge_a=3, edge_b=4, loop="l1", s_a="s31", s_b="s41"),
        Pairing(edge_a=3, edge_b=4, loop="l2", s_a="s32", s_b="s42"),
    ]
    result = compute_lorentz_factor(pairings)
    # Both (3,4) → z1, so product = z1 * z1 = z1^2
    assert result == "z1^2"


def test_compute_lorentz_factor_mixed():
    """Mixed adjacent and opposite pairings."""
    pairings = [
        Pairing(edge_a=3, edge_b=4, loop="l1", s_a="s31", s_b="s41"),
        Pairing(edge_a=1, edge_b=3, loop="l2", s_a="s11", s_b="s32"),
    ]
    result = compute_lorentz_factor(pairings)
    # (3,4) → z1, (1,3) → -(1+z1)
    assert "z1" in result
    assert "(1+z1)" in result


def test_compute_j_4e01():
    """Full J_n computation for 4e-01."""
    orbit = Orbit(
        orbit_id="4e-01",
        emission_type="4e",
        n_gluons=2,
        edge_distribution=EdgeDistribution(
            counts={3: 2, 4: 2},
            s_vars={3: ["s31", "s32"], 4: ["s41", "s42"]},
        ),
        pairings=[
            Pairing(edge_a=3, edge_b=4, loop="l1", s_a="s31", s_b="s41"),
            Pairing(edge_a=3, edge_b=4, loop="l2", s_a="s32", s_b="s42"),
        ],
        bounds=[
            Bound(var="s31", low="0", high="s32"),
            Bound(var="s41", low="0", high="s42"),
            Bound(var="s32", low="0", high="1"),
            Bound(var="s42", low="0", high="1"),
        ],
        vertices=[],
        propagators=[],
    )
    hc = HypercubeResult(
        s_to_x={"s32": "x1", "s42": "x2", "s31": "x1*x3", "s41": "x2*x4"},
        jacobian="x1*x2",
        x_vars=["x1", "x2", "x3", "x4"],
    )

    j = compute_j(orbit, hc)

    # Lorentz factor: z1 * z1 = z1^2
    assert j.lorentz_z == "z1^2"
    assert j.lorentz_v14_power == "2"

    # Jacobian passed through
    assert j.jacobian == "x1*x2"

    # f1 = (1-s31)*s41*z1 → (1-x1*x3)*x2*x4*z1
    assert "f1" in j.f_polynomials
    assert "f2" in j.f_polynomials

    # f_product should contain both f1 and f2 with (-1+eps) exponent
    assert "(-1+eps)" in j.f_product
    assert "f1" in j.f_product
    assert "f2" in j.f_product


from montopar.j_factor import compute_j_expanded


def test_compute_j_expanded_single():
    """单个 f_p 被正确展开。"""
    j = JFactors(
        lorentz_z="z1",
        jacobian="1",
        f_product="f1^(-1+eps)",
        f_polynomials={"f1": "x2*z1 - x1*x2*z1"},
        lorentz_v14_power="1",
    )
    result = compute_j_expanded(j)
    assert result == "(x2*z1 - x1*x2*z1)^(-1+eps)"
    assert "f1" not in result


def test_compute_j_expanded_multiple():
    """多个 f_p 全部展开，f_name 按长度降序替换避免部分匹配。"""
    j = JFactors(
        lorentz_z="z1^2",
        jacobian="x1*x2",
        f_product="f1^(-1+eps) * f2^(-1+eps)",
        f_polynomials={
            "f1": "x2*x4*z1 - x1*x2*x3*x4*z1",
            "f2": "x2*z1 - x1*x2*z1",
        },
        lorentz_v14_power="2",
    )
    result = compute_j_expanded(j)
    assert "f1" not in result
    assert "f2" not in result
    assert "(x2*x4*z1 - x1*x2*x3*x4*z1)^(-1+eps)" in result
    assert "(x2*z1 - x1*x2*z1)^(-1+eps)" in result


def test_compute_j_expanded_4e01_golden():
    """4e-01 golden 值：展开后与已知结果一致。"""
    j = JFactors(
        lorentz_z="z1^2",
        jacobian="x1*x2",
        f_product="f1^(-1+eps) * f2^(-1+eps)",
        f_polynomials={
            "f1": "x2*x4*z1 - x1*x2*x3*x4*z1",
            "f2": "x2*z1 - x1*x2*z1",
        },
        lorentz_v14_power="2",
    )
    result = compute_j_expanded(j)
    expected = "(x2*x4*z1 - x1*x2*x3*x4*z1)^(-1+eps) * (x2*z1 - x1*x2*z1)^(-1+eps)"
    assert result == expected


class TestComputeJTGV:
    def test_3e3g_01_j(self):
        """TGV J for 3e-3g-01 should produce Pnum and S in x-variables."""
        orbit = Orbit(
            orbit_id="3e-3g-01",
            emission_type="3e-3g",
            n_gluons=3,
            dispatch="tgv",
            edge_distribution=EdgeDistribution(
                counts={3: 1, 4: 2},
                s_vars={3: ["s31"], 4: ["s41", "s42"]},
            ),
            tgv_legs=[
                TGVLeg(edge=3, s_var="s31", loop="l1"),
                TGVLeg(edge=4, s_var="s41", loop="l2"),
                TGVLeg(edge=4, s_var="s42", loop="l3"),
            ],
            pairings=[], bounds=[], vertices=[], propagators=[], source={},
        )
        hc = HypercubeResult(
            s_to_x={"s31": "x1", "s42": "x2", "s41": "x2*x3"},
            jacobian="x2",
            x_vars=["x1", "x2", "x3"],
        )
        j = compute_j_tgv(orbit, hc)
        assert isinstance(j, TGVJFactors)
        assert j.jacobian == "x2"
        assert j.alpha_monomial == "(alpha1*alpha2*alpha3)^(-eps)"
        assert j.denom_power == "2 - 2*eps"
        assert j.x_vars == ["x1", "x2", "x3"]
        # Pnum should contain alpha and x variables (not s)
        assert "alpha1" in j.pnum
        assert "s31" not in j.pnum  # s-vars replaced by x
        assert "s41" not in j.pnum
        # S should contain alpha and x variables
        assert "alpha1" in j.S_polynomial
        assert "s31" not in j.S_polynomial


class TestComputeJFGV:
    def test_4e4g_02_j(self):
        """FGV J for 4e-4g-02 should produce Q' and lorentzZ."""
        orbit = Orbit(
            orbit_id="4e-4g-02",
            emission_type="4e-4g",
            n_gluons=4,
            dispatch="fourv",
            edge_distribution=EdgeDistribution(
                counts={3: 2, 4: 2},
                s_vars={3: ["s31", "s32"], 4: ["s41", "s42"]},
            ),
            fgv_legs=[
                FGVLeg(edge=3, s_var="s31", loop="l1"),
                FGVLeg(edge=3, s_var="s32", loop="l2"),
                FGVLeg(edge=4, s_var="s41", loop="l3"),
                FGVLeg(edge=4, s_var="s42", loop="l4"),
            ],
            pairings=[], bounds=[], vertices=[], propagators=[], source={},
        )
        hc = HypercubeResult(
            s_to_x={
                "s32": "x1", "s42": "x2",
                "s31": "x1*x3", "s41": "x2*x4",
            },
            jacobian="x1*x2",
            x_vars=["x1", "x2", "x3", "x4"],
        )
        j = compute_j_fgv(orbit, hc)
        assert isinstance(j, FGVJFactors)
        assert j.jacobian == "x1*x2"
        assert j.alpha_monomial == "(alpha1*alpha2*alpha3*alpha4)^(-eps)"
        assert j.denom_power == "2 - 3*eps"
        assert j.x_vars == ["x1", "x2", "x3", "x4"]
        assert j.physical_zero is False
        # lorentzZ should be a z polynomial
        assert "z1" in j.lorentzZ
        # Q' should contain alpha and x variables
        assert "alpha1" in j.Q_prime_polynomial
        assert "s31" not in j.Q_prime_polynomial


class TestComputeJMixed:
    def test_5e3g_01_j(self):
        """5e-3g-01: pure(3,4) + TGV(3,4,4)."""
        import json
        from pathlib import Path
        from montopar.parser import parse_orbit
        from montopar.hypercube import hypercube_transform

        fixture = Path(__file__).parent / "fixtures" / "5e-3g-01_orbit.json"
        with open(fixture) as f:
            orbit = parse_orbit(json.load(f))

        hc = hypercube_transform(orbit.bounds)
        tgv_loops = {leg.loop for leg in orbit.tgv_legs}
        pure = [p for p in orbit.pairings if p.loop not in tgv_loops]

        j = compute_j_mixed(orbit, hc, pure)

        # Pure part
        assert j.lorentz_z == "z1"
        assert j.n_pure == 1
        assert "f1" in j.f_polynomials
        assert "f1^(-1+eps)" in j.f_product

        # TGV part
        assert "alpha1" in j.pnum
        assert j.alpha_monomial == "(alpha1*alpha2*alpha3)^(-eps)"
        assert j.denom_power == "2 - 2*eps"

        # Shared
        assert len(j.x_vars) == 5  # 3+2*1 = 5
        assert "alpha2" in j.S_polynomial


def test_mixed_pure_part_matches_2e():
    """5e-3g-01 pure part f_p matches 2e-01 (up to variable renaming)."""
    import json
    from pathlib import Path
    from montopar.parser import parse_orbit

    fixtures = Path(__file__).parent / "fixtures"

    # 5e-3g-01 pure part
    with open(fixtures / "5e-3g-01_orbit.json") as f:
        o5 = parse_orbit(json.load(f))
    tgv_loops = {leg.loop for leg in o5.tgv_legs}
    pure_5 = [p for p in o5.pairings if p.loop not in tgv_loops][0]
    f_5 = compute_f_p(pure_5.edge_a, pure_5.edge_b, pure_5.s_a, pure_5.s_b)

    # 2e-01 uses s31, s41; 5e-3g uses s32, s43
    f_5_renamed = f_5.replace("s32", "s31").replace("s43", "s41")
    expected = "(1-s31)*s41*z1"
    assert f_5_renamed == expected, f"Got {f_5_renamed}, expected {expected}"


def test_mixed_tgv_part_matches_3e3g():
    """5e-3g-01 TGV S and Pnum match 3e-3g-01."""
    import json
    from pathlib import Path
    from montopar.parser import parse_orbit
    from montopar.tgv import compute_S, compute_Pnum

    fixtures = Path(__file__).parent / "fixtures"

    with open(fixtures / "5e-3g-01_orbit.json") as f:
        o5 = parse_orbit(json.load(f))
    with open(fixtures / "3e-3g-01_orbit.json") as f:
        o3 = parse_orbit(json.load(f))

    assert compute_S(o5.tgv_legs) == compute_S(o3.tgv_legs)
    assert compute_Pnum(o5.tgv_legs) == compute_Pnum(o3.tgv_legs)
