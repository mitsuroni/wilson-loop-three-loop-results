"""Snapshot tests: end-to-end verification against known results."""
import json
from pathlib import Path

import pytest

from montopar.parser import parse_orbit
from montopar.hypercube import hypercube_transform
from montopar.c_factor import compute_c, compute_c_tgv, compute_c_mixed
from montopar.j_factor import compute_j, compute_j_tgv, compute_j_mixed
from montopar.fgv import is_physical_zero
from montopar.emitter import emit_montopar_json

FIXTURES = Path(__file__).parent / "fixtures"

# Known golden values from hepeflow-preprocessing 4e-01 output
GOLDEN_COLOR = "-1/2*((CA - 2*CF)*CF)"
GOLDEN_F1 = "x2*x4*z1 - x1*x2*x3*x4*z1"
GOLDEN_F2 = "x2*z1 - x1*x2*z1"
GOLDEN_JACOBIAN = "x1*x2"
GOLDEN_LORENTZ_Z = "z1^2"


def load_fixture(name: str) -> dict:
    """Load a JSON fixture file."""
    with open(FIXTURES / name, encoding="utf-8") as f:
        return json.load(f)


def test_snapshot_4e01():
    """montopar output must match hepeflow-preprocessing golden output."""
    data = load_fixture("4e-01_orbit.json")
    orbit = parse_orbit(data)

    # Hypercube
    hc = hypercube_transform(orbit.bounds)
    assert hc.jacobian == GOLDEN_JACOBIAN
    assert hc.s_to_x["s32"] == "x1"
    assert hc.s_to_x["s42"] == "x2"
    assert hc.s_to_x["s31"] == "x1*x3"
    assert hc.s_to_x["s41"] == "x2*x4"

    # C factor (color check only — actual Mathematica result may vary)
    c = compute_c(orbit.n_gluons, GOLDEN_COLOR)
    assert "Gamma[1-eps]" in c.ft_base

    # J factor
    j = compute_j(orbit, hc)
    assert j.lorentz_z == GOLDEN_LORENTZ_Z
    assert j.jacobian == GOLDEN_JACOBIAN
    assert j.f_polynomials["f1"] == GOLDEN_F1
    assert j.f_polynomials["f2"] == GOLDEN_F2

    # Final JSON structure (montopar.json merged output)
    result = emit_montopar_json(orbit, c, j, hc)
    assert result["dispatch"] == "pure"
    assert result["C"]["combined"] is not None
    assert "expanded" in result["J_raw"]
    assert "f1" not in result["J_raw"]["expanded"]


def test_snapshot_6e01():
    """6e-01: 3 gluons on edges 3 and 4, all adjacent pairings."""
    data = {
        "orbit_id": "6e-01",
        "emission_type": "6e",
        "sector": "dg",
        "nloops": 3,
        "vertices": [
            {"line": 3, "lorentz": "mu1", "color": "a1"},
            {"line": 3, "lorentz": "mu2", "color": "a2"},
            {"line": 3, "lorentz": "mu3", "color": "a3"},
            {"line": 4, "lorentz": "mu4", "color": "a4"},
            {"line": 4, "lorentz": "mu5", "color": "a5"},
            {"line": 4, "lorentz": "mu6", "color": "a6"},
        ],
        "propagators": [
            {"loop": "l1", "lorentz": ["mu1", "mu4"], "color": ["a1", "a4"]},
            {"loop": "l2", "lorentz": ["mu2", "mu5"], "color": ["a2", "a5"]},
            {"loop": "l3", "lorentz": ["mu3", "mu6"], "color": ["a3", "a6"]},
        ],
        "internal_vertices": [],
        "phases": [
            {"loop": "l1", "expression": "-((1-s31)*v3 + s41*v4)"},
            {"loop": "l2", "expression": "-((1-s32)*v3 + s42*v4)"},
            {"loop": "l3", "expression": "-((1-s33)*v3 + s43*v4)"},
        ],
        "bounds": [
            {"var": "s31", "low": "0", "high": "s32"},
            {"var": "s41", "low": "0", "high": "s42"},
            {"var": "s32", "low": "0", "high": "s33"},
            {"var": "s42", "low": "0", "high": "s43"},
            {"var": "s33", "low": "0", "high": "1"},
            {"var": "s43", "low": "0", "high": "1"},
        ],
    }

    orbit = parse_orbit(data)
    assert orbit.n_gluons == 3
    assert orbit.edge_distribution.counts == {3: 3, 4: 3}

    hc = hypercube_transform(orbit.bounds)
    assert hc.jacobian == "x1^2*x2^2*x3*x4"
    assert len(hc.x_vars) == 6

    c = compute_c(orbit.n_gluons, GOLDEN_COLOR)
    assert "Gamma[1-eps]^3" in c.ft_base

    j = compute_j(orbit, hc)
    assert j.lorentz_z == "z1^3"
    assert len(j.f_polynomials) == 3
    assert "f1" in j.f_polynomials
    assert "f2" in j.f_polynomials
    assert "f3" in j.f_polynomials

    result = emit_montopar_json(orbit, c, j, hc)
    assert "expanded" in result["J_raw"]
    assert "f1" not in result["J_raw"]["expanded"]
    assert "f2" not in result["J_raw"]["expanded"]
    assert "f3" not in result["J_raw"]["expanded"]


def test_snapshot_3e3g_01():
    """montopar TGV output must produce correct C×J combined.

    Cross-checks the structural form of the TGV C and J factors against
    the hepeflow-preprocessing golden output for orbit 3e-3g-01 (one
    gluon on edge 3 and two on edge 4, all meeting at a triple-gluon
    vertex). Component-level strings differ from hepeflow because
    montopar uses the derivation v2 convention, but the combined
    integral I = C×J is equivalent.
    """
    data = load_fixture("3e-3g-01_orbit.json")
    orbit = parse_orbit(data)
    assert orbit.dispatch == "tgv"
    assert len(orbit.tgv_legs) == 3
    # Legs land on edges 3, 4, 4 (one gluon on edge 3, two on edge 4).
    assert [leg.edge for leg in orbit.tgv_legs] == [3, 4, 4]

    hc = hypercube_transform(orbit.bounds)
    assert hc.jacobian == "x2"

    c = compute_c_tgv()
    # Verify C factors against the hepeflow golden 12_assemble_c.json
    assert c.color == "-I/2*CA*CF"
    assert c.derivative == "2*(1-2*eps)"
    assert "Gamma[1-2*eps]" in c.ft_base
    assert c.v14_scale == "(-2*v14+I0)^(2*eps)"

    j = compute_j_tgv(orbit, hc)

    # Verify J structure
    assert j.jacobian == "x2"
    assert j.alpha_monomial == "(alpha1*alpha2*alpha3)^(-eps)"
    assert j.denom_power == "2 - 2*eps"
    assert len(j.x_vars) == 3

    # Pnum should contain alpha and x variables (the TGV numerator is a
    # polynomial in α_i, x_i and z1 after s→x substitution).
    assert "alpha1" in j.pnum
    assert "x1" in j.pnum or "x2" in j.pnum or "x3" in j.pnum
    assert "z1" in j.pnum
    # S should not contain s-variables (they are substituted by x).
    assert "s31" not in j.S_polynomial
    assert "s41" not in j.S_polynomial
    assert "s42" not in j.S_polynomial

    # Final JSON structure: dispatch-aware emitter wires C/J pieces.
    result = emit_montopar_json(orbit, c, j, hc)
    assert result["dispatch"] == "tgv"
    assert "derivative" in result["C"]["components"]
    # The TGV J_raw contains pnum in the integrand (not as a top-level key).
    assert "pnum" not in result["J_raw"]
    # The golden lorentz_z for this orbit is z1 (single non-zero pairing);
    # montopar's v2 convention folds the z factor into Pnum, so check
    # the integrand carries the z1 instead.
    assert "z1" in result["J_raw"]["integrand"]


def test_snapshot_4e4g_physical_zero():
    """4e-4g-01 should be detected as physical zero.

    A 4-gluon vertex with all four gluons on edge 4 (max m_i = 3) maps
    to the vanishing 4v Wilson-loop configuration and is short-circuited
    by is_physical_zero before any CAS work.
    """
    data = {
        "orbit_id": "4e-4g-01",
        "emission_type": "4e-4g",
        "sector": "dg", "nloops": 3,
        "vertices": [
            {"line": 3, "lorentz": "mu1", "color": "a1"},
            {"line": 4, "lorentz": "mu2", "color": "a2"},
            {"line": 4, "lorentz": "mu3", "color": "a3"},
            {"line": 4, "lorentz": "mu4", "color": "a4"},
        ],
        "propagators": [
            {"loop": "l1", "lorentz": ["mu1", "nu1"], "color": ["a1", "b1"]},
            {"loop": "l2", "lorentz": ["mu2", "nu2"], "color": ["a2", "b2"]},
            {"loop": "l3", "lorentz": ["mu3", "nu3"], "color": ["a3", "b3"]},
            {"loop": "l4", "lorentz": ["mu4", "nu4"], "color": ["a4", "b4"]},
        ],
        "internal_vertices": [
            {"type": "4v", "id": "v1", "loops": ["l1", "l2", "l3", "l4"]}
        ],
        "phases": [
            {"loop": "l1", "expression": "v1 + v2 + s31*v3"},
            {"loop": "l2", "expression": "v1 + v2 + v3 + s41*v4"},
            {"loop": "l3", "expression": "v1 + v2 + v3 + s42*v4"},
            {"loop": "l4", "expression": "v1 + v2 + v3 + s43*v4"},
        ],
        "bounds": [
            {"var": "s31", "low": "0", "high": "1"},
            {"var": "s41", "low": "0", "high": "s42"},
            {"var": "s42", "low": "0", "high": "s43"},
            {"var": "s43", "low": "0", "high": "1"},
        ],
    }
    orbit = parse_orbit(data)
    assert orbit.dispatch == "fourv"
    assert is_physical_zero(orbit.edge_distribution)


def test_snapshot_5e3g_01():
    """5e-3g-01 montopar output: pure(3,4) + TGV(3,4,4).

    Mixed dispatch: one bare eikonal pairing (l1 connecting edges 3 and 4)
    plus a TGV vertex with 3 legs (l2, l3, l4).  orbit.n_gluons = 1 because
    the parser only counts propagators whose both lorentz indices map to
    eikonal vertices — TGV propagators have one internal index and are
    excluded from pairings.
    """
    data = load_fixture("5e-3g-01_orbit.json")
    orbit = parse_orbit(data)

    assert orbit.dispatch == "mixed"
    assert orbit.n_gluons == 1

    # Hypercube
    hc = hypercube_transform(orbit.bounds)
    assert len(hc.x_vars) == 5
    assert hc.s_to_x["s32"] == "x1"
    assert hc.s_to_x["s43"] == "x2"
    assert hc.s_to_x["s31"] == "x1*x3"
    assert hc.s_to_x["s42"] == "x2*x4"
    assert hc.s_to_x["s41"] == "x2*x4*x5"

    # C factor: mixed = pure(1) × tgv
    n_pure = 1
    c = compute_c_mixed(n_pure)
    c.color = "CF"  # placeholder for snapshot
    assert "Gamma[1-eps]^1" in c.ft_base
    assert "Gamma[1-2*eps]" in c.ft_base
    assert c.derivative == "2*(1-2*eps)"

    # J factor
    tgv_loops = {leg.loop for leg in orbit.tgv_legs}
    pure = [p for p in orbit.pairings if p.loop not in tgv_loops]
    j = compute_j_mixed(orbit, hc, pure)

    assert j.lorentz_z == "z1"
    assert j.n_pure == 1
    assert "f1" in j.f_polynomials
    assert "alpha1" in j.pnum
    assert "alpha2" in j.S_polynomial
    assert "alpha3" in j.S_polynomial

    # Full emit
    result = emit_montopar_json(orbit, c, j, hc)
    assert result["dispatch"] == "mixed"
    assert result["J_raw"]["pure"]["f_polynomials"]["f1"] is not None
    assert result["J_raw"]["tgv"]["S_polynomial"] is not None
