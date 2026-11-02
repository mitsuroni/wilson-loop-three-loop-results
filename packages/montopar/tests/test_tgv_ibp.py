"""Tests for pure 3e-3g contour-IBP generation."""
from montopar.c_factor import compute_c_tgv
from montopar.models import Bound, EdgeDistribution, Orbit, TGVLeg
from montopar.rotation import rotate_orbit_full
from montopar.tgv_ibp import (
    classify_tgv_ibp_case,
    compute_tgv_ibp,
    find_ordered_pair,
)


def _orbit(orbit_id, edges, variables, bounds):
    return Orbit(
        orbit_id=orbit_id,
        emission_type="3e-3g",
        n_gluons=3,
        edge_distribution=EdgeDistribution(
            counts={edge: edges.count(edge) for edge in set(edges)},
            s_vars={},
        ),
        pairings=[],
        bounds=[Bound(*bound) for bound in bounds],
        vertices=[],
        propagators=[],
        dispatch="tgv",
        tgv_legs=[
            TGVLeg(edge=edge, s_var=var, loop=f"l{i}")
            for i, (edge, var) in enumerate(zip(edges, variables), start=1)
        ],
    )


def _orbits():
    return {
        "01": _orbit(
            "3e-3g-01", [3, 4, 4], ["s31", "s41", "s42"],
            [("s31", "0", "1"), ("s41", "0", "s42"), ("s42", "0", "1")],
        ),
        "02": _orbit(
            "3e-3g-02", [3, 3, 4], ["s31", "s32", "s41"],
            [("s31", "0", "s32"), ("s32", "0", "1"), ("s41", "0", "1")],
        ),
        "03": _orbit(
            "3e-3g-03", [2, 4, 4], ["s21", "s41", "s42"],
            [("s21", "0", "1"), ("s41", "0", "s42"), ("s42", "0", "1")],
        ),
        "04": _orbit(
            "3e-3g-04", [2, 3, 4], ["s21", "s31", "s41"],
            [("s21", "0", "1"), ("s31", "0", "1"), ("s41", "0", "1")],
        ),
    }


def test_case_classification():
    orbits = _orbits()
    assert classify_tgv_ibp_case(orbits["01"].tgv_legs) == "degenerate_bc"
    assert classify_tgv_ibp_case(orbits["02"].tgv_legs) == "degenerate_ab"
    assert classify_tgv_ibp_case(orbits["04"].tgv_legs) == "nondegenerate"


def test_find_ordered_pair():
    orbit = _orbits()["02"]
    assert find_ordered_pair(orbit.bounds, {"s31", "s32"}) == ("s31", "s32")


def test_all_degenerate_orbits_have_two_cusps_and_one_contact():
    expected = {
        "01": ["z1", "z1", "-2*z1"],
        "02": ["z1", "z1", "-2*z1"],
        "03": ["-1 - z1", "-1 - z1", "2 + 2*z1"],
    }
    for key, coefficients in expected.items():
        result = compute_tgv_ibp(_orbits()[key], compute_c_tgv())
        assert result.status == "complete"
        assert result.certificates["pnum_identity"] is True
        assert len(result.terms) == 3
        assert [term.kind for term in result.terms] == [
            "cusp_boundary", "cusp_boundary", "contact_boundary",
        ]
        assert [term.coefficient for term in result.terms] == coefficients
        assert all(len(term.variables) == 2 for term in result.terms)
        assert len(result.terms[-1].identified_parameters) == 1


def test_degenerate_rotated_coefficients():
    expected = {
        "01": ["1", "1", "-2"],
        "02": ["1", "1", "-2"],
        "03": ["-1 - z1", "-1 - z1", "2 + 2*z1"],
    }
    for key, coefficients in expected.items():
        rotated = rotate_orbit_full(_orbits()[key], 1)
        result = compute_tgv_ibp(rotated, compute_c_tgv())
        assert result.status == "complete"
        assert [term.coefficient for term in result.terms] == coefficients


def test_nondegenerate_a_has_six_faces_and_recoil():
    result = compute_tgv_ibp(_orbits()["04"], compute_c_tgv())
    assert result.status == "complete"
    assert result.derivative_vector == ["z1", "-1 - z1", "1"]
    assert len(result.terms) == 7
    assert [term.coefficient for term in result.terms[:-1]] == [
        "-z1", "z1", "1 + z1", "-1 - z1", "-1", "1",
    ]
    recoil = result.terms[-1]
    assert recoil.kind == "recoil"
    assert recoil.c_ref == "C_raw"
    assert len(recoil.variables) == 3
    assert result.certificates["pnum_identity"] is True
    for fragment in (
        "-1/2*alpha1*alpha2",
        "1/2*alpha1*alpha3*z1",
        "1/2*x3*alpha1*alpha3*z1^2",
    ):
        assert fragment in recoil.integrand


def test_nondegenerate_b_is_recomputed_after_rotation():
    rotated = rotate_orbit_full(_orbits()["04"], 1)
    result = compute_tgv_ibp(rotated, compute_c_tgv())
    assert result.status == "complete"
    assert result.derivative_vector == ["1", "-1 - z1", "z1"]
    assert len(result.terms) == 7
    assert result.certificates == {
        "operator_identity": True,
        "boundary_orientation": True,
        "pnum_identity": True,
    }


def test_bad_degenerate_bounds_are_explicitly_unsupported():
    orbit = _orbits()["02"]
    orbit.bounds = [
        Bound("s31", "0", "1"), Bound("s32", "0", "1"), Bound("s41", "0", "1"),
    ]
    result = compute_tgv_ibp(orbit, compute_c_tgv())
    assert result.status == "unsupported"
    assert result.reason == "unrecognized_ordered_pair"
    assert result.terms == []
