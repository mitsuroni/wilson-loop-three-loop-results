"""Tests for type-level initial integral packages."""
import json

from montopar.package_builder import (
    build_initial_package,
    c4_symmetry,
    direct_opposite_edge_links,
    integral_structure,
)


def _raw(orbit_id="4e-01", emission_type="4e", vertices=None, propagators=None, internal=None):
    return {
        "orbit_id": orbit_id,
        "emission_type": emission_type,
        "nloops": 2,
        "vertices": vertices or [],
        "propagators": propagators or [],
        "internal_vertices": internal or [],
        "phases": [],
        "bounds": [],
    }


def _vertex(line, name):
    return {"line": line, "lorentz": name, "color": name.replace("mu", "a")}


def _prop(loop, left, right):
    return {"loop": loop, "lorentz": [left, right], "color": ["a", "b"]}


def test_c4_multiplicity_matches_known_4e_examples():
    paired = _raw(vertices=[_vertex(3, "mu1"), _vertex(3, "mu2"), _vertex(4, "mu3"), _vertex(4, "mu4")],
                  propagators=[_prop("l1", "mu1", "mu3"), _prop("l2", "mu2", "mu4")])
    alternating = _raw(vertices=[_vertex(2, "mu1"), _vertex(2, "mu2"), _vertex(4, "mu3"), _vertex(4, "mu4")],
                       propagators=[_prop("l1", "mu1", "mu3"), _prop("l2", "mu2", "mu4")])
    all_edges = _raw(vertices=[_vertex(i, f"mu{i}") for i in range(1, 5)],
                     propagators=[_prop("l1", "mu1", "mu3"), _prop("l2", "mu2", "mu4")])
    assert c4_symmetry(paired)["orbit_multiplicity"] == 4
    assert c4_symmetry(alternating)["orbit_multiplicity"] == 2
    assert c4_symmetry(all_edges)["orbit_multiplicity"] == 1


def test_opposite_links_distinguish_squared_channel_from_tensor_product():
    same = _raw(vertices=[_vertex(2, "mu1"), _vertex(2, "mu2"), _vertex(4, "mu3"), _vertex(4, "mu4")],
                propagators=[_prop("l1", "mu1", "mu3"), _prop("l2", "mu2", "mu4")])
    split = _raw(vertices=[_vertex(i, f"mu{i}") for i in range(1, 5)],
                 propagators=[_prop("l1", "mu1", "mu3"), _prop("l2", "mu2", "mu4")])
    assert direct_opposite_edge_links(same)["structure"] == "L24^2"
    assert direct_opposite_edge_links(split)["structure"] == "L13 × L24"


def test_mixed_tgv_core_and_direct_link_are_combined():
    raw = _raw(
        orbit_id="5e-3g-11", emission_type="5e-3g",
        vertices=[_vertex(3, "mu1"), _vertex(4, "mu2"), _vertex(4, "mu3"), _vertex(2, "mu4"), _vertex(4, "mu5")],
        propagators=[
            _prop("bare", "mu4", "mu5"),
            _prop("t1", "mu1", "nu1"), _prop("t2", "mu2", "nu2"), _prop("t3", "mu3", "nu3"),
        ],
        internal=[{"type": "tgv", "loops": ["t1", "t2", "t3"]}],
    )
    structure = integral_structure(raw)
    assert structure["core"]["topology_class"] == "A"
    assert structure["direct_opposite_edge_links"]["structure"] == "L24"
    assert structure["symbol"] == "A × L24"


def test_package_keeps_one_c_and_j_without_rotated_variants(tmp_path):
    raw = _raw(vertices=[_vertex(1, "mu1"), _vertex(3, "mu2")], propagators=[_prop("l1", "mu1", "mu2")])
    source = {"orbits": [raw]}
    target = tmp_path / "2loop" / "4e" / "4e-01"
    target.mkdir(parents=True)
    (target / "montopar.json").write_text(json.dumps({
        "dispatch": "pure", "geometry": {"n_lines": 4},
        "C": {"combined": "C"}, "J_raw": {"integrand": "J"},
        "C_rotated": {"combined": "bad"}, "J_raw_rotated": {"integrand": "bad"},
    }), encoding="utf-8")
    package, _ = build_initial_package(source, tmp_path, "4e")
    record = package["orbits"][0]
    assert record["C"] == {"combined": "C"}
    assert record["J"] == {"integrand": "J"}
    assert "C_rotated" not in record and "J_raw_rotated" not in record


def test_package_passes_through_ibp_representation(tmp_path):
    raw = _raw(vertices=[_vertex(1, "mu1"), _vertex(3, "mu2")], propagators=[_prop("l1", "mu1", "mu2")])
    source = {"orbits": [raw]}
    target = tmp_path / "2loop" / "4e" / "4e-01"
    target.mkdir(parents=True)
    ibp = {"status": "complete", "terms": [{"term_id": "b1"}]}
    (target / "montopar.json").write_text(json.dumps({
        "dispatch": "tgv", "geometry": {"n_lines": 4}, "schema_version": "1.1",
        "representative": {"canonical_topology": {"fixture": True}},
        "J_raw_s": {"alpha_measure": {"kind": "simplex"}},
        "normalization": {"euler_compensation_count": 0},
        "C": {"combined": "C"}, "J_raw": {"integrand": "J"},
        "J_ibp": ibp,
    }), encoding="utf-8")

    package, _ = build_initial_package(source, tmp_path, "4e")
    assert package["orbits"][0]["J_ibp"] == ibp


def test_package_rejects_legacy_tgv_without_canonical_semantics(tmp_path):
    import pytest
    raw = _raw(vertices=[_vertex(1, "mu1"), _vertex(3, "mu2")], propagators=[_prop("l1", "mu1", "mu2")])
    target = tmp_path / "2loop" / "4e" / "4e-01"
    target.mkdir(parents=True)
    (target / "montopar.json").write_text(json.dumps({
        "dispatch": "mixed", "C": {"combined": "old"}, "J_raw": {"integrand": "J"}
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="regenerate"):
        build_initial_package({"orbits": [raw]}, tmp_path, "4e")
