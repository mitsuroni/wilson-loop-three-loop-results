"""Tests for montopar.parser."""
import json
from pathlib import Path

from montopar.parser import parse_orbit, parse_all_orbits, filter_pure_ngluon
from montopar.models import TGVLeg, FGVLeg

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture(name: str) -> dict:
    with open(FIXTURES / name, encoding="utf-8") as f:
        return json.load(f)


def test_parse_orbit_4e01():
    data = load_fixture("4e-01_orbit.json")
    orbit = parse_orbit(data)

    assert orbit.orbit_id == "4e-01"
    assert orbit.emission_type == "4e"
    assert orbit.n_gluons == 2

    # Edge distribution
    assert orbit.edge_distribution.counts == {3: 2, 4: 2}

    # Pairings
    assert len(orbit.pairings) == 2
    p1 = orbit.pairings[0]
    assert p1.edge_a == 3
    assert p1.edge_b == 4
    assert p1.s_a == "s31"
    assert p1.s_b == "s41"

    p2 = orbit.pairings[1]
    assert p2.loop == "l2"
    assert p2.s_a == "s32"
    assert p2.s_b == "s42"

    # Bounds
    assert len(orbit.bounds) == 4

    # Vertices preserved
    assert len(orbit.vertices) == 4
    assert orbit.vertices[0].line == 3

    # Propagators preserved
    assert len(orbit.propagators) == 2


def test_parse_orbit_edge_distribution():
    """Test edge distribution with 3 attachment points on one edge."""
    data = {
        "orbit_id": "test-01",
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
    assert len(orbit.pairings) == 3


def test_parse_all_orbits():
    data = load_fixture("4e-01_orbit.json")
    source = {"n_loops": 2, "n_orbits": 1, "orbits": [data]}
    orbits = parse_all_orbits(source)

    assert len(orbits) == 1
    assert orbits[0].orbit_id == "4e-01"


def test_filter_pure_ngluon():
    from montopar.models import Orbit, EdgeDistribution

    pure = Orbit(
        orbit_id="4e-01", emission_type="4e", n_gluons=2,
        edge_distribution=EdgeDistribution(counts={}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[],
    )
    tgv = Orbit(
        orbit_id="3e-3g-01", emission_type="3e-3g", n_gluons=1,
        edge_distribution=EdgeDistribution(counts={}, s_vars={}),
        pairings=[], bounds=[], vertices=[], propagators=[],
    )

    filtered = filter_pure_ngluon([pure, tgv])
    assert len(filtered) == 1
    assert filtered[0].emission_type == "4e"


class TestParseTGVOrbit:
    def test_parse_3e3g_01(self):
        """Parse a 3e-3g orbit with TGV internal vertex."""
        data = {
            "orbit_id": "3e-3g-01",
            "emission_type": "3e-3g",
            "sector": "dg",
            "nloops": 2,
            "vertices": [
                {"line": 3, "lorentz": "mu1", "color": "a1"},
                {"line": 4, "lorentz": "mu2", "color": "a2"},
                {"line": 4, "lorentz": "mu3", "color": "a3"},
            ],
            "propagators": [
                {"loop": "l1", "lorentz": ["mu1", "nu1"], "color": ["a1", "b1"]},
                {"loop": "l2", "lorentz": ["mu2", "nu2"], "color": ["a2", "b2"]},
                {"loop": "l3", "lorentz": ["mu3", "nu3"], "color": ["a3", "b3"]},
            ],
            "internal_vertices": [
                {"type": "tgv", "id": "v1", "loops": ["l1", "l2", "l3"]}
            ],
            "phases": [
                {"loop": "l1", "expression": "-(v1 + v2 + s31*v3)"},
                {"loop": "l2", "expression": "-(v1 + v2 + v3 + s41*v4)"},
                {"loop": "l3", "expression": "-(v1 + v2 + v3 + s42*v4)"},
            ],
            "bounds": [
                {"var": "s31", "low": "0", "high": "1"},
                {"var": "s41", "low": "0", "high": "s42"},
                {"var": "s42", "low": "0", "high": "1"},
            ],
        }
        orbit = parse_orbit(data)
        assert orbit.dispatch == "tgv"
        assert len(orbit.tgv_legs) == 3
        assert orbit.tgv_legs[0].edge == 3
        assert orbit.tgv_legs[0].s_var == "s31"
        assert orbit.tgv_legs[1].edge == 4
        assert orbit.tgv_legs[1].s_var == "s41"
        assert orbit.tgv_legs[2].edge == 4
        assert orbit.tgv_legs[2].s_var == "s42"
        assert orbit.pairings == []  # no bare gluons in pure TGV

    def test_parse_4e4g_01(self):
        """Parse a 4e-4g orbit with FGV internal vertex."""
        data = {
            "orbit_id": "4e-4g-01",
            "emission_type": "4e-4g",
            "sector": "dg",
            "nloops": 3,
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
        assert len(orbit.fgv_legs) == 4
        assert orbit.fgv_legs[0].edge == 3
        assert orbit.fgv_legs[3].edge == 4
        assert orbit.pairings == []  # no bare gluons in pure FGV
