"""Tests for C4 rotation of full Orbit objects."""
import pytest
from montopar.models import (
    Vertex, Propagator, Bound, Pairing, EdgeDistribution,
    Orbit, TGVLeg, FGVLeg,
)
from montopar.rotation import rotate_edge, rotate_s_var, rotate_orbit_full


class TestRotateEdge:
    def test_rotate_by_0_is_identity(self):
        assert rotate_edge(1, 0) == 1
        assert rotate_edge(2, 0) == 2
        assert rotate_edge(3, 0) == 3
        assert rotate_edge(4, 0) == 4

    def test_rotate_by_1_shifts_forward(self):
        assert rotate_edge(1, 1) == 2
        assert rotate_edge(2, 1) == 3
        assert rotate_edge(3, 1) == 4
        assert rotate_edge(4, 1) == 1

    def test_rotate_by_2_gives_diagonal(self):
        assert rotate_edge(1, 2) == 3
        assert rotate_edge(2, 2) == 4
        assert rotate_edge(3, 2) == 1
        assert rotate_edge(4, 2) == 2

    def test_rotate_by_3_gives_reverse_shift(self):
        assert rotate_edge(1, 3) == 4
        assert rotate_edge(2, 3) == 1
        assert rotate_edge(3, 3) == 2
        assert rotate_edge(4, 3) == 3


class TestRotateSVar:
    def test_rotate_s_var_by_1(self):
        assert rotate_s_var("s31", 1) == "s41"
        assert rotate_s_var("s41", 1) == "s11"
        assert rotate_s_var("s21", 1) == "s31"

    def test_rotate_s_var_by_2(self):
        assert rotate_s_var("s31", 2) == "s11"
        assert rotate_s_var("s41", 2) == "s21"

    def test_non_s_var_unchanged(self):
        assert rotate_s_var("x1", 1) == "x1"
        assert rotate_s_var("alpha1", 2) == "alpha1"
        assert rotate_s_var("1", 3) == "1"


class TestRotateOrbitFull:
    """Tests for rotate_orbit_full — C4 rotation of entire Orbit."""

    def _make_tgv_orbit(self) -> Orbit:
        """Create a minimal 3e-3g-04 style orbit: TGV on (2,3,4)."""
        return Orbit(
            orbit_id="3e-3g-04",
            emission_type="3e-3g",
            n_gluons=3,
            edge_distribution=EdgeDistribution(
                counts={2: 1, 3: 1, 4: 1},
                s_vars={2: ["s21"], 3: ["s31"], 4: ["s41"]},
            ),
            pairings=[
                Pairing(edge_a=2, edge_b=3, loop="l1", s_a="s21", s_b="s31"),
                Pairing(edge_a=3, edge_b=4, loop="l2", s_a="s31", s_b="s41"),
                Pairing(edge_a=2, edge_b=4, loop="l3", s_a="s21", s_b="s41"),
            ],
            bounds=[
                Bound(var="s21", low="0", high="1"),
                Bound(var="s31", low="0", high="1"),
                Bound(var="s41", low="0", high="1"),
            ],
            vertices=[
                Vertex(line=2, lorentz="mu1", color="a1"),
                Vertex(line=3, lorentz="mu2", color="a2"),
                Vertex(line=4, lorentz="mu3", color="a3"),
            ],
            propagators=[
                Propagator(loop="l1", lorentz=("mu1", "nu1"), color=("a1", "b1")),
                Propagator(loop="l2", lorentz=("mu2", "nu2"), color=("a2", "b2")),
                Propagator(loop="l3", lorentz=("mu3", "nu3"), color=("a3", "b3")),
            ],
            dispatch="tgv",
            tgv_legs=[
                TGVLeg(edge=2, s_var="s21", loop="l1"),
                TGVLeg(edge=3, s_var="s31", loop="l2"),
                TGVLeg(edge=4, s_var="s41", loop="l3"),
            ],
        )

    def test_rotate_by_0_returns_equal_orbit(self):
        orbit = self._make_tgv_orbit()
        rotated = rotate_orbit_full(orbit, 0)
        assert rotated.orbit_id == orbit.orbit_id
        assert rotated.edge_distribution.counts == orbit.edge_distribution.counts
        for i, v in enumerate(rotated.vertices):
            assert v.line == orbit.vertices[i].line

    def test_rotate_by_1_shifts_all_edges(self):
        """C4+1: (2,3,4) → (3,4,1)."""
        orbit = self._make_tgv_orbit()
        rotated = rotate_orbit_full(orbit, 1)

        # Vertices
        assert rotated.vertices[0].line == 3  # was 2
        assert rotated.vertices[1].line == 4  # was 3
        assert rotated.vertices[2].line == 1  # was 4

        # Pairings
        assert rotated.pairings[0].edge_a == 3  # was 2
        assert rotated.pairings[0].edge_b == 4  # was 3
        assert rotated.pairings[0].s_a == "s31"  # was s21
        assert rotated.pairings[0].s_b == "s41"  # was s31

        # TGV legs
        assert rotated.tgv_legs[0].edge == 3  # was 2
        assert rotated.tgv_legs[0].s_var == "s31"  # was s21
        assert rotated.tgv_legs[1].edge == 4  # was 3
        assert rotated.tgv_legs[1].s_var == "s41"  # was s31
        assert rotated.tgv_legs[2].edge == 1  # was 4
        assert rotated.tgv_legs[2].s_var == "s11"  # was s41

        # Bounds
        s_vars_rotated = {b.var for b in rotated.bounds}
        assert s_vars_rotated == {"s31", "s41", "s11"}

        # Edge distribution
        assert rotated.edge_distribution.counts == {3: 1, 4: 1, 1: 1}
        assert rotated.edge_distribution.s_vars == {3: ["s31"], 4: ["s41"], 1: ["s11"]}

    def test_rotate_preserves_dispatch(self):
        orbit = self._make_tgv_orbit()
        rotated = rotate_orbit_full(orbit, 1)
        assert rotated.dispatch == "tgv"
        assert rotated.emission_type == "3e-3g"

    def test_rotate_bounds_with_dependencies(self):
        """Bounds with nested references: s31 in [0, s32], s32 in [0, 1]."""
        orbit = Orbit(
            orbit_id="test",
            emission_type="3e-3g",
            n_gluons=3,
            edge_distribution=EdgeDistribution(
                counts={3: 2, 4: 1},
                s_vars={3: ["s31", "s32"], 4: ["s41"]},
            ),
            pairings=[],
            bounds=[
                Bound(var="s31", low="0", high="s32"),
                Bound(var="s32", low="0", high="1"),
                Bound(var="s41", low="0", high="1"),
            ],
            vertices=[],
            propagators=[],
            dispatch="tgv",
            tgv_legs=[],
        )
        rotated = rotate_orbit_full(orbit, 1)

        # s31 → s41, s32 → s42, s41 → s11
        by_var = {b.var: b for b in rotated.bounds}
        assert set(by_var.keys()) == {"s41", "s42", "s11"}
        assert by_var["s41"].high == "s42"  # was s31 → s32
        assert by_var["s42"].high == "1"  # was s32 → 1
        assert by_var["s11"].high == "1"  # was s41 → 1

    def test_rotate_fgv_orbit(self):
        """4e-4g FGV orbit rotation."""
        orbit = Orbit(
            orbit_id="4e-4g-01",
            emission_type="4e-4g",
            n_gluons=4,
            edge_distribution=EdgeDistribution(
                counts={3: 1, 4: 3},
                s_vars={3: ["s31"], 4: ["s41", "s42", "s43"]},
            ),
            pairings=[],
            bounds=[
                Bound(var="s31", low="0", high="1"),
                Bound(var="s41", low="0", high="s42"),
                Bound(var="s42", low="0", high="s43"),
                Bound(var="s43", low="0", high="1"),
            ],
            vertices=[
                Vertex(line=3, lorentz="mu1", color="a1"),
                Vertex(line=4, lorentz="mu2", color="a2"),
                Vertex(line=4, lorentz="mu3", color="a3"),
                Vertex(line=4, lorentz="mu4", color="a4"),
            ],
            propagators=[],
            dispatch="fourv",
            fgv_legs=[
                FGVLeg(edge=3, s_var="s31", loop="l1"),
                FGVLeg(edge=4, s_var="s41", loop="l2"),
                FGVLeg(edge=4, s_var="s42", loop="l3"),
                FGVLeg(edge=4, s_var="s43", loop="l4"),
            ],
        )
        rotated = rotate_orbit_full(orbit, 1)

        # Vertices: (3,4,4,4) → (4,1,1,1)
        assert rotated.vertices[0].line == 4
        assert rotated.vertices[1].line == 1
        assert rotated.vertices[2].line == 1
        assert rotated.vertices[3].line == 1

        # FGV legs
        assert rotated.fgv_legs[0].edge == 4
        assert rotated.fgv_legs[0].s_var == "s41"
        assert rotated.fgv_legs[1].edge == 1
        assert rotated.fgv_legs[1].s_var == "s11"

        # Edge distribution
        assert rotated.edge_distribution.counts == {4: 1, 1: 3}

        # Bounds
        by_var = {b.var: b for b in rotated.bounds}
        assert set(by_var.keys()) == {"s41", "s11", "s12", "s13"}
        assert by_var["s11"].high == "s12"  # was s41 → s42

    def test_rotate_mixed_orbit(self):
        """5e-3g mixed dispatch: TGV + pure gluon."""
        orbit = Orbit(
            orbit_id="5e-3g-14",
            emission_type="5e-3g",
            n_gluons=2,  # pure gluons only
            edge_distribution=EdgeDistribution(
                counts={2: 1, 3: 2, 4: 2},
                s_vars={2: ["s21"], 3: ["s31", "s32"], 4: ["s41", "s42"]},
            ),
            pairings=[
                Pairing(edge_a=3, edge_b=4, loop="l1", s_a="s32", s_b="s42"),
            ],
            bounds=[
                Bound(var="s21", low="0", high="1"),
                Bound(var="s31", low="0", high="1"),
                Bound(var="s32", low="0", high="1"),
                Bound(var="s41", low="0", high="1"),
                Bound(var="s42", low="0", high="1"),
            ],
            vertices=[
                Vertex(line=2, lorentz="mu1", color="a1"),
                Vertex(line=3, lorentz="mu2", color="a2"),
                Vertex(line=4, lorentz="mu3", color="a3"),
                Vertex(line=4, lorentz="mu4", color="a4"),
                Vertex(line=4, lorentz="mu5", color="a5"),
            ],
            propagators=[],
            dispatch="mixed",
            tgv_legs=[
                TGVLeg(edge=3, s_var="s31", loop="l2"),
                TGVLeg(edge=4, s_var="s41", loop="l3"),
                TGVLeg(edge=4, s_var="s42", loop="l4"),
            ],
        )
        rotated = rotate_orbit_full(orbit, 1)

        # Vertices: (2,3,4,4,4) → (3,4,1,1,1)
        assert rotated.vertices[0].line == 3
        assert rotated.vertices[1].line == 4
        assert rotated.vertices[2].line == 1
        assert rotated.vertices[3].line == 1
        assert rotated.vertices[4].line == 1

        # TGV legs: (3,4,4) → (4,1,1)
        assert rotated.tgv_legs[0].edge == 4
        assert rotated.tgv_legs[0].s_var == "s41"  # was s31
        assert rotated.tgv_legs[1].edge == 1
        assert rotated.tgv_legs[1].s_var == "s11"  # was s41
        assert rotated.tgv_legs[2].edge == 1
        assert rotated.tgv_legs[2].s_var == "s12"  # was s42

        # Pure pairing: (3,4) → (4,1)
        assert rotated.pairings[0].edge_a == 4
        assert rotated.pairings[0].edge_b == 1
        assert rotated.pairings[0].s_a == "s42"  # was s32
        assert rotated.pairings[0].s_b == "s12"  # was s42

        # Dispatch preserved
        assert rotated.dispatch == "mixed"
        assert rotated.n_gluons == 2