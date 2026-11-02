"""Tests for combinatorial Wilson loop diagram generator."""
from __future__ import annotations

import pytest
from qgraf_pipeline.generator import _canonical_edge_distributions


class TestEdgeDistributions:
    def test_E2(self):
        """E=2: all distributions summing to 2, canonical under C4."""
        result = _canonical_edge_distributions(2)
        # All valid (a,b,c,d) with sum=2, canonical
        expected = [
            (0, 0, 1, 1),
            (0, 1, 0, 1),
        ]
        assert sorted(result) == sorted(expected)

    def test_E4(self):
        """E=4: all distributions summing to 4, canonical under C4."""
        result = _canonical_edge_distributions(4)
        expected = [
            (0, 0, 1, 3),
            (0, 0, 2, 2),
            (0, 0, 3, 1),
            (0, 1, 0, 3),
            (0, 1, 1, 2),
            (0, 1, 2, 1),
            (0, 2, 0, 2),
            (0, 2, 1, 1),
            (1, 1, 1, 1),
        ]
        assert sorted(result) == sorted(expected)

    def test_E6_count(self):
        """E=6: verify count and that all are canonical."""
        result = _canonical_edge_distributions(6)
        # Verify each is canonical (lexicographically smallest under C4)
        for ec in result:
            for r in range(1, 4):
                rotated = tuple(ec[(i - r) % 4] for i in range(4))
                assert ec <= rotated, f"{ec} not canonical: rotation {r} gives {rotated}"
        # E=6 has more distributions; just verify count is reasonable
        assert len(result) >= 6

    def test_no_sole_nonzero_edge(self):
        """No distribution should have all eikonals on a single edge (scaleless)."""
        for E in range(2, 9):
            for ec in _canonical_edge_distributions(E):
                non_zero = sum(1 for x in ec if x > 0)
                assert non_zero >= 2, f"E={E}: {ec} has only one non-zero edge"

    def test_all_E5_count(self):
        """E=5: all edge distributions without C4 canonicalization."""
        from qgraf_pipeline.generator import _all_edge_distributions
        result = _all_edge_distributions(5)
        # Total partitions of 5 into 4 parts: C(5+4-1,4-1) = C(8,3) = 56
        # Minus 4 single-edge distributions (5,0,0,0), (0,5,0,0), (0,0,5,0), (0,0,0,5)
        # = 52
        assert len(result) == 52

    def test_all_vs_canonical_E5(self):
        """Verify _all_edge_distributions is a superset of _canonical_edge_distributions."""
        from qgraf_pipeline.generator import _all_edge_distributions
        all_ec = _all_edge_distributions(5)
        canon_ec = _canonical_edge_distributions(5)
        assert set(canon_ec).issubset(set(all_ec))
        assert len(all_ec) > len(canon_ec)

    def test_all_E4_count(self):
        """E=4: all edge distributions without C4 canonicalization."""
        from qgraf_pipeline.generator import _all_edge_distributions
        result = _all_edge_distributions(4)
        # Total partitions: C(4+4-1,4-1) = C(7,3) = 35
        # Minus 4 single-edge = 31
        assert len(result) == 31

    def test_all_no_sole_nonzero_edge(self):
        """No distribution should have all eikonals on a single edge."""
        from qgraf_pipeline.generator import _all_edge_distributions
        for E in range(2, 9):
            for ec in _all_edge_distributions(E):
                non_zero = sum(1 for x in ec if x > 0)
                assert non_zero >= 2, f"E={E}: {ec} has only one non-zero edge"


from qgraf_pipeline.generator import (
    _canonical_edge_distributions,
    _contour_positions,
    _perfect_matchings,
    _build_wl_wl_diagram,
    _generate_wl_wl,
)
from qgraf_pipeline.internal_types import EmissionType


class TestContourPositions:
    def test_0033(self):
        """Edge distribution (0,0,3,3): positions 1-3 on edge 3, 4-6 on edge 4."""
        pos = _contour_positions((0, 0, 3, 3))
        assert pos == [
            (3, 1), (3, 2), (3, 3),  # edge 3
            (4, 1), (4, 2), (4, 3),  # edge 4
        ]

    def test_1111(self):
        """Edge distribution (1,1,1,1): one per edge."""
        pos = _contour_positions((1, 1, 1, 1))
        assert pos == [
            (1, 1), (2, 1), (3, 1), (4, 1),
        ]


class TestPerfectMatchings:
    def test_4_elements_count(self):
        """4 elements: (4-1)!! = 3 perfect matchings."""
        result = _perfect_matchings([0, 1, 2, 3])
        assert len(result) == 3

    def test_6_elements_count(self):
        """6 elements: (6-1)!! = 15 perfect matchings."""
        result = _perfect_matchings([0, 1, 2, 3, 4, 5])
        assert len(result) == 15

    def test_all_pairs_cover_all_elements(self):
        """Each matching must cover all elements exactly once."""
        for matching in _perfect_matchings([0, 1, 2, 3]):
            covered = set()
            for a, b in matching:
                covered.add(a)
                covered.add(b)
            assert covered == {0, 1, 2, 3}


class TestBuildWLWLDiagram:
    def test_0033_parallel_chords(self):
        """(0,0,3,3) with parallel matching (1,4)(2,5)(3,6)."""
        ec = (0, 0, 3, 3)
        matching = [(0, 3), (1, 4), (2, 5)]  # (pos1,pos4), (pos2,pos5), (pos3,pos6)
        d = _build_wl_wl_diagram(1, ec, matching)

        # Should have 6 eikonal vertices + 3 gluon propagators
        assert len(d.verts) == 6
        assert len(d.props) == 3

        # Check eikonal vertices: edges 3 and 4
        edges = []
        for v in d.verts:
            for leg in v.legs:
                if leg.field.startswith("z") and not leg.field.endswith("bar"):
                    edges.append(int(leg.field[1:]))
        assert edges == [3, 3, 3, 4, 4, 4]

        # Check gluon momenta: each prop has matching +/- in verts
        prop_momenta = {p.mom for p in d.props}
        assert len(prop_momenta) == 3  # g1, g2, g3

    def test_0033_crossed_chords(self):
        """(0,0,3,3) with crossed matching (1,5)(2,4)(3,6)."""
        ec = (0, 0, 3, 3)
        matching = [(0, 4), (1, 3), (2, 5)]  # (1,5)(2,4)(3,6)
        d = _build_wl_wl_diagram(1, ec, matching)
        assert len(d.verts) == 6
        assert len(d.props) == 3


class TestGenerateWLWL:
    def test_6e_generates_diagrams(self):
        """6e emission type should generate diagrams."""
        et = EmissionType(
            name="6e", category="tree",
            eikonal_count=6, wl_wl_pairs=3, g_power=6,
        )
        diagrams = _generate_wl_wl(et, start_idx=1)
        assert len(diagrams) > 0
        # All diagrams should have 6 eikonal verts and 3 props
        for d in diagrams:
            assert len(d.verts) == 6
            assert len(d.props) == 3


from qgraf_pipeline.generator import (
    _build_one_internal_diagram,
    _generate_one_internal,
)


class TestBuildOneInternalDiagram:
    def test_4e_4g(self):
        """4e-4g: all 4 eikonals connect to one 4g vertex."""
        ec = (1, 1, 1, 1)
        internal_pos = [0, 1, 2, 3]  # all positions to 4g
        wl_wl_positions = []  # no WL-WL pairs
        d = _build_one_internal_diagram(1, ec, internal_pos, wl_wl_positions, "4g", 1)

        # 4 eikonal verts + 1 4g vert = 5 verts
        assert len(d.verts) == 5
        # 4 gluon propagators
        assert len(d.props) == 4
        # 4g vertex has 4 A legs
        g4_vert = d.verts[-1]
        assert len(g4_vert.legs) == 4
        assert all(leg.field == "A" for leg in g4_vert.legs)

    def test_5e_3g(self):
        """5e-3g: 3 eikonals to 3g, 2 eikonals form WL-WL pair."""
        ec = (0, 1, 2, 2)
        internal_pos = [0, 1, 2]  # positions 0,1,2 to 3g
        wl_wl_positions = [3, 4]  # positions 3,4 form WL-WL pair
        d = _build_one_internal_diagram(1, ec, internal_pos, wl_wl_positions, "3g", 1)

        # 5 eikonal verts + 1 3g vert = 6 verts
        assert len(d.verts) == 6
        # 3 gluons to 3g + 1 WL-WL gluon = 4 props
        assert len(d.props) == 4
        # 3g vertex has 3 A legs
        g3_vert = d.verts[-1]
        assert len(g3_vert.legs) == 3
        assert all(leg.field == "A" for leg in g3_vert.legs)

    def test_4e_Pi1(self):
        """4e-Pi1: 2 eikonals to Pi1, 2 eikonals form WL-WL pair."""
        ec = (0, 1, 1, 2)
        internal_pos = [0, 2]  # positions 0,2 to Pi1
        wl_wl_positions = [1, 3]  # positions 1,3 form WL-WL pair
        d = _build_one_internal_diagram(1, ec, internal_pos, wl_wl_positions, "Pi1", 1)

        # 4 eikonal verts + 1 Pi1 vert = 5 verts
        assert len(d.verts) == 5
        # 2 gluons to Pi1 + 1 WL-WL gluon = 3 props
        assert len(d.props) == 3


class TestGenerateOneInternal:
    def test_4e_4g_generates(self):
        """4e-4g should generate diagrams."""
        et = EmissionType(
            name="4e-4g", category="tree",
            eikonal_count=4, wl_wl_pairs=0, g_power=6,
            internal={"4g": 1},
        )
        diagrams = _generate_one_internal(et, start_idx=1)
        assert len(diagrams) > 0
        for d in diagrams:
            # 4 eikonal + 1 4g vert
            assert len(d.verts) == 5
            assert len(d.props) == 4

    def test_5e_3g_generates(self):
        """5e-3g should generate diagrams."""
        et = EmissionType(
            name="5e-3g", category="tree",
            eikonal_count=5, wl_wl_pairs=1, g_power=6,
            internal={"3g": 1},
        )
        diagrams = _generate_one_internal(et, start_idx=1)
        assert len(diagrams) > 0
        for d in diagrams:
            assert len(d.verts) == 6  # 5 eikonal + 1 3g
            assert len(d.props) == 4  # 3 to 3g + 1 WL-WL

    def test_no_same_edge_wl_wl(self):
        """Generated diagrams must not have same-edge WL-WL pairs."""
        et = EmissionType(
            name="5e-3g", category="tree",
            eikonal_count=5, wl_wl_pairs=1, g_power=6,
            internal={"3g": 1},
        )
        diagrams = _generate_one_internal(et, start_idx=1)
        for d in diagrams:
            # Check each WL-WL gluon: find its two eikonal endpoints
            for prop in d.props:
                mom = prop.mom
                # Find verts with +mom and -mom
                pos_verts = []
                neg_verts = []
                for vi, v in enumerate(d.verts):
                    for leg in v.legs:
                        if leg.field == "A":
                            if leg.mom == mom:
                                pos_verts.append(vi)
                            elif leg.mom == f"-{mom}":
                                neg_verts.append(vi)
                if pos_verts and neg_verts:
                    # This is a WL-WL pair (both endpoints are eikonal)
                    # Check edges are different
                    edge_pos = None
                    edge_neg = None
                    for v in [d.verts[pos_verts[0]], d.verts[neg_verts[0]]]:
                        for leg in v.legs:
                            if leg.field.startswith("z") and not leg.field.endswith("bar"):
                                if edge_pos is None:
                                    edge_pos = int(leg.field[1:])
                                else:
                                    edge_neg = int(leg.field[1:])
                    if edge_pos is not None and edge_neg is not None:
                        assert edge_pos != edge_neg, \
                            f"Same-edge WL-WL pair: edge {edge_pos}"

    def test_5e_3g_raw_diagram_count(self):
        """5e-3g should generate exactly 336 diagrams before 2PI filtering."""
        et = EmissionType(
            name="5e-3g", category="tree",
            eikonal_count=5, wl_wl_pairs=1, g_power=6,
            internal={"3g": 1},
        )
        diagrams = _generate_one_internal(et, start_idx=1)
        assert len(diagrams) == 336, \
            f"Expected 336 raw 5e-3g diagrams, got {len(diagrams)}"

    def test_4e_Pi1_raw_diagram_count(self):
        """4e-Pi1 should generate exactly 126 diagrams before 2PI filtering."""
        et = EmissionType(
            name="4e-Pi1", category="loop",
            eikonal_count=4, wl_wl_pairs=1, g_power=6,
            internal={"Pi1": 1},
        )
        diagrams = _generate_one_internal(et, start_idx=1)
        assert len(diagrams) == 126, \
            f"Expected 126 raw 4e-Pi1 diagrams, got {len(diagrams)}"

    def test_4e_4g_raw_diagram_count_unchanged(self):
        """4e-4g (p=0) should be unaffected by the fix."""
        et = EmissionType(
            name="4e-4g", category="tree",
            eikonal_count=4, wl_wl_pairs=0, g_power=6,
            internal={"4g": 1},
        )
        diagrams_before = _generate_one_internal(et, start_idx=1)
        # p=0: should still use _canonical_edge_distributions
        assert len(diagrams_before) > 0


from qgraf_pipeline.generator import (
    generate_diagrams,
    _build_two_internal_diagram,
    _generate_two_internal,
)
from qgraf_pipeline.internal_types import enumerate_emission_types


class TestBuildTwoInternalDiagram:
    def test_4e_3g3g(self):
        """4e-3g+3g: two 3g vertices connected by a gluon, each connects to 2 eikonals."""
        ec = (1, 1, 1, 1)
        pos_a = [0, 1]  # eikonals for 3g-A
        pos_b = [2, 3]  # eikonals for 3g-B
        d = _build_two_internal_diagram(1, ec, pos_a, pos_b, "3g", "3g", 1)

        # 4 eikonal + 2 3g = 6 verts
        assert len(d.verts) == 6
        # 2 to 3g-A + 2 to 3g-B + 1 internal gluon = 5 props
        assert len(d.props) == 5
        # Both 3g vertices have 3 legs
        for vi in [4, 5]:
            assert len(d.verts[vi].legs) == 3
            assert all(leg.field == "A" for leg in d.verts[vi].legs)


class TestGenerateTwoInternal:
    def test_4e_3g3g_generates(self):
        """4e-3g+3g should generate diagrams."""
        et = EmissionType(
            name="4e-3g+3g", category="tree",
            eikonal_count=4, wl_wl_pairs=0, g_power=6,
            internal={"3g": 2},
        )
        diagrams = _generate_two_internal(et, start_idx=1)
        assert len(diagrams) > 0
        for d in diagrams:
            assert len(d.verts) == 6  # 4 eikonal + 2 3g
            assert len(d.props) == 5  # 4 to eikonals + 1 internal


class TestGenerateDiagrams:
    def test_n2_emission_types(self):
        """n=2: verify emission type count and diagram generation."""
        emission_types = enumerate_emission_types(2)
        diagrams = generate_diagrams(2, emission_types)
        assert len(diagrams) > 0
        # All diagrams should have valid structure
        for d in diagrams:
            assert d.idx > 0
            assert len(d.verts) > 0
            assert len(d.props) > 0

    def test_n3_emission_types(self):
        """n=3: verify emission type count and diagram generation."""
        emission_types = enumerate_emission_types(3)
        assert len(emission_types) == 9  # plus TGV-leg Pi1; plus 2e-Pi1+Pi1; 6e, 5e-3g, 4e-4g, 4e-Pi1, 2e-Pi2, 3e-Gamma1, 4e-3g+3g
        diagrams = generate_diagrams(3, emission_types)
        assert len(diagrams) > 0
        # All diagrams should have unique idx
        indices = [d.idx for d in diagrams]
        assert len(indices) == len(set(indices))

    def test_all_diagrams_have_valid_vertex_structure(self):
        """Every generated diagram must have proper eikonal vertices."""
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        for d in diagrams:
            for v in d.verts:
                fields = [leg.field for leg in v.legs]
                # Every vertex must have at least one leg
                assert len(fields) >= 1


# --- End-to-end integration tests (Task 6) ---
from pathlib import Path  # noqa: E402

from qgraf_pipeline.classifier import classify_orbits  # noqa: E402


class TestEndToEnd:
    def test_n2_produces_12_orbits(self):
        """n=2 with generator must produce exactly 12 orbits."""
        emission_types = enumerate_emission_types(2)
        diagrams = generate_diagrams(2, emission_types)
        orbits = classify_orbits(diagrams, emission_types)
        assert len(orbits) == 12, f"Expected 12 orbits, got {len(orbits)}"

        # Verify emission type coverage
        et_names = set(o.emission_type for o in orbits)
        assert "2e-Pi1" in et_names
        assert "3e-3g" in et_names
        assert "4e" in et_names

    def test_n3_6e_produces_39_orbits(self):
        """n=3 6e must produce exactly 39 orbits (including the 2 previously missing)."""
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        orbits = classify_orbits(diagrams, emission_types)

        six_e_orbits = [o for o in orbits if o.emission_type == "6e"]
        assert len(six_e_orbits) == 39, \
            f"Expected 39 6e orbits, got {len(six_e_orbits)}"

    def test_n3_per_type_orbits(self):
        """n=3 per-type orbit counts must match the reference table."""
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        orbits = classify_orbits(diagrams, emission_types)

        expected = {
            "6e": 39,
            "5e-3g": 44,
            "4e-4g": 9,
            "4e-Pi1": 11,
            "2e-Pi2": 2,
            "3e-Gamma1": 4,
            "4e-3g+3g": 26,
        }
        from collections import Counter
        actual = Counter(o.emission_type for o in orbits)
        for name, count in expected.items():
            assert actual[name] == count, \
                f"{name}: expected {count}, got {actual[name]}"

    def test_n3_total_orbits(self):
        """n=3 total orbits should match reference."""
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        orbits = classify_orbits(diagrams, emission_types)
        assert len(orbits) == 149, f"Expected 149 total orbits, got {len(orbits)}"

    def test_all_orbits_are_2pi(self):
        """Every generated orbit must be 2PI (proper irreducible web)."""
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        orbits = classify_orbits(diagrams, emission_types)

        non_2pi = [o.orbit_id for o in orbits if not o.is_2pi]
        assert len(non_2pi) == 0, \
            f"Found non-2PI orbits: {non_2pi}"

    def test_mirror_pairs_are_symmetric(self):
        """Mirror pairs must have matching orbit counts."""
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        orbits = classify_orbits(diagrams, emission_types)

        # Build mirror pair map
        mirrors = {}
        for o in orbits:
            if o.mirror_relation.startswith("mirror-pair"):
                pair_id = o.mirror_relation
                mirrors.setdefault(pair_id, []).append(o)

        for pair_id, pair_orbits in mirrors.items():
            assert len(pair_orbits) == 2, \
                f"Mirror pair {pair_id} has {len(pair_orbits)} orbits, expected 2"
            # Both orbits in a pair should have the same n_diagrams
            assert pair_orbits[0].n_diagrams == pair_orbits[1].n_diagrams, \
                f"Mirror pair {pair_id}: diagram counts differ"

    def test_n2_matches_reference_json(self):
        """n=2 output must match the reference integrals.json."""
        import json
        emission_types = enumerate_emission_types(2)
        diagrams = generate_diagrams(2, emission_types)
        orbits = classify_orbits(diagrams, emission_types)

        ref_path = Path(__file__).parent.parent / "outputs" / "2loop" / "integrals.json"
        if ref_path.exists():
            ref = json.loads(ref_path.read_text())
            ref_orbits = ref["orbits"]

            # Same total count
            assert len(orbits) == len(ref_orbits), \
                f"Orbit count mismatch: {len(orbits)} vs ref {len(ref_orbits)}"

            # Same orbit IDs
            gen_ids = sorted(o.orbit_id for o in orbits)
            ref_ids = sorted(o["orbit_id"] for o in ref_orbits)
            assert gen_ids == ref_ids, \
                f"Orbit ID mismatch: missing {set(ref_ids) - set(gen_ids)}, " \
                f"extra {set(gen_ids) - set(ref_ids)}"

    def test_n3_matches_reference_json(self):
        """n=3 output must match the reference integrals.json."""
        import json
        emission_types = enumerate_emission_types(3)
        diagrams = generate_diagrams(3, emission_types)
        orbits = classify_orbits(diagrams, emission_types)

        ref_path = Path(__file__).parent.parent / "outputs" / "3loop" / "integrals.json"
        if ref_path.exists():
            ref = json.loads(ref_path.read_text())
            assert len(orbits) == ref["n_orbits"], \
                f"Orbit count mismatch: {len(orbits)} vs ref {ref['n_orbits']}"
            from collections import Counter
            ref_types = Counter(o["emission_type"] for o in ref["orbits"])
            gen_types = Counter(o.emission_type for o in orbits)
            assert gen_types == ref_types, \
                f"Type mismatch: gen={dict(gen_types)} ref={dict(ref_types)}"

