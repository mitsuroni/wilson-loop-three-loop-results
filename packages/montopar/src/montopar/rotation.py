# -*- coding: utf-8 -*-
"""
C4 rotation for TGV edges and full Orbit objects.

Canonical bases (by cyclic structure, not specific edge numbers):
  (3,4,4): b=c, non-degenerate edge adjacent to degenerate edges
  (3,3,4): a=b, third edge adjacent to degenerate edges
  (1,3,3): b=c, non-degenerate edge opposite to degenerate edges
  (1,3,4): non-degenerate, cyclic pattern (adjacent, adjacent, opposite)
"""
import re
from typing import Tuple, Optional

from montopar.models import (
    Orbit, Vertex, Propagator, Bound, Pairing, EdgeDistribution,
    TGVLeg, FGVLeg, DoubleTGVLeg,
)

_S_VAR_RE = re.compile(r"^s(\d)(\d+)$")


def rotate_edge(edge: int, k: int) -> int:
    """Apply C4 rotation: edge -> (edge + k - 1) % 4 + 1."""
    return ((edge + k - 1) % 4) + 1


def rotate_s_var(s_var: str, k: int) -> str:
    """Rotate the edge part of a Schwinger parameter name.

    ``s31`` with k=1 -> ``s41``.  Non-s variables are returned unchanged.
    """
    m = _S_VAR_RE.match(s_var)
    if not m:
        return s_var
    edge = int(m.group(1))
    idx = m.group(2)
    new_edge = rotate_edge(edge, k)
    return f"s{new_edge}{idx}"


def cyclic_dist(a: int, b: int) -> int:
    """Distance between two edges in cyclic order: (b - a) mod 4, in {1,2,3}."""
    d = (b - a) % 4
    return 4 if d == 0 else d


def classify_structure(a: int, b: int, c: int) -> str:
    """Classify the cyclic structure of three edges (a,b,c) in order."""
    if a == b and b == c:
        return "triple"
    if a == b:
        d = cyclic_dist(a, c)
        return "degenerate_ab_adjacent" if d in (1, 3) else "degenerate_ab_opposite"
    if b == c:
        d = cyclic_dist(a, b)
        return "degenerate_bc_adjacent" if d in (1, 3) else "degenerate_bc_opposite"
    if a == c:
        d = cyclic_dist(a, b)
        return "degenerate_ac_adjacent" if d in (1, 3) else "degenerate_ac_opposite"
    d_ab = cyclic_dist(a, b)
    d_bc = cyclic_dist(b, c)
    if d_ab == 1 and d_bc == 1:
        return "nondegenerate_consecutive"
    return "nondegenerate_mixed"


def find_canonical_rotation(a: int, b: int, c: int) -> Optional[Tuple[int, str]]:
    """Find the C4 rotation that maps (a,b,c) to the canonical form.

    Returns (rotation_k, structure_type) or None.
    """
    struct = classify_structure(a, b, c)
    if struct == "triple":
        return None
    if struct in ("degenerate_bc_adjacent", "degenerate_bc_opposite"):
        return 0, struct
    if struct in ("degenerate_ab_adjacent", "degenerate_ab_opposite"):
        return 0, struct
    if struct in ("degenerate_ac_adjacent", "degenerate_ac_opposite"):
        return 0, struct.replace("ac", "bc")
    if struct == "nondegenerate_consecutive":
        return 0, struct
    if struct == "nondegenerate_mixed":
        return 0, struct
    return None


def rotate_orbit_edges(tgv_edges: Tuple[int, int, int], bare_edges: list) -> dict:
    """Rotate all edges to canonical form.

    Returns dict with rotation_k, canonical_tgv, canonical_bare, structure_type.
    """
    a, b, c = tgv_edges
    result = find_canonical_rotation(a, b, c)
    if result is None:
        raise ValueError(f"Cannot classify TGV edges {tgv_edges}")

    k, struct = result
    rotated_tgv = tuple(rotate_edge(e, k) for e in tgv_edges)
    rotated_bare = [(rotate_edge(be[0], k), rotate_edge(be[1], k)) for be in bare_edges]

    return {
        "rotation_k": k,
        "canonical_tgv": rotated_tgv,
        "canonical_bare": rotated_bare,
        "structure": struct,
    }


def rotate_orbit_full(orbit: Orbit, k: int) -> Orbit:
    """Apply C4 rotation by k to all edge-indexed fields of an Orbit.

    Rotates vertices, bounds, edge_distribution, pairings, tgv_legs,
    and fgv_legs.  Returns a new Orbit; the original is unchanged.
    """
    # Vertices
    rotated_vertices = [
        Vertex(line=rotate_edge(v.line, k), lorentz=v.lorentz, color=v.color)
        for v in orbit.vertices
    ]

    # Propagators (unchanged - lorentz/color labels don't depend on edge numbers)
    rotated_propagators = [
        Propagator(loop=p.loop, lorentz=p.lorentz, color=p.color)
        for p in orbit.propagators
    ]

    # Bounds: rotate var names and any references in low/high
    rotated_bounds = []
    for b in orbit.bounds:
        new_var = rotate_s_var(b.var, k)
        new_low = rotate_s_var(b.low, k) if b.low not in ("0", "1") else b.low
        new_high = rotate_s_var(b.high, k) if b.high not in ("0", "1") else b.high
        rotated_bounds.append(Bound(var=new_var, low=new_low, high=new_high))

    # Edge distribution
    rotated_counts: dict[int, int] = {}
    rotated_s_vars: dict[int, list[str]] = {}
    for edge, cnt in orbit.edge_distribution.counts.items():
        new_edge = rotate_edge(edge, k)
        rotated_counts[new_edge] = cnt
    for edge, svars in orbit.edge_distribution.s_vars.items():
        new_edge = rotate_edge(edge, k)
        rotated_s_vars[new_edge] = [rotate_s_var(s, k) for s in svars]

    rotated_edge_dist = EdgeDistribution(
        counts=rotated_counts,
        s_vars=rotated_s_vars,
    )

    # Pairings
    rotated_pairings = []
    for p in orbit.pairings:
        rotated_pairings.append(Pairing(
            edge_a=rotate_edge(p.edge_a, k),
            edge_b=rotate_edge(p.edge_b, k),
            loop=p.loop,
            s_a=rotate_s_var(p.s_a, k),
            s_b=rotate_s_var(p.s_b, k),
        ))

    # TGV legs
    rotated_tgv_legs = [
        TGVLeg(edge=rotate_edge(leg.edge, k), s_var=rotate_s_var(leg.s_var, k), loop=leg.loop)
        for leg in orbit.tgv_legs
    ]

    # FGV legs
    rotated_fgv_legs = [
        FGVLeg(edge=rotate_edge(leg.edge, k), s_var=rotate_s_var(leg.s_var, k), loop=leg.loop)
        for leg in orbit.fgv_legs
    ]

    # Double TGV legs
    rotated_double_tgv_legs = [
        DoubleTGVLeg(edge=rotate_edge(leg.edge, k), s_var=rotate_s_var(leg.s_var, k),
                     loop=leg.loop, tgv_id=leg.tgv_id)
        for leg in orbit.double_tgv_legs
    ]

    return Orbit(
        orbit_id=orbit.orbit_id,
        emission_type=orbit.emission_type,
        n_gluons=orbit.n_gluons,
        edge_distribution=rotated_edge_dist,
        pairings=rotated_pairings,
        bounds=rotated_bounds,
        vertices=rotated_vertices,
        propagators=rotated_propagators,
        source=orbit.source,
        dispatch=orbit.dispatch,
        tgv_legs=rotated_tgv_legs,
        fgv_legs=rotated_fgv_legs,
        double_tgv_legs=rotated_double_tgv_legs,
    )


# Test
if __name__ == "__main__":
    tests = [
        ((3, 4, 4), [(3, 4)]),
        ((2, 3, 4), [(3, 4)]),
        ((2, 4, 4), [(3, 4)]),
        ((1, 3, 3), [(3, 4)]),
        ((3, 3, 4), [(3, 4)]),
    ]
    for tgv, bare in tests:
        r = rotate_orbit_edges(tgv, bare)
        print(f"TGV {tgv}: struct={r['structure']}, rotated={r['canonical_tgv']}, bare={r['canonical_bare']}")
