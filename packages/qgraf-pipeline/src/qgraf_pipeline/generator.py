"""Combinatorial Wilson loop diagram generator.

Replaces QGRAF by enumerating all valid eikonal-vertex placements
and gluon-connection patterns directly from the emission type
constraints.
"""
from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Tuple

from .parser import Diagram, Prop, Vert, Leg
from .internal_types import EmissionType, _name_to_mn


def _canonical_edge_distributions(E: int) -> List[Tuple[int, int, int, int]]:
    """Enumerate all (a,b,c,d) summing to E, canonical under C4 rotation.

    Returns each distribution in its lexicographically smallest form.
    """
    seen: set = set()
    result: List[Tuple[int, int, int, int]] = []

    for a in range(E + 1):
        for b in range(E + 1 - a):
            for c in range(E + 1 - a - b):
                d = E - a - b - c
                ec = (a, b, c, d)

                # Skip: all eikonals on a single edge (scaleless type B)
                non_zero = sum(1 for x in ec if x > 0)
                if non_zero < 2:
                    continue

                # Canonicalize: find lexicographically smallest rotation
                best = ec
                for r in range(1, 4):
                    rotated = tuple(ec[(i - r) % 4] for i in range(4))
                    if rotated < best:
                        best = rotated

                if best not in seen:
                    seen.add(best)
                    result.append(best)

    return sorted(result)


def _all_edge_distributions(E: int) -> List[Tuple[int, int, int, int]]:
    """Enumerate all (a,b,c,d) summing to E, without C4 canonicalization.

    Only filters out single-edge distributions (scaleless type B).
    Used for mixed dispatch types (p > 0) where C4 rotation of edge
    distributions is not valid because TGV/Pi legs and WL-WL gluon
    legs are not interchangeable.
    """
    result: List[Tuple[int, int, int, int]] = []

    for a in range(E + 1):
        for b in range(E + 1 - a):
            for c in range(E + 1 - a - b):
                d = E - a - b - c
                ec = (a, b, c, d)

                # Skip: all eikonals on a single edge (scaleless type B)
                non_zero = sum(1 for x in ec if x > 0)
                if non_zero < 2:
                    continue

                result.append(ec)

    return sorted(result)


def _contour_positions(
    edge_counts: Tuple[int, int, int, int]
) -> List[Tuple[int, int]]:
    """Return list of (edge, slot) for positions 1..E along the contour.

    Contour order: edge 1 slots 1..a, edge 2 slots 1..b,
    edge 3 slots 1..c, edge 4 slots 1..d.
    """
    positions: List[Tuple[int, int]] = []
    for edge in (1, 2, 3, 4):
        n = edge_counts[edge - 1]
        for slot in range(1, n + 1):
            positions.append((edge, slot))
    return positions


def _perfect_matchings(elements: list) -> List[List[Tuple[int, int]]]:
    """Generate all perfect matchings of a list of elements.

    A perfect matching pairs all elements into len(elements)/2 disjoint pairs.
    Yields lists of (a, b) tuples where a < b.
    """
    if not elements:
        return [[]]

    first = elements[0]
    result = []
    for i in range(1, len(elements)):
        second = elements[i]
        rest = elements[1:i] + elements[i + 1:]
        for submatch in _perfect_matchings(rest):
            pair = (first, second) if first < second else (second, first)
            result.append([pair] + submatch)
    return result


def _build_wl_wl_diagram(
    idx: int,
    edge_counts: Tuple[int, int, int, int],
    matching: List[Tuple[int, int]],
    et_name: str = "",
) -> Diagram:
    """Build a synthetic Diagram for a pure WL-WL configuration.

    Args:
        idx: Diagram index (unique across all generated diagrams)
        edge_counts: (a, b, c, d) eikonal counts per edge
        matching: List of (pos_i, pos_j) pairs, 0-indexed contour positions
    """
    positions = _contour_positions(edge_counts)
    E = len(positions)

    # Build eikonal vertices in contour order
    verts = []
    pos_to_mom: Dict[int, str] = {}  # contour position -> gluon momentum string

    for gluon_idx, (pa, pb) in enumerate(matching, 1):
        mom = f"g{gluon_idx}"
        pos_to_mom[pa] = mom
        pos_to_mom[pb] = f"-{mom}"

    for pos_idx in range(E):
        edge, slot = positions[pos_idx]
        mom = pos_to_mom[pos_idx]
        verts.append(Vert([
            Leg(f"z{edge}", "0"),
            Leg(f"z{edge}bar", "0"),
            Leg("A", mom),
        ]))

    # Build gluon propagators
    props = []
    for gluon_idx in range(1, len(matching) + 1):
        props.append(Prop("A", "A", f"g{gluon_idx}"))

    return Diagram(
        idx=idx,
        sign=1,
        sym_factor="1",
        props=props,
        verts=verts,
        emission_type=et_name,
    )


def _generate_wl_wl(
    et: EmissionType,
    start_idx: int,
) -> List[Diagram]:
    """Generate all pure WL-WL diagrams for an emission type.

    Enumerates all edge distributions, perfect matchings, and filters
    out same-edge pairs (scaleless). The full 2PI/WEB check is done
    later by classify_orbits().
    """
    E = et.eikonal_count
    p = et.wl_wl_pairs  # E = 2p for pure WL-WL

    diagrams = []
    idx = start_idx

    for ec in _canonical_edge_distributions(E):
        positions = _contour_positions(ec)
        all_matchings = _perfect_matchings(list(range(E)))

        for matching in all_matchings:
            # Filter: no same-edge pair (scaleless type A)
            valid = True
            for pa, pb in matching:
                edge_a = positions[pa][0]
                edge_b = positions[pb][0]
                if edge_a == edge_b:
                    valid = False
                    break
            if not valid:
                continue

            diagrams.append(_build_wl_wl_diagram(idx, ec, matching, et.name))
            idx += 1

    return diagrams


def _build_one_internal_diagram(
    idx: int,
    edge_counts: Tuple[int, int, int, int],
    internal_positions: List[int],
    wl_wl_positions: List[int],
    internal_name: str,
    gluon_start: int,
    et_name: str = "",
) -> Diagram:
    """Build a synthetic Diagram with one internal vertex.

    Args:
        idx: Diagram index
        edge_counts: (a, b, c, d) eikonal counts per edge
        internal_positions: Contour positions (0-indexed) connecting to internal vertex
        wl_wl_positions: Contour positions (0-indexed) forming WL-WL pairs
        internal_name: Name of internal vertex type (e.g. "3g", "4g", "Pi1")
        gluon_start: Starting index for gluon momentum labels
    """
    positions = _contour_positions(edge_counts)
    E = len(positions)
    N = len(internal_positions)
    p = len(wl_wl_positions) // 2

    verts = []
    pos_to_mom: Dict[int, str] = {}

    # Assign momenta: gluons gluon_start .. gluon_start+N-1 go to internal vertex
    gluon_idx = gluon_start
    for pi in internal_positions:
        mom = f"g{gluon_idx}"
        pos_to_mom[pi] = mom
        gluon_idx += 1

    # Assign momenta for WL-WL pairs
    for pair_i in range(p):
        pa = wl_wl_positions[2 * pair_i]
        pb = wl_wl_positions[2 * pair_i + 1]
        mom = f"g{gluon_idx}"
        pos_to_mom[pa] = mom
        pos_to_mom[pb] = f"-{mom}"
        gluon_idx += 1

    # Build eikonal vertices in contour order
    for pos_idx in range(E):
        edge, slot = positions[pos_idx]
        mom = pos_to_mom[pos_idx]
        verts.append(Vert([
            Leg(f"z{edge}", "0"),
            Leg(f"z{edge}bar", "0"),
            Leg("A", mom),
        ]))

    # Build internal vertex
    internal_legs = []
    for gi in range(gluon_start, gluon_start + N):
        internal_legs.append(Leg("A", f"-g{gi}"))
    verts.append(Vert(internal_legs))

    # Build propagators
    props = []
    for gi in range(gluon_start, gluon_idx):
        props.append(Prop("A", "A", f"g{gi}"))

    return Diagram(
        idx=idx,
        sign=1,
        sym_factor="1",
        props=props,
        verts=verts,
        emission_type=et_name,
    )


def _generate_one_internal(
    et: EmissionType,
    start_idx: int,
) -> List[Diagram]:
    """Generate all diagrams for types with one internal vertex.

    Covers: 5e-3g, 4e-4g, 4e-Pi1, 2e-Pi2, 3e-Gamma1.
    """
    E = et.eikonal_count
    p = et.wl_wl_pairs
    internal_name = list(et.internal.keys())[0]
    _, N = _name_to_mn(internal_name)

    diagrams = []
    idx = start_idx

    if p > 0:
        edge_distributions = _all_edge_distributions(E)
    else:
        edge_distributions = _canonical_edge_distributions(E)

    for ec in edge_distributions:
        positions = _contour_positions(ec)
        all_pos = list(range(E))

        # Choose which positions connect to the internal vertex
        for internal_pos in combinations(all_pos, N):
            internal_set = set(internal_pos)
            wl_wl_positions = [pos for pos in all_pos if pos not in internal_set]

            if p > 0:
                # Generate WL-WL pairings for the remaining positions
                for matching in _perfect_matchings(wl_wl_positions):
                    # Flatten matching into ordered list
                    wl_wl_flat = []
                    valid = True
                    for pa, pb in matching:
                        edge_a = positions[pa][0]
                        edge_b = positions[pb][0]
                        if edge_a == edge_b:
                            valid = False
                            break
                        wl_wl_flat.extend([pa, pb])
                    if not valid:
                        continue

                    diagrams.append(_build_one_internal_diagram(
                        idx, ec, list(internal_pos), wl_wl_flat,
                        internal_name, gluon_start=1,
                        et_name=et.name,
                    ))
                    idx += 1
            else:
                # No WL-WL pairs: all eikonals connect to internal vertex
                diagrams.append(_build_one_internal_diagram(
                    idx, ec, list(internal_pos), [],
                    internal_name, gluon_start=1,
                    et_name=et.name,
                ))
                idx += 1

    return diagrams


def _build_two_internal_diagram(
    idx: int,
    edge_counts: Tuple[int, int, int, int],
    pos_a: List[int],
    pos_b: List[int],
    name_a: str,
    name_b: str,
    gluon_start: int,
    et_name: str = "",
    wl_wl_positions: List[int] | None = None,
) -> Diagram:
    """Build a synthetic Diagram with two internal vertices connected by a gluon.

    Args:
        idx: Diagram index
        edge_counts: (a, b, c, d) eikonal counts per edge
        pos_a: Contour positions (0-indexed) connecting to vertex A
        pos_b: Contour positions (0-indexed) connecting to vertex B
        name_a, name_b: Internal vertex type names
        gluon_start: Starting index for gluon momentum labels
        wl_wl_positions: Optional list of positions forming WL-WL pairs,
            in order [p1a, p1b, p2a, p2b, ...]. Each pair gets +gN/-gN.
    """
    positions = _contour_positions(edge_counts)
    E = len(positions)
    N_a = len(pos_a)
    N_b = len(pos_b)

    verts = []
    pos_to_mom: Dict[int, str] = {}

    # Gluon g{gluon_start} is the internal gluon connecting the two vertices
    internal_gluon = gluon_start
    gluon_idx = gluon_start + 1

    # Assign momenta for eikonals -> vertex A
    for pi in pos_a:
        mom = f"g{gluon_idx}"
        pos_to_mom[pi] = mom
        gluon_idx += 1

    # Assign momenta for eikonals -> vertex B
    for pi in pos_b:
        mom = f"g{gluon_idx}"
        pos_to_mom[pi] = mom
        gluon_idx += 1

    # Assign momenta for WL-WL pairs (if any)
    if wl_wl_positions:
        for pair_i in range(len(wl_wl_positions) // 2):
            pa = wl_wl_positions[2 * pair_i]
            pb = wl_wl_positions[2 * pair_i + 1]
            mom = f"g{gluon_idx}"
            pos_to_mom[pa] = mom
            pos_to_mom[pb] = f"-{mom}"
            gluon_idx += 1

    # Build eikonal vertices in contour order
    for pos_idx in range(E):
        edge, slot = positions[pos_idx]
        mom = pos_to_mom.get(pos_idx, "0")
        verts.append(Vert([
            Leg(f"z{edge}", "0"),
            Leg(f"z{edge}bar", "0"),
            Leg("A", mom),
        ]))

    # Build vertex A: internal gluon leg + eikonal legs
    legs_a = [Leg("A", f"g{internal_gluon}")]  # +g0 on vertex A
    for gi in range(gluon_start + 1, gluon_start + 1 + N_a):
        legs_a.append(Leg("A", f"-g{gi}"))
    verts.append(Vert(legs_a))

    # Build vertex B: internal gluon leg + eikonal legs
    legs_b = [Leg("A", f"-g{internal_gluon}")]  # -g0 on vertex B (opposite sign)
    for gi in range(gluon_start + 1 + N_a, gluon_idx):
        legs_b.append(Leg("A", f"-g{gi}"))
    verts.append(Vert(legs_b))

    # Build propagators: internal gluon + eikonal gluons
    props = []
    for gi in range(gluon_start, gluon_idx):
        props.append(Prop("A", "A", f"g{gi}"))

    return Diagram(
        idx=idx,
        sign=1,
        sym_factor="1",
        props=props,
        verts=verts,
        emission_type=et_name,
    )


def _generate_two_internal(
    et: EmissionType,
    start_idx: int,
) -> List[Diagram]:
    """Generate all diagrams for types with two internal vertices coupled by a gluon.

    Covers: 4e-3g+3g, 6e-3g+3g (and future types like 3e-3g+4g, etc.)
    """
    E = et.eikonal_count
    p = et.wl_wl_pairs
    internal_names = list(et.internal.keys())
    name_a = internal_names[0]
    name_b = internal_names[1] if len(internal_names) > 1 else internal_names[0]

    _, N_total_a = _name_to_mn(name_a)
    _, N_total_b = _name_to_mn(name_b)

    # Available legs: total legs minus 1 (consumed by internal gluon)
    N_a = N_total_a - 1
    N_b = N_total_b - 1

    diagrams = []
    idx = start_idx

    if p > 0:
        edge_distributions = _all_edge_distributions(E)
    else:
        edge_distributions = _canonical_edge_distributions(E)

    for ec in edge_distributions:
        positions = _contour_positions(ec)
        all_pos = list(range(E))

        # Choose which positions connect to vertex A
        for pos_a in combinations(all_pos, N_a):
            remaining_a = [p for p in all_pos if p not in set(pos_a)]

            if p == 0:
                # All remaining positions go to vertex B
                diagrams.append(_build_two_internal_diagram(
                    idx, ec, list(pos_a), remaining_a,
                    name_a, name_b, gluon_start=1, et_name=et.name,
                ))
                idx += 1
            else:
                # Choose N_b positions for vertex B from remaining
                for pos_b in combinations(remaining_a, N_b):
                    wl_wl_positions = [p for p in remaining_a if p not in set(pos_b)]

                    # Enumerate WL-WL pairings of the remaining positions
                    for matching in _perfect_matchings(wl_wl_positions):
                        wl_wl_flat = []
                        valid = True
                        for pa, pb in matching:
                            edge_a = positions[pa][0]
                            edge_b = positions[pb][0]
                            if edge_a == edge_b:
                                valid = False
                                break
                            wl_wl_flat.extend([pa, pb])
                        if not valid:
                            continue

                        diagrams.append(_build_two_internal_diagram(
                            idx, ec, list(pos_a), list(pos_b),
                            name_a, name_b, gluon_start=1, et_name=et.name,
                            wl_wl_positions=wl_wl_flat,
                        ))
                        idx += 1

    return diagrams


def generate_diagrams(
    n: int,
    emission_types: List[EmissionType],
) -> List[Diagram]:
    """Generate all valid Wilson loop diagrams for n-loop.

    Args:
        n: Number of loops (used for logging/validation only)
        emission_types: Pre-derived emission types from enumerate_emission_types(n)

    Returns:
        List of synthetic Diagram objects, each with a unique idx.
        These can be fed directly into classify_orbits().
    """
    all_diagrams: List[Diagram] = []
    idx = 1

    for et in emission_types:
        if et.name == "2e-Pi1+Pi1":
            from .self_energy_chain import generate_pi1_chain
            diagrams = generate_pi1_chain(et, idx)
        elif not et.internal:
            # Pure WL-WL: e.g. 6e, 4e
            diagrams = _generate_wl_wl(et, idx)
        elif sum(et.internal.values()) == 1:
            # One internal vertex: e.g. 5e-3g, 4e-4g, 4e-Pi1, 2e-Pi2, 3e-Gamma1
            diagrams = _generate_one_internal(et, idx)
        elif sum(et.internal.values()) == 2:
            # Two internal vertices coupled by a gluon: e.g. 4e-3g+3g
            diagrams = _generate_two_internal(et, idx)
        else:
            # Future: more complex internal structures
            diagrams = []

        all_diagrams.extend(diagrams)
        idx += len(diagrams)

    return all_diagrams
