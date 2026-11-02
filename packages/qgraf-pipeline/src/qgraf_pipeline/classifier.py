"""Group diagrams into orbits: C4 rotation + mirror + 2PI screening."""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from typing import Dict, List, Optional, Set, Tuple

from .parser import Diagram
from .topology import TopologyInfo, classify_topology
from .internal_types import EmissionType


@dataclass
class ClassifiedOrbit:
    orbit_id: str
    emission_type: str
    category: str
    edge_counts: Tuple[int, int, int, int]
    eikonal_count: int
    g_power: int
    n_diagrams: int
    wl_wl_pairs: int = 0
    internal: Dict[str, int] = field(default_factory=dict)
    representative: Optional[Diagram] = None
    mirror_relation: str = "self-mirror"
    mirror_of: Optional[str] = None
    is_2pi: bool = True
    reducible_to: List[str] = field(default_factory=list)
    diagram_indices: List[int] = field(default_factory=list)
    canon_chords: frozenset = field(default_factory=frozenset)
    canon_internal: frozenset = field(default_factory=frozenset)


def rotate_edge_counts(
    ec: Tuple[int, int, int, int], r: int
) -> Tuple[int, int, int, int]:
    return tuple(ec[(i - r) % 4] for i in range(4))  # type: ignore[return-value]


def canonical_edge_counts(
    ec: Tuple[int, int, int, int]
) -> Tuple[int, int, int, int]:
    best = ec
    for r in range(1, 4):
        rotated = rotate_edge_counts(ec, r)
        if rotated < best:
            best = rotated
    return best


def mirror_edge_counts(
    ec: Tuple[int, int, int, int]
) -> Tuple[int, int, int, int]:
    return (ec[0], ec[3], ec[2], ec[1])


def _contour_order(edge: int, slot: int, edge_counts: Tuple[int, int, int, int]) -> int:
    """Return the 1-indexed contour position of (edge, slot)."""
    pos = 0
    for e in range(1, edge):
        pos += edge_counts[e - 1]
    pos += slot
    return pos


def _connection_sets(info: TopologyInfo) -> List[frozenset]:
    """Build the list of contour-position connection sets.

    Each WL-WL pair contributes a 2-element chord; each internal
    (non-eikonal) gluon vertex contributes the set of eikonal positions
    it attaches to.  Two internal vertices joined by a dedicated
    internal gluon line (detected in ``info.internal_gluon_links``) are
    merged into a single connection set, because the internal gluon line
    makes them one connected sub-structure spanning the union of their
    eikonal endpoints — the 2PI cut test must not be able to separate
    them.
    """
    ec = info.edge_counts
    conns: List[frozenset] = []
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        a = _contour_order(edge_a, slot_a, ec)
        b = _contour_order(edge_b, slot_b, ec)
        conns.append(frozenset([a, b]))

    internal_sets = [frozenset(c) for c in info.internal_connections]

    # Union-Find over internal sets, merging those linked by an
    # internal gluon line (and, defensively, those that happen to share
    # a contour position).
    parent = list(range(len(internal_sets)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(internal_sets)):
        for j in range(i + 1, len(internal_sets)):
            if internal_sets[i] & internal_sets[j]:
                parent[find(i)] = find(j)
    for a, b in info.internal_gluon_links:
        if a < len(internal_sets) and b < len(internal_sets):
            parent[find(a)] = find(b)

    groups: Dict[int, set] = {}
    for i, s in enumerate(internal_sets):
        groups.setdefault(find(i), set()).update(s)

    conns.extend(frozenset(g) for g in groups.values())
    return conns


def _has_one_sided_wrap(connections: List[frozenset]) -> bool:
    """Detect a one-sided wrapping factorisation.

    A connection A "one-sided wraps" the rest of the diagram when EVERY
    other connection lies strictly inside A's contour interval AND no
    position of A lies strictly inside any other connection's interval.
    This is the 1-loop × (n-1)-loop web factorisation: the outer
    connection is a 1-loop sub-web whose interval contains the entire
    (n-1)-loop sub-web, with no mutual crossing.

    Partial wraps (A contains some but not all other connections) do
    NOT factorise the whole diagram — the remaining connections tie the
    outer web to the inner one — and are kept as genuine 2PI webs.
    Mutual wraps (where A has a position inside another connection's
    interval) are also kept.
    """
    for a in connections:
        amin, amax = min(a), max(a)
        wraps_all = True
        for b in connections:
            if a is b or a == b:
                continue
            bmin, bmax = min(b), max(b)
            # B must be strictly inside A's interval.
            if not (bmin > amin and bmax < amax):
                wraps_all = False
                break
            # Must be one-sided (no A position inside B's interval).
            if any(bmin < p < bmax for p in a):
                wraps_all = False
                break
        if wraps_all and len(connections) > 1:
            return True
    return False


def _check_2pi(info: TopologyInfo) -> Tuple[bool, List[str]]:
    """Check if a diagram is a proper irreducible (2PI) web.

    A web is reducible (1PR — factorises into independent sub-webs) iff
    there exists a 2-cut of the contour that separates the eikonal
    positions into two non-empty contiguous arcs with NO connection
    crossing the cut.  A connection (WL-WL chord or internal-vertex
    hyperchord) crosses the cut iff it has endpoints in both arcs.

    For pure WL-WL chord diagrams this is equivalent to the chord
    (overlap) graph being connected, and also subsumes the "wrapping
    gluon" factorisation case.

    An additional one-sided wrapping check rejects diagrams where one
    connection's interval strictly contains another with no mutual
    crossing — the 1-loop × (n-1)-loop factorisation that survives the
    cut test for mixed diagrams (e.g. a 3g vertex spanning the whole
    contour around an adjacent WL-WL pair).
    """
    E = info.eikonal_count
    if E <= 1:
        return True, []

    connections = _connection_sets(info)
    if not connections:
        return True, []

    # Enumerate all contiguous 2-cuts: arc A = positions (s, s+1, ...,
    # s+L-1) mod E, for start s in 0..E-1 and length L in 1..E-1.
    for s in range(E):
        for L in range(1, E):
            A = {((s + i) % E) + 1 for i in range(L)}
            crossed = False
            for conn in connections:
                in_a = any(p in A for p in conn)
                out_a = any(p not in A for p in conn)
                if in_a and out_a:
                    crossed = True
                    break
            if not crossed:
                return False, ["factorizable-cut"]

    if _has_one_sided_wrap(connections):
        return False, ["one-sided-wrap"]

    return True, []


def _build_type_map(
    emission_types: List[EmissionType],
) -> Dict[tuple, EmissionType]:
    """Build lookup: (E, has_loop, internal_key) -> EmissionType.

    internal_key = frozenset of internal structure items, distinguishing
    e.g. 4e-3g ({"3g":2}) from 4e-4g ({"4g":1}) which both have E=4, tree.
    """
    type_map: Dict[tuple, EmissionType] = {}
    for et in emission_types:
        if et.name == "3e-3g+Pi1":
            continue  # explicit mixed two-point/three-point structure below
        has_loop = et.category == "loop"
        # For loop types, use empty key (g-power constraint gives unique match)
        # For tree types, use internal structure to distinguish 4e-3g vs 4e-4g
        if has_loop:
            ikey = frozenset()
        else:
            ikey = frozenset(et.internal.items()) if et.internal else frozenset()
        type_map[(et.eikonal_count, has_loop, ikey)] = et
    return type_map


def _diagram_internal_key(info: TopologyInfo) -> frozenset:
    """Compute internal structure key from diagram vertex counts.

    For tree diagrams: use vertex counts directly (3g, 4g).
    For loop diagrams: use empty key — (E, has_loop) uniquely determines
    the type via g-power constraint.
    """
    if info.sector == "tree":
        internal: Dict[str, int] = {}
        if info.n_three_gluon > 0:
            internal["3g"] = info.n_three_gluon
        if info.n_four_gluon > 0:
            internal["4g"] = info.n_four_gluon
        return frozenset(internal.items())
    else:
        return frozenset()


def _canonical_chord_pairing(
    info, canon_ec: Tuple[int, int, int, int]
) -> frozenset:
    """Compute canonical chord pairing under C4 rotation.

    Returns a frozenset of frozensets, each inner frozenset being
    {pos_a, pos_b} for a chord.  Empty frozenset if no WL-WL pairs.
    """
    wl_wl_pairs = info.wl_wl_pairs
    if not wl_wl_pairs:
        return frozenset()

    ec = info.edge_counts
    best_chords = None
    for r in range(4):
        rotated_ec = rotate_edge_counts(ec, r)
        if rotated_ec != canon_ec:
            continue

        chords = []
        for edge_a, slot_a, edge_b, slot_b in wl_wl_pairs:
            new_edge_a = ((edge_a - 1 + r) % 4) + 1
            new_edge_b = ((edge_b - 1 + r) % 4) + 1
            pos_a = _contour_order(new_edge_a, slot_a, canon_ec)
            pos_b = _contour_order(new_edge_b, slot_b, canon_ec)
            chords.append(frozenset([pos_a, pos_b]))

        chords_key = frozenset(chords)
        if best_chords is None:
            best_chords = chords_key
        else:
            best_sorted = sorted(sorted(c) for c in best_chords)
            curr_sorted = sorted(sorted(c) for c in chords_key)
            if curr_sorted < best_sorted:
                best_chords = chords_key

    return best_chords  # type: ignore[return-value]


def _canonical_internal_pairing(
    info, canon_ec: Tuple[int, int, int, int]
) -> frozenset:
    """Compute canonical internal-vertex connectivity under C4 rotation.

    Returns a frozenset of frozensets, each inner frozenset being the set
    of contour positions attached to one internal vertex.  Empty
    frozenset if the diagram has no internal vertices.
    """
    if not info.internal_connections:
        return frozenset()

    ec = info.edge_counts
    # Build per-vertex (edge, slot) endpoint lists so we can rotate them.
    # Re-derive from the diagram: for each internal vertex, the list of
    # (edge, slot) of its eikonal endpoints, in the original frame.
    # internal_connections already stores contour positions in the
    # original (ec) frame; we rotate those positions by mapping each
    # position back to (edge, slot), rotating the edge, and recomputing
    # the contour position in canon_ec.
    def _pos_to_edge_slot(pos: int, ec_tuple: Tuple[int, int, int, int]
                         ) -> Tuple[int, int]:
        cum = 0
        for edge in range(1, 5):
            cnt = ec_tuple[edge - 1]
            if pos <= cum + cnt:
                return edge, pos - cum
            cum += cnt
        return 1, 1

    best = None
    for r in range(4):
        rotated_ec = rotate_edge_counts(ec, r)
        if rotated_ec != canon_ec:
            continue

        vertex_sets = []
        for conn in info.internal_connections:
            new_positions = []
            for pos in conn:
                edge, slot = _pos_to_edge_slot(pos, ec)
                new_edge = ((edge - 1 + r) % 4) + 1
                new_positions.append(_contour_order(new_edge, slot, canon_ec))
            vertex_sets.append(frozenset(new_positions))

        key = frozenset(vertex_sets)
        if best is None:
            best = key
        else:
            best_sorted = sorted(sorted(c) for c in best)
            curr_sorted = sorted(sorted(c) for c in key)
            if curr_sorted < best_sorted:
                best = key

    return best  # type: ignore[return-value]


def classify_orbits(
    diagrams: List[Diagram],
    emission_types: List[EmissionType],
    *,
    include_reducible: bool = False,
) -> List[ClassifiedOrbit]:
    infos: List[Tuple[Diagram, TopologyInfo]] = []
    for d in diagrams:
        info = classify_topology(d)
        if info.is_scaleless:
            continue
        # Discard 1PR diagrams (can be factorized into lower-order products)
        is_2pi, _ = _check_2pi(info)
        if not is_2pi and not include_reducible:
            continue
        infos.append((d, info))

    type_map = _build_type_map(emission_types)
    by_name: Dict[str, EmissionType] = {et.name: et for et in emission_types}

    groups: Dict[Tuple[str, Tuple[int, int, int, int], frozenset, frozenset],
                 List[Diagram]] = {}
    for d, info in infos:
        # Use the generator's emission-type tag when present (synthetic
        # diagrams); otherwise fall back to structural reconstruction.
        if d.emission_type and d.emission_type in by_name:
            et = by_name[d.emission_type]
            et_name = et.name
        elif (info.eikonal_count == 3 and info.n_two_gluon == 1
              and info.n_three_gluon == 1 and "3e-3g+Pi1" in by_name):
            et = by_name["3e-3g+Pi1"]
            et_name = et.name
        else:
            key = (info.eikonal_count, info.sector == "loop",
                   _diagram_internal_key(info))
            et = type_map.get(key)
            et_name = et.name if et else f"{info.eikonal_count}e-unknown"

        # Skip diagrams where the momentum matching is ambiguous:
        # the emission type expects WL-WL pairs but the diagram has fewer.
        # Only applies to pure-WL-WL types (no internal structure) where
        # the chord pairing IS the entire structure.
        if (et is not None and et.wl_wl_pairs > 0
                and not et.internal
                and info.n_wl_wl_gluons < et.wl_wl_pairs):
            continue

        canon_ec = canonical_edge_counts(info.edge_counts)
        canon_chords = _canonical_chord_pairing(info, canon_ec)
        canon_internal = _canonical_internal_pairing(info, canon_ec)
        group_key = (et_name, canon_ec, canon_chords, canon_internal)
        groups.setdefault(group_key, []).append(d)

    orbits: List[ClassifiedOrbit] = []
    orbit_counter: Dict[str, int] = {}

    for (et_name, canon_ec, canon_chords, canon_internal), diags in sorted(
            groups.items()):
        # Audit-only graphs must never renumber existing production orbits.
        group_is_2pi, _ = _check_2pi(classify_topology(diags[0]))
        if group_is_2pi:
            orbit_counter[et_name] = orbit_counter.get(et_name, 0) + 1
            oid = f"{et_name}-{orbit_counter[et_name]:02d}"
        else:
            identity = [et_name, list(canon_ec),
                        sorted(sorted(c) for c in canon_chords),
                        sorted(sorted(c) for c in canon_internal)]
            digest = hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:16]
            oid = f"{et_name}-aux-{digest}"

        # Choose representative: prefer diagrams with parallel (slot-matching)
        # chord pairing over crossed pairing for the same edge distribution.
        def _pairing_score(diag):
            info = classify_topology(diag)
            if not info.wl_wl_pairs:
                return 0
            # Count how many pairs have matching slots (parallel)
            matching = sum(1 for a, sa, b, sb in info.wl_wl_pairs if sa == sb)
            return matching  # higher = more parallel = preferred

        diags_sorted = sorted(diags, key=_pairing_score, reverse=True)
        rep_diag = diags_sorted[0]

        rep_info = classify_topology(rep_diag)
        is_2pi, reducible = _check_2pi(rep_info)

        # Find the matching emission type for metadata
        et_meta = by_name.get(et_name)
        if et_meta is None:
            g_power = 2 * 2
            category = rep_info.sector
        else:
            g_power = et_meta.g_power
            # For synthetic diagrams the structural sector detection is
            # unreliable (loop vertices carry no loop-momentum symbols),
            # so prefer the emission type's declared category.
            category = et_meta.category

        orbits.append(ClassifiedOrbit(
            orbit_id=oid,
            emission_type=et_name,
            category=category,
            edge_counts=canon_ec,
            eikonal_count=rep_info.eikonal_count,
            g_power=g_power,
            n_diagrams=len(diags),
            wl_wl_pairs=rep_info.n_wl_wl_gluons,
            internal=et_meta.internal if et_meta else {},
            representative=rep_diag,
            is_2pi=is_2pi,
            reducible_to=reducible,
            diagram_indices=[d.idx for d in diags],
            canon_chords=canon_chords,
            canon_internal=canon_internal,
        ))

    # Preserve the legacy global mirror labels when appending a new family.
    legacy = [o for o in orbits if o.is_2pi and o.emission_type != "3e-3g+Pi1"]
    added = [o for o in orbits if o.is_2pi and o.emission_type == "3e-3g+Pi1"]
    _assign_mirrors(legacy)
    _assign_mirrors(added)
    for orbit in added:
        if orbit.mirror_of is not None:
            orbit.mirror_relation = "tgv-pi1-" + orbit.mirror_relation
    auxiliary = [o for o in orbits if not o.is_2pi]
    _assign_mirrors(auxiliary)
    for orbit in auxiliary:
        if orbit.mirror_of is not None:
            orbit.mirror_relation = "audit-" + orbit.mirror_relation
    return orbits


def _mirror_positions(
    positions: frozenset, canon_ec: Tuple[int, int, int, int]
) -> frozenset:
    """Mirror a set of contour positions under the Wilson-loop mirror.

    Mirror the edge distribution, canonicalise back via C4 rotation,
    and transform the (edge, slot) of each position accordingly.
    """
    if not positions:
        return frozenset()

    mir_ec = mirror_edge_counts(canon_ec)
    mir_canon_ec = canonical_edge_counts(mir_ec)

    # Find rotation r such that rotate_edge_counts(mir_ec, r) == mir_canon_ec
    r = 0
    for rr in range(4):
        if rotate_edge_counts(mir_ec, rr) == mir_canon_ec:
            r = rr
            break

    # Mirror edge map: 1→1, 2→4, 3→3, 4→2
    mir_edge_map = {1: 1, 2: 4, 3: 3, 4: 2}

    def _pos_to_edge_slot(pos, ec):
        cum = 0
        for edge in range(1, 5):
            cnt = ec[edge - 1]
            if pos <= cum + cnt:
                return edge, pos - cum
            cum += cnt
        return 1, 1  # fallback

    result = []
    for pos in positions:
        edge, slot = _pos_to_edge_slot(pos, canon_ec)
        mir_edge = mir_edge_map[edge]
        # Reverse slot ordering on mirrored edge
        mir_slot = mir_ec[mir_edge - 1] + 1 - slot
        # Rotate to canonical frame
        rot_edge = ((mir_edge - 1 + r) % 4) + 1
        new_pos = _contour_order(rot_edge, mir_slot, mir_canon_ec)
        result.append(new_pos)

    return frozenset(result)


def _mirror_canon_chords(
    canon_chords: frozenset, canon_ec: Tuple[int, int, int, int]
) -> frozenset:
    """Compute the mirror image of a canonical chord pairing."""
    if not canon_chords:
        return frozenset()
    return frozenset(_mirror_positions(chord, canon_ec) for chord in canon_chords)


def _mirror_canon_internal(
    canon_internal: frozenset, canon_ec: Tuple[int, int, int, int]
) -> frozenset:
    """Compute the mirror image of the canonical internal connectivity."""
    if not canon_internal:
        return frozenset()
    return frozenset(
        _mirror_positions(vertex_set, canon_ec)
        for vertex_set in canon_internal
    )


def _assign_mirrors(orbits: List[ClassifiedOrbit]) -> None:
    """Assign mirror relations between orbits.

    Two orbits are a mirror pair if their canonical edge distributions,
    chord pairings, and internal connectivities are all mirror images.
    """
    # Index orbits by (edge_counts, canon_chords, canon_internal) for lookup
    by_key: Dict[Tuple[Tuple[int, int, int, int], frozenset, frozenset], int] = {}
    for i, o in enumerate(orbits):
        by_key[(o.edge_counts, o.canon_chords, o.canon_internal)] = i

    mirror_pair = 1
    assigned: Set[int] = set()

    for i, o in enumerate(orbits):
        if i in assigned:
            continue

        mir_ec = canonical_edge_counts(mirror_edge_counts(o.edge_counts))
        mir_chords = _mirror_canon_chords(o.canon_chords, o.edge_counts)
        mir_internal = _mirror_canon_internal(o.canon_internal, o.edge_counts)

        if (mir_ec == o.edge_counts
                and mir_chords == o.canon_chords
                and mir_internal == o.canon_internal):
            # Mirror maps to itself → self-mirror
            o.mirror_relation = "self-mirror"
            assigned.add(i)
        else:
            # Look for the mirror partner
            j = by_key.get((mir_ec, mir_chords, mir_internal))
            if j is not None and j != i and j not in assigned:
                o.mirror_relation = f"mirror-pair-{mirror_pair}"
                o.mirror_of = orbits[j].orbit_id
                orbits[j].mirror_relation = f"mirror-pair-{mirror_pair}"
                orbits[j].mirror_of = o.orbit_id
                assigned.add(i)
                assigned.add(j)
                mirror_pair += 1
            else:
                # No mirror partner found (shouldn't normally happen)
                o.mirror_relation = "self-mirror"
                assigned.add(i)
