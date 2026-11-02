"""Classify a parsed Diagram: discard scaleless, classify tree/loop,
extract vertex counts and edge distribution.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict

from .parser import Diagram, Vert, mom_dict, mom_eq, mom_neg


WILSON_FIELDS = {f"z{i}" for i in (1, 2, 3, 4)}
WILSON_BAR_FIELDS = {f"z{i}bar" for i in (1, 2, 3, 4)}
CORNER_FIELDS = {f"phi{i}" for i in (1, 2, 3, 4)}


def edge_of_wilson(field: str) -> int:
    s = field
    if s.endswith("bar"):
        s = s[:-3]
    assert s.startswith("z"), s
    return int(s[1:])


def vert_kind(v: Vert) -> str:
    fields = [l.field for l in v.legs]
    if any(f in CORNER_FIELDS for f in fields):
        return "corner"
    if len(fields) == 4 and all(f == "A" for f in fields):
        return "four_gluon"
    if len(fields) == 3 and all(f == "A" for f in fields):
        return "three_gluon"
    if len(fields) == 2 and all(f == "A" for f in fields):
        # 2-point A-A vertex: a self-energy insertion (Pi-type loop).
        # Tree-level T2 would be a bare propagator (handled as WL-WL),
        # so any actual 2-gluon vertex is a loop self-energy.
        return "two_gluon"
    if any(f in WILSON_BAR_FIELDS for f in fields) and any(
        f in WILSON_FIELDS for f in fields):
        return "eikonal"
    if any(f.startswith("c") for f in fields) and "A" in fields:
        return "ghost_gluon"
    if any(f.startswith("q") for f in fields) and "A" in fields:
        return "quark_gluon"
    return "other"


@dataclass
class TopologyInfo:
    sector: str = ""              # "tree" or "loop"
    is_scaleless: bool = False
    scaleless_reason: str = ""    # "A" or "B"
    eikonal_count: int = 0
    edge_counts: Tuple[int, int, int, int] = (0, 0, 0, 0)
    n_three_gluon: int = 0
    n_four_gluon: int = 0
    n_ghost_gluon: int = 0
    n_quark_gluon: int = 0
    n_wl_wl_gluons: int = 0
    n_complex_gluons: int = 0
    n_two_gluon: int = 0
    wl_wl_pairs: List[Tuple[int, int, int, int]] = field(default_factory=list)
    # (edge_a, slot_a, edge_b, slot_b) — direct WL-WL exchanges only
    gluon_connections: List[Tuple[int, List[Tuple[int, int]]]] = field(default_factory=list)
    # [(gluon_prop_idx, [(edge, slot), ...]), ...] — ALL gluon-to-eikonal connections
    internal_connections: List[Tuple[int, ...]] = field(default_factory=list)
    # For each internal (non-eikonal) vertex, the sorted tuple of contour
    # positions (1-indexed, in `edge_counts` frame) of the eikonal vertices
    # it attaches to.  Used to canonicalise the internal connectivity.
    internal_gluon_links: List[Tuple[int, int]] = field(default_factory=list)
    # Pairs of indices into `internal_connections` joined by a dedicated
    # internal gluon line (no eikonal endpoint), e.g. the gluon coupling
    # the two 3g vertices in 4e-3g+3g.  Such linked vertices form a single
    # connected sub-structure for the 2PI cut test.
    loop_momenta: Set[str] = field(default_factory=set)


def _get_eikonal_edge(v: Vert) -> Optional[int]:
    for leg in v.legs:
        if leg.field in WILSON_FIELDS | WILSON_BAR_FIELDS:
            return edge_of_wilson(leg.field)
    return None


def _get_gluon_momentum(v: Vert) -> Optional[str]:
    for leg in v.legs:
        if leg.field == "A":
            return leg.mom
    return None


def classify_topology(d: Diagram) -> TopologyInfo:
    info = TopologyInfo()

    eikonal_verts: List[Tuple[int, Vert, int]] = []
    for vi, v in enumerate(d.verts):
        kind = vert_kind(v)
        if kind == "eikonal":
            edge = _get_eikonal_edge(v)
            if edge is not None:
                eikonal_verts.append((edge, v, vi))
        elif kind == "three_gluon":
            info.n_three_gluon += 1
        elif kind == "four_gluon":
            info.n_four_gluon += 1
        elif kind == "two_gluon":
            info.n_two_gluon += 1
        elif kind == "ghost_gluon":
            info.n_ghost_gluon += 1
        elif kind == "quark_gluon":
            info.n_quark_gluon += 1

    edge_counts = [0, 0, 0, 0]
    for edge, _, _ in eikonal_verts:
        edge_counts[edge - 1] += 1
    info.edge_counts = tuple(edge_counts)  # type: ignore[assignment]
    info.eikonal_count = len(eikonal_verts)

    # Gluon classification
    gluon_props = [(i, p) for i, p in enumerate(d.props)
                   if p.field1 == "A" and p.field2 == "A"]

    # Type B: all eikonal vertices live on the same Wilson-line edge.
    # Such diagrams are scaleless (no net color flow between distinct
    # Wilson lines => vacuum-like subdiagram), regardless of whether
    # gluon propagators connect the eikonal vertices on that single edge.
    non_zero = [e for e in edge_counts if e > 0]
    if len(non_zero) == 1:
        info.is_scaleless = True
        info.scaleless_reason = "B"
        return info

    eikonal_gluon_momenta: Dict[int, str] = {}
    for edge, v, vi in eikonal_verts:
        mom = _get_gluon_momentum(v)
        if mom:
            eikonal_gluon_momenta[vi] = mom

    # Compute slot numbers: for each edge, sort eikonal vertices by
    # vertex index (QGRAF emits in contour order, so vertex index
    # correlates with path-ordering position).
    edge_slots: Dict[int, Dict[int, int]] = {1: {}, 2: {}, 3: {}, 4: {}}
    for edge in (1, 2, 3, 4):
        verts_on_edge = sorted(
            [(vi, v) for vi, v in enumerate(d.verts)
             if vert_kind(v) == "eikonal" and _get_eikonal_edge(v) == edge],
            key=lambda x: x[0])
        for slot, (vi, _) in enumerate(verts_on_edge, 1):
            edge_slots[edge][vi] = slot

    # A WL-WL gluon connects exactly ONE eikonal vertex whose gluon
    # momentum equals the propagator momentum and exactly ONE whose
    # gluon momentum is the negative.  Using mom_eq OR mom_neg for
    # both ends fails when multiple gluons carry the same momentum.
    used_vis: Set[int] = set()
    wl_wl_pairs = []
    ambiguous = False
    for _, prop in gluon_props:
        eq_vis = [vi for vi, mom in eikonal_gluon_momenta.items()
                  if mom_eq(prop.mom, mom)]
        neg_vis = [vi for vi, mom in eikonal_gluon_momenta.items()
                   if mom_neg(prop.mom, mom)]
        if len(eq_vis) == 1 and len(neg_vis) == 1:
            vi_a, vi_b = eq_vis[0], neg_vis[0]
            edge_a = _get_eikonal_edge(d.verts[vi_a])
            edge_b = _get_eikonal_edge(d.verts[vi_b])
            if edge_a is not None and edge_b is not None:
                slot_a = edge_slots[edge_a].get(vi_a, 1)
                slot_b = edge_slots[edge_b].get(vi_b, 1)
                used_vis.add(vi_a)
                used_vis.add(vi_b)
                wl_wl_pairs.append((edge_a, slot_a, edge_b, slot_b))
        else:
            ambiguous = True

    # Remove duplicate pairs: if two gluons connect the exact same
    # pair of (edge, slot) endpoints, it is a blob insertion
    # (e.g. 4g-tadpole), not a direct WL-WL exchange.
    seen_pairs: Dict[frozenset, int] = {}
    for pair in wl_wl_pairs:
        key = frozenset([(pair[0], pair[1]), (pair[2], pair[3])])
        seen_pairs[key] = seen_pairs.get(key, 0) + 1
    wl_wl_pairs = [p for p in wl_wl_pairs
                   if seen_pairs.get(frozenset([(p[0], p[1]), (p[2], p[3])]), 0) == 1]

    # Type A fallback: when momentum matching is ambiguous (multiple
    # gluons carry the same momenta), use the pigeonhole principle:
    # if any edge has more eikonal vertices than there are gluon
    # propagators, at least one gluon must be a same-edge pair.
    if ambiguous and not wl_wl_pairs:
        max_on_edge = max(edge_counts)
        if max_on_edge > len(gluon_props):
            info.is_scaleless = True
            info.scaleless_reason = "A"
            return info

    # Type A scaleless: gluon connecting two eikonal vertices on the SAME edge
    # => v_i^2 = 0 => scaleless
    for (edge_a, slot_a, edge_b, slot_b) in wl_wl_pairs:
        if edge_a == edge_b:
            info.is_scaleless = True
            info.scaleless_reason = "A"
            return info

    info.wl_wl_pairs = wl_wl_pairs  # type: ignore[assignment]
    info.n_wl_wl_gluons = len(wl_wl_pairs)
    info.n_complex_gluons = len(gluon_props) - info.n_wl_wl_gluons

    # Populate gluon_connections: for EVERY gluon propagator, record
    # which eikonal vertices it connects to (edge, slot).
    # This is needed by the emitter to write coordinate differences
    # for ALL gluon legs, not just WL-WL pairs.
    gluon_connections = []
    for pi, prop in gluon_props:
        endpoints = []
        for vi, mom in eikonal_gluon_momenta.items():
            if mom_eq(prop.mom, mom) or mom_neg(prop.mom, mom):
                edge = _get_eikonal_edge(d.verts[vi])
                if edge is not None:
                    slot = edge_slots[edge].get(vi, 1)
                    endpoints.append((edge, slot))
        if endpoints:
            gluon_connections.append((pi, endpoints))
    info.gluon_connections = gluon_connections

    # Sector: "loop" if any non-Wilson propagator (gluon, ghost, or quark)
    # carries a loop momentum not present in Wilson line propagators
    # (genuine internal loop — gluon/ghost/quark self-energy or vertex);
    # "tree" otherwise (loop closes entirely through Wilson lines).
    wilson_k: Set[str] = set()
    internal_k: Set[str] = set()
    for prop in d.props:
        km = {sym for sym in mom_dict(prop.mom) if sym.startswith("k")}
        if prop.field1.startswith("z") and prop.field2.startswith("z"):
            wilson_k |= km
        elif (prop.field1 == "A" and prop.field2 == "A"
              or prop.field1.startswith("c") and prop.field2.startswith("c")
              or prop.field1.startswith("q") and prop.field2.startswith("q")):
            internal_k |= km
    extra_internal_k = internal_k - wilson_k
    info.sector = "loop" if extra_internal_k else "tree"

    for prop in d.props:
        for sym in mom_dict(prop.mom):
            if sym.startswith("k"):
                info.loop_momenta.add(sym)

    # Extract per-internal-vertex eikonal endpoint positions.
    # For each non-eikonal gluon vertex (3g/4g/Pi), find which eikonal
    # vertices its gluon legs attach to (by matching ±momentum) and
    # record their contour positions.  This captures the partition of
    # eikonal positions among internal vertices, which the classifier
    # canonicalises to distinguish orbits that share the same edge
    # distribution and WL-WL chord pairing.
    internal_conns: List[Tuple[int, ...]] = []
    internal_vert_indices: List[int] = []
    internal_vert_edges: List[Set[int]] = []
    for vi, v in enumerate(d.verts):
        kind = vert_kind(v)
        if kind not in ("three_gluon", "four_gluon", "two_gluon"):
            continue
        positions: Set[int] = set()
        edges: Set[int] = set()
        for leg in v.legs:
            if leg.field != "A":
                continue
            for ek_vi, ek_mom in eikonal_gluon_momenta.items():
                if mom_eq(leg.mom, ek_mom) or mom_neg(leg.mom, ek_mom):
                    ek_edge = _get_eikonal_edge(d.verts[ek_vi])
                    if ek_edge is not None:
                        slot = edge_slots[ek_edge].get(ek_vi, 1)
                        pos = 0
                        for e in range(1, ek_edge):
                            pos += edge_counts[e - 1]
                        pos += slot
                        positions.add(pos)
                        edges.add(ek_edge)
        if positions:
            internal_conns.append(tuple(sorted(positions)))
            internal_vert_indices.append(vi)
            internal_vert_edges.append(edges)

    info.internal_connections = internal_conns

    # Detect internal gluon lines that couple two internal vertices
    # (e.g. the gluon joining the two 3g vertices in 4e-3g+3g).  Such a
    # gluon propagator carries a momentum that appears on two internal
    # vertices (with opposite signs) and on NO eikonal vertex.  The
    # classifier merges the two internal-vertex connection sets joined
    # by an internal gluon, since they form a single connected
    # sub-structure for the 2PI cut test.
    internal_gluon_links: List[Tuple[int, int]] = []
    ek_mom_values = list(eikonal_gluon_momenta.values())
    for _pi, prop in gluon_props:
        m = prop.mom
        # Skip gluons that touch an eikonal vertex (WL-WL or eikonal-to-vertex).
        if any(mom_eq(m, em) or mom_neg(m, em) for em in ek_mom_values):
            continue
        matching = []
        for idx, vi in enumerate(internal_vert_indices):
            v = d.verts[vi]
            if any(mom_eq(lg.mom, m) or mom_neg(lg.mom, m)
                   for lg in v.legs if lg.field == "A"):
                matching.append(idx)
        for a in range(len(matching)):
            for b in range(a + 1, len(matching)):
                link = (matching[a], matching[b])
                if link not in internal_gluon_links:
                    internal_gluon_links.append(link)
    info.internal_gluon_links = internal_gluon_links

    # Scaleless type D: an internal (non-eikonal) gluon vertex whose
    # eikonal endpoints all lie on a single Wilson-line edge.  Such a
    # sub-diagram is a single-line correction (self-energy / vertex
    # renormalisation on one Wilson line) with no inter-line colour
    # flow, hence not a genuine web exchange.
    #
    # Internal vertices connected by a dedicated internal gluon line
    # (e.g. the two 3g vertices in 4e-3g+3g) form a single connected
    # sub-structure — their edge sets must be merged before checking
    # whether the combined structure spans multiple Wilson lines.
    if internal_vert_edges:
        # Union-find: merge edge sets of vertices linked by internal gluons
        n_iv = len(internal_vert_edges)
        parent = list(range(n_iv))

        def _find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, b in internal_gluon_links:
            if a < n_iv and b < n_iv:
                parent[_find(a)] = _find(b)

        merged_edges: Dict[int, set] = {}
        for i, edges in enumerate(internal_vert_edges):
            root = _find(i)
            merged_edges.setdefault(root, set()).update(edges)

        for edges in merged_edges.values():
            if len(edges) <= 1:
                info.is_scaleless = True
                info.scaleless_reason = "D"
                return info

    return info
