"""Parse orbit_definitions.json into Orbit dataclasses."""
import re

from montopar.models import (
    Vertex, Propagator, Bound, Pairing,
    EdgeDistribution, Orbit, TGVLeg, FGVLeg, DoubleTGVLeg,
)


def parse_orbit(orbit_json: dict) -> Orbit:
    """Parse a single orbit definition into an Orbit dataclass.

    Args:
        orbit_json: A dict from orbit_definitions.json["orbits"][i].

    Returns:
        Orbit with parsed edge distribution, pairings, bounds, etc.
    """
    if ("Gamma" in orbit_json.get("emission_type", "")
            or "vertex_correction" in orbit_json
            or any(v.get("type") == "gamma1" for v in orbit_json.get("internal_vertices", []))):
        from montopar.gamma1_source import parse_gamma1_orbit
        return parse_gamma1_orbit(orbit_json)

    if orbit_json.get("emission_type") == "3e-3g+Pi1":
        from montopar.tgv_pi1 import parse_tgv_pi1_orbit
        return parse_tgv_pi1_orbit(orbit_json)
    if orbit_json.get("emission_type") == "2e-Pi2":
        from montopar.self_energy_pi2 import parse_pi2_orbit
        return parse_pi2_orbit(orbit_json)
    if orbit_json.get("emission_type") == "2e-Pi1+Pi1":
        from montopar.self_energy_chain import parse_chain_orbit
        return parse_chain_orbit(orbit_json)
    if ("-Pi" in orbit_json.get("emission_type", "")
            or any(v.get("type") == "se" for v in orbit_json.get("internal_vertices", []))):
        from montopar.self_energy import parse_self_energy_orbit
        return parse_self_energy_orbit(orbit_json)
    if orbit_json.get("emission_type") == "4e-3g+3g":
        _validate_double_tgv_topology(orbit_json)

    orbit_id = orbit_json["orbit_id"]
    emission_type = orbit_json["emission_type"]

    # Parse vertices
    vertices = [
        Vertex(line=v["line"], lorentz=v["lorentz"], color=v["color"])
        for v in orbit_json["vertices"]
    ]

    # Parse propagators
    propagators = [
        Propagator(
            loop=p["loop"],
            lorentz=(p["lorentz"][0], p["lorentz"][1]),
            color=(p["color"][0], p["color"][1]),
        )
        for p in orbit_json["propagators"]
    ]

    # Build lorentz → edge lookup from vertices
    lorentz_to_edge: dict[str, int] = {v.lorentz: v.line for v in vertices}

    # Build phases lookup: loop → expression
    phases: dict[str, str] = {
        p["loop"]: p["expression"] for p in orbit_json["phases"]
    }

    # Parse bounds
    bounds = [
        Bound(var=b["var"], low=str(b["low"]), high=str(b["high"]))
        for b in orbit_json["bounds"]
    ]

    # Build edge distribution and s_var mapping
    # s variable naming convention: s{edge}{index} → s31 = edge 3, index 1
    s_var_to_edge: dict[str, int] = {}
    for b in bounds:
        m = re.match(r"s(\d)(\d)", b.var)
        if m:
            edge = int(m.group(1))
            s_var_to_edge[b.var] = edge

    # Group s_vars by edge, preserving nesting order from bounds
    edge_s_vars: dict[int, list[str]] = {}
    for b in bounds:
        edge = s_var_to_edge.get(b.var)
        if edge is not None:
            edge_s_vars.setdefault(edge, []).append(b.var)

    counts = {edge: len(vars) for edge, vars in edge_s_vars.items()}
    edge_distribution = EdgeDistribution(counts=counts, s_vars=edge_s_vars)

    # Build pairings: for each propagator, determine which edges it connects.
    # Skip propagators with internal lorentz indices (e.g. nu1 in SE blobs)
    # that don't map to any eikonal vertex.
    pairings = []
    for prop in propagators:
        try:
            edge_a = lorentz_to_edge[prop.lorentz[0]]
            edge_b = lorentz_to_edge[prop.lorentz[1]]
        except KeyError:
            continue

        # Find which s variables on edges a and b belong to this loop
        phase_expr = phases[prop.loop]
        s_a = _find_s_var_for_edge(phase_expr, edge_a, edge_s_vars)
        s_b = _find_s_var_for_edge(phase_expr, edge_b, edge_s_vars)

        pairings.append(Pairing(
            edge_a=edge_a, edge_b=edge_b,
            loop=prop.loop, s_a=s_a, s_b=s_b,
        ))

    # --- Detect dispatch and extract TGV/FGV legs ---
    dispatch = _detect_dispatch(
        orbit_json.get("internal_vertices", []),
        len(orbit_json["vertices"]),
    )

    tgv_legs = []
    fgv_legs = []
    double_tgv_legs = []
    if dispatch in ("tgv", "mixed"):
        tgv_legs = _extract_tgv_legs(orbit_json, edge_s_vars)
    elif dispatch == "fourv":
        fgv_legs = _extract_fgv_legs(orbit_json, edge_s_vars)
    elif dispatch == "double_tgv":
        double_tgv_legs = _extract_double_tgv_legs(orbit_json, edge_s_vars)

    return Orbit(
        orbit_id=orbit_id,
        emission_type=emission_type,
        n_gluons=len(pairings),
        edge_distribution=edge_distribution,
        pairings=pairings,
        bounds=bounds,
        vertices=vertices,
        propagators=propagators,
        source=orbit_json,
        dispatch=dispatch,
        tgv_legs=tgv_legs,
        fgv_legs=fgv_legs,
        double_tgv_legs=double_tgv_legs,
    )


def _detect_dispatch(internal_vertices: list[dict], num_vertices: int) -> str:
    """Classify diagram type from internal_vertices and vertex count."""
    if not internal_vertices:
        return "pure"
    # Count TGV vertices
    tgv_count = sum(1 for iv in internal_vertices if iv.get("type") == "tgv")
    if tgv_count == 1:
        iv = internal_vertices[0]
        n_loops = len(iv.get("loops", []))
        if n_loops == 4:
            # 4e-3g+3g: effective vertex with 4 external legs
            # (two TGV vertices connected by internal gluon)
            return "double_tgv"
        if num_vertices > 3:
            return "mixed"  # TGV + pure gluons
        return "tgv"       # pure TGV (3 legs only)
    if tgv_count == 2:
        return "double_tgv"  # explicit 2 TGV vertices
    iv_type = internal_vertices[0]["type"]
    if iv_type in ("4v", "fgv"):
        return "fourv"
    return "unsupported"


def _extract_tgv_legs(
    orbit_json: dict, edge_s_vars: dict[int, list[str]]
) -> list[TGVLeg]:
    """Extract the 3 TGV legs from propagators connected to the TGV."""
    tgv_loops = set()
    for iv in orbit_json["internal_vertices"]:
        if iv["type"] == "tgv":
            tgv_loops.update(iv["loops"])

    lorentz_to_edge = {v["lorentz"]: v["line"] for v in orbit_json["vertices"]}
    phases = {p["loop"]: p["expression"] for p in orbit_json["phases"]}

    legs = []
    for prop in orbit_json["propagators"]:
        if prop["loop"] not in tgv_loops:
            continue
        mu_index = prop["lorentz"][0]
        edge = lorentz_to_edge[mu_index]
        s_var = _find_s_var_for_edge(phases[prop["loop"]], edge, edge_s_vars)
        legs.append(TGVLeg(edge=edge, s_var=s_var, loop=prop["loop"]))

    if orbit_json.get("emission_type") == "3e-3g":
        assert len(legs) == 3, f"Pure TGV (3e-3g) must have exactly 3 legs, got {len(legs)}"
    return legs


def _extract_fgv_legs(
    orbit_json: dict, edge_s_vars: dict[int, list[str]]
) -> list[FGVLeg]:
    """Extract the 4 FGV legs from propagators connected to the FGV."""
    fgv_loops = set()
    for iv in orbit_json["internal_vertices"]:
        if iv["type"] in ("4v", "fgv"):
            fgv_loops.update(iv["loops"])

    lorentz_to_edge = {v["lorentz"]: v["line"] for v in orbit_json["vertices"]}
    phases = {p["loop"]: p["expression"] for p in orbit_json["phases"]}

    legs = []
    for prop in orbit_json["propagators"]:
        if prop["loop"] not in fgv_loops:
            continue
        mu_index = prop["lorentz"][0]
        edge = lorentz_to_edge[mu_index]
        s_var = _find_s_var_for_edge(phases[prop["loop"]], edge, edge_s_vars)
        legs.append(FGVLeg(edge=edge, s_var=s_var, loop=prop["loop"]))

    if orbit_json.get("emission_type") == "4e-4g":
        assert len(legs) == 4, f"Pure FGV (4e-4g) must have exactly 4 legs, got {len(legs)}"
    return legs


def _validate_double_tgv_topology(raw):
    """Closed four-external-leg contract; the internal connector is implicit.

    Accept one effective four-leg TGV, or two cyclic three-leg loop lists
    sharing one connector name. Explicit connector propagators and vertex
    index annotations are unsupported and rejected, never silently ignored.
    """
    def fail(reason):
        raise ValueError("double_tgv_" + reason)

    def label(value):
        return isinstance(value, str) and bool(value.strip())

    vertices = raw.get("vertices", [])
    props = raw.get("propagators", [])
    ivs = raw.get("internal_vertices", [])
    phases = raw.get("phases", [])
    if len(vertices) != 4 or len(props) != 4:
        fail("requires_four_vertices_and_four_external_propagators")
    if any(not isinstance(v, dict) or set(v) != {"line", "lorentz", "color"}
           or type(v["line"]) is not int or v["line"] not in (1, 2, 3, 4)
           or not label(v["lorentz"]) or not label(v["color"]) for v in vertices):
        fail("unsupported_wilson_vertex_encoding")
    external_mu = {v["lorentz"] for v in vertices}
    external_color = {v["color"] for v in vertices}
    if len(external_mu) != 4 or len(external_color) != 4:
        fail("requires_independent_wilson_indices")
    if any(not isinstance(p, dict) or set(p) != {"loop", "lorentz", "color"}
           or not label(p["loop"])
           or any(not isinstance(p[k], (list, tuple)) or len(p[k]) != 2
                  or not all(label(x) for x in p[k]) for k in ("lorentz", "color"))
           for p in props):
        fail("unsupported_propagator_encoding")
    loops = {p["loop"] for p in props}
    if len(loops) != 4:
        fail("duplicate_external_loop")
    if any(not isinstance(v, dict) or not {"type", "loops"} <= set(v)
           or set(v) - {"type", "id", "loops"} or v["type"] != "tgv"
           or not isinstance(v["loops"], (list, tuple))
           or not all(label(x) for x in v["loops"])
           or len(set(v["loops"])) != len(v["loops"]) for v in ivs):
        fail("unsupported_internal_vertex_encoding")
    if len(ivs) == 1 and len(ivs[0]["loops"]) == 4:
        order = list(ivs[0]["loops"])
    elif len(ivs) == 2 and all(len(v["loops"]) == 3 for v in ivs):
        common = set(ivs[0]["loops"]) & set(ivs[1]["loops"])
        if len(common) != 1:
            fail("requires_one_implicit_connector")
        connector = next(iter(common))
        if connector in loops:
            fail("connector_cannot_be_external")
        order = []
        for v in ivs:
            seq = list(v["loops"]); i = seq.index(connector)
            order.extend((seq[i+1:] + seq[:i+1])[:2])
    else:
        fail("unsupported_internal_vertex_count_or_arity")
    if len(set(order)) != 4 or set(order) != loops:
        fail("vertex_loop_set_mismatch")
    if (len(phases) != 4 or any(not isinstance(p, dict)
            or set(p) != {"loop", "expression"} or not label(p["loop"])
            or not label(p["expression"]) for p in phases)
            or {p["loop"] for p in phases} != loops):
        fail("phase_loop_set_mismatch")
    by_mu = {v["lorentz"]: v for v in vertices}
    used = []; internal_mu = []; internal_color = []
    for p in props:
        ends = [i for i, mu in enumerate(p["lorentz"]) if mu in external_mu]
        if len(ends) != 1:
            fail("leg_requires_one_wilson_endpoint")
        i = ends[0]; mu = p["lorentz"][i]
        if p["color"][i] != by_mu[mu]["color"]:
            fail("color_endpoint_mismatch")
        used.append(mu)
        internal_mu.append(p["lorentz"][1-i])
        internal_color.append(p["color"][1-i])
    if set(used) != external_mu or len(set(used)) != 4:
        fail("wilson_endpoint_bijection_required")
    if (len(set(internal_mu)) != 4 or len(set(internal_color)) != 4
            or set(internal_mu) & external_mu or set(internal_color) & external_color):
        fail("internal_indices_must_be_independent")
    return order


def _extract_double_tgv_legs(orbit_json, edge_s_vars):
    """Bind legs to the oriented vertex loop order, independently of array order."""
    order = _validate_double_tgv_topology(orbit_json)
    vertices = {v["lorentz"]: v for v in orbit_json["vertices"]}
    props = {p["loop"]: p for p in orbit_json["propagators"]}
    phases = {p["loop"]: p["expression"] for p in orbit_json["phases"]}
    if len(phases) != len(orbit_json["phases"]):
        raise ValueError("duplicate_double_tgv_phase")
    legs = []
    for i, loop in enumerate(order):
        if loop not in props or loop not in phases:
            raise ValueError("missing_double_tgv_external_leg")
        ends = [mu for mu in props[loop]["lorentz"] if mu in vertices]
        if len(ends) != 1:
            raise ValueError("double_tgv_leg_requires_one_eikonal_endpoint")
        edge = vertices[ends[0]]["line"]
        candidates = [v for v in edge_s_vars.get(edge, [])
                      if re.search(r"\b" + re.escape(v) + r"\b", phases[loop])]
        if len(candidates) != 1:
            raise ValueError("ambiguous_double_tgv_parameter")
        legs.append(DoubleTGVLeg(edge, candidates[0], loop, "A" if i < 2 else "B"))
    if len({l.s_var for l in legs}) != 4:
        raise ValueError("repeated_double_tgv_endpoint")
    return legs


def _find_s_var_for_edge(
    phase_expr: str, edge: int, edge_s_vars: dict[int, list[str]]
) -> str:
    """Find which s variable on the given edge appears in the phase expression.

    The phase expression contains terms like "(1-s31)*v3" or "s41*v4".
    We match s variables on the given edge against the expression.
    """
    candidates = edge_s_vars.get(edge, [])
    for s_var in candidates:
        if s_var in phase_expr:
            return s_var
    raise ValueError(
        f"No s variable for edge {edge} found in phase: {phase_expr}"
    )


def parse_all_orbits(source_json: dict) -> list[Orbit]:
    """Parse all orbits from orbit_definitions.json format.

    Args:
        source_json: Full orbit_definitions.json dict with "orbits" key.

    Returns:
        List of Orbit objects.
    """
    return [parse_orbit(o) for o in source_json["orbits"]]


def filter_pure_ngluon(orbits: list[Orbit]) -> list[Orbit]:
    """Filter orbits to pure n-gluon types.

    Pure n-gluon types have emission_type matching the pattern
    '{2n}e' where n is the number of gluon propagators, and no
    internal vertices.
    """
    return [
        o for o in orbits
        if o.emission_type in _PURE_NGLUON_TYPES
    ]


# Emission types for pure n-gluon (no internal vertices)
_PURE_NGLUON_TYPES = {"2e", "4e", "6e", "8e", "10e", "12e"}
