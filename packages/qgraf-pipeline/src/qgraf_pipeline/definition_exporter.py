"""Export structured orbit definitions for hepeflow-preprocessing consumption."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, List, Tuple

from .sign_rules import attach_contract
from .sign_rules import attach_contract
from .classifier import ClassifiedOrbit
from .topology import classify_topology, TopologyInfo


def _coord_to_origin_expr(edge: int, slot: int) -> str:
    """X vector from origin to contour position (edge, slot)."""
    if edge == 1:
        return f"s{edge}{slot}*v1"
    elif edge == 2:
        return f"v1 + s{edge}{slot}*v2"
    elif edge == 3:
        return f"v1 + v2 + s{edge}{slot}*v3"
    else:
        return f"v1 + v2 + v3 + s{edge}{slot}*v4"


def _coord_diff_expr(edge_a: int, slot_a: int,
                     edge_b: int, slot_b: int) -> str:
    """X vector difference between two contour positions."""
    if edge_a == edge_b:
        return f"(s{edge_a}{slot_a} - s{edge_b}{slot_b})*v{edge_a}"
    forward = (edge_a - edge_b) % 4
    backward = (edge_b - edge_a) % 4
    if forward <= backward:
        terms = [f"(1-s{edge_b}{slot_b})*v{edge_b}"]
        cur = edge_b
        for _ in range(forward - 1):
            cur = cur % 4 + 1
            terms.append(f"v{cur}")
        terms.append(f"s{edge_a}{slot_a}*v{edge_a}")
        return " + ".join(terms)
    else:
        terms = [f"(1-s{edge_a}{slot_a})*v{edge_a}"]
        cur = edge_a
        for _ in range(backward - 1):
            cur = cur % 4 + 1
            terms.append(f"v{cur}")
        terms.append(f"s{edge_b}{slot_b}*v{edge_b}")
        return f"-({(' + '.join(terms))})"


def _negate_phase(expr: str) -> str:
    """Negate a phase expression: Φ → -Φ.

    The phase convention is e^{iΦ·ℓ}.  The standard Feynman rule is
    e^{-iℓ·X}.  For these to agree we need Φ = -X, so the phase
    expression must be the negative of the coordinate (difference).
    """
    expr = expr.strip()
    if expr.startswith('-(') and expr.endswith(')'):
        # -(X) → X
        return expr[2:-1].strip()
    elif expr.startswith('-'):
        # -X → X
        return expr[1:].strip()
    else:
        # X → -(X)
        return f"-({expr})"


def _extract_vars(expr: str) -> List[str]:
    """Extract unique s-variable names from a phase expression."""
    return sorted(set(re.findall(r's\d+', expr)))


def _make_bound(var: str, edge: int, slot: int,
                ec: Tuple[int, int, int, int]) -> dict:
    """Determine bound range for a Schwinger parameter variable."""
    m = re.match(r's(\d)(\d+)', var)
    if not m:
        return {"var": var, "low": "0", "high": "1"}
    e, s = int(m.group(1)), int(m.group(2))
    n = ec[e - 1]
    if s == n:
        return {"var": var, "low": "0", "high": "1"}
    else:
        return {"var": var, "low": "0", "high": f"s{e}{s+1}"}


def _nloops_from_orbit(orbit: ClassifiedOrbit) -> int:
    """Derive nloops from g_power = 2 * nloops."""
    return orbit.g_power // 2


def _build_orbit_definition(orbit: ClassifiedOrbit) -> dict | None:
    """Extract structured definition for one orbit's representative diagram."""
    if orbit.representative is None:
        return None

    d = orbit.representative
    info: TopologyInfo = classify_topology(d)
    ec = info.edge_counts
    from .self_energy_chain import is_pi1_chain, build_chain_definition
    if is_pi1_chain(orbit):
        return attach_contract(build_chain_definition(orbit, info))

    from .tgv_self_energy import is_tgv_pi1, build_definition
    if is_tgv_pi1(orbit):
        return attach_contract(build_definition(orbit, info))

    if orbit.emission_type == '3e-Gamma1':
        from .central_vertex import build_definition
        return attach_contract(build_definition(orbit, info))

    mu_map: Dict[Tuple[int, int], int] = {}
    k = 0
    for edge in (1, 2, 3, 4):
        for slot in range(1, ec[edge - 1] + 1):
            k += 1
            mu_map[(edge, slot)] = k

    vertices = []
    for edge in (1, 2, 3, 4):
        for slot in range(1, ec[edge - 1] + 1):
            k = mu_map[(edge, slot)]
            vertices.append({
                "line": edge,
                "lorentz": f"mu{k}",
                "color": f"a{k}",
            })

    wl_wl_set = set()
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        wl_wl_set.add((edge_a, slot_a))
        wl_wl_set.add((edge_b, slot_b))

    propagators = []
    phases = []
    internal_vertices = []
    bounds = []
    seen_vars = set()
    loop_idx = 0

    # WL-WL gluons
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        loop_idx += 1
        loop_name = f"l{loop_idx}"
        ka = mu_map[(edge_a, slot_a)]
        kb = mu_map[(edge_b, slot_b)]
        propagators.append({
            "loop": loop_name,
            "lorentz": [f"mu{ka}", f"mu{kb}"],
            "color": [f"a{ka}", f"a{kb}"],
        })
        dX = _coord_diff_expr(edge_a, slot_a, edge_b, slot_b)
        phases.append({"loop": loop_name, "expression": _negate_phase(dX)})
        for var in _extract_vars(dX):
            if var not in seen_vars:
                seen_vars.add(var)
                bounds.append(_make_bound(var, edge_a, slot_a, ec))

    # Internal-vertex gluons
    complex_gluons = []
    seen_complex = set()
    for _prop_idx, endpoints in info.gluon_connections:
        for ep in endpoints:
            if ep not in wl_wl_set and ep not in seen_complex:
                seen_complex.add(ep)
                complex_gluons.append(ep)

    if complex_gluons:
        itype = list(orbit.internal.keys())[0]
        if itype == "3g":
            tgv_loops = []
            for ep in complex_gluons:
                loop_idx += 1
                loop_name = f"l{loop_idx}"
                tgv_loops.append(loop_name)
                edge, slot = ep
                k = mu_map[(edge, slot)]
                propagators.append({
                    "loop": loop_name,
                    "lorentz": [f"mu{k}", f"nu{loop_idx}"],
                    "color": [f"a{k}", f"b{loop_idx}"],
                })
                dX = _coord_to_origin_expr(edge, slot)
                phases.append({"loop": loop_name, "expression": _negate_phase(dX)})
                for var in _extract_vars(dX):
                    if var not in seen_vars:
                        seen_vars.add(var)
                        bounds.append(_make_bound(var, edge, slot, ec))
            internal_vertices.append({
                "type": "tgv", "id": "v1", "loops": tgv_loops,
            })
        elif itype == "4g":
            fgv_loops = []
            for ep in complex_gluons:
                loop_idx += 1
                loop_name = f"l{loop_idx}"
                fgv_loops.append(loop_name)
                edge, slot = ep
                k = mu_map[(edge, slot)]
                propagators.append({
                    "loop": loop_name,
                    "lorentz": [f"mu{k}", f"nu{loop_idx}"],
                    "color": [f"a{k}", f"b{loop_idx}"],
                })
                dX = _coord_to_origin_expr(edge, slot)
                phases.append({"loop": loop_name, "expression": _negate_phase(dX)})
                for var in _extract_vars(dX):
                    if var not in seen_vars:
                        seen_vars.add(var)
                        bounds.append(_make_bound(var, edge, slot, ec))
            internal_vertices.append({
                "type": "4v", "id": "v1", "loops": fgv_loops,
            })
        elif itype.startswith("Pi"):
            (edge_a, slot_a) = complex_gluons[0]
            (edge_b, slot_b) = complex_gluons[1]
            loop_idx += 1
            loop_name = f"l{loop_idx}"
            ka = mu_map[(edge_a, slot_a)]
            kb = mu_map[(edge_b, slot_b)]
            propagators.append({
                "loop": loop_name,
                "lorentz": [f"mu{ka}", f"nu1"],
                "color": [f"a{ka}", f"b1"],
            })
            propagators.append({
                "loop": loop_name,
                "lorentz": [f"nu2", f"mu{kb}"],
                "color": [f"b2", f"a{kb}"],
            })
            dX = _coord_diff_expr(edge_a, slot_a, edge_b, slot_b)
            phases.append({"loop": loop_name, "expression": _negate_phase(dX)})
            for var in _extract_vars(dX):
                if var not in seen_vars:
                    seen_vars.add(var)
                    bounds.append(_make_bound(var, edge_a, slot_a, ec))
            m_val = int(itype[2:]) if len(itype) > 2 else 1
            se_type = "vacuum-pol" if m_val == 1 else "vacuum-pol-2loop"
            internal_vertices.append({
                "type": "se", "propagator": loop_name, "se_type": se_type,
            })

    return attach_contract({
        "orbit_id": orbit.orbit_id,
        "emission_type": orbit.emission_type,
        "sector": "dg",
        "nloops": _nloops_from_orbit(orbit),
        "vertices": vertices,
        "propagators": propagators,
        "internal_vertices": internal_vertices,
        "phases": phases,
        "bounds": bounds,
    })


def export_orbit_definitions(
    orbits: List[ClassifiedOrbit],
    n: int,
    output_dir: Path,
) -> Path:
    """Export orbit_definitions.json for n-loop."""
    definitions = []
    for o in orbits:
        d = _build_orbit_definition(o)
        if d is not None:
            definitions.append(d)

    data = {
        "n_loops": n,
        "n_orbits": len(definitions),
        "orbits": definitions,
    }

    loop_dir = output_dir / f"{n}loop"
    loop_dir.mkdir(parents=True, exist_ok=True)
    out_path = loop_dir / "orbit_definitions.json"
    out_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return out_path
