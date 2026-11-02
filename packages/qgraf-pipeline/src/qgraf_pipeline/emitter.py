"""Render LaTeX momentum-space integral for a classified orbit."""
from __future__ import annotations

from typing import Dict, List, Tuple

from .classifier import ClassifiedOrbit
from .topology import classify_topology
from .internal_types import _name_to_mn, get_vertex_info


def _emit_measure(edge_counts: Tuple[int, int, int, int]) -> str:
    pieces = []
    for edge in (1, 2, 3, 4):
        n = edge_counts[edge - 1]
        if n == 0:
            continue
        if n == 1:
            pieces.append(f"\\int_0^1 ds_{{{edge},1}}")
        else:
            for k in range(n, 0, -1):
                ub = "1" if k == n else f"s_{{{edge},{k+1}}}"
                pieces.append(f"\\int_0^{{{ub}}} ds_{{{edge},{k}}}")
    return " ".join(pieces)


def _emit_eikonal_chain(
    edge_counts: Tuple[int, int, int, int]
) -> Tuple[str, Dict[Tuple[int, int], int]]:
    parts = []
    mu_map = {}
    k = 0
    for edge in (1, 2, 3, 4):
        for slot in range(1, edge_counts[edge - 1] + 1):
            k += 1
            mu_map[(edge, slot)] = k
            parts.append(
                f"(ig\\,v_{{{edge}}}^{{\\mu_{{{k}}}}} t^{{a_{{{k}}}}})")
    return " ".join(reversed(parts)), mu_map


def _emit_coord_diff(
    edge_a: int, slot_a: int, edge_b: int, slot_b: int
) -> str:
    if edge_a == edge_b:
        return (f"(s_{{{edge_a},{slot_a}}} - s_{{{edge_b},{slot_b}}})"
                f"\\,v_{{{edge_a}}}")
    forward = (edge_a - edge_b) % 4
    backward = (edge_b - edge_a) % 4
    if forward <= backward:
        terms = [f"(1-s_{{{edge_b},{slot_b}}})\\,v_{{{edge_b}}}"]
        cur = edge_b
        for _ in range(forward - 1):
            cur = cur % 4 + 1
            terms.append(f"v_{{{cur}}}")
        terms.append(f"s_{{{edge_a},{slot_a}}}\\,v_{{{edge_a}}}")
        return " + ".join(terms)
    else:
        terms = [f"(1-s_{{{edge_a},{slot_a}}})\\,v_{{{edge_a}}}"]
        cur = edge_a
        for _ in range(backward - 1):
            cur = cur % 4 + 1
            terms.append(f"v_{{{cur}}}")
        terms.append(f"s_{{{edge_b},{slot_b}}}\\,v_{{{edge_b}}}")
        body = " + ".join(terms)
        return f"-\\bigl[{body}\\bigr]"


def emit_integral(orbit: ClassifiedOrbit) -> str:
    if orbit.representative is None:
        return f"% [No representative for {orbit.orbit_id}]"

    d = orbit.representative
    info = classify_topology(d)
    ec = info.edge_counts  # use the representative's edge counts, not canonical
    from .self_energy_chain import is_pi1_chain, emit_chain_integral
    if is_pi1_chain(orbit):
        return emit_chain_integral(orbit, info)

    from .tgv_self_energy import is_tgv_pi1, emit_integral as emit_tgv_pi1
    if is_tgv_pi1(orbit):
        return emit_tgv_pi1(orbit, info)

    if orbit.emission_type == '3e-Gamma1':
        from .central_vertex import emit_integral as emit_gamma1
        return emit_gamma1(orbit, info)
    if orbit.emission_type == '4e-3g+3g':
        from .double_vertex import emit_integral as emit_double_vertex
        return emit_double_vertex(orbit,info)

    measure = _emit_measure(ec)
    chain, mu_map = _emit_eikonal_chain(ec)

    lines = [f"I_{{({orbit.orbit_id})}} = {measure}\\;\\frac{{1}}{{N_c}}\\operatorname{{Tr}}\\bigl[{chain}\\bigr]"]

    # Separate gluons into WL-WL and complex.
    # A WL-WL gluon is identified by matching wl_wl_pairs.
    # Complex gluons are the rest (part of internal structure).
    wl_wl_endpoints = set()
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        wl_wl_endpoints.add((edge_a, slot_a))
        wl_wl_endpoints.add((edge_b, slot_b))

    complex_gluons = []
    seen_complex = set()
    for _prop_idx, endpoints in info.gluon_connections:
        for ep in endpoints:
            if ep not in wl_wl_endpoints and ep not in seen_complex:
                seen_complex.add(ep)
                complex_gluons.append(ep)

    # ---- WL-WL gluons: each has independent loop momentum ----
    line_id = 0
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        line_id += 1
        ka = mu_map[(edge_a, slot_a)]
        kb = mu_map[(edge_b, slot_b)]
        dX = _emit_coord_diff(edge_a, slot_a, edge_b, slot_b)
        lines.append(
            f"\\times \\int \\frac{{d^D\\ell_{{{line_id}}}}}{{(2\\pi)^D}}\\,"
            f"\\frac{{-ig_{{\\mu_{{{ka}}}\\mu_{{{kb}}}}}"
            f"\\delta^{{a_{{{ka}}}a_{{{kb}}}}}}}"
            f"{{\\ell_{{{line_id}}}^2+i0}}\\,"
            f"e^{{-i\\ell_{{{line_id}}}\\cdot\\bigl({dX}\\bigr)}}"
        )

    # ---- Complex gluons: data-driven emission from internal type registry ----
    if complex_gluons:
        et = orbit.emission_type
        # For single internal structure, determine (M, N) and emit
        if orbit.internal:
            # Get the first (and typically only) internal structure type
            itype = list(orbit.internal.keys())[0]
            M, N = _name_to_mn(itype)
            if M == 0 and N == 0:
                # Fallback: try to infer from emission type name
                pass
            vi = get_vertex_info(M, N)

            # Apply overall sign
            if vi["sign"] == -1:
                lines[0] = lines[0].replace(" = ", " = -")

            line_id = 0
            if vi["shared"]:
                # SE-type: all legs share the same loop momentum.
                # The N legs are actually two halves of N/2 gluon lines
                # that go through the effective vertex.
                # For N=2 (Pi): one gluon line with the blob in the middle.
                if N == 2:
                    (edge_a, slot_a) = complex_gluons[0]
                    (edge_b, slot_b) = complex_gluons[1]
                    ka = mu_map[(edge_a, slot_a)]
                    kb = mu_map[(edge_b, slot_b)]
                    dX = _emit_coord_diff(edge_a, slot_a, edge_b, slot_b)
                    lines.append(
                        "\\times \\int \\frac{d^D\\ell}{(2\\pi)^D}\\,"
                        "\\frac{-ig_{\\mu_{" + str(ka) + "}\\nu_1}\\delta^{a_{" + str(ka) + "}b_1}}{\\ell^2+i0}\\,"
                        + vi["latex"] + "\\,"
                        "\\frac{-ig_{\\nu_2\\mu_{" + str(kb) + "}}\\delta^{b_2a_{" + str(kb) + "}}}{\\ell^2+i0}\\,"
                        "e^{-i\\ell\\cdot\\bigl(" + dX + "\\bigr)}"
                    )
                    lines.append(
                        "\\quad\\text{(SE: ghost+quark+gluon-bubble+4g-tadpole)}"
                    )
                else:
                    # General shared case (N>2): emit as two halves
                    mid = len(complex_gluons) // 2
                    for ep in complex_gluons[:mid]:
                        line_id += 1
                        edge, slot = ep
                        k = mu_map[(edge, slot)]
                        dX = _emit_coord_to_origin(edge, slot)
                        lines.append(
                            "\\times \\int \\frac{d^D\\ell}{(2\\pi)^D}\\," +
                            "\\frac{-ig_{\\mu_{" + str(k) + "}\\nu_{" + str(line_id) + "}}}{\\ell^2+i0}\\," +
                            "e^{-i\\ell\\cdot " + dX + "}"
                        )
                    lines.append("\\times " + vi["latex"])
                    for ep in complex_gluons[mid:]:
                        line_id += 1
                        edge, slot = ep
                        k = mu_map[(edge, slot)]
                        dX = _emit_coord_to_origin(edge, slot)
                        lines.append(
                            "\\times " +
                            "\\frac{-ig_{\\nu_{" + str(line_id) + "}\\mu_{" + str(k) + "}}}{\\ell^2+i0}\\," +
                            "e^{-i\\ell\\cdot " + dX + "}"
                        )
            else:
                # Each leg has independent loop momentum
                # Emit vertex first, then each leg
                lines.append("\\times " + vi["latex"])
                lines.append(r'\times (2\pi)^D\delta^{(D)}('+ '+'.join(r'\ell_{'+str(j)+'}' for j in range(1,len(complex_gluons)+1))+')')
                for i, ep in enumerate(complex_gluons, 1):
                    edge, slot = ep
                    k = mu_map[(edge, slot)]
                    dX = _emit_coord_to_origin(edge, slot)
                    lines.append(
                        "\\times \\int \\frac{d^D\\ell_{" + str(i) + "}}{(2\\pi)^D}\\," +
                        "\\frac{-ig_{\\mu_{" + str(k) + "}\\nu_{" + str(i) + "}}" +
                        "\\delta^{a_{" + str(k) + "}b_{" + str(i) + "}}}" +
                        "{\\ell_{" + str(i) + "}^2+i0}\\," +
                        "e^{-i\\ell_{" + str(i) + "}\\cdot " + dX + "}"
                    )

    return "\n".join(lines)


def _emit_coord_to_origin(edge: int, slot: int) -> str:
    """X_i(s) = corner_{i-1} + s v_i, expanded in v-basis."""
    if edge == 1:
        return "s_{" + str(edge) + "," + str(slot) + "}\\,v_1"
    elif edge == 2:
        return "v_1 + s_{" + str(edge) + "," + str(slot) + "}\\,v_2"
    elif edge == 3:
        return "v_1 + v_2 + s_{" + str(edge) + "," + str(slot) + "}\\,v_3"
    else:
        return "v_1 + v_2 + v_3 + s_{" + str(edge) + "," + str(slot) + "}\\,v_4"


def emit_md_document(orbits: List[ClassifiedOrbit], n: int) -> str:
    lines = [
        f"# {n}-Loop Wilson Loop — Formal Integrals\n",
        "Auto-generated by `qgraf-pipeline`.\n",
        "## Conventions\n",
        "- **Wilson loop**: Square contour, edges $z_1, z_2, z_3, z_4$.",
        "- **Eikonal vertex**: $(i g v_i^{\\mu_k} t^{a_k})$ on edge $i$.",
        "- **Gluon propagator**: $\\frac{-i g_{\\mu\\nu} \\delta^{ab}}{\\ell^2 + i0}$",
        "- **Dimension**: $D = 4 - 2\\varepsilon$, Feynman gauge.\n",
        "---\n",
    ]

    for o in orbits:
        lines.append(f"## {o.orbit_id}\n")
        lines.append(f"- **Type**: `{o.emission_type}`")
        lines.append(f"- **Category**: `{o.category}`")
        lines.append(f"- **Edge counts**: `{list(o.edge_counts)}`")
        lines.append(f"- **Eikonal count**: {o.eikonal_count}")
        lines.append(f"- **Diagrams**: {o.n_diagrams}")
        lines.append(f"- **Mirror**: {o.mirror_relation}")
        if o.mirror_of:
            lines.append(f"  (mirror of `{o.mirror_of}`)")
        lines.append(f"- **2PI**: {o.is_2pi}")
        if o.reducible_to:
            lines.append(f"  (reducible to: {', '.join(o.reducible_to)})")
        lines.append("")
        latex = emit_integral(o)
        lines.append(f"$$\n{latex}\n$$\n")
        lines.append("---\n")

    return "\n".join(lines)
