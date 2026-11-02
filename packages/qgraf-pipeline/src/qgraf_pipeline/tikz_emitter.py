"""Generate TikZ Feynman diagrams matching reference style.

Reference: wilsonloop-langgraph-agent/docs/two_loop_initial_expressions/
Style: standalone, square contour, corner labels 1-4, attachment dots,
gluon lines, internal vertex dots (3g/4g) or open blobs (SE).
"""
from __future__ import annotations

from typing import List, Tuple

from .classifier import ClassifiedOrbit
from .topology import classify_topology


# Square contour: corners 1=(0,0), 2=(3,0), 3=(3,3), 4=(0,3)
# Edge 1: bottom (0,0)->(3,0), edge 2: right (3,0)->(3,3)
# Edge 3: top (3,3)->(0,3), edge 4: left (0,3)->(0,0)
CORNERS = [(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)]
CENTER = (1.5, 1.5)


def _attachment_xy(edge: int, slot: int, total: int) -> Tuple[float, float]:
    """Position of attachment (edge, slot) along the square contour."""
    t = slot / (total + 1)
    if edge == 1:    # bottom: corner 1 -> corner 2
        return (t * 3.0, 0.0)
    elif edge == 2:  # right: corner 2 -> corner 3
        return (3.0, t * 3.0)
    elif edge == 3:  # top: corner 3 -> corner 4 (right to left)
        return ((1 - t) * 3.0, 3.0)
    else:            # edge 4: left: corner 4 -> corner 1 (top to bottom)
        return (0.0, (1 - t) * 3.0)


def _push_inward(p: Tuple[float, float], dist: float) -> Tuple[float, float]:
    """Push point p toward the square center by dist."""
    dx, dy = CENTER[0] - p[0], CENTER[1] - p[1]
    n = (dx * dx + dy * dy) ** 0.5
    if n < 1e-6:
        return p
    return (p[0] + dist * dx / n, p[1] + dist * dy / n)


def _emit_one_diagram(orbit: ClassifiedOrbit) -> str:
    """Generate TikZ for a single orbit."""
    d = orbit.representative
    info = classify_topology(d)
    ec = info.edge_counts

    lines = []
    lines.append("\\begin{tikzpicture}[scale=0.85, line width=0.5pt]")

    # Wilson contour (4 edges)
    for e in range(1, 5):
        s = CORNERS[(e - 1) % 4]
        p = CORNERS[e % 4]
        lines.append(
            f"  \\draw ({s[0]:.2f},{s[1]:.2f}) -- "
            f"({p[0]:.2f},{p[1]:.2f});")

    # Corner labels
    offsets = [(-0.25, -0.25), (0.25, -0.25), (0.25, 0.25), (-0.25, 0.25)]
    for i in range(4):
        x, y = CORNERS[i]
        ox, oy = offsets[i]
        lines.append(
            f"  \\node[font=\\scriptsize] at "
            f"({x+ox:.2f},{y+oy:.2f}) {{{i + 1}}};")

    # Attachment positions
    pos = {}
    for edge in (1, 2, 3, 4):
        n = ec[edge - 1]
        for slot in range(1, n + 1):
            px, py = _attachment_xy(edge, slot, n)
            pos[(edge, slot)] = (px, py)
            lines.append(f"  \\filldraw ({px:.2f},{py:.2f}) circle (1.2pt);")

    # Separate gluons into WL-WL and complex
    wl_wl_endpoints = set()
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        wl_wl_endpoints.add((edge_a, slot_a))
        wl_wl_endpoints.add((edge_b, slot_b))

    complex_eps = []
    seen = set()
    for _prop_idx, endpoints in info.gluon_connections:
        for ep in endpoints:
            if ep not in wl_wl_endpoints and ep not in seen:
                seen.add(ep)
                complex_eps.append(ep)

    # ---- Draw WL-WL gluon lines ----
    for edge_a, slot_a, edge_b, slot_b in info.wl_wl_pairs:
        a = pos[(edge_a, slot_a)]
        b = pos[(edge_b, slot_b)]
        same_edge = (edge_a == edge_b)
        if same_edge:
            mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            ctrl = _push_inward(mid, 0.65)
            lines.append(
                f"  \\draw ({a[0]:.2f},{a[1]:.2f}) .. controls "
                f"({ctrl[0]:.2f},{ctrl[1]:.2f}) .. "
                f"({b[0]:.2f},{b[1]:.2f});")
        else:
            lines.append(
                f"  \\draw ({a[0]:.2f},{a[1]:.2f}) -- "
                f"({b[0]:.2f},{b[1]:.2f});")

    # ---- Draw complex structure ----
    et = orbit.emission_type
    if complex_eps:
        # Two-internal-vertex case (e.g. 4e-3g+3g):
        # draw each vertex separately with internal gluon line.
        if et == "3e-3g+Pi1":
            from .tgv_self_energy import emit_drawing
            lines.extend(emit_drawing(orbit, info, pos))
        elif (info.internal_connections and len(info.internal_connections) >= 2
                and "3g" in et):
            # Map contour positions to (edge, slot) for the current frame
            ec = info.edge_counts
            _pos_to_xy = {}
            for edge in (1, 2, 3, 4):
                n = ec[edge - 1]
                for slot in range(1, n + 1):
                    # Compute contour position
                    cpos = 0
                    for e in range(1, edge):
                        cpos += ec[e - 1]
                    cpos += slot
                    _pos_to_xy[cpos] = pos[(edge, slot)]

            # Draw each internal vertex
            vertex_centers = []
            for conn_set in info.internal_connections:
                eps_xy = [_pos_to_xy[p] for p in conn_set]
                cx = sum(p[0] for p in eps_xy) / len(eps_xy)
                cy = sum(p[1] for p in eps_xy) / len(eps_xy)
                vx, vy = _push_inward((cx, cy), 0.30)
                vertex_centers.append((vx, vy))
                for ex, ey in eps_xy:
                    lines.append(
                        f"  \\draw ({ex:.2f},{ey:.2f}) -- "
                        f"({vx:.2f},{vy:.2f});")
                lines.append(
                    f"  \\filldraw ({vx:.2f},{vy:.2f}) circle (1.5pt);")

            # Draw internal gluon lines between linked vertices
            for a, b in info.internal_gluon_links:
                if a < len(vertex_centers) and b < len(vertex_centers):
                    ax, ay = vertex_centers[a]
                    bx, by = vertex_centers[b]
                    lines.append(
                        f"  \\draw ({ax:.2f},{ay:.2f}) -- "
                        f"({bx:.2f},{by:.2f});")

        elif et == "2e-Pi1+Pi1":
            # Three visible propagator segments separated by two distinct blobs.
            a, b = pos[complex_eps[0]], pos[complex_eps[1]]
            lines.append(f"  \\draw ({a[0]:.2f},{a[1]:.2f}) -- ({b[0]:.2f},{b[1]:.2f});")
            for fraction in (1/3, 2/3):
                cx = a[0] + fraction*(b[0]-a[0])
                cy = a[1] + fraction*(b[1]-a[1])
                lines.append(f"  \\draw[fill=white] ({cx:.2f},{cy:.2f}) circle (0.12cm);")

        elif "Pi" in et:
            # SE: gluon line with open blob in the middle
            a = pos[complex_eps[0]]
            b = pos[complex_eps[1]]
            mid = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            lines.append(
                f"  \\draw ({a[0]:.2f},{a[1]:.2f}) -- "
                f"({mid[0]:.2f},{mid[1]:.2f});")
            lines.append(
                f"  \\draw ({mid[0]:.2f},{mid[1]:.2f}) -- "
                f"({b[0]:.2f},{b[1]:.2f});")
            lines.append(
                f"  \\draw[fill=white] ({mid[0]:.2f},{mid[1]:.2f}) "
                f"circle (0.12cm);")

        elif "3g" in et:
            # 3g vertex: lines from attachments to central filled dot
            cx = sum(pos[ep][0] for ep in complex_eps) / len(complex_eps)
            cy = sum(pos[ep][1] for ep in complex_eps) / len(complex_eps)
            vx, vy = _push_inward((cx, cy), 0.30)
            for ep in complex_eps:
                p = pos[ep]
                lines.append(
                    f"  \\draw ({p[0]:.2f},{p[1]:.2f}) -- "
                    f"({vx:.2f},{vy:.2f});")
            lines.append(f"  \\filldraw ({vx:.2f},{vy:.2f}) circle (1.5pt);")

        elif "Gamma" in et:
            # Gamma vertex: like 3g but with open blob (loop correction)
            cx = sum(pos[ep][0] for ep in complex_eps) / len(complex_eps)
            cy = sum(pos[ep][1] for ep in complex_eps) / len(complex_eps)
            vx, vy = _push_inward((cx, cy), 0.30)
            for ep in complex_eps:
                p = pos[ep]
                lines.append(
                    f"  \\draw ({p[0]:.2f},{p[1]:.2f}) -- "
                    f"({vx:.2f},{vy:.2f});")
            lines.append(
                f"  \\draw[fill=white] ({vx:.2f},{vy:.2f}) "
                f"circle (0.12cm);")

        elif "4g" in et:
            # 4g vertex: lines from attachments to central filled dot
            cx = sum(pos[ep][0] for ep in complex_eps) / len(complex_eps)
            cy = sum(pos[ep][1] for ep in complex_eps) / len(complex_eps)
            vx, vy = _push_inward((cx, cy), 0.30)
            for ep in complex_eps:
                p = pos[ep]
                lines.append(
                    f"  \\draw ({p[0]:.2f},{p[1]:.2f}) -- "
                    f"({vx:.2f},{vy:.2f});")
            lines.append(f"  \\filldraw ({vx:.2f},{vy:.2f}) circle (1.5pt);")

    # Diagram label
    lines.append(
        f"  \\node[font=\\small] at (1.5,-0.55) {{{orbit.orbit_id}}};")
    lines.append("\\end{tikzpicture}")

    return "\n".join(lines)


def emit_tikz_document(orbits: List[ClassifiedOrbit]) -> str:
    """Generate a complete LaTeX document with one diagram per orbit."""
    lines = [
        "\\documentclass{article}",
        "\\usepackage{tikz}",
        "\\usetikzlibrary{decorations.pathmorphing, "
        "decorations.markings, decorations.pathreplacing}",
        "\\usepackage[margin=0.5cm]{geometry}",
        "\\begin{document}",
        "\\centering",
    ]

    for o in orbits:
        lines.append(_emit_one_diagram(o))
        lines.append("\\vspace{0.5cm}")
        lines.append("\\newpage")
        lines.append("")

    lines.append("\\end{document}")
    return "\n".join(lines)