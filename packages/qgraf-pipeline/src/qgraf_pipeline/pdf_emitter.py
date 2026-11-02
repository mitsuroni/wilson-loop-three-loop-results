"""Generate per-emission-type PDFs with classification and mirror annotations."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List

from .classifier import ClassifiedOrbit
from .topology import classify_topology
from .tikz_emitter import _emit_one_diagram, _attachment_xy, _push_inward, CORNERS, CENTER


def _emit_type_document(
    orbits: List[ClassifiedOrbit],
    type_name: str,
    n: int,
) -> str:
    """Generate a complete LaTeX document for one emission type.

    Orbits are grouped by canonical edge distribution, with mirror
    annotations on each diagram.
    """
    # Group by edge distribution
    by_edges: Dict[tuple, List[ClassifiedOrbit]] = {}
    for o in orbits:
        by_edges.setdefault(o.edge_counts, []).append(o)

    lines = [
        r"\documentclass[a4paper,12pt]{article}",
        r"\usepackage{tikz}",
        r"\usepackage[margin=0.8cm]{geometry}",
        r"\usepackage{hyperref}",
        r"\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}",
        r"\usetikzlibrary{decorations.pathmorphing,decorations.markings}",
        r"\begin{document}",
        r"\centering",
        r"",
        r"\title{" + type_name + r" -- " + str(len(orbits)) + r" Orbits}",
        r"\author{qgraf-pipeline (n=" + str(n) + r"-loop)}",
        r"\maketitle",
        r"",
        # Summary table
        r"\section*{Summary}",
        r"\begin{tabular}{lll}",
        r"\hline",
        r"Edge Distribution & Orbits & Mirror Structure \\",
        r"\hline",
    ]

    for ec in sorted(by_edges.keys()):
        ec_orbits = by_edges[ec]
        self_count = sum(1 for o in ec_orbits if o.mirror_relation == "self-mirror")
        pair_count = len(ec_orbits) - self_count
        mirror_desc = f"{self_count} self-mirror"
        if type_name == "3e-3g+Pi1":
            local_ids = {o.orbit_id for o in ec_orbits}
            within = sum(o.mirror_of in local_ids for o in ec_orbits) // 2
            outside = sum(o.mirror_of is not None and o.mirror_of not in local_ids for o in ec_orbits)
            if within:
                mirror_desc += f", {within} mirror pairs"
            if outside:
                mirror_desc += f", {outside} cross-group partners"
        elif pair_count > 0:
            mirror_desc += f", {pair_count // 2} mirror pairs"
        lines.append(
            rf"${list(ec)}$ & {len(ec_orbits)} & {mirror_desc} \\"
        )

    lines.extend([
        r"\hline",
        r"\end{tabular}",
        r"\newpage",
    ])

    # TOC-style quick links
    lines.append(r"\section*{Edge Distributions}")
    lines.append(r"\begin{description}")
    for ec in sorted(by_edges.keys()):
        label = f"ec-{ec[0]}{ec[1]}{ec[2]}{ec[3]}"
        lines.append(
            rf"  \item[{list(ec)}] "
            rf"\hyperref[{{{label}}}]{{{len(by_edges[ec])} orbits}}"
        )
    lines.append(r"\end{description}")
    lines.append(r"\newpage")

    # Detailed sections per edge distribution
    for ec in sorted(by_edges.keys()):
        ec_orbits = by_edges[ec]
        label = f"ec-{ec[0]}{ec[1]}{ec[2]}{ec[3]}"
        lines.append(rf"\section{{Edge Distribution ${list(ec)}$}}")
        lines.append(rf"\label{{{label}}}")
        lines.append(rf"{len(ec_orbits)} orbits")
        lines.append(r"\vspace{0.5cm}")
        lines.append("")

        # Sort: self-mirror first, then mirror pairs together
        self_mirrors = [o for o in ec_orbits if o.mirror_relation == "self-mirror"]
        pair_orbits = [o for o in ec_orbits if o.mirror_relation != "self-mirror"]

        # Group mirror pairs
        pair_groups: Dict[str, List[ClassifiedOrbit]] = {}
        for o in pair_orbits:
            pair_groups.setdefault(o.mirror_relation, []).append(o)

        orbit_idx = 0
        for o in self_mirrors:
            orbit_idx += 1
            lines.append(rf"\subsection*{{{orbit_idx}. {o.orbit_id}}}")
            lines.append(r"{\small \textbf{Mirror:} self-mirror}")
            lines.append(r"\vspace{0.2cm}")
            lines.append("")
            lines.append(_emit_one_diagram(o))
            lines.append(r"\vspace{0.3cm}")
            lines.append("")

        for pair_id, pair_list in sorted(pair_groups.items()):
            orbit_idx += 1
            if len(pair_list) == 2:
                a, b = pair_list[0], pair_list[1]
                lines.append(
                    rf"\subsection*{{{orbit_idx}. {a.orbit_id} $\leftrightarrow$ {b.orbit_id}}}"
                )
                lines.append(rf"{{\small \textbf{{Mirror:}} {pair_id}}}")
                lines.append(r"\vspace{0.2cm}")
                lines.append("")
                lines.append(r"\begin{minipage}[t]{0.48\textwidth}")
                lines.append(r"\centering")
                lines.append(_emit_one_diagram(a))
                lines.append(rf"\vspace{{0.1cm}}")
                lines.append(rf"{{\footnotesize {a.orbit_id}}}")
                lines.append(r"\end{minipage}")
                lines.append(r"\hfill")
                lines.append(r"\begin{minipage}[t]{0.48\textwidth}")
                lines.append(r"\centering")
                lines.append(_emit_one_diagram(b))
                lines.append(rf"\vspace{{0.1cm}}")
                lines.append(rf"{{\footnotesize {b.orbit_id}}}")
                lines.append(r"\end{minipage}")
            else:
                # Mirror partner in a different edge distribution
                o = pair_list[0]
                lines.append(rf"\subsection*{{{orbit_idx}. {o.orbit_id}}}")
                lines.append(rf"{{\small \textbf{{Mirror:}} {pair_id} (partner: {o.mirror_of})}}")
                lines.append(r"\vspace{0.2cm}")
                lines.append("")
                lines.append(_emit_one_diagram(o))
            lines.append(r"\vspace{0.3cm}")
            lines.append("")

        lines.append(r"\newpage")

    lines.append(r"\end{document}")
    return "\n".join(lines)


def emit_per_type_pdfs(orbits: List[ClassifiedOrbit], n: int, out_dir: Path) -> None:
    """Generate per-emission-type PDFs with classification and mirror labels.

    Output: out_dir/pdfs/{type_name}.pdf
    """
    # Group orbits by emission type
    by_type: Dict[str, List[ClassifiedOrbit]] = {}
    for o in orbits:
        by_type.setdefault(o.emission_type, []).append(o)

    pdf_dir = out_dir / "pdfs"
    pdf_dir.mkdir(parents=True, exist_ok=True)

    for type_name in sorted(by_type.keys()):
        type_orbits = by_type[type_name]
        tex_content = _emit_type_document(type_orbits, type_name, n)

        tex_path = pdf_dir / f"{type_name}.tex"
        tex_path.write_text(tex_content, encoding="utf-8")

        # Compile with pdflatex (twice for TOC/cross-refs)
        for _ in range(2):
            subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", tex_path.name],
                capture_output=True,
                timeout=30,
                cwd=str(pdf_dir.resolve()),
            )

        # Clean up aux/log files
        for ext in [".aux", ".log", ".out"]:
            aux = pdf_dir / f"{type_name}{ext}"
            if aux.exists():
                aux.unlink()

        pdf_path = pdf_dir / f"{type_name}.pdf"
        if pdf_path.exists():
            print(f"  {type_name}: {pdf_path} ({len(type_orbits)} orbits)")
        else:
            print(f"  {type_name}: FAILED (check {tex_path})")