"""CLI entry point: qgraf-pipeline <n>."""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from .parser import parse_file  # keep for --out-file backward compat
from .internal_types import enumerate_emission_types
from .classifier import classify_orbits
from .emitter import emit_md_document
from .json_emitter import write_json
from .tikz_emitter import emit_tikz_document
from .pdf_emitter import emit_per_type_pdfs
from .generator import generate_diagrams
from .definition_exporter import export_orbit_definitions


def main() -> int:
    parser = argparse.ArgumentParser(
        description="n-loop Wilson loop diagram pipeline")
    parser.add_argument("n", type=int, help="Number of loops")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory")
    parser.add_argument("--out-file", type=str, default=None,
                        help="Use existing QGRAF .out file (skip generator)")
    parser.add_argument("--audit-recursion", action="store_true",
                        help="Export pre-2PI recursion audit only (loops 1-3; separate filenames)")
    parser.add_argument("--no-pdf", action="store_true",
                        help="Skip optional pdflatex rendering; preserve all mathematical exports")
    args = parser.parse_args()

    if args.audit_recursion:
        if args.out_file or args.n not in (1, 2, 3):
            parser.error('--audit-recursion requires the generator and loop order 1, 2 or 3')
        from .recursion_audit import write_audit
        audit_dir = Path(args.output_dir) if args.output_dir else Path(f'outputs/{args.n}loop/audit')
        audit = write_audit(args.n, audit_dir)
        print(f"Recursion audit: {len(audit['orbits'])} orbits -> {audit_dir}")
        return 0

    n = args.n
    t0 = time.time()

    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        out_dir = Path(f"outputs/{n}loop")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== {n}-loop Wilson Loop Pipeline ===\n")

    print("Step 1: Deriving emission types from g-power table...")
    emission_types = enumerate_emission_types(n)
    print(f"  Found {len(emission_types)} emission types:")
    for et in emission_types:
        print(f"    {et.name}  (E={et.eikonal_count}, "
              f"wl_wl={et.wl_wl_pairs}, cat={et.category})")

    if args.out_file:
        out_path = Path(args.out_file)
        print(f"\nStep 2: Using existing .out file: {out_path}")
        diagrams = parse_file(out_path)
        print(f"  Parsed {len(diagrams)} diagrams")
    else:
        print(f"\nStep 2: Generating diagrams combinatorially...")
        diagrams = generate_diagrams(n, emission_types)
        print(f"  Generated {len(diagrams)} diagrams")

    print(f"\nStep 3: Classifying diagrams into orbits...")
    orbits = classify_orbits(diagrams, emission_types)
    print(f"  {len(orbits)} orbits:")
    for o in orbits:
        print(f"    {o.orbit_id}  edges={o.edge_counts}  "
              f"n={o.n_diagrams}  mirror={o.mirror_relation}  "
              f"2pi={o.is_2pi}")

    print(f"\nStep 4: Emitting LaTeX...")
    md_path = out_dir / "integrals.md"
    md_content = emit_md_document(orbits, n)
    md_path.write_text(md_content, encoding="utf-8")
    print(f"  MD: {md_path} ({len(md_content):,} chars)")

    print(f"\nStep 5: Writing JSON...")
    json_path = out_dir / "integrals.json"
    write_json(orbits, n, json_path)
    print(f"  JSON: {json_path}")

    print(f"\nStep 5b: Exporting orbit definitions for hepeflow...")
    def_path = export_orbit_definitions(orbits, n, out_dir.parent)
    print(f"  orbit_definitions.json -> {def_path}")

    print(f"\nStep 7: Generating TikZ...")
    tex_path = out_dir / "diagrams.tex"
    tex_content = emit_tikz_document(orbits)
    tex_path.write_text(tex_content, encoding="utf-8")
    print(f"  TeX: {tex_path}")

    if args.no_pdf:
        print("\nStep 8: PDF rendering skipped (--no-pdf).")
    else:
        print(f"\nStep 8: Generating per-type PDFs...")
        emit_per_type_pdfs(orbits, n, out_dir)

    elapsed = time.time() - t0
    print(f"\n=== Done in {elapsed:.1f}s ===")
    print(f"Output: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
