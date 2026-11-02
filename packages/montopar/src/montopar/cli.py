"""montopar CLI entry point.

Usage:
    montopar --source orbit_definitions.json --orbit-id 4e-01
    montopar --source orbit_definitions.json --all
    montopar --source orbit_definitions.json --type 6e
    montopar --source orbit_definitions.json --type 3e-3g
    montopar --source orbit_definitions.json --type 4e-4g
"""
import argparse
import json
import sys
from pathlib import Path

from montopar.parser import parse_all_orbits
from montopar.hypercube import hypercube_transform
from montopar.color import compute_color
from montopar.c_factor import compute_c, compute_c_tgv, compute_c_fgv, compute_c_mixed, compute_c_double_tgv
from montopar.j_factor import compute_j, compute_j_tgv, compute_j_fgv, compute_j_mixed
from montopar.j_factor import compute_f_p, compute_lorentz_factor
from montopar.fgv import is_physical_zero
from montopar.double_tgv import compute_j_double_tgv
from montopar.models import CFactors, FGVJFactors, DoubleTGVJFactors, Orbit, HypercubeResult
from montopar.rotation import (
    rotate_orbit_full, classify_structure, find_canonical_rotation,
)
from montopar.tgv import compute_S, compute_Pnum
from montopar.emitter import (
    emit_montopar_json,
    write_json,
    _build_jraw_s_tgv,
    _build_jraw_s_mixed,
)

# Canonical bases for TGV normalization (used by 5e-3g mixed dispatch).
# Keys are structure types from classify_structure().
_CANONICAL_BASES: dict[str, tuple[int, int, int]] = {
    "degenerate_bc_adjacent": (3, 4, 4),
    "degenerate_ab_adjacent": (3, 3, 4),
    "degenerate_bc_opposite": (1, 3, 3),
    "nondegenerate_consecutive": (1, 2, 3),
}


def _find_repo_root() -> Path:
    """Find the wilsonloop repository root."""
    current = Path(__file__).resolve().parent
    for _ in range(10):
        if (current / ".git").exists():
            return current
        current = current.parent
    return Path(__file__).resolve().parent.parent


def _reorder_tgv_legs(orbit: Orbit) -> Orbit:
    """Reorder TGV legs to canonical order.

    Only reorders when the leg labeling hides the true physical type:
    - degenerate_ab_opposite → degenerate_bc_opposite (same edge set)
    - degenerate_ac_opposite → degenerate_bc_opposite
    - nondegenerate_mixed → nondegenerate_consecutive (same edge set)

    Does NOT reorder degenerate_ab_adjacent → degenerate_bc_adjacent
    (these are genuinely different C4 orbits).
    """
    if not orbit.tgv_legs or len(orbit.tgv_legs) != 3:
        return orbit

    from montopar.rotation import cyclic_dist

    a = orbit.tgv_legs[0].edge
    b = orbit.tgv_legs[1].edge
    c = orbit.tgv_legs[2].edge
    legs = list(orbit.tgv_legs)

    # Degenerate: two legs on same edge, unique leg on opposite edge
    if a == b and a != c and cyclic_dist(a, c) == 2:
        # degenerate_ab_opposite → put unique (c) first
        legs = [orbit.tgv_legs[2], orbit.tgv_legs[0], orbit.tgv_legs[1]]
    elif a == c and a != b and cyclic_dist(a, b) == 2:
        # degenerate_ac_opposite → put unique (b) first
        legs = [orbit.tgv_legs[1], orbit.tgv_legs[0], orbit.tgv_legs[2]]

    # Nondegenerate: reorder along contour if not consecutive
    elif a != b and b != c and a != c:
        if cyclic_dist(a, b) != 1 or cyclic_dist(b, c) != 1:
            # nondegenerate_mixed → find the cyclic shift that makes consecutive
            for shift in range(1, 3):
                rotated = legs[shift:] + legs[:shift]
                ra, rb, rc = rotated[0].edge, rotated[1].edge, rotated[2].edge
                if cyclic_dist(ra, rb) == 1 and cyclic_dist(rb, rc) == 1:
                    legs = rotated
                    break

    if legs == list(orbit.tgv_legs):
        return orbit

    return Orbit(
        orbit_id=orbit.orbit_id,
        emission_type=orbit.emission_type,
        n_gluons=orbit.n_gluons,
        edge_distribution=orbit.edge_distribution,
        pairings=orbit.pairings,
        bounds=orbit.bounds,
        vertices=orbit.vertices,
        propagators=orbit.propagators,
        source=orbit.source,
        dispatch=orbit.dispatch,
        tgv_legs=legs,
        fgv_legs=orbit.fgv_legs,
    )


def _normalize_tgv_to_canonical(orbit: Orbit) -> Orbit:
    """Rotate orbit so that its TGV legs are on the canonical base.

    First reorders TGV legs to canonical order (unique leg first for
    degenerate, along contour for nondegenerate), then applies C4
    rotation to the canonical base.
    """
    if not orbit.tgv_legs:
        return orbit

    # Step 1: reorder legs to canonical order
    orbit = _reorder_tgv_legs(orbit)

    a, b, c = orbit.tgv_legs[0].edge, orbit.tgv_legs[1].edge, orbit.tgv_legs[2].edge
    struct = classify_structure(a, b, c)

    if struct == "triple":
        return orbit

    canonical = _CANONICAL_BASES.get(struct)
    if canonical is None:
        return orbit

    # Step 2: find C4 rotation to canonical base
    for k in range(4):
        from montopar.rotation import rotate_edge
        if tuple(rotate_edge(e, k) for e in (a, b, c)) == canonical:
            return rotate_orbit_full(orbit, k)

    return orbit


def _compute_one(
    orbit: Orbit,
    out_dir: Path,
    suffix: str = "",
    tgv_ibp_enabled: bool = False,
    double_tgv_forms_enabled: bool = False,
    double_tgv_operator_enabled: bool = False,
) -> None:
    """Compute and emit C × J; TGV/mixed use one canonical representative."""
    if double_tgv_forms_enabled and orbit.dispatch != "double_tgv":
        raise ValueError("double_tgv_forms_requires_double_tgv_source")
    if double_tgv_operator_enabled and orbit.dispatch != "double_tgv":
        raise ValueError("double_tgv_operator_requires_double_tgv_source")
    from montopar.sign_rules import validate_contract
    validate_contract(orbit.source,required=True)
    if orbit.dispatch == "gamma1":
        from montopar.gamma1 import build_result
        result = build_result(orbit, ibp_enabled=tgv_ibp_enabled)
        directory = out_dir / f"{orbit.source['nloops']}loop" / orbit.emission_type / orbit.orbit_id
        write_json(result, directory / (f"montopar{suffix}.json" if suffix else "montopar.json"))
        return
    if orbit.dispatch == "tgv_pi1":
        from montopar.tgv_pi1 import build_result
        result = build_result(orbit, ibp_enabled=tgv_ibp_enabled)
        directory = out_dir / f"{orbit.source['nloops']}loop" / orbit.emission_type / orbit.orbit_id
        write_json(result, directory / (f"montopar{suffix}.json" if suffix else "montopar.json"))
        return
    if orbit.dispatch in ("self_energy", "self_energy_chain", "self_energy_pi2"):
        if tgv_ibp_enabled:
            raise ValueError("Self energy currently exports complete ordered raw C/J; no contour IBP")
        if orbit.dispatch == "self_energy_pi2":
            from montopar.self_energy_pi2 import build_result
        elif orbit.dispatch == "self_energy_chain":
            from montopar.self_energy_chain import build_result
        else:
            from montopar.self_energy import build_result
        result = build_result(orbit)
        directory = out_dir / f"{orbit.source['nloops']}loop" / orbit.emission_type / orbit.orbit_id
        write_json(result, directory / (f"montopar{suffix}.json" if suffix else "montopar.json"))
        return
    original_orbit = orbit
    hc = hypercube_transform(orbit.bounds)
    c_rot = None
    j_rot = None
    hc_rot = None
    j_raw_s = None
    j_raw_s_rot = None
    ibp = None
    ibp_rot = None

    if orbit.dispatch == "pure":
        color = compute_color(orbit.vertices, orbit.propagators)
        c = compute_c(orbit.n_gluons, color)
        j = compute_j(orbit, hc)

        # Rotated variant (C4+1)
        orbit_rot = rotate_orbit_full(orbit, 1)
        hc_rot = hypercube_transform(orbit_rot.bounds)
        color_rot = compute_color(orbit_rot.vertices, orbit_rot.propagators)
        c_rot = compute_c(orbit_rot.n_gluons, color_rot)
        j_rot = compute_j(orbit_rot, hc_rot)

    elif orbit.dispatch == "tgv":
        c = compute_c_tgv()
        j = compute_j_tgv(orbit, hc)
        if tgv_ibp_enabled:
            from montopar.tgv_ibp import compute_tgv_ibp
            ibp = compute_tgv_ibp(orbit, c)

        # s-variable version (before hypercube)
        S_s = compute_S(orbit.tgv_legs)
        Pnum_s = compute_Pnum(orbit.tgv_legs)
        s_vars = [b.var for b in orbit.bounds]
        j_raw_s = _build_jraw_s_tgv(Pnum_s, S_s, orbit.bounds, s_vars)

    elif orbit.dispatch == "fourv":
        # Keep one contour representative. Even zero graphs retain the common C;
        # their vanishing is explicitly a Lorentz zero in J.
        c = compute_c_fgv()
        j = compute_j_fgv(orbit, hc)

    elif orbit.dispatch == "mixed":
        # Normalize TGV to canonical form first
        orbit_canon = _normalize_tgv_to_canonical(orbit)

        # --- A group: canonical form ---
        tgv_loops = {leg.loop for leg in orbit_canon.tgv_legs}
        tgv_color_map: dict[str, str] = {}
        for prop in orbit_canon.propagators:
            if prop.loop in tgv_loops and prop.color:
                tgv_color_map[prop.loop] = prop.color[0]
        tgv_colors = [
            tgv_color_map[leg.loop]
            for leg in orbit_canon.tgv_legs
            if leg.loop in tgv_color_map
        ]
        color = compute_color(
            orbit_canon.vertices, orbit_canon.propagators,
            tgv_color_indices=tgv_colors if tgv_colors else None,
        )
        n_pure = orbit_canon.n_gluons
        pure = [p for p in orbit_canon.pairings if p.loop not in tgv_loops]

        c = compute_c_mixed(n_pure)
        c.color = color
        hc = hypercube_transform(orbit_canon.bounds)
        j = compute_j_mixed(orbit_canon, hc, pure)

        # s-variable version
        S_s = compute_S(orbit_canon.tgv_legs)
        Pnum_s = compute_Pnum(orbit_canon.tgv_legs)
        lorentz_z_s = compute_lorentz_factor(pure)
        f_polys_s: dict[str, str] = {}
        f_factors_s: list[str] = []
        for i, p in enumerate(pure, start=1):
            f_s = compute_f_p(p.edge_a, p.edge_b, p.s_a, p.s_b)
            f_polys_s[f"f{i}"] = f_s
            f_factors_s.append(f"f{i}^(-1+eps)")
        f_product_s = " * ".join(f_factors_s)
        s_vars = [b.var for b in orbit_canon.bounds]
        j_raw_s = _build_jraw_s_mixed(
            lorentz_z_s, f_product_s, Pnum_s, S_s, orbit_canon.bounds, s_vars,
            edge_pairs=[(p.edge_a, p.edge_b) for p in pure],
            f_polynomials=f_polys_s,
        )

        # Emit metadata from the same canonical frame used for C/J.
        orbit = orbit_canon

    elif orbit.dispatch == "double_tgv":
        from montopar.double_tgv_output import build_double_tgv_result
        result = build_double_tgv_result(orbit, ibp_enabled=tgv_ibp_enabled,
                                         forms_enabled=double_tgv_forms_enabled,
                                         operator_enabled=double_tgv_operator_enabled)
        directory = out_dir / f"{orbit.source['nloops']}loop" / orbit.emission_type / orbit.orbit_id
        filename = f"montopar{suffix}.json" if suffix else "montopar.json"
        write_json(result, directory / filename)
        return

    else:
        raise ValueError(f"unsupported dispatch: {orbit.dispatch}")

    # Emit
    nloops = orbit.source.get("nloops", "?")
    out_dir = out_dir / f"{nloops}loop" / orbit.emission_type / orbit.orbit_id
    montopar = emit_montopar_json(
        orbit, c, j, hc,
        c_rotated=c_rot,
        j_rotated=j_rot,
        hc_rotated=hc_rot,
        j_raw_s=j_raw_s,
        j_raw_s_rotated=j_raw_s_rot,
        ibp=ibp,
        ibp_rotated=ibp_rot,
    )
    if orbit.dispatch in ("tgv", "mixed"):
        from montopar.mixed_ibp import compute_mixed_ibp, finalize_single_representative
        if orbit.dispatch == "mixed" and tgv_ibp_enabled:
            montopar["J_ibp"] = compute_mixed_ibp(orbit, montopar["C"]["combined"])
        montopar = finalize_single_representative(montopar, original_orbit, orbit, hc)
    if orbit.dispatch == "fourv":
        from montopar.fgv_output import finalize_fgv
        montopar = finalize_fgv(montopar, orbit, hc)
    filename = f"montopar{suffix}.json" if suffix else "montopar.json"
    write_json(montopar, out_dir / filename)


def main():
    parser = argparse.ArgumentParser(
        description="montopar — direct parameter-space integral generation"
    )
    parser.add_argument(
        "--source", required=True, type=Path,
        help="Path to orbit_definitions.json",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=_find_repo_root() / "workspaces",
        help="Output directory (default: <repo_root>/workspaces/)",
    )
    parser.add_argument(
        "--orbit-id", type=str, default=None,
        help="Process a specific orbit by ID",
    )
    parser.add_argument(
        "--type", type=str, default=None,
        help="Process all orbits of a given emission type (e.g. 4e, 6e, 3e-3g, 4e-4g)",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Process all orbits (pure, tgv, fourv)",
    )
    parser.add_argument(
        "--suffix", type=str, default="",
        help="Suffix for output filename (e.g. '_new' writes montopar_new.json)",
    )
    parser.add_argument(
        "--tgv-ibp", "--gamma1-ibp", dest="tgv_ibp", action="store_true",
        help="Add contour-IBP sums for TGV, mixed, double-TGV, TGV-Pi1, and Gamma1 orbits",
    )
    parser.add_argument(
        "--package-initial", action="store_true",
        help="Write one montopar-initial.json package for the selected emission type(s)",
    )
    parser.add_argument(
        "--double-tgv-forms", action="store_true",
        help="Add compact and analytically reduced double-TGV C/J forms; retain all raw fields",
    )
    parser.add_argument(
        "--double-tgv-operator", action="store_true",
        help="Add an unapplied parent-position double-TGV operator; retain all scalar fields",
    )
    parser.add_argument(
        "--package-only", action="store_true",
        help="Build the initial package from existing montopar.json files without regeneration",
    )

    parser.add_argument("--fgv-compat-out-dir", type=Path, default=None,
        help="Also regenerate FGV final/state/export.dimreduce files under this output root")
    args = parser.parse_args()

    # Load source
    with open(args.source, encoding="utf-8") as f:
        source = json.load(f)

    # Filter before parsing so unrelated unsupported operators do not block a selected family.
    selected = source["orbits"]
    if args.orbit_id:
        selected = [o for o in selected if o["orbit_id"] == args.orbit_id]
    elif args.type:
        selected = [o for o in selected if o["emission_type"] == args.type]
    all_orbits = parse_all_orbits({"orbits": selected})

    # Filter targets
    if args.orbit_id:
        targets = [o for o in all_orbits if o.orbit_id == args.orbit_id]
        if not targets:
            print(f"Orbit {args.orbit_id} not found")
            sys.exit(1)
    elif args.type:
        targets = [o for o in all_orbits if o.emission_type == args.type]
        if not targets:
            print(f"No orbits of type {args.type} found")
            sys.exit(1)
    elif args.all:
        targets = all_orbits
    else:
        parser.print_help()
        sys.exit(1)

    if args.fgv_compat_out_dir and (any(o.dispatch != "fourv" for o in targets) or args.suffix):
        parser.error("--fgv-compat-out-dir requires only FGV targets and no --suffix")
    if args.package_only and not args.package_initial:
        parser.error("--package-only requires --package-initial")
    if args.double_tgv_forms and (args.package_only or any(o.dispatch != "double_tgv" for o in targets)):
        parser.error("--double-tgv-forms requires double-TGV targets and regeneration (no --package-only)")
    if args.double_tgv_operator and (args.package_only or any(o.dispatch != "double_tgv" for o in targets)):
        parser.error("--double-tgv-operator requires double-TGV targets and regeneration (no --package-only)")

    if not args.package_only:
        print(f"Processing {len(targets)} orbit(s)...")

        failures = []
        for orbit in targets:
            print(f"  {orbit.orbit_id}...", end=" ", flush=True)

            try:
                _compute_one(
                    orbit, args.out_dir, args.suffix,
                    tgv_ibp_enabled=args.tgv_ibp,
                    double_tgv_forms_enabled=args.double_tgv_forms,
                    double_tgv_operator_enabled=args.double_tgv_operator,
                )
                print("OK")
            except Exception as e:
                print(f"FAILED: {e}")
                failures.append(orbit.orbit_id)
                continue
        if failures:
            raise SystemExit("Generation failed; package not written: " + ", ".join(failures))

    if args.fgv_compat_out_dir:
        from montopar.fgv_output import write_compatibility_exports
        for orbit in targets:
            relative = Path(f"{orbit.source['nloops']}loop") / orbit.emission_type / orbit.orbit_id
            data = json.loads((args.out_dir / relative / "montopar.json").read_text(encoding="utf-8"))
            write_compatibility_exports(data, args.fgv_compat_out_dir / relative)

    if args.package_initial:
        from montopar.package_builder import write_initial_package
        types = sorted({orbit.emission_type for orbit in targets})
        for emission_type in types:
            path = write_initial_package(source, args.out_dir, emission_type)
            print(f"Initial package: {path}")

    print("Done.")


if __name__ == "__main__":
    main()
