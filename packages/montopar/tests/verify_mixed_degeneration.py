"""Verify (3+2n)e-3g mixed formula: degeneration checks.

Verification 1: 5e-3g-01 pure part ≡ 2e-01 (remove TGV)
  - Compare f_p, Lorentz factor, Jacobian in s-variables

Verification 2: 5e-3g-01 TGV part ≡ 3e-3g-01 (remove pure)
  - Compare S, N (Pnum), Delta'_ij in s-variables
"""
import json
import sys
sys.path.insert(0, "src")
from montopar.parser import parse_all_orbits
from montopar.hypercube import hypercube_transform
from montopar.j_factor import (
    compute_f_p, compute_lorentz_factor, substitute_s_to_x, expand_polynomial,
)
from montopar.tgv import compute_S, compute_Pnum
from montopar.geometry import compute_Delta_prime


def load_orbits(path):
    with open(path) as f:
        return parse_all_orbits(json.load(f))


def find(orbits, oid):
    for o in orbits:
        if o.orbit_id == oid:
            return o
    raise KeyError(oid)


def main():
    # Load data
    orbits_2loop = load_orbits("../qgraf-pipeline/outputs/2loop/orbit_definitions.json")
    orbits_1loop = load_orbits("../qgraf-pipeline/outputs/1loop/orbit_definitions.json")
    orbits_3loop = load_orbits("../qgraf-pipeline/outputs/3loop/orbit_definitions.json")

    o5e3g = find(orbits_3loop, "5e-3g-01")
    o2e = find(orbits_1loop, "2e-01")
    o3e3g = find(orbits_2loop, "3e-3g-01")

    print("=" * 70)
    print("VERIFICATION 1: 5e-3g-01 pure part vs 2e-01")
    print("=" * 70)

    # 5e-3g-01 pure part
    print("\n--- 5e-3g-01 ---")
    print(f"Edge distribution: {o5e3g.edge_distribution.counts}")
    print(f"TGV legs: {[(leg.edge, leg.s_var) for leg in o5e3g.tgv_legs]}")

    pure_pairings = [p for p in o5e3g.pairings
                     if p.loop not in {leg.loop for leg in o5e3g.tgv_legs}]
    print(f"Pure pairings: {[(p.edge_a, p.edge_b, p.s_a, p.s_b) for p in pure_pairings]}")

    # Compute pure part for 5e-3g-01
    for p in pure_pairings:
        f_s = compute_f_p(p.edge_a, p.edge_b, p.s_a, p.s_b)
        print(f"  Pure f_p (s): {f_s}")
        print(f"  Lorentz: {compute_lorentz_factor([p])}")

    # 2e-01
    print("\n--- 2e-01 ---")
    print(f"Edge distribution: {o2e.edge_distribution.counts}")
    print(f"Pairings: {[(p.edge_a, p.edge_b, p.s_a, p.s_b) for p in o2e.pairings]}")
    for p in o2e.pairings:
        f_s = compute_f_p(p.edge_a, p.edge_b, p.s_a, p.s_b)
        print(f"  f_p (s): {f_s}")
        print(f"  Lorentz: {compute_lorentz_factor([p])}")

    # Compare functional forms (up to variable renaming)
    f_5e3g = compute_f_p(pure_pairings[0].edge_a, pure_pairings[0].edge_b,
                          pure_pairings[0].s_a, pure_pairings[0].s_b)
    f_2e = compute_f_p(o2e.pairings[0].edge_a, o2e.pairings[0].edge_b,
                        o2e.pairings[0].s_a, o2e.pairings[0].s_b)
    # Rename: 5e-3g uses s32,s43; 2e uses s31,s41
    f_5e3g_renamed = f_5e3g.replace("s32", "s31").replace("s43", "s41")
    print(f"\n  f_p(5e-3g renamed): {f_5e3g_renamed}")
    print(f"  f_p(2e-01):         {f_2e}")
    print(f"  MATCH: {f_5e3g_renamed == f_2e}")

    print("\n" + "=" * 70)
    print("VERIFICATION 2: 5e-3g-01 TGV part vs 3e-3g-01")
    print("=" * 70)

    # 5e-3g-01 TGV part
    print("\n--- 5e-3g-01 TGV legs ---")
    for leg in o5e3g.tgv_legs:
        print(f"  edge={leg.edge}, s_var={leg.s_var}, loop={leg.loop}")

    S_5e3g = compute_S(o5e3g.tgv_legs)
    Pnum_5e3g = compute_Pnum(o5e3g.tgv_legs)
    print(f"\n  S (s): {S_5e3g}")
    print(f"  Pnum (s): {Pnum_5e3g}")

    # Delta' for 5e-3g-01 TGV
    for i, j in [(0,1), (1,2), (0,2)]:
        D = compute_Delta_prime(
            o5e3g.tgv_legs[i].edge, o5e3g.tgv_legs[j].edge,
            o5e3g.tgv_legs[i].s_var, o5e3g.tgv_legs[j].s_var,
        )
        print(f"  Delta'_{i+1}{j+1}: {D}")

    # 3e-3g-01 TGV part
    print("\n--- 3e-3g-01 TGV legs ---")
    for leg in o3e3g.tgv_legs:
        print(f"  edge={leg.edge}, s_var={leg.s_var}, loop={leg.loop}")

    S_3e3g = compute_S(o3e3g.tgv_legs)
    Pnum_3e3g = compute_Pnum(o3e3g.tgv_legs)
    print(f"\n  S (s): {S_3e3g}")
    print(f"  Pnum (s): {Pnum_3e3g}")

    for i, j in [(0,1), (1,2), (0,2)]:
        D = compute_Delta_prime(
            o3e3g.tgv_legs[i].edge, o3e3g.tgv_legs[j].edge,
            o3e3g.tgv_legs[i].s_var, o3e3g.tgv_legs[j].s_var,
        )
        print(f"  Delta'_{i+1}{j+1}: {D}")

    # Compare S and Pnum (same s-variable names for TGV legs!)
    print(f"\n  S MATCH: {S_5e3g == S_3e3g}")
    print(f"  Pnum MATCH: {Pnum_5e3g == Pnum_3e3g}")

    print("\n" + "=" * 70)
    print("VERIFICATION 3: C factor decomposition")
    print("=" * 70)

    n_pure = len(pure_pairings)
    print(f"\n  5e-3g-01 n_pure = {n_pure}")
    print(f"  Expected C structure:")
    print(f"    C = C_pure({n_pure}) x C_tgv  (with unified color)")
    print(f"    C_pure: color x (4pi a_s)^{n_pure} x Gamma(1-eps)^{n_pure}/(4pi^(2-eps))^{n_pure}")
    print(f"            x v_ref^({n_pure}*eps) x (-2)^(-{n_pure})")
    print(f"    C_tgv:  i/2 CA CF x i^3 x (4pi a_s)^2 x Gamma(1-2eps)/(64pi^(4-2eps))")
    print(f"            x 2(1-2eps) x v_ref^(2eps)")
    print(f"    Total v_ref power: ({n_pure}+2)*eps = {n_pure+2}*eps")


if __name__ == "__main__":
    main()