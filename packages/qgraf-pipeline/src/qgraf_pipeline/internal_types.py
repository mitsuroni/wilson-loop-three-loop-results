"""Internal structure registry: (M, N) → type name, g-power, weight.
Emission type enumeration for given n-loop.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations_with_replacement
from typing import Dict, List


@dataclass
class EmissionType:
    """A pre-derived emission type for a given n-loop order."""
    name: str                          # e.g. "4e", "3e-3g", "4e-Pi1"
    category: str                      # "tree" or "loop"
    eikonal_count: int                 # E = total eikonal vertices
    internal: Dict[str, int] = field(default_factory=dict)
    wl_wl_pairs: int = 0               # p = number of WL-WL gluon pairs
    g_power: int = 0


def get_internal_type(M: int, N: int) -> dict:
    """Given (M, N), return internal structure type info.

    M = loop count, N = external leg count.
    g_power = 2M + N - 2
    weight = M + N - 1  (used in enumeration constraint)
    """
    weight = M + N - 1
    g_power = 2 * M + N - 2

    if M == 0:
        if N == 3:
            name = "3g"
        elif N == 4:
            name = "4g"
        else:
            name = f"T{N}"
    elif N == 2:
        name = "Pi1" if M == 1 else f"Pi{M}"
    elif N == 3:
        name = "Gamma1" if M == 1 else f"Gamma{M}"
    else:
        name = f"L{M}N{N}"

    return {"name": name, "M": M, "N": N, "g_power": g_power, "weight": weight}


def get_vertex_info(M: int, N: int) -> dict:
    """Return LaTeX emission info for an N-point M-loop effective vertex.

    Returns:
        legs: number of external gluon legs connecting to eikonal vertices
        shared: True if all legs share the same loop momentum (M>0)
        sign: external insertion coefficient (+1); internal Grassmann signs belong to the operator
        latex: LaTeX symbol for the effective vertex
    """
    shared = M > 0
    sign = 1  # Included-i 1PI objects use positive Dyson chaining.

    if N == 2:
        base = "\\Pi"
        indices = "_{b_1 b_2}^{\\nu_1\\nu_2}(\\ell)"
    elif N == 3:
        base = "\\Gamma"
        indices = "_{\\nu_1\\nu_2\\nu_3}^{b_1 b_2 b_3}(-\\ell_1,-\\ell_2,-\\ell_3)"
    elif N == 4:
        base = "V"
        indices = "_{\\nu_1\\nu_2\\nu_3\\nu_4}^{b_1 b_2 b_3 b_4}"
    else:
        base = "\\mathcal{I}"
        indices = f"_{{\\nu_1\\ldots\\nu_{{{N}}}}}^{{b_1\\ldots b_{{{N}}}}}"

    if M >= 2:
        latex = f"{base}^{{({M})}}{indices}"
    elif M == 1 and N >= 3:
        latex = f"{base}^{{(1)}}{indices}"
    else:
        latex = f"{base}{indices}"

    return {"legs": N, "shared": shared, "sign": sign, "latex": latex}


def _name_to_mn(name: str) -> tuple:
    """Reverse lookup: internal type name -> (M, N)."""
    if name == "3g":
        return (0, 3)
    if name == "4g":
        return (0, 4)
    if name.startswith("Pi"):
        m = int(name[2:]) if len(name) > 2 else 1
        return (m, 2)
    if name.startswith("Gamma"):
        m = int(name[5:]) if len(name) > 5 else 1
        return (m, 3)
    # Generic: L{M}N{N}
    if name.startswith("L"):
        parts = name[1:].split("N")
        return (int(parts[0]), int(parts[1]))
    return (0, 0)


def enumerate_emission_types(n: int) -> List[EmissionType]:
    """Enumerate all emission types for n-loop from g-power constraint.

    Constraint: Σ (M_i + N_i - 1) + p = n
    where p = number of WL-WL gluon pairs.
    E = Σ N_i + 2p
    """
    total_g = 2 * n
    types: List[EmissionType] = []
    seen: set = set()

    # Collect all possible internal structures with weight ≤ n
    possible_internal = []
    for N in (2, 3, 4):
        for M in range(0, n + 1):
            # Bug fix: skip tree-level N=2 (T2 = gluon propagator, covered by WL-WL pairs)
            if M == 0 and N == 2:
                continue
            info = get_internal_type(M, N)
            if info["weight"] <= n:
                possible_internal.append(info)

    # Case: no internal structures (pure WL-WL)
    p = n
    E = 2 * p
    name = f"{E}e"
    seen.add(name)
    types.append(EmissionType(
        name=name, category="tree",
        eikonal_count=E, wl_wl_pairs=p, g_power=total_g,
    ))

    # Case: 1 internal structure
    for info in possible_internal:
        w = info["weight"]
        if w > n:
            continue
        p = n - w
        E = info["N"] + 2 * p
        name = f"{E}e-{info['name']}"
        if name not in seen:
            seen.add(name)
            cat = "tree" if info["M"] == 0 else "loop"
            types.append(EmissionType(
                name=name, category=cat,
                eikonal_count=E,
                internal={info["name"]: 1},
                wl_wl_pairs=p, g_power=total_g,
            ))

    # Case: 2 internal structures connected by an internal gluon line.
    # The generic branch below handles vertex-type structures (N >= 3).
    # The N=2 serial insertion topology is registered separately below,
    # retaining two distinct 1PI blobs and their middle propagator.
    for a, b in combinations_with_replacement(possible_internal, 2):
        if a["N"] < 3 or b["N"] < 3:
            continue
        w = a["weight"] + b["weight"] - 1
        if w > n:
            continue
        p = n - w
        E = a["N"] + b["N"] - 2 + 2 * p

        internal: Dict[str, int] = {}
        for info in (a, b):
            internal[info["name"]] = internal.get(info["name"], 0) + 1

        parts = []
        for tname in sorted(internal.keys()):
            # Repeat by multiplicity so {"3g": 2} -> ["3g", "3g"]
            # giving "4e-3g+3g" rather than the ambiguous "4e-3g".
            for _ in range(internal[tname]):
                parts.append(tname)
        name = f"{E}e-" + "+".join(parts)

        if name not in seen:
            seen.add(name)
            cat = "tree" if a["M"] == 0 and b["M"] == 0 else "loop"
            types.append(EmissionType(
                name=name, category=cat,
                eikonal_count=E,
                internal=internal,
                wl_wl_pairs=p, g_power=total_g,
            ))

    # The existing weight formula for two joined components also admits
    # N=2 + N=2: 2+2-1=3, E=2. Keep this reducible propagator chain
    # distinct from Pi2 (a single two-loop 1PI blob). Current scope: n=3.
    if n == 3:
        types.append(EmissionType(name="2e-Pi1+Pi1", category="loop",
            eikonal_count=2, internal={"Pi1": 2}, wl_wl_pairs=0, g_power=6))

    # One Pi1 on a bare TGV leg; distinct from the central 1PI Gamma1.
    if n == 3:
        types.append(EmissionType(name="3e-3g+Pi1", category="loop",
            eikonal_count=3, internal={"3g": 1, "Pi1": 1}, wl_wl_pairs=0, g_power=6))

    return types


if __name__ == "__main__":
    for n_test in (2, 3):
        print(f"\nn={n_test}:")
        for et in enumerate_emission_types(n_test):
            print(f"  {et.name}  E={et.eikonal_count}  "
                  f"wl_wl={et.wl_wl_pairs}  internal={et.internal}  "
                  f"cat={et.category}")
