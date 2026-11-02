"""N=4 Wilson loop geometry utilities.

Provides dot products, vector differences, and Delta'_ij construction
for TGV and FGV Pnum/Q computations.
"""
from montopar.models import TGVLeg


def dot_product_factor(edge_a: int, edge_b: int) -> str:
    """v_a·v_b / v_ref for N=4 quadrilateral.

    Kinematics:
      v_14 = v_23 = v_ref
      v_12 = v_34 = z * v_ref
      v_13 = v_24 = -(1+z) * v_ref
      v_ii = 0  (lightlike edges)
    """
    if edge_a == edge_b:
        return "0"
    a, b = (edge_a, edge_b) if edge_a < edge_b else (edge_b, edge_a)
    if (a, b) in ((1, 2), (3, 4)):
        return "z1"
    if (a, b) in ((2, 3), (1, 4)):
        return "1"
    if (a, b) in ((1, 3), (2, 4)):
        return "-(1+z1)"
    raise ValueError(f"Unknown edge pair: ({edge_a}, {edge_b})")


# ---------------------------------------------------------------------------
# Internal: (constant, z1_coefficient) representation for linear-in-z1 exprs
# ---------------------------------------------------------------------------

def _parse_dot(dot: str) -> tuple[int, int]:
    """Parse a dot_product_factor result into (constant, z1_coefficient)."""
    if dot == "0":
        return (0, 0)
    if dot == "1":
        return (1, 0)
    if dot == "z1":
        return (0, 1)
    if dot == "-(1+z1)":
        return (-1, -1)
    raise ValueError(f"Unknown dot product expression: {dot!r}")


def _format_linear(c: int, z: int) -> str:
    """Format c + z*z1 as a string."""
    if c == 0 and z == 0:
        return "0"
    if z == 0:
        return str(c)
    if c == 0:
        if z == 1:
            return "z1"
        if z == -1:
            return "-z1"
        return f"{z}*z1"
    # Both nonzero
    if z == 1:
        z_part = "z1"
    elif z == -1:
        z_part = "-z1"
    else:
        z_part = f"{z}*z1"

    if c > 0 and z > 0:
        return f"{c} + {z_part}"
    if c > 0 and z < 0:
        # e.g. 2 + (-z1) → 2 - z1
        if z == -1:
            return f"{c} - z1"
        return f"{c} - {-z}*z1"
    if c < 0 and z > 0:
        return f"{c} + {z_part}"
    # c < 0 and z < 0
    return f"{c} + {z_part}"


def v_dot_V_diff(t: int, k: int, j: int) -> str:
    """Compute v_t · (V_k - V_j) / v_ref.

    V_1 = 0, V_i = Σ_{m=1}^{i-1} v_m for N=4.

    Args:
        t: Edge index of the vector v_t.
        k: Edge index of the first point.
        j: Edge index of the second point.

    Returns:
        String expression in z1 (or "1", "0").
    """
    if k == j:
        return "0"
    lo, hi = (j, k) if j < k else (k, j)
    sign = 1 if j < k else -1
    total_c, total_z = 0, 0
    for m in range(lo, hi):
        dot = dot_product_factor(t, m)
        c, z = _parse_dot(dot)
        total_c += c
        total_z += z
    total_c *= sign
    total_z *= sign
    return _format_linear(total_c, total_z)


def compute_Delta_prime(e_i: int, e_j: int, s_i: str, s_j: str) -> str:
    """Compute Δ'_ij = -(X_i - X_j)^2 / A_ref for a pair of gluon legs.

    Args:
        e_i, e_j: Edge indices (1-4) for the two legs.
        s_i, s_j: Schwinger parameter variable names.

    Returns:
        String expression in s_i, s_j, and z1. "0" if same edge.
    """
    if e_i == e_j:
        return "0"

    diff = abs(e_i - e_j) % 4
    if diff == 1 or diff == 3:
        # Adjacent edges: Δ'_ij = (1-s_early) * s_late * (v_early·v_late / v_ref)
        # where "early" is the edge that comes first in the contour 1→2→3→4→1.
        # diff=1 → early=min, late=max; diff=3 → early=max, late=min (wrap-around).
        dot = dot_product_factor(e_i, e_j)
        if diff == 1:
            e_early, e_late = (e_i, e_j) if e_i < e_j else (e_j, e_i)
        else:
            e_early, e_late = (e_i, e_j) if e_i > e_j else (e_j, e_i)
        s_early = s_i if e_early == e_i else s_j
        s_late = s_i if e_late == e_i else s_j
        return f"(1-{s_early})*{s_late}*{dot}"

    # Opposite edges (diff == 2)
    a, b = (e_i, e_j) if e_i < e_j else (e_j, e_i)
    if (a, b) == (1, 3):
        return f"z1*(1-{s_i})*(1-{s_j}) + {s_i}*{s_j}"
    else:  # (2, 4)
        return f"(1-{s_i})*(1-{s_j}) + z1*{s_i}*{s_j}"


def v_dot_X_diff(t: int, k_leg: TGVLeg, j_leg: TGVLeg) -> str:
    """Compute v_t · (X_k - X_j) / v_ref.

    X_k = V_{e_k} + s_k v_{e_k}, so:
      v_t · (X_k - X_j) / v_ref = v_t·(V_k-V_j)/v_ref + s_k(v_t·v_k)/v_ref - s_j(v_t·v_j)/v_ref

    Args:
        t: Edge index of the vector v_t.
        k_leg: The TGV leg at position k.
        j_leg: The TGV leg at position j.

    Returns:
        String expression in Schwinger parameters and z1.
    """
    V_diff = v_dot_V_diff(t, k_leg.edge, j_leg.edge)
    dot_t_k = dot_product_factor(t, k_leg.edge)
    dot_t_j = dot_product_factor(t, j_leg.edge)

    parts = []
    if V_diff != "0":
        parts.append(V_diff)

    if dot_t_k == "1":
        parts.append(k_leg.s_var)
    elif dot_t_k != "0":
        parts.append(f"{k_leg.s_var}*{dot_t_k}")

    if dot_t_j == "1":
        parts.append(f"-{j_leg.s_var}")
    elif dot_t_j != "0":
        parts.append(f"-{j_leg.s_var}*{dot_t_j}")

    return " + ".join(parts) if parts else "0"
