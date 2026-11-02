"""TGV Pnum and S computation.

Implements the O(Q) → Pnum algorithm from TGV_通用因子分解推导_v2.md §6-§7.

The O operator acts on Q to produce Pnum = N:
  O(Q^(3-D)) = -(D-3) * Q^(2-D) * O(Q)
  O(Q) = -2i * A_ref^2 * N
  N = 1/4 * Σ α_i α_j W_ij
"""
from montopar.models import TGVLeg
from montopar.geometry import compute_Delta_prime, v_dot_X_diff, dot_product_factor


def compute_S(legs: list[TGVLeg]) -> str:
    """Compute S = α_1 α_2 Δ'_12 + α_2 α_3 Δ'_23 + α_3 α_1 Δ'_31.

    S is in s-variables; s→x substitution is done downstream.

    Args:
        legs: 3 TGV legs in order (leg 1, leg 2, leg 3).

    Returns:
        String expression in α_i, s_i, z1.
    """
    D12 = compute_Delta_prime(
        legs[0].edge, legs[1].edge, legs[0].s_var, legs[1].s_var
    )
    D23 = compute_Delta_prime(
        legs[1].edge, legs[2].edge, legs[1].s_var, legs[2].s_var
    )
    D31 = compute_Delta_prime(
        legs[2].edge, legs[0].edge, legs[2].s_var, legs[0].s_var
    )

    terms = []
    if D12 != "0":
        terms.append(f"alpha1*alpha2*({D12})")
    if D23 != "0":
        terms.append(f"alpha2*alpha3*({D23})")
    if D31 != "0":
        terms.append(f"alpha3*alpha1*({D31})")

    return " + ".join(terms) if terms else "0"


def compute_Pnum(legs: list[TGVLeg]) -> str:
    """Compute Pnum = N = 1/4 Σ α_i α_j W_ij.

    Constructs the O operator's 6 terms, applies each to Q via
    v_t·∂_{X_k}Q = -2 α_k Σ_{j≠k} α_j · v_t·(X_k-X_j),
    collects α_i α_j coefficients, and assembles W_ij.

    Pnum is in s-variables; s→x substitution is done downstream.

    Args:
        legs: 3 TGV legs in order (leg 1, leg 2, leg 3).

    Returns:
        String expression in α_i, s_i, z1.
    """
    a, b, c = legs[0].edge, legs[1].edge, legs[2].edge

    # O operator terms: (sign, v_ab_factor, t_vector_edge, k_index)
    # sign = ±1 from ±i, v_ab_factor = v_ab/v_ref, t = edge of v_c
    O_terms = [
        (+1, dot_product_factor(a, b), c, 1),   # +i·v_ab · v_c·∂_{X_2}
        (-1, dot_product_factor(a, b), c, 0),   # -i·v_ab · v_c·∂_{X_1}
        (+1, dot_product_factor(b, c), a, 2),   # +i·v_bc · v_a·∂_{X_3}
        (-1, dot_product_factor(b, c), a, 1),   # -i·v_bc · v_a·∂_{X_2}
        (+1, dot_product_factor(c, a), b, 0),   # +i·v_ca · v_b·∂_{X_1}
        (-1, dot_product_factor(c, a), b, 2),   # -i·v_ca · v_b·∂_{X_3}
    ]

    # Collect α_i α_j coefficients → W[pair] = list of term strings
    W: dict[tuple[int, int], list[str]] = {
        (0, 1): [], (1, 2): [], (0, 2): [],
    }

    for sign, v_ab, t, k in O_terms:
        if v_ab == "0":
            continue  # degenerate: v_ab = 0 → term vanishes
        for j in range(3):
            if j == k:
                continue
            pair = tuple(sorted([k, j]))
            v_X_diff = v_dot_X_diff(t, legs[k], legs[j])
            if v_X_diff == "0":
                continue
            # Term: sign * i * (-2) * v_ab * v_t·(X_k - X_j)
            # O(Q) = -2i * A_ref^2 * N, so dividing by -2i A_ref^2 cancels
            # the i*(-2) prefactor: contribution to N is
            #   sign * i * (-2) * v_ab * v_t·(X_k-X_j) / (-2i * A_ref^2)
            #   = sign * (v_ab/v_ref) * (v_t·(X_k-X_j)/v_ref)
            #
            # sign=+1 → +v_ab * v_X_diff; sign=-1 → -v_ab * v_X_diff
            inner = f"{v_ab}*({v_X_diff})"
            if sign > 0:
                coeff = inner  # +i term → +v_ab*v_X_diff
            else:
                coeff = f"-{inner}"  # -i term → -v_ab*v_X_diff
            W[pair].append(coeff)

    # Assemble Pnum = 1/4 Σ α_i α_j W_ij
    alpha_pairs = {
        (0, 1): ("alpha1", "alpha2"),
        (1, 2): ("alpha2", "alpha3"),
        (0, 2): ("alpha1", "alpha3"),
    }
    terms = []
    for pair, (a1, a2) in alpha_pairs.items():
        if W[pair]:
            coeffs = " + ".join(W[pair])
            terms.append(f"{a1}*{a2}*({coeffs})")

    if not terms:
        return "0"
    return f"1/4*({' + '.join(terms)})"
