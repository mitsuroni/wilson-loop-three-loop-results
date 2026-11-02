"""FGV Q' and lorentzZ computation.

Implements the Q' and lorentzZ algorithms from 4GV-general-derivation.md.

The FGV J integrand has no Pnum (no partial_X derivative), but has:
  - lorentzZ: pure z polynomial from color x lorentz channel sum
  - Q': denominator polynomial from 4-point Delta'_ij sum
"""
from montopar.models import FGVLeg, EdgeDistribution
from montopar.geometry import compute_Delta_prime, dot_product_factor


def is_physical_zero(edge_dist: EdgeDistribution) -> bool:
    """Check if max(m_i) >= 3 -> all Lorentz channels vanish -> I_4v = 0."""
    return max(edge_dist.counts.values(), default=0) >= 3


def compute_Q_prime(legs: list[FGVLeg]) -> str:
    """Compute Q' = Sum_{i<j} alpha_i alpha_j Delta'_ij.

    Q' is in s-variables; s->x substitution is done downstream.

    Args:
        legs: 4 FGV legs in order.

    Returns:
        String expression in alpha_i, s_i, z1.
    """
    terms = []
    for i in range(4):
        for j in range(i + 1, 4):
            D = compute_Delta_prime(
                legs[i].edge, legs[j].edge,
                legs[i].s_var, legs[j].s_var,
            )
            if D != "0":
                terms.append(f"alpha{i+1}*alpha{j+1}*({D})")
    return " + ".join(terms) if terms else "0"


def compute_lorentzZ(legs: list[FGVLeg]) -> str:
    """Compute lorentzZ = 1/4 * (-lorentz_1 + lorentz_3) / v_ref^2.

    Channels (color_2 = 0 always):
      lorentz_1 = v_e1.e3 * v_e2.e4 - v_e1.e4 * v_e2.e3
      lorentz_3 = v_e1.e2 * v_e3.e4 - v_e1.e3 * v_e2.e4

    Returns a pure z polynomial (string expression in z1).
    """
    e1, e2, e3, e4 = [leg.edge for leg in legs]

    dot_13 = dot_product_factor(e1, e3)
    dot_24 = dot_product_factor(e2, e4)
    dot_14 = dot_product_factor(e1, e4)
    dot_23 = dot_product_factor(e2, e3)
    dot_12 = dot_product_factor(e1, e2)
    dot_34 = dot_product_factor(e3, e4)

    # lorentz_1 = v_e1.e3 * v_e2.e4 - v_e1.e4 * v_e2.e3
    lorentz_1 = f"({dot_13})*({dot_24}) - ({dot_14})*({dot_23})"

    # lorentz_3 = v_e1.e2 * v_e3.e4 - v_e1.e3 * v_e2.e4
    lorentz_3 = f"({dot_12})*({dot_34}) - ({dot_13})*({dot_24})"

    # lorentzZ = 1/4 * (-lorentz_1 + lorentz_3)
    return f"1/4*(-({lorentz_1}) + ({lorentz_3}))"
