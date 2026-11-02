"""C_n prefactor assembly: 5-factor decomposition.

From the derivation doc §6b:
  C_n = color × (4πα_s)^n × Γ(1-ε)^n/(4π^{2-ε})^n × (-2v_ref+i0)^{nε} × (-2)^{-n}

All five factors depend only on n and the color trace result.
None depend on the specific edge distribution.

TGV and FGV have their own closed-form C factorizations:
  TGV: C = color × coupling × ft_base × v14_scale × twoscale × derivative
  FGV: C = color × coupling × ft_base × v14_scale × twoscale  (no derivative)
  Mixed: C = pure(n) × tgv = color × coupling × ft_base × v14_scale × twoscale × derivative
"""
from montopar.models import CFactors


def compute_c(n: int, color: str, v_ref: str = "v14") -> CFactors:
    """Compute the C_n prefactor 5-factor decomposition.

    Args:
        n: Number of gluon propagators.
        color: Color factor string from SUNSimplify (e.g. "-1/2*((CA - 2*CF)*CF)").
        v_ref: Reference scale variable name (default "v14").

    Returns:
        CFactors with all 5 components.
    """
    return CFactors(
        color=color,
        coupling=f"(4*Pi*alphaS)^{n}",
        ft_base=f"Gamma[1-eps]^{n}/(4*Pi^(2-eps))^{n}",
        v14_scale=f"(-2*{v_ref}+I0)^({n}*eps)",
        twoscale=f"(-2)^(-{n})",
    )


def compute_c_tgv(v_ref: str = "v14") -> CFactors:
    """TGV C: 6 factors, all independent of edge distribution.

    From TGV_通用因子分解推导_v2.md §7a:
      C = ① color × ② vertex × ③ coupling × ④ ft_base × ⑤ derivative × ⑥ v14_scale
    """
    return CFactors(
        color="-I/2*CA*CF",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-2*eps]/(64*Pi^(4-2*eps))",
        v14_scale=f"(-2*{v_ref}+I0)^(2*eps)",
        twoscale="1",
        derivative="2*(1-2*eps)",
    )


def compute_c_star(v_ref: str = "v14") -> CFactors:
    """Undifferentiated three-point star-kernel prefactor for TGV IBP terms.

    It is related to the existing fully differentiated TGV prefactor by
    C_star = C_raw / (4 * (1 - 2 eps)).
    """
    return CFactors(
        color="CA*CF",
        coupling="(4*Pi*alphaS)^2",
        ft_base="Gamma[1-2*eps]/(256*Pi^(4-2*eps))",
        v14_scale=f"(-2*{v_ref}+I0)^(2*eps)",
        twoscale="1",
    )


def compute_c_fgv(v_ref: str = "v14") -> CFactors:
    """FGV C: 6 factors. lorentzZ is in J, not here.

    From 4GV-general-derivation.md, section 8:
      C = ① color × ② vertex × ③ coupling × ④ ft_base × ⑤ v14_scale × ⑥ twoscale
    """
    return CFactors(
        color="CA^2*CF",
        coupling="(4*Pi*alphaS)^3",
        ft_base="Gamma[2-3*eps]/(256*Pi^(6-3*eps))",
        v14_scale=f"(-2*{v_ref}+I0)^(3*eps)",
        twoscale="1/4",
    )


def compute_c_mixed(n_pure: int, v_ref: str = "v14") -> CFactors:
    """Mixed C: 7 factors = pure(n) × tgv.

    Args:
        n_pure: Number of pure gluon pairings.
        v_ref: Reference scale variable name.

    Returns:
        CFactors with color placeholder (replaced by CLI with
        wolframscript result).
    """
    return CFactors(
        color="MIXED_COLOR_PLACEHOLDER",
        coupling=f"(4*Pi*alphaS)^{n_pure + 2}",
        ft_base=(
            f"Gamma[1-eps]^{n_pure}*Gamma[1-2*eps]/"
            f"((4*Pi^(2-eps))^{n_pure}*64*Pi^(4-2*eps))"
        ),
        v14_scale=f"(-2*{v_ref}+I0)^({n_pure + 2}*eps)",
        twoscale=f"(-2)^(-{n_pure})",
        derivative="2*(1-2*eps)",
    )


def compute_c_double_tgv(color: str, v_ref: str = "v14") -> CFactors:
    """Raw factor from the corrected general derivation; no Euler preparation."""
    return CFactors(
        color=color,
        coupling="(alphaS/Pi)^3",
        ft_base="Pi^(3*eps)*Gamma[2-3*eps]/16",
        v14_scale=f"(-2*{v_ref}+I0)^(3*eps)",
        twoscale="1",
        derivative="1",
    )
