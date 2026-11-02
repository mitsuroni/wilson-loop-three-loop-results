"""Contour integration-by-parts representations for pure 3e-3g orbits."""
from __future__ import annotations

from montopar.c_factor import compute_c_star
from montopar.sign_rules import coordinate_phase
from montopar.geometry import dot_product_factor
from montopar.hypercube import hypercube_transform
from montopar.j_factor import (
    differentiate_polynomial,
    expand_polynomial,
    linear_combination_polynomials,
    substitute_s_to_x,
)
from montopar.models import (
    Bound,
    CFactors,
    Orbit,
    TGVIBPResult,
    TGVIBPTerm,
    TGVLeg,
)
from montopar.tgv import compute_Pnum, compute_S


_ALPHA_MONOMIAL = "(alpha1*alpha2*alpha3)^(-eps)"
_ALPHA_VARS = ["alpha1", "alpha2", "alpha3"]


def _negative(expr: str) -> str:
    return linear_combination_polynomials([("-1", expr)])


def _sum_expr(left: str, right: str, right_sign: str = "1") -> str:
    return linear_combination_polynomials([("1", left), (right_sign, right)])


def _combine_c_raw(c: CFactors) -> str:
    factors = [
        coordinate_phase(3,3,1,derivatives=1,chain_negatives=1), c.derivative, c.color, c.coupling,
        c.ft_base, c.v14_scale, c.twoscale,
    ]
    return " * ".join(f"({factor})" for factor in factors if factor)


def _combine_c_star(c: CFactors) -> str:
    factors = [c.color, c.coupling, c.ft_base, c.v14_scale, c.twoscale]
    return " * ".join(f"({factor})" for factor in factors)


def classify_tgv_ibp_case(legs: list[TGVLeg]) -> str:
    """Classify equality among the three ordered TGV legs."""
    if len(legs) != 3:
        return "unsupported"
    a, b, c = (leg.edge for leg in legs)
    if a == b == c:
        return "triple"
    if a == b:
        return "degenerate_ab"
    if b == c:
        return "degenerate_bc"
    if c == a:
        return "degenerate_ca"
    return "nondegenerate"


def degenerate_derivative_vector(legs: list[TGVLeg]) -> list[str]:
    """Return coefficients relative to the positive common C_star."""
    case = classify_tgv_ibp_case(legs)
    a, b, c = (leg.edge for leg in legs)
    if case == "degenerate_ab":
        rho = dot_product_factor(a, c)
        return [_negative(rho), expand_polynomial(rho), "0"]
    if case == "degenerate_bc":
        rho = dot_product_factor(a, b)
        return ["0", _negative(rho), expand_polynomial(rho)]
    if case == "degenerate_ca":
        rho = dot_product_factor(a, b)
        return [expand_polynomial(rho), "0", _negative(rho)]
    raise ValueError(f"not a two-leg degenerate TGV case: {case}")


def nondegenerate_derivative_vector(legs: list[TGVLeg]) -> list[str]:
    """Return (rho_bc, rho_ca, rho_ab) for a nondegenerate TGV."""
    if classify_tgv_ibp_case(legs) != "nondegenerate":
        raise ValueError("nondegenerate derivative vector requires distinct edges")
    a, b, c = (leg.edge for leg in legs)
    return [
        expand_polynomial(dot_product_factor(b, c)),
        expand_polynomial(dot_product_factor(c, a)),
        expand_polynomial(dot_product_factor(a, b)),
    ]


def find_ordered_pair(
    bounds: list[Bound], repeated_vars: set[str],
) -> tuple[str, str]:
    """Find p, q for a unique 0 <= p <= q <= 1 repeated-edge domain."""
    candidates = [
        (bound.var, bound.high)
        for bound in bounds
        if bound.var in repeated_vars
        and bound.low == "0"
        and bound.high in repeated_vars
        and bound.high != bound.var
    ]
    if len(candidates) != 1:
        raise ValueError("unrecognized_ordered_pair")
    inner, outer = candidates[0]
    outer_bounds = [bound for bound in bounds if bound.var == outer]
    if len(outer_bounds) != 1 or outer_bounds[0].low != "0" or outer_bounds[0].high != "1":
        raise ValueError("unrecognized_ordered_pair")
    return inner, outer


def _substitute_parameters(
    expr: str,
    fixed: dict[str, str],
    identified: dict[str, str],
) -> str:
    substitutions = dict(identified)
    substitutions.update(fixed)
    return substitute_s_to_x(expr, substitutions)


def _localize_boundary(
    S_s: str,
    bounds: list[Bound],
    *,
    fixed: dict[str, str] | None = None,
    identified: dict[str, str] | None = None,
) -> tuple[str, list[str], dict[str, str], dict[str, list[int]]]:
    """Apply a boundary constraint, then assign fresh local x variables."""
    fixed = fixed or {}
    identified = identified or {}
    constrained = _substitute_parameters(S_s, fixed, identified)
    removed = set(fixed) | set(identified)
    remaining = [bound.var for bound in bounds if bound.var not in removed]
    s_to_x = {var: f"x{i}" for i, var in enumerate(remaining, start=1)}
    S_x = expand_polynomial(substitute_s_to_x(constrained, s_to_x))
    x_vars = list(s_to_x.values())
    x_bounds = {var: [0, 1] for var in x_vars}
    return S_x, x_vars, s_to_x, x_bounds


def _boundary_integrand(S_x: str) -> str:
    return f"({_ALPHA_MONOMIAL}) * ({S_x})^(-(1 - 2*eps))"


def _raw_integrand(jacobian: str, pnum: str, S_x: str) -> str:
    return (
        f"({jacobian}) * ({pnum}) * ({_ALPHA_MONOMIAL})"
        f" * ({S_x})^(-(2 - 2*eps))"
    )


def _make_boundary_term(
    *,
    term_id: str,
    kind: str,
    coefficient: str,
    source_derivative: str,
    S_s: str,
    bounds: list[Bound],
    c_star_combined: str,
    fixed: dict[str, str] | None = None,
    identified: dict[str, str] | None = None,
) -> TGVIBPTerm:
    coefficient = expand_polynomial(coefficient)
    S_x, x_vars, s_to_x, x_bounds = _localize_boundary(
        S_s, bounds, fixed=fixed, identified=identified,
    )
    return TGVIBPTerm(
        term_id=term_id,
        kind=kind,
        coefficient=coefficient,
        c_ref="C_star",
        combined_C=f"({coefficient}) * ({c_star_combined})",
        integrand=_boundary_integrand(S_x),
        variables=x_vars,
        alpha_vars=list(_ALPHA_VARS),
        bounds=x_bounds,
        delta_constraints={},
        denom_power="1 - 2*eps",
        source_derivative=source_derivative,
        fixed_parameters=fixed or {},
        identified_parameters=identified or {},
        s_to_x=s_to_x,
    )


def build_triangle_boundary_terms(
    orbit: Orbit,
    S_s: str,
    derivative_vector: list[str],
    c_star_combined: str,
) -> list[TGVIBPTerm]:
    """Apply the ordered-triangle flux formula for a repeated edge."""
    repeated = {}
    for index, leg in enumerate(orbit.tgv_legs):
        repeated.setdefault(leg.edge, []).append((index, leg.s_var))
    repeated_entries = [entries for entries in repeated.values() if len(entries) == 2]
    if len(repeated_entries) != 1:
        raise ValueError("unrecognized_repeated_edge")
    entries = repeated_entries[0]
    index_by_var = {var: index for index, var in entries}
    inner, outer = find_ordered_pair(orbit.bounds, set(index_by_var))
    c_inner = derivative_vector[index_by_var[inner]]
    c_outer = derivative_vector[index_by_var[outer]]
    lower_coeff = _negative(c_inner)
    upper_coeff = expand_polynomial(c_outer)
    diagonal_coeff = _sum_expr(c_inner, c_outer, "-1")

    return [
        _make_boundary_term(
            term_id=f"boundary_{inner}_0",
            kind="cusp_boundary",
            coefficient=lower_coeff,
            source_derivative=f"{c_inner}*d/d{inner}",
            S_s=S_s,
            bounds=orbit.bounds,
            c_star_combined=c_star_combined,
            fixed={inner: "0"},
        ),
        _make_boundary_term(
            term_id=f"boundary_{outer}_1",
            kind="cusp_boundary",
            coefficient=upper_coeff,
            source_derivative=f"{c_outer}*d/d{outer}",
            S_s=S_s,
            bounds=orbit.bounds,
            c_star_combined=c_star_combined,
            fixed={outer: "1"},
        ),
        _make_boundary_term(
            term_id=f"boundary_{inner}_eq_{outer}",
            kind="contact_boundary",
            coefficient=diagonal_coeff,
            source_derivative=f"{c_inner}*d/d{inner} + {c_outer}*d/d{outer}",
            S_s=S_s,
            bounds=orbit.bounds,
            c_star_combined=c_star_combined,
            identified={inner: outer},
        ),
    ]


def _validate_cube(bounds: list[Bound], leg_vars: set[str]) -> None:
    relevant = [bound for bound in bounds if bound.var in leg_vars]
    if len(relevant) != 3 or any(bound.low != "0" or bound.high != "1" for bound in relevant):
        raise ValueError("nondegenerate_domain_is_not_unit_cube")


def build_cube_boundary_terms(
    orbit: Orbit,
    S_s: str,
    derivative_vector: list[str],
    c_star_combined: str,
) -> list[TGVIBPTerm]:
    """Build the six oriented faces of a three-dimensional unit cube."""
    leg_index = {leg.s_var: index for index, leg in enumerate(orbit.tgv_legs)}
    _validate_cube(orbit.bounds, set(leg_index))
    terms: list[TGVIBPTerm] = []
    for bound in orbit.bounds:
        if bound.var not in leg_index:
            continue
        coefficient = derivative_vector[leg_index[bound.var]]
        for value, oriented in (("0", _negative(coefficient)), ("1", coefficient)):
            terms.append(_make_boundary_term(
                term_id=f"boundary_{bound.var}_{value}",
                kind="cusp_boundary",
                coefficient=oriented,
                source_derivative=f"{coefficient}*d/d{bound.var}",
                S_s=S_s,
                bounds=orbit.bounds,
                c_star_combined=c_star_combined,
                fixed={bound.var: value},
            ))
    return terms


def _pnum_boundary(S_s: str, legs: list[TGVLeg], coefficients: list[str]) -> str:
    terms = []
    for coefficient, leg in zip(coefficients, legs):
        derivative = differentiate_polynomial(S_s, leg.s_var)
        terms.append((f"-1/4*({coefficient})", derivative))
    return linear_combination_polynomials(terms)


def build_recoil_term(
    orbit: Orbit,
    S_s: str,
    pnum_recoil_s: str,
    c_raw_combined: str,
) -> TGVIBPTerm:
    """Build the full-dimensional recoil integral in the orbit hypercube."""
    hc = hypercube_transform(orbit.bounds)
    S_x = expand_polynomial(substitute_s_to_x(S_s, hc.s_to_x))
    pnum_x = expand_polynomial(substitute_s_to_x(pnum_recoil_s, hc.s_to_x))
    return TGVIBPTerm(
        term_id="recoil",
        kind="recoil",
        coefficient="1",
        c_ref="C_raw",
        combined_C=c_raw_combined,
        integrand=_raw_integrand(hc.jacobian, pnum_x, S_x),
        variables=list(hc.x_vars),
        alpha_vars=list(_ALPHA_VARS),
        bounds={var: [0, 1] for var in hc.x_vars},
        delta_constraints={},
        denom_power="2 - 2*eps",
        s_to_x=dict(hc.s_to_x),
    )


def _unsupported(
    case: str,
    c_star: CFactors,
    reason: str,
    derivative_vector: list[str] | None = None,
) -> TGVIBPResult:
    return TGVIBPResult(
        method="contour_ibp",
        status="unsupported",
        case=case,
        C_star=c_star,
        terms=[],
        sum_expression="0",
        certificates={},
        derivative_vector=derivative_vector or [],
        reason=reason,
    )


def compute_tgv_ibp(
    orbit: Orbit,
    c_raw: CFactors,
    *,
    verify: bool = True,
) -> TGVIBPResult:
    """Compute the contour-IBP term sum for one pure TGV orbit variant."""
    c_star = compute_c_star()
    case = classify_tgv_ibp_case(orbit.tgv_legs)
    if orbit.dispatch != "tgv":
        return _unsupported(case, c_star, "dispatch_is_not_tgv")
    if case == "unsupported":
        return _unsupported(case, c_star, "expected_three_tgv_legs")
    if case == "triple":
        return TGVIBPResult(
            method="contour_ibp", status="complete", case=case,
            C_star=c_star, terms=[], sum_expression="0",
            certificates={"physical_zero": True}, derivative_vector=["0", "0", "0"],
        )

    S_s = compute_S(orbit.tgv_legs)
    pnum_s = compute_Pnum(orbit.tgv_legs)
    c_star_combined = _combine_c_star(c_star)
    c_raw_combined = _combine_c_raw(c_raw)

    try:
        if case.startswith("degenerate_"):
            derivative_vector = degenerate_derivative_vector(orbit.tgv_legs)
            pnum_boundary = _pnum_boundary(S_s, orbit.tgv_legs, derivative_vector)
            certificate = linear_combination_polynomials([
                ("1", pnum_s), ("-1", pnum_boundary),
            ]) == "0"
            if verify and not certificate:
                return _unsupported(case, c_star, "pnum_identity_failed", derivative_vector)
            terms = build_triangle_boundary_terms(
                orbit, S_s, derivative_vector, c_star_combined,
            )
            certificates = {
                "operator_identity": True,
                "boundary_orientation": True,
                "pnum_identity": certificate,
            }
        else:
            derivative_vector = nondegenerate_derivative_vector(orbit.tgv_legs)
            pnum_boundary = _pnum_boundary(S_s, orbit.tgv_legs, derivative_vector)
            pnum_recoil = linear_combination_polynomials([
                ("1", pnum_s), ("-1", pnum_boundary),
            ])
            certificate = linear_combination_polynomials([
                ("1", pnum_s),
                ("-1", pnum_boundary),
                ("-1", pnum_recoil),
            ]) == "0"
            if verify and not certificate:
                return _unsupported(case, c_star, "pnum_identity_failed", derivative_vector)
            terms = build_cube_boundary_terms(
                orbit, S_s, derivative_vector, c_star_combined,
            )
            terms.append(build_recoil_term(
                orbit, S_s, pnum_recoil, c_raw_combined,
            ))
            certificates = {
                "operator_identity": True,
                "boundary_orientation": True,
                "pnum_identity": certificate,
            }
    except ValueError as exc:
        return _unsupported(case, c_star, str(exc), locals().get("derivative_vector", []))

    sum_expression = " + ".join(
        f"term({term.term_id})" for term in terms
    ) if terms else "0"
    return TGVIBPResult(
        method="contour_ibp",
        status="complete",
        case=case,
        C_star=c_star,
        terms=terms,
        sum_expression=sum_expression,
        certificates=certificates,
        derivative_vector=derivative_vector,
    )
