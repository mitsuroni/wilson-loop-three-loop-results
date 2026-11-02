from montopar.sign_rules import coordinate_phase, annotate_result
"""Output emitter: generate montopar.json for workspace."""
import json
from pathlib import Path

from montopar.models import (
    Orbit, CFactors, JFactors, HypercubeResult,
    TGVJFactors, FGVJFactors, MixedJFactors, DoubleTGVJFactors,
    TGVIBPResult,
)
from montopar.j_factor import compute_j_expanded


def build_kinematics_map() -> dict[str, str]:
    """N=4 kinematics map for v_i.v_j -> f_ij(z1) * v14."""
    return {
        "SP[v1,v4]": "v14",
        "SP[v3,v4]": "z1*v14",
        "SP[v1,v2]": "z1*v14",
        "SP[v1,v3]": "(-1 - z1)*v14",
        "SP[v2,v3]": "v14",
        "SP[v2,v4]": "(-1 - z1)*v14",
    }


def build_c_components(c: CFactors, n: int) -> dict[str, str]:
    """Build C components dict in hepeflow-compatible format.

    The hepeflow format has:
      - ft_prefactor: combined ft_base x twoscale (the pure FT prefactor)
      - ft_v14_scale: v14 scale from FT
      - v14_scale: Lorentz scalar v14 power
    """
    # ft_prefactor = ft_base with the Pi/Gamma part
    ft_prefactor = c.ft_base

    # ft_v14_scale: 2^(n*(-1+eps)) * (-v14)^(n*(-1+eps))
    ft_v14_scale = f"2^({n}*(-1+eps))*(-v14)^({n}*(-1+eps))"

    # v14_scale: the remaining v14 power from Lorentz scalar
    v14_scale = f"v14^{n}"

    return {
        "i_factor": f"I^{2*n}",
        "color": c.color,
        "coupling": c.coupling,
        "ft_prefactor": ft_prefactor,
        "ft_v14_scale": ft_v14_scale,
        "v14_scale": v14_scale,
    }


def combine_c(c: CFactors) -> str:
    """Combine all C factors into a single Mathematica expression."""
    factors = [c.color, c.coupling, c.ft_base, c.v14_scale, c.twoscale]
    return " * ".join(f"({f})" for f in factors)


def _combine_c_tgv(c: CFactors) -> str:
    """Combine TGV C factors including i_factor and derivative (7 factors).

    Derivation: ① color = I/2·CA·CF, ② vertex = I³ = -I.
    ① × ② = 1/2·CA·CF (real).  The i_factor I³ must be included
    in the combined expression so downstream consumers can verify
    I-cancellation.
    """
    factors = [coordinate_phase(3,3,1,derivatives=1,chain_negatives=1), c.derivative, c.color, c.coupling, c.ft_base, c.v14_scale, c.twoscale]
    return " * ".join(f"({f})" for f in factors if f)


def build_j_components(j: JFactors) -> dict[str, str]:
    """Build J components dict in hepeflow-compatible format."""
    return {
        "lorentz_z": j.lorentz_z,
        "jacobian": j.jacobian,
        "f_product": j.f_product,
    }


def build_j_integrand(j: JFactors) -> str:
    """Build the full J integrand string."""
    parts = [f"({j.lorentz_z})", f"({j.jacobian})", j.f_product]
    return " * ".join(parts)


def build_geometry() -> dict:
    """Build the geometry metadata block shared by all dispatch types."""
    return {
        "n_lines": 4,
        "line_vectors": ["v1", "v2", "v3", "v4"],
        "lightlike": True,
        "cross_ratios": [{"name": "z1", "definition": "v34/v14"}],
        "reference_scale": "v14",
    }


def _build_c_pure(c: CFactors, n: int) -> dict:
    """Build C block for pure n-gluon dispatch."""
    components={
            "color": c.color,
            "coupling": c.coupling,
            "ft_prefactor": c.ft_base,
            "v14_scale": c.v14_scale,
            "twoscale": c.twoscale,
            "i_factor": coordinate_phase(2*n,n,0),
    }
    return {'components':components,'component_semantics':'multiplicative','factor_order':list(components),
        'combined':' * '.join(f'({v})' for v in components.values()),
        'phase_provenance':{'wilson_vertices':f'I^{2*n}','coordinate_propagators':f'(-1)^{n}',
                            'product':coordinate_phase(2*n,n,0)}}


def _build_c_tgv_block(c: CFactors) -> dict:
    """Build C block for TGV dispatch."""
    return {
        "components": {
            "color": c.color,
            "coupling": c.coupling,
            "ft_prefactor": c.ft_base,
            "derivative": c.derivative,
            "v14_scale": c.v14_scale,
            "twoscale": c.twoscale,
            "i_factor": coordinate_phase(3,3,1,derivatives=1,chain_negatives=1),
        },
        "component_semantics": "multiplicative",
        "factor_order": ["i_factor","derivative","color","coupling","ft_prefactor","v14_scale","twoscale"],
        "phase_provenance": {"wilson_vertices":"I^3","coordinate_propagators":"(-1)^3",
            "coordinate_centres":"-I","momentum_derivative":"-I","chain_rule":"-1",
            "product":coordinate_phase(3,3,1,derivatives=1,chain_negatives=1)},
        "combined": _combine_c_tgv(c),
    }


def _build_c_fgv_block(c: CFactors) -> dict:
    """Multiplicative factors; phase provenance is kept outside their product."""
    components = {"color": c.color, "coupling": c.coupling,
        "ft_prefactor": c.ft_base, "v14_scale": c.v14_scale,
        "twoscale": c.twoscale, "i_factor": coordinate_phase(4,4,1,vertex4=1)}
    return {"component_semantics": "multiplicative",
        "components": components, "factor_order": list(components),
        "combined": "0" if c.color == "0" else " * ".join(f"({v})" for v in components.values()),
        "phase_provenance": {"wilson_vertices": "I^4", "four_gluon_vertex": "-I",
            "four_coordinate_propagators": "(-1)^4", "star_Wick_rotation": "-I", "product": coordinate_phase(4,4,1,vertex4=1)}}


def _build_jraw_pure(j: JFactors, hc: HypercubeResult) -> dict:
    """Build J_raw block for pure n-gluon dispatch."""
    return {
        "integrand": build_j_integrand(j),
        "expanded": (
            f"({j.lorentz_z}) * ({j.jacobian}) * {compute_j_expanded(j)}"
        ),
        "variables": hc.x_vars,
        "bounds": {x: [0, 1] for x in hc.x_vars},
        "delta_constraints": {},
    }


def _build_jraw_tgv(j: TGVJFactors, hc: HypercubeResult) -> dict:
    """Build J_raw block for TGV dispatch."""
    integrand = (
        f"({j.jacobian}) * ({j.pnum}) * ({j.alpha_monomial})"
        f" * ({j.S_polynomial})^(-({j.denom_power}))"
    )
    return {
        "integrand": integrand,
        "expanded": None,
        "variables": hc.x_vars,
        "alpha_vars": ["alpha1", "alpha2", "alpha3"],
        "bounds": {x: [0, 1] for x in hc.x_vars},
        "delta_constraints": {},
    }


def _build_jraw_s_tgv(
    pnum_s: str, S_s: str, bounds: list, s_vars: list[str],
) -> dict:
    """Build J_raw_s block for TGV dispatch (before hypercube transform).

    J = Pnum(s) * (α₁α₂α₃)^(-ε) * S(s)^(-(2-2ε))
    with Schwinger-parameter bounds.
    """
    integrand = (
        f"({pnum_s}) * ((alpha1*alpha2*alpha3)^(-eps))"
        f" * ({S_s})^(-(2 - 2*eps))"
    )
    return {
        "integrand": integrand,
        "variables": s_vars,
        "alpha_vars": ["alpha1", "alpha2", "alpha3"],
        "bounds": _format_bounds_s(bounds),
        "delta_constraints": {},
    }


def _build_jibp_tgv(result: TGVIBPResult, c_raw: CFactors) -> dict:
    """Serialize a contour-IBP TGV result without doing algebra in emitter."""
    c_star_combined = combine_c(result.C_star)
    payload = {
        "schema_version": "1.0",
        "method": result.method,
        "status": result.status,
        "case": result.case,
        "prefactors": {
            "C_raw": _combine_c_tgv(c_raw),
            "C_star": c_star_combined,
            "relation": "C_star = C_raw/(4*(1-2*eps))",
        },
        "operator": {
            "derivative_vector": result.derivative_vector,
        },
        "terms": [
            {
                "term_id": term.term_id,
                "kind": term.kind,
                "coefficient": term.coefficient,
                "C_ref": term.c_ref,
                "combined_C": term.combined_C,
                "source_derivative": term.source_derivative,
                "fixed_parameters": term.fixed_parameters,
                "identified_parameters": term.identified_parameters,
                "s_to_x": term.s_to_x,
                "integrand": term.integrand,
                "variables": term.variables,
                "alpha_vars": term.alpha_vars,
                "bounds": term.bounds,
                "delta_constraints": term.delta_constraints,
                "denom_power": term.denom_power,
            }
            for term in result.terms
        ],
        "sum_expression": result.sum_expression,
        "certificates": result.certificates,
    }
    if result.reason is not None:
        payload["reason"] = result.reason
    return payload


def _build_jraw_s_mixed(
    lorentz_z: str, f_product: str,
    pnum_s: str, S_s: str, bounds: list, s_vars: list[str],
    edge_pairs: list[tuple[int, int]] | None = None,
    f_polynomials: dict[str, str] | None = None,
) -> dict:
    """Build J_raw_s block for mixed dispatch (before hypercube transform).

    J = lorentzZ * f_product(s) * (α₁α₂α₃)^(-ε) * Pnum(s) * S(s)^(-(2-2ε))
    """
    result: dict = {
        "integrand": (
            f"({lorentz_z}) * {f_product} * ((alpha1*alpha2*alpha3)^(-eps))"
            f" * ({pnum_s}) * ({S_s})^(-(2 - 2*eps))"
        ),
        "variables": s_vars,
        "alpha_vars": ["alpha1", "alpha2", "alpha3"],
        "bounds": _format_bounds_s(bounds),
        "delta_constraints": {},
    }
    if not f_polynomials:
        raise ValueError("original_s_f_polynomials_required")
    result["f_polynomials"] = dict(f_polynomials)
    result["pure"] = {"lorentz_z": lorentz_z, "f_product": f_product,
                      "f_polynomials": dict(f_polynomials)}
    result["tgv"] = {"pnum": pnum_s, "S_polynomial": S_s,
                    "alpha_monomial": "(alpha1*alpha2*alpha3)^(-eps)",
                    "denom_power": "2-2*eps"}
    if edge_pairs:
        result["edge_pairs"] = [[a, b] for a, b in edge_pairs]
    return result


def _format_bounds_s(bounds: list) -> dict:
    """Format Schwinger-parameter bounds as {var: [low, high]}."""
    result = {}
    for b in bounds:
        result[b.var] = [b.low, b.high]
    return result


def _build_jraw_fgv(j: FGVJFactors, hc: HypercubeResult) -> dict:
    """Cube kernel with an explicit four-alpha simplex, including zero graphs."""
    from montopar.fgv_output import attach_measure
    integrand = "0" if j.physical_zero else (
        f"({j.lorentzZ}) * ({j.jacobian}) * ({j.alpha_monomial})"
        f" * ({j.Q_prime_polynomial})^(-({j.denom_power}))")
    return attach_measure({"integrand": integrand, "expanded": None,
        "lorentzZ": j.lorentzZ, "jacobian": j.jacobian,
        "Q_prime": j.Q_prime_polynomial, "alpha_monomial": j.alpha_monomial,
        "denom_power": j.denom_power, "physical_zero": j.physical_zero,
        "variables": hc.x_vars, "bounds": {x: [0, 1] for x in hc.x_vars}})


def _build_c_mixed_block(c: CFactors, n_pure: int) -> dict:
    """Build C block for mixed dispatch."""
    i_factor = coordinate_phase(2*n_pure+3,n_pure+3,1,derivatives=1,chain_negatives=1)
    pure_ft_sign = f"(-1)^{n_pure}"  # one minus sign per pure Fourier propagator
    return {
        "components": {
            "color": c.color,
            "coupling": c.coupling,
            "ft_prefactor": c.ft_base,
            "derivative": c.derivative,
            "v14_scale": c.v14_scale,
            "i_factor": i_factor,
            "twoscale": c.twoscale,
        },
        "component_semantics": "multiplicative",
        "factor_order": ["i_factor", "derivative", "color", "coupling", "ft_prefactor", "v14_scale", "twoscale"],
        "phase_provenance": {"wilson_vertices": f"I^{2*n_pure+3}",
            "coordinate_propagators": f"(-1)^{n_pure+3}", "coordinate_centres": "-I",
            "momentum_derivative": "-I", "chain_rule": "-1", "product": i_factor},
        "combined": (
            f"({i_factor}) * ({c.derivative}) * ({c.color})"
            f" * ({c.coupling}) * ({c.ft_base})"
            f" * ({c.v14_scale}) * ({c.twoscale})"
        ),
    }


def _build_jraw_mixed(j: MixedJFactors, hc: HypercubeResult) -> dict:
    """Build J_raw block for mixed dispatch."""
    integrand = (
        f"({j.lorentz_z}) * ({j.jacobian}) * {j.f_product}"
        f" * ({j.alpha_monomial})"
        f" * ({j.pnum})"
        f" * ({j.S_polynomial})^(-({j.denom_power}))"
    )
    return {
        "pure": {
            "lorentz_z": j.lorentz_z,
            "f_product": j.f_product,
            "f_polynomials": j.f_polynomials,
            "edge_pairs": [[a, b] for a, b in j.edge_pairs],
        },
        "tgv": {
            "pnum": j.pnum,
            "alpha_monomial": j.alpha_monomial,
            "S_polynomial": j.S_polynomial,
            "denom_power": j.denom_power,
        },
        "integrand": integrand,
        "variables": hc.x_vars,
        "alpha_vars": ["alpha1", "alpha2", "alpha3"],
        "bounds": {x: [0, 1] for x in hc.x_vars},
        "delta_constraints": {},
    }


def _build_c_double_tgv_block(c: CFactors) -> dict:
    from montopar.double_tgv_output import build_c_block
    return build_c_block(c)


def _build_jraw_double_tgv(j: DoubleTGVJFactors, hc: HypercubeResult) -> dict:
    raise ValueError("double_tgv_export_requires_orbit_source_and_ordered_domain")


def emit_montopar_json(
    orbit: Orbit,
    c: CFactors,
    j: JFactors | TGVJFactors | FGVJFactors | MixedJFactors | DoubleTGVJFactors,
    hc: HypercubeResult,
    c_rotated: CFactors | None = None,
    j_rotated: JFactors | TGVJFactors | FGVJFactors | MixedJFactors | DoubleTGVJFactors | None = None,
    hc_rotated: HypercubeResult | None = None,
    j_raw_s: dict | None = None,
    j_raw_s_rotated: dict | None = None,
    ibp: TGVIBPResult | None = None,
    ibp_rotated: TGVIBPResult | None = None,
) -> dict:
    """Generate montopar.json — merged C×J result for workspace.

    Optionally includes a second C4-rotated variant (C_rotated, J_raw_rotated)
    and s-variable versions (J_raw_s, J_raw_s_rotated) before hypercube.
    """
    result: dict = {
        "orbit_id": orbit.orbit_id,
        "dispatch": orbit.dispatch,
        "geometry": build_geometry(),
    }

    if orbit.dispatch == "tgv":
        result["C"] = _build_c_tgv_block(c)
        result["J_raw"] = _build_jraw_tgv(j, hc)  # type: ignore[arg-type]
        if j_raw_s is not None:
            result["J_raw_s"] = j_raw_s
        if ibp is not None:
            result["J_ibp"] = _build_jibp_tgv(ibp, c)
        if c_rotated is not None and j_rotated is not None and hc_rotated is not None:
            result["C_rotated"] = _build_c_tgv_block(c_rotated)
            result["J_raw_rotated"] = _build_jraw_tgv(j_rotated, hc_rotated)  # type: ignore[arg-type]
            if j_raw_s_rotated is not None:
                result["J_raw_s_rotated"] = j_raw_s_rotated
            if ibp_rotated is not None:
                result["J_ibp_rotated"] = _build_jibp_tgv(ibp_rotated, c_rotated)
    elif orbit.dispatch == "fourv":
        result["C"] = _build_c_fgv_block(c)
        result["J_raw"] = _build_jraw_fgv(j, hc)  # type: ignore[arg-type]
        if c_rotated is not None and j_rotated is not None and hc_rotated is not None:
            result["C_rotated"] = _build_c_fgv_block(c_rotated)
            result["J_raw_rotated"] = _build_jraw_fgv(j_rotated, hc_rotated)  # type: ignore[arg-type]
    elif orbit.dispatch == "mixed":
        result["n_pure"] = j.n_pure  # type: ignore[union-attr]
        result["C"] = _build_c_mixed_block(c, j.n_pure)  # type: ignore[union-attr]
        result["J_raw"] = _build_jraw_mixed(j, hc)  # type: ignore[arg-type]
        if j_raw_s is not None:
            result["J_raw_s"] = j_raw_s
        if c_rotated is not None and j_rotated is not None and hc_rotated is not None:
            result["C_rotated"] = _build_c_mixed_block(c_rotated, j_rotated.n_pure)  # type: ignore[union-attr]
            result["J_raw_rotated"] = _build_jraw_mixed(j_rotated, hc_rotated)  # type: ignore[arg-type]
            if j_raw_s_rotated is not None:
                result["J_raw_s_rotated"] = j_raw_s_rotated
    elif orbit.dispatch == "double_tgv":
        from montopar.double_tgv_output import build_double_tgv_result, build_c_block
        from montopar.j_factor import expand_polynomial
        corrected = build_double_tgv_result(orbit)
        if build_c_block(c) != corrected["C"]:
            raise ValueError("double_tgv_C_must_use_corrected_normalization")
        for field, value in (("B5", j.B5_polynomial), ("U", j.U_polynomial),
                             ("pnum1", j.pnum1), ("pnum2", j.pnum2)):
            if expand_polynomial(value) != corrected["J_raw"][field]:
                raise ValueError("double_tgv_legacy_numerator_mismatch_" + field)
        return corrected

    else:
        result["n_gluons"] = orbit.n_gluons
        result["C"] = _build_c_pure(c, orbit.n_gluons)
        result["J_raw"] = _build_jraw_pure(j, hc)  # type: ignore[arg-type]
        if c_rotated is not None and j_rotated is not None and hc_rotated is not None:
            result["C_rotated"] = _build_c_pure(c_rotated, orbit.n_gluons)
            result["J_raw_rotated"] = _build_jraw_pure(j_rotated, hc_rotated)  # type: ignore[arg-type]

    return annotate_result(result, orbit.source)


def write_json(data: dict, path: str | Path) -> None:
    """Write a dict as JSON to the given path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
