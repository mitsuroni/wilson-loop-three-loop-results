"""Validated, single-representative double-TGV raw C/J exports."""
from copy import deepcopy
from dataclasses import asdict
from fractions import Fraction
from montopar.models import Orbit
from montopar.parser import parse_orbit
from montopar.hypercube import hypercube_transform
from montopar.mixed_ibp import ordered_chains
from montopar.double_tgv import compute_j_double_tgv, ALPHAS
from montopar.j_factor import expand_polynomial, substitute_s_to_x
from montopar.c_factor import compute_c_double_tgv

SCHEMA = "1.3"
CONSTRAINT = "+".join(ALPHAS)
MONOMIAL = "(" + "*".join(ALPHAS) + ")^(-eps)"

def alpha_simplex5():
    return {"kind": "simplex", "variables": list(ALPHAS),
            "constraint": CONSTRAINT + "=1", "nonnegative_variables": list(ALPHAS),
            "independent_variables": ALPHAS[:4],
            "substitution": {"alpha5": "1-alpha1-alpha2-alpha3-alpha4"},
            "bounds": {a: ["0", "1" + "".join("-" + b for b in ALPHAS[:i])]
                       for i, a in enumerate(ALPHAS[:4])},
            "iterated_order_outer_to_inner": ALPHAS[:4],
            "jacobian": "1", "weight_location": "integrand"}

def branch_prescription():
    return {"propagator": "-I/(p^2+I0)",
            "geometric_B5": "-Sum[i<j,alpha_i*alpha_j*c_ij*(Xi-Xj)^2]",
            "regulated_B5_before_simplex": "B5_dimensionful+I*eta*U*Sum[alpha_i]",
            "regulated_B5_on_simplex": "B5_dimensionful+I*eta*U",
            "eta": "positive, taken to zero from above",
            "dimensionless_kernel_in_starting_region": "B5+I0*U",
            "reference_power": "(-2*v14+I0)^(3*eps)",
            "parameter_starting_region": {"z1": ">0", "v14": "<0"},
            "power_convention": "boundary value inherited from common Feynman prescription; joint analytic continuation",
            "complex_continuation_verified": False}

def validate_source(orbit):
    if (orbit.dispatch != "double_tgv" or orbit.emission_type != "4e-3g+3g"
            or orbit.source.get("nloops") != 3):
        raise ValueError("expected_three_loop_double_tgv_source")
    legs = orbit.double_tgv_legs
    if len(legs) != 4 or len(orbit.vertices) != 4 or orbit.pairings:
        raise ValueError("expected_four_external_double_tgv_legs")
    if len({l.loop for l in legs}) != 4 or len({l.s_var for l in legs}) != 4:
        raise ValueError("duplicate_double_tgv_leg")
    vertices = {v.lorentz: v for v in orbit.vertices}
    external = [p for p in orbit.propagators if sum(mu in vertices for mu in p.lorentz) == 1]
    if len(external) != 4 or {p.loop for p in external} != {l.loop for l in legs}:
        raise ValueError("double_tgv_external_connectivity_mismatch")
    if len({v.color for v in orbit.vertices}) != 4:
        raise ValueError("double_tgv_external_colors_must_be_independent")
    for prop in external:
        index = next(i for i, mu in enumerate(prop.lorentz) if mu in vertices)
        if prop.color[index] != vertices[prop.lorentz[index]].color:
            raise ValueError("double_tgv_color_endpoint_mismatch")
    # Re-bind from the source to catch stale/reordered in-memory leg metadata.
    parsed = parse_orbit(orbit.source)
    if [asdict(l) for l in legs] != [asdict(l) for l in parsed.double_tgv_legs]:
        raise ValueError("double_tgv_leg_metadata_disagrees_with_source")
    if (orbit.orbit_id != parsed.orbit_id or orbit.emission_type != parsed.emission_type
            or orbit.dispatch != parsed.dispatch or orbit.bounds != parsed.bounds
            or orbit.vertices != parsed.vertices or orbit.propagators != parsed.propagators):
        raise ValueError("double_tgv_metadata_disagrees_with_source")
    chains = ordered_chains(orbit.bounds)
    edge = {l.s_var: l.edge for l in legs}
    if set(edge) != {b.var for b in orbit.bounds}:
        raise ValueError("double_tgv_endpoint_domain_mismatch")
    if any(len({edge[v] for v in ch}) != 1 for ch in chains):
        raise ValueError("double_tgv_order_crosses_edges")
    if len(chains) != len(set(edge.values())):
        raise ValueError("double_tgv_same_edge_requires_one_ordered_chain")
    chains.sort(key=lambda ch: edge[ch[0]])
    phases = {p["loop"]: p["expression"] for p in orbit.source["phases"]}
    for leg in legs:
        expected = "- (" + " + ".join([f"v{k}" for k in range(1, leg.edge)]
                                        + [f"{leg.s_var}*v{leg.edge}"]) + ")"
        if expand_polynomial(f"({phases[leg.loop]}) - ({expected})") != "0":
            raise ValueError("double_tgv_phase_geometry_mismatch")
    return chains

def color_data(orbit, chains=None):
    """Four-generator Fierz result, preserving both vertex orientations."""
    if chains is None:
        chains = validate_source(orbit)
    contour = [v for ch in chains for v in ch]
    legs = orbit.double_tgv_legs
    ranks = [contour.index(l.s_var) for l in legs]
    sign = 1
    pairs = []
    for a, b in (ranks[:2], ranks[2:]):
        if a > b:
            sign = -sign
        pairs.append(tuple(sorted((a, b))))
    key = tuple(sorted(pairs))
    ratios = {((0, 1), (2, 3)): Fraction(-1, 4),
              ((0, 2), (1, 3)): Fraction(0),
              ((0, 3), (1, 2)): Fraction(1, 4)}
    ratio = sign * ratios[key]
    factor = "0" if not ratio else f"({ratio})*CA^2*CF"
    return {"factor": factor, "ratio_to_CA2_CF": str(ratio),
            "contour_leg_order": [next(i+1 for i, l in enumerate(legs) if l.s_var == v) for v in contour],
            "oriented_vertex_pairs": [[1, 2], [3, 4]],
            "method": "exact fundamental SU(N) four-generator Fierz contraction",
            "normalization": "Tr[identity]/N=1; Tr[Ta*Tb]=delta_ab/2"}

def double_tgv_c4_symmetry(orbit, chains=None):
    if chains is None:
        chains = validate_source(orbit)
    slots = {v: i+1 for ch in chains for i, v in enumerate(ch)}
    ends = [(l.edge, slots[l.s_var]) for l in orbit.double_tgv_legs]
    def fingerprint(shift):
        points = [((edge-1+shift) % 4+1, slot) for edge, slot in ends]
        return tuple(sorted(tuple(sorted(pair)) for pair in (points[:2], points[2:])))
    base = fingerprint(0)
    stab = [k for k in range(4) if fingerprint(k) == base]
    return {"orbit_multiplicity": 4 // len(stab), "stabilizer_order": len(stab),
            "stabilizer_rotations": stab,
            "definition": "unoriented endpoint pairs and exchange of the two identical internal vertices; same-edge order retained"}

def build_c_block(c):
    components = {"color": c.color, "coupling": c.coupling, "ft_prefactor": c.ft_base,
                  "derivative": c.derivative, "twoscale": c.twoscale, "v14_scale": c.v14_scale}
    order = list(components)
    return {"components": components, "factor_order": order,
            "component_semantics": "multiplicative",
            "combined": " * ".join(f"({components[k]})" for k in order),
            "phase_provenance": {"five_coordinate_propagators": "(-1)^5",
                                 "two_Wick_rotations": "(-I)^2", "product": "1",
                                 "two_momentum_derivatives": "minus sign already in T1_hat and T2_hat"},
            "Gamma_recurrence": "(1-3*eps)*Gamma[1-3*eps]=Gamma[2-3*eps]"}

def _measure(block):
    block.update({"alpha_vars": list(ALPHAS), "alpha_measure": alpha_simplex5(),
                  "delta_constraints": {CONSTRAINT: "1"},
                  "integration_order": "path_parameters_outer_alpha_simplex_inner"})
    return block

def _block(core, variables, bounds, order, jacobian, zero, substitution=None):
    sub = substitution or {}
    poly = lambda value: expand_polynomial(substitute_s_to_x(value, sub)) if sub else expand_polynomial(value)
    fields = {k: poly(core[k]) for k in ("B5", "U", "pnum1", "pnum2", "T1_hat", "T2_hat", "K_plus", "K_minus")}
    bracket = (f"({fields['pnum1']}) * ({fields['B5']})^(-(3-3*eps))"
               f" + ({fields['pnum2']}) * ({fields['B5']})^(-(2-3*eps))")
    kernel = f"({jacobian}) * ({MONOMIAL}) * ({fields['U']})^(-(1+2*eps)) * ({bracket})"
    return _measure({**fields, "integrand": "0" if zero else kernel,
                     "geometric_integrand_before_color_screening": kernel,
                     "variables": list(variables), "bounds": bounds,
                     "iterated_order_outer_to_inner": list(order), "jacobian": jacobian,
                     "alpha_monomial": MONOMIAL, "denom_power_1": "3 - 3*eps",
                     "denom_power_2": "2 - 3*eps", "physical_zero": zero,
                     "numerator_convention": "pnum1=-(2-3eps)*T1_hat/8; pnum2=-T2_hat/4; D=4-2eps included"})

def build_double_tgv_result(orbit, *, ibp_enabled=False, forms_enabled=False, operator_enabled=False):
    chains = validate_source(orbit)
    color = color_data(orbit, chains)
    hc = hypercube_transform(orbit.bounds)
    core = compute_j_double_tgv(orbit.double_tgv_legs, tuple(l.edge for l in orbit.double_tgv_legs), hc.jacobian, hc.x_vars)
    zero = color["factor"] == "0"
    ordered = _block(core, [b.var for b in orbit.bounds],
                     {b.var: [b.low, b.high] for b in orbit.bounds},
                     [v for ch in chains for v in reversed(ch)], "1", zero)
    ordered.update({"representation": "ordered_s", "ordered_chains": chains})
    cube = _block(core, hc.x_vars, {v: [0, 1] for v in hc.x_vars}, hc.x_vars,
                  hc.jacobian, zero, hc.s_to_x)
    cube.update({"representation": "unit_cube", "s_to_x": dict(hc.s_to_x)})
    result = {"schema_version": SCHEMA, "orbit_id": orbit.orbit_id,
            "dispatch": "double_tgv", "emission_type": "4e-3g+3g", "loop_order": 3,
            "geometry": {"n_lines": 4, "z_definition": "z1=v12/v14", "reference_scale": "v14",
                         "lightlike_edges": True},
            "representative": {"scheme": "single_source_oriented_vertex_representative", "rotation": 0,
                               "leg_permutation": [1, 2, 3, 4],
                               "legs": [asdict(l) for l in orbit.double_tgv_legs],
                               "canonical_topology": deepcopy(orbit.source),
                               "source_topology": deepcopy(orbit.source),
                               "s_to_x": dict(hc.s_to_x), "jacobian": hc.jacobian},
            "color_contraction": color, "c4_symmetry": double_tgv_c4_symmetry(orbit, chains),
            "normalization": {"scheme": "montopar-raw-CJ", "loop_order": 3,
                              "color_projection": "none", "euler_compensation_count": 0,
                              "orbit_multiplicity_included": False, "crossing_assembled": False},
            "C": build_c_block(compute_c_double_tgv(color["factor"])),
            "J_raw_s": ordered, "J_raw": cube,
            "branch_prescription": branch_prescription(),
            "evaluation_status": "algebraic_zero" if zero else "not_evaluated",
            "zero_reason": "normalized_fundamental_color_trace_zero" if zero else None,
            "zero_policy": "effective J is zero after color screening; geometric kernel retained without evaluation",
            "formal": "not_run"}
    if ibp_enabled:
        from montopar.double_tgv_ibp import build_double_tgv_ibp
        result["J_ibp"] = build_double_tgv_ibp(orbit, result)
    if forms_enabled:
        from montopar.double_tgv_forms import build_compact, build_reduced
        result["J_compact"] = build_compact(orbit, result)
        result["J_reduced"] = build_reduced(orbit, result)
    if operator_enabled:
        from montopar.double_tgv_operator import build_operator
        result["J_operator"] = build_operator(orbit, result)
    from montopar.sign_rules import annotate_result
    return annotate_result(result, orbit.source)

def validate_double_tgv_result(result):
    if result.get("schema_version") != SCHEMA or result.get("dispatch") != "double_tgv":
        raise ValueError("regenerate_corrected_double_tgv_source")
    if any(k.endswith("_rotated") for k in result):
        raise ValueError("double_tgv_requires_single_representative")
    try:
        orbit = parse_orbit(result["representative"]["canonical_topology"])
        expected = build_double_tgv_result(orbit, ibp_enabled="J_ibp" in result,
                                          forms_enabled="J_compact" in result or "J_reduced" in result,
                                          operator_enabled="J_operator" in result)
    except (KeyError, TypeError) as error:
        raise ValueError("incomplete_double_tgv_source_semantics") from error
    for key in expected:
        if result.get(key) != expected[key]:
            raise ValueError("double_tgv_result_mismatch_" + key)
    if "J_operator" in result:
        from montopar.double_tgv_operator import to_scalar
        if to_scalar(result["J_operator"]) != result["J_raw_s"]:
            raise ValueError("double_tgv_operator_scalar_conversion_mismatch")
    return True
