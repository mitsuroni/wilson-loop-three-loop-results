"""Mixed TGV contour IBP on full ordered endpoint domains.

Only polynomial algebra is performed here. "complete" means representation
generated, not integration performed or formal proof admitted.
"""
from __future__ import annotations
from dataclasses import asdict
from copy import deepcopy
from montopar.models import Bound
from montopar.hypercube import hypercube_transform
from montopar.j_factor import (
    compute_f_p, compute_lorentz_factor, differentiate_polynomial,
    expand_polynomial, linear_combination_polynomials, substitute_s_to_x,
)
from montopar.tgv import compute_S, compute_Pnum
from montopar.tgv_ibp import (
    classify_tgv_ibp_case, degenerate_derivative_vector,
    nondegenerate_derivative_vector,
)

ALPHA_MONOMIAL = "(alpha1*alpha2*alpha3)^(-eps)"


def alpha_simplex():
    return {
        "kind": "simplex", "variables": ["alpha1", "alpha2", "alpha3"],
        "constraint": "alpha1+alpha2+alpha3=1",
        "independent_variables": ["alpha1", "alpha2"],
        "substitution": {"alpha3": "1-alpha1-alpha2"},
        "bounds": {"alpha1": ["0", "1"], "alpha2": ["0", "1-alpha1"]},
        "iterated_order_outer_to_inner": ["alpha1", "alpha2"],
        "jacobian": "1",
    }


def ordered_chains(bounds):
    """Accept only disjoint unit ordered chains; reject cycles and branches."""
    by = {b.var: b for b in bounds}
    if len(by) != len(bounds) or any(b.low != "0" for b in bounds):
        raise ValueError("expected_unique_zero_lower_ordered_bounds")
    incoming = {}
    for b in bounds:
        if b.high == "1":
            continue
        if b.high not in by or b.high in incoming or b.high == b.var:
            raise ValueError("invalid_ordered_chain")
        incoming[b.high] = b.var
    chains, seen = [], set()
    for first in by:
        if first in incoming:
            continue
        chain, var = [], first
        while var != "1":
            if var in seen:
                raise ValueError("cycle_or_shared_endpoint")
            seen.add(var); chain.append(var); var = by[var].high
        chains.append(chain)
    if seen != set(by):
        raise ValueError("cyclic_ordered_bounds")
    return chains


def chain_bounds(chains):
    return [Bound(v, "0", chain[i+1] if i+1 < len(chain) else "1")
            for chain in chains for i, v in enumerate(chain)]


def _lin(terms):
    return linear_combination_polynomials(terms)


def oriented_faces(bounds, vector):
    """Flux coefficients in inherited parameter measure (no sqrt(2))."""
    chains = ordered_chains(bounds)
    if set(vector) != {b.var for b in bounds}:
        raise ValueError("derivative_vector_domain_mismatch")
    result = []
    for index, chain in enumerate(chains):
        specs = [
            ("cusp_boundary", {chain[0]: "0"}, {}, _lin([("-1", vector[chain[0]])])),
            ("cusp_boundary", {chain[-1]: "1"}, {}, expand_polynomial(vector[chain[-1]])),
        ]
        specs += [
            ("contact_boundary", {}, {left: right},
             _lin([("1", vector[left]), ("-1", vector[right])]))
            for left, right in zip(chain, chain[1:])
        ]
        for kind, fixed, identified, coefficient in specs:
            if coefficient == "0":
                continue
            removed = set(fixed) | set(identified)
            remaining = [[v for v in c if v not in removed] for c in chains]
            remaining = [c for c in remaining if c]
            b = chain_bounds(remaining)
            label = next(iter(fixed.items()), None) or next(iter(identified.items()))
            result.append({
                "term_id": f"boundary_{label[0]}_{label[1]}",
                "kind": kind, "coefficient": coefficient,
                "fixed_parameters": fixed, "identified_parameters": identified,
                "ordered_chains": remaining,
                "ordered_bounds": {x.var: [x.low, x.high] for x in b},
                "iterated_order_outer_to_inner": [v for c in remaining for v in reversed(c)],
                "origin_chain": index,
            })
    return result


def _polynomial(expr, sub):
    return expand_polynomial(substitute_s_to_x(expr, sub)) if sub else expand_polynomial(expr)


def _local_term(orbit, face, S, numerator, pure, c_raw):
    constraints = {**face["identified_parameters"], **face["fixed_parameters"]}
    bounds = [Bound(v, *b) for v, b in face["ordered_bounds"].items()]
    hc = hypercube_transform(bounds)
    S_s = _polynomial(S, constraints)
    P_s = _polynomial(numerator, constraints)
    fs_s = {f"f{i}": _polynomial(compute_f_p(p.edge_a, p.edge_b, p.s_a, p.s_b), constraints)
            for i, p in enumerate(pure, 1)}
    S_x, P_x = _polynomial(S_s, hc.s_to_x), _polynomial(P_s, hc.s_to_x)
    fs_x = {k: _polynomial(v, hc.s_to_x) for k, v in fs_s.items()}
    is_boundary = face["kind"] != "recoil"
    power = "1-2*eps" if is_boundary else "2-2*eps"
    ratio = f"({face['coefficient']})/(4*(1-2*eps))" if is_boundary else "1"
    lorentz = compute_lorentz_factor(pure)
    def integrand(S0, P0, fs, jac):
        fp = " * ".join(f"({v})^(-1+eps)" for v in fs.values()) or "1"
        return f"({jac}) * ({lorentz}) * ({fp}) * ({ALPHA_MONOMIAL}) * ({P0}) * ({S0})^(-({power}))"
    return {
        **face, "C_ref": "C_star" if is_boundary else "C_raw",
        "C_ratio_to_raw": ratio, "combined_C": f"({ratio}) * ({c_raw})",
        "integrand": integrand(S_x, P_x, fs_x, hc.jacobian),
        "integrand_s": integrand(S_s, P_s, fs_s, "1"),
        "variables": hc.x_vars, "bounds": {v: [0, 1] for v in hc.x_vars},
        "s_to_x": hc.s_to_x, "jacobian": hc.jacobian,
        "alpha_vars": ["alpha1", "alpha2", "alpha3"],
        "alpha_measure": alpha_simplex(),
        "delta_constraints": {"alpha1+alpha2+alpha3": "1"},
        "denom_power": power, "S_polynomial": S_x, "pnum": P_x,
        "f_polynomials": fs_x, "S_polynomial_s": S_s,
        "pnum_s": P_s, "f_polynomials_s": fs_s, "lorentz_z": lorentz,
    }


def compute_mixed_ibp(orbit, c_raw_combined):
    if orbit.dispatch != "mixed":
        raise ValueError("expected_mixed_dispatch")
    case = classify_tgv_ibp_case(orbit.tgv_legs)
    if case in ("triple", "unsupported"):
        raise ValueError(f"unsupported_tgv_case:{case}")
    chains = ordered_chains(orbit.bounds)
    coeffs = (nondegenerate_derivative_vector(orbit.tgv_legs)
              if case == "nondegenerate" else degenerate_derivative_vector(orbit.tgv_legs))
    vector = {b.var: "0" for b in orbit.bounds}
    for leg, c in zip(orbit.tgv_legs, coeffs):
        if leg.s_var not in vector:
            raise ValueError("missing_tgv_endpoint")
        vector[leg.s_var] = c
    tgv_loops = {leg.loop for leg in orbit.tgv_legs}
    pure = [p for p in orbit.pairings if p.loop not in tgv_loops]
    endpoints = [t.s_var for t in orbit.tgv_legs] + [v for p in pure for v in (p.s_a, p.s_b)]
    if len(set(endpoints)) != len(endpoints) or set(endpoints) != set(vector):
        raise ValueError("endpoints_not_distinct_or_incomplete")
    endpoint_edges = {t.s_var: t.edge for t in orbit.tgv_legs}
    endpoint_edges.update({v: edge for p in pure for v, edge in ((p.s_a,p.edge_a),(p.s_b,p.edge_b))})
    if any(len({endpoint_edges[v] for v in ch}) != 1 for ch in chains):
        raise ValueError("ordered_chain_crosses_edges")
    if len(chains) != len({endpoint_edges[ch[0]] for ch in chains}):
        raise ValueError("same_edge_split_into_independent_chains")
    S, P = compute_S(orbit.tgv_legs), compute_Pnum(orbit.tgv_legs)
    Pb = _lin([(f"-1/4*({c})", differentiate_polynomial(S, v)) for v,c in vector.items()])
    Pr = _lin([("1", P), ("-1", Pb)])
    if case.startswith("degenerate") and Pr != "0":
        raise ValueError("degenerate_numerator_identity_failed")
    # The first original-s IBP differentiates only distinct TGV endpoints.
    for p in pure:
        f = compute_f_p(p.edge_a,p.edge_b,p.s_a,p.s_b)
        if _lin([(c,differentiate_polynomial(f,v)) for v,c in vector.items()]) != "0":
            raise ValueError("pure_weight_has_nonzero_directional_derivative")
    faces = oriented_faces(orbit.bounds, vector)
    terms = [_local_term(orbit, f, S, "1", pure, c_raw_combined) for f in faces]
    if Pr != "0":
        face = {
            "term_id": "recoil", "kind": "recoil", "coefficient": "1",
            "fixed_parameters": {}, "identified_parameters": {},
            "ordered_chains": chains,
            "ordered_bounds": {b.var:[b.low,b.high] for b in orbit.bounds},
            "iterated_order_outer_to_inner": [v for ch in chains for v in reversed(ch)],
        }
        terms.append(_local_term(orbit,face,S,Pr,pure,c_raw_combined))
    return {
        "schema_version": "1.1", "method": "contour_ibp_full_ordered_domain",
        "status": "complete", "evaluation_status": "not_evaluated", "formal": "not_run",
        "case": case, "prefactors": {
            "C_raw": c_raw_combined, "C_star": f"({c_raw_combined})/(4*(1-2*eps))",
            "relation": "C_star=C_raw/(4*(1-2*eps))", "euler_compensation_count": 0,
        },
        "operator": {"derivative_vector_by_parameter": vector, "ordered_chains": chains,
                     "pnum_boundary": Pb, "pnum_recoil": Pr},
        "terms": terms, "sum_expression": " + ".join(f"term({t['term_id']})" for t in terms),
        "validity": {
            "initial_real_domain": "z1>0 and 1/2<Re(eps)<1",
            "continuation": "meromorphic continuation from common absolute-convergence strip",
            "complex_physical_branch_proved": False,
        },
        "algebra_checks": {"degenerate_recoil_zero": Pr=="0",
                          "pure_weight_directional_derivative_zero": True},
    }


def canonical_metadata(original, canonical, hc):
    """Bind exported expressions to actual topology, not stale source topology."""
    from montopar.rotation import rotate_edge, rotate_s_var
    ks = [k for k in range(4)
          if [(rotate_edge(v.line,k),v.lorentz,v.color) for v in original.vertices]
          == [(v.line,v.lorentz,v.color) for v in canonical.vertices]]
    if not ks:
        raise ValueError("canonical_rotation_not_reconstructible")
    k = ks[0]
    old_loops = [t.loop for t in original.tgv_legs]
    permutation = [old_loops.index(t.loop)+1 for t in canonical.tgv_legs]
    return {
        "policy": "single_canonical_representative",
        "canonical_topology": {
            "tgv_legs": [asdict(t) for t in canonical.tgv_legs],
            "pure_pairings": [asdict(p) for p in canonical.pairings],
            "vertices": [asdict(v) for v in canonical.vertices],
            "propagators": [asdict(p) for p in canonical.propagators],
            "bounds": {b.var:[b.low,b.high] for b in canonical.bounds},
        },
        "transformation": {
            "c4_rotation": k, "s_renaming": {b.var:rotate_s_var(b.var,k) for b in original.bounds},
            "tgv_leg_permutation_new_to_old": permutation,
            "alpha_renaming_old_to_new": {f"alpha{old}":f"alpha{i}" for i,old in enumerate(permutation,1)},
            "color_recomputed_after_permutation": True,
            "z_relation": "z_canonical=1/z_source" if k%2 else "z_canonical=z_source",
            "reference_scale_relation": "v14_canonical=z_source*v14_source" if k%2 else "v14_canonical=v14_source",
        },
        "hypercube": {"s_to_x":dict(hc.s_to_x),"jacobian":hc.jacobian},
    }


def finalize_single_representative(result, original, canonical, hc):
    """Explicit source semantics; no Euler preparation and no duplicated B."""
    result = deepcopy(result)
    for key in list(result):
        if key.endswith("_rotated"):
            del result[key]
    result["schema_version"] = "1.1"
    result["representative"] = canonical_metadata(original,canonical,hc)
    result["normalization"] = {"stage":"raw_montopar", "euler_compensation_count":0}
    for key in ("J_raw","J_raw_s"):
        if key in result:
            result[key]["alpha_measure"] = alpha_simplex()
            result[key]["delta_constraints"] = {"alpha1+alpha2+alpha3":"1"}
    if "J_raw_s" in result:
        ch = ordered_chains(canonical.bounds)
        result["J_raw_s"]["iterated_order_outer_to_inner"] = [v for c in ch for v in reversed(c)]
    if "J_ibp" in result:
        for term in result["J_ibp"]["terms"]:
            term["alpha_measure"] = alpha_simplex()
            term["delta_constraints"] = {"alpha1+alpha2+alpha3":"1"}
    return result
