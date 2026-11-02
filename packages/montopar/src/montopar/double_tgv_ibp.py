"""Two-layer contour IBP, with ordered contact-face normal corrections.

Implements theory-r003 sections 8--10 using exact polynomial arithmetic.
No CAS runtime, alpha division, or inverse external Gram matrix is used.
Equality is for the grouped regulated integrals, not pointwise integrands.
"""
from copy import deepcopy
from montopar.double_tgv import (
    ALPHAS, _add, _mul, _scale, _vadd, _dot, _tangent, _gradients,
)
from montopar.j_factor import expand_polynomial, substitute_s_to_x, differentiate_polynomial
from montopar.mixed_ibp import chain_bounds
from montopar.hypercube import hypercube_transform

ZERO = ["0"] * 4


def restrict(poly, mapping):
    """Simultaneous parent-to-chart substitution, never sequential identification."""
    return expand_polynomial(substitute_s_to_x(poly, mapping))


def ordered_faces(chains):
    """All codimension-one faces with inherited coordinate measure (Jacobian 1)."""
    variables = [v for chain in chains for v in chain]
    faces = []
    for chain in chains:
        specs = [("lower", chain[0], "0"), ("upper", chain[-1], "1")]
        specs += [("contact", a, b) for a, b in zip(chain, chain[1:])]
        for kind, removed, target in specs:
            normal = {v: "0" for v in variables}
            normal[removed] = "-1" if kind == "lower" else "1"
            if kind == "contact":
                normal[target] = "-1"
            remaining = [[v for v in c if v != removed] for c in chains]
            remaining = [c for c in remaining if c]
            faces.append({"face_id": f"{kind}_{removed}_{target}", "kind": kind,
                          "removed": removed, "target": target, "normal": normal,
                          "parent_to_face": {v: target if v == removed else v for v in variables},
                          "ordered_chains": remaining, "measure_jacobian": "1"})
    return faces


def recoil_vector(u, v, p, q):
    """R B0 / v14 in the four-edge vector basis; dot() means dot/v14."""
    total = _vadd(p, q)
    return _vadd(_scale(_dot(u, v), _vadd(q, _scale("-1", p))),
                 _scale(_mul("2", _dot(v, total)), u),
                 _scale(_mul("-2", _dot(u, total)), v))


def operator_data(legs, B):
    u = [_tangent(l) for l in legs]
    p = _gradients(legs)
    left = recoil_vector(u[0], u[1], p[0], p[1])
    right = recoil_vector(u[2], u[3], p[2], p[3])
    # GL=R_L Bhat = -left/2; the edge-basis dot supplies one v14.
    gl, gr = _scale("-1/2", left), _scale("-1/2", right)
    pair = _mul(_dot(u[0], u[1]), _dot(u[2], u[3]))
    cross = _add(_mul(_dot(u[0], u[2]), _dot(u[1], u[3])),
                 _mul("-1", _dot(u[0], u[3]), _dot(u[1], u[2])))
    h = _mul("1/2", ALPHAS[4], _add(
        _mul("4-2*eps", pair, "alpha1-alpha2", "alpha3-alpha4"),
        _mul("8", cross, "alpha1+alpha2", "alpha3+alpha4")))
    return {"G_L_edge_basis": gl, "G_R_edge_basis": gr,
            "G_L_dot_G_R_over_A0": _mul("-1/2", _dot(gl, gr)),
            "H_RR_over_A0": h,
            "dB_ds": {l.s_var: differentiate_polynomial(B, l.s_var) for l in legs},
            "left_recoil_zero": left == ZERO, "right_recoil_zero": right == ZERO}


def _term(kind, mapping, chains, numerator, power, ratio, origins, raw):
    from montopar.double_tgv_output import alpha_simplex5, MONOMIAL, CONSTRAINT
    bounds = chain_bounds(chains)
    hc = hypercube_transform(bounds)
    ordered = {b.var: [b.low, b.high] for b in bounds}
    B = restrict(raw["J_raw_s"]["B5"], mapping)
    U = raw["J_raw_s"]["U"]
    p = restrict(numerator, mapping)
    bx, px = restrict(B, hc.s_to_x), restrict(p, hc.s_to_x)

    def kernel(b, n, jac, regulated=False):
        denominator = f"({b})+I*eta*({U})" if regulated else b
        return (f"({jac}) * ({MONOMIAL}) * ({U})^(-1-2*eps)"
                f" * ({n}) * ({denominator})^(-({power}))")

    zero = raw["J_raw_s"]["physical_zero"]
    return {"kind": kind, "C_ref": "C_raw", "C_ratio_to_raw": ratio,
            "combined_C": f"({ratio}) * ({raw['C']['combined']})",
            "integrand": "0" if zero else kernel(bx, px, hc.jacobian),
            "integrand_s": "0" if zero else kernel(B, p, "1"),
            "geometric_integrand_s_before_color_screening": kernel(B, p, "1"),
            "regulated_integrand_s_before_color_screening": kernel(B, p, "1", True),
            "variables": hc.x_vars, "bounds": {v: [0, 1] for v in hc.x_vars},
            "iterated_order_outer_to_inner": hc.x_vars,
            "ordered_variables": list(ordered), "ordered_bounds": ordered,
            "ordered_chains": deepcopy(chains),
            "ordered_iterated_order_outer_to_inner": [v for c in chains for v in reversed(c)],
            "parent_to_face": mapping, "s_to_x": hc.s_to_x, "jacobian": hc.jacobian,
            "B5_s": B, "B5": bx, "U": U, "pnum_s": p, "pnum": px,
            "denom_power": power, "alpha_monomial": MONOMIAL,
            "alpha_vars": list(ALPHAS), "alpha_measure": alpha_simplex5(),
            "delta_constraints": {CONSTRAINT: "1"},
            "integration_order": "path_parameters_outer_alpha_simplex_inner",
            "origins": origins, "physical_zero": zero,
            "regulator_group": raw["orbit_id"] + ":all_ibp_terms"}


def build_double_tgv_ibp(orbit, raw):
    """Generate scalar term C/J pairs, plus exact parent-space provenance."""
    legs = orbit.double_tgv_legs
    chains = raw["J_raw_s"]["ordered_chains"]
    svars = [l.s_var for l in legs]
    data = operator_data(legs, raw["J_raw_s"]["B5"])
    u = [_tangent(l) for l in legs]
    avec = dict(zip(svars, [u[1], _scale("-1", u[0]), ZERO, ZERO]))
    bvec = dict(zip(svars, [ZERO, ZERO, u[3], _scale("-1", u[2])]))
    doubles, singles, face_records = {}, [], []
    for face in ordered_faces(chains):
        n = face["normal"]
        wl = _vadd(*(_scale(n[v], avec[v]) for v in svars))
        wr = _vadd(*(_scale(n[v], bvec[v]) for v in svars))
        # b is divided by v14. Projection retains parent components on free coordinates.
        b = {v: _dot(wl, bvec[v]) for v in svars}
        removed, target = face["removed"], face["target"]
        beta = {v: b[v] for v in svars if v != removed}
        chi = _add(b[removed], _mul("-1", b[target])) if face["kind"] == "contact" else b[removed]
        tr = _mul("-1/2", _dot(wl, data["G_R_edge_basis"]))
        rt = _mul("-1/2", _dot(wr, data["G_L_edge_basis"]))
        normal = _mul("-1/2", chi, data["dB_ds"][removed])
        numerator = _add(tr, rt, normal)
        record = {**face, "omega_L_edge_basis": wl, "omega_R_edge_basis": wr,
                  "b_over_v14": b, "beta_over_v14": beta, "chi_over_v14": chi,
                  "parent_numerators_over_A0": {"TR": tr, "RT": rt, "TT_normal": normal},
                  "combined_parent_numerator_over_A0": numerator}
        face_records.append(record)
        if restrict(numerator, face["parent_to_face"]) != "0":
            singles.append(_term("single_boundary", face["parent_to_face"], face["ordered_chains"],
                                 numerator, "2-3*eps", "1", [record], raw))
        for subface in ordered_faces(face["ordered_chains"]):
            k = _add(*(_mul(subface["normal"][v], beta[v]) for v in beta))
            if k == "0":
                continue
            mapping = {v: restrict(value, subface["parent_to_face"])
                       for v, value in face["parent_to_face"].items()}
            key = (tuple(sorted(mapping.items())), tuple(tuple(c) for c in subface["ordered_chains"]))
            origin = {"first_face": face["face_id"], "second_face": subface,
                      "k_over_v14": k, "k_over_A0": _mul("-1/2", k)}
            if key not in doubles:
                doubles[key] = {"mapping": mapping, "chains": subface["ordered_chains"],
                                "coefficient": "0", "origins": []}
            entry = doubles[key]
            entry["coefficient"] = _add(entry["coefficient"], k)
            entry["origins"].append(origin)
    terms = []
    for entry in doubles.values():
        if entry["coefficient"] == "0":
            continue
        ratio = f"({_mul('1/2', entry['coefficient'])})/(1-3*eps)"
        terms.append(_term("double_boundary", entry["mapping"], entry["chains"], "1",
                           "1-3*eps", ratio, entry["origins"], raw))
    terms.extend(singles)
    for kind, key, power, ratio in [
        ("double_recoil_gradient", "G_L_dot_G_R_over_A0", "3-3*eps", "-(2-3*eps)"),
        ("double_recoil_hessian", "H_RR_over_A0", "2-3*eps", "1"),
    ]:
        if data[key] != "0":
            terms.append(_term(kind, {v: v for v in svars}, chains, data[key], power, ratio,
                               [{"operator_block": "RR", "group": "double_recoil_body"}], raw))
    for i, term in enumerate(terms, 1):
        term["term_id"] = f"ibp_{i:03d}"
    return {"schema_version": "double-tgv-ibp-1.0", "representation": "ordered_two_layer_contour_ibp",
            "theory": "theory-r003 sections 8-10", "C_raw": deepcopy(raw["C"]),
            "relation": "C_raw*J_raw_s = Sum[term.combined_C*Integral[term.integrand_s]]",
            "operator_data": data, "faces": face_records, "terms": terms,
            "counts": {kind: sum(t["kind"] == kind for t in terms) for kind in
                       ("double_boundary", "single_boundary", "double_recoil_gradient", "double_recoil_hessian")},
            "normalization": deepcopy(raw["normalization"]),
            "regulator": {"dimensionless_eta": ">0", "kernel": "B5+I*eta*U",
                          "relation_to_dimensionful_kappa": "eta=kappa/A0; A0=-2*v14",
                          "identity_domain": "fixed interior alpha, finite eta, inherited ordered face measure",
                          "alpha_cutoff": "alpha_i>=delta>0 before removal",
                          "removal": "sum all terms first; common Feynman boundary and meromorphic continuation",
                          "termwise_convergence_established": False},
            "vector_convention": "four-edge basis; dot_product_factor=ui.uj/v14; A0=-2*v14",
            "differentiation_order": "differentiate parent B5 before face restriction",
            "physical_zero": raw["J_raw_s"]["physical_zero"],
            "evaluation_status": raw["evaluation_status"], "formal": "not_run"}
