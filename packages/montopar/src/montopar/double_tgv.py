"""Double TGV parameterization from the corrected two-vertex contraction.

Raw convention: B5=B_dimensionful/(-2*v14), T1=T1_dimensionful/v14^3,
T2=T2_dimensionful/v14^2.  Pnum1=-(2-3eps)*T1/8, Pnum2=-T2/4.
D=4-2eps is retained in the mixed Hessian.  No CAS runtime is required.
"""
from montopar.models import DoubleTGVLeg
from montopar.geometry import dot_product_factor, compute_Delta_prime
from montopar.j_factor import expand_polynomial

ALPHAS = [f"alpha{i}" for i in range(1, 6)]

def _add(*terms):
    values = [f"({t})" for t in terms if t != "0"]
    return expand_polynomial(" + ".join(values)) if values else "0"

def _mul(*terms):
    if "0" in terms:
        return "0"
    return expand_polynomial(" * ".join(f"({t})" for t in terms))

def _scale(c, v):
    return [_mul(c, x) for x in v]

def _vadd(*vectors):
    return [_add(*(v[i] for v in vectors)) for i in range(4)]

def _dot(v, w):
    return _add(*(_mul(v[i], dot_product_factor(i+1, j+1), w[j])
                  for i in range(4) for j in range(4)))

def _point(leg):
    return ["1" if k < leg.edge-1 else leg.s_var if k == leg.edge-1 else "0"
            for k in range(4)]

def _tangent(leg):
    return ["1" if k == leg.edge-1 else "0" for k in range(4)]

def _c(i, j, av=ALPHAS):
    if i < 2 and j < 2:
        return _add(av[2], av[3], av[4])
    if i >= 2 and j >= 2:
        return _add(av[0], av[1], av[4])
    return av[4]

def _validate_legs(legs):
    if len(legs) != 4 or len({l.loop for l in legs}) != 4 or len({l.s_var for l in legs}) != 4:
        raise ValueError("double_tgv_requires_four_distinct_external_legs")
    if any(l.edge not in (1, 2, 3, 4) for l in legs):
        raise ValueError("double_tgv_requires_quadrilateral_edges")

def compute_U():
    return "(alpha1+alpha2+alpha5)*(alpha3+alpha4+alpha5) - alpha5^2"

def compute_B5(legs):
    _validate_legs(legs)
    return _add(*(_mul(ALPHAS[i], ALPHAS[j], _c(i, j),
                       compute_Delta_prime(legs[i].edge, legs[j].edge,
                                           legs[i].s_var, legs[j].s_var))
                  for i in range(4) for j in range(i+1, 4)))

def _gradients(legs, av=ALPHAS):
    positions = [_point(l) for l in legs]
    return [_scale(_mul("-2", av[i]), _vadd(*(
        _scale(_mul(av[j], _c(i, j, av)),
               _vadd(positions[i], _scale("-1", positions[j])))
        for j in range(4) if j != i))) for i in range(4)]

def v_dot_X_diff_4leg(p_edge, leg_i, leg_j):
    unit = ["1" if k == p_edge-1 else "0" for k in range(4)]
    return _dot(unit, _vadd(_point(leg_i), _scale("-1", _point(leg_j))))

def compute_F_i_p(i, p_edge, legs, alpha_vars):
    """Dimensional directional derivative divided by v_ref (legacy API)."""
    _validate_legs(legs)
    return _dot(["1" if k == p_edge-1 else "0" for k in range(4)],
                _gradients(legs, alpha_vars)[i-1])

def compute_H_ij_pq(i, j, p_edge, q_edge, alpha_vars):
    """Dimensional directional Hessian divided by v_ref, not d^2 Bhat."""
    i -= 1
    j -= 1
    dot = dot_product_factor(p_edge, q_edge)
    if i != j:
        return _mul("2", alpha_vars[i], alpha_vars[j], _c(i, j, alpha_vars), dot)
    return _mul("-2", alpha_vars[i], dot,
                _add(*(_mul(alpha_vars[k], _c(i, k, alpha_vars))
                       for k in range(4) if k != i)))

def _vertex(u, v, p, q):
    return _vadd(_scale(_dot(u, v), _vadd(q, _scale("-1", p))),
                 _scale(_add(_mul("2", _dot(v, p)), _dot(v, q)), u),
                 _scale(_mul("-1", _add(_dot(u, p), _mul("2", _dot(u, q)))), v))

def compute_numerators(legs, alpha_vars=None):
    """Complete T1/T2 via the vertex product and the D-dimensional Hessian."""
    _validate_legs(legs)
    av = list(alpha_vars) if alpha_vars is not None else ALPHAS
    if len(av) != 5 or len(set(av)) != 5:
        raise ValueError("double_tgv_requires_five_alpha_parameters")
    g = _gradients(legs, av)
    u = [_tangent(l) for l in legs]
    left = _vertex(u[0], u[1], g[0], g[1])
    right = _vertex(u[2], u[3], g[2], g[3])
    t1 = _mul("-1", _dot(left, right))
    pair = _mul(_dot(u[0], u[1]), _dot(u[2], u[3]))
    cross = _mul(_dot(u[0], u[2]), _dot(u[1], u[3]))
    nested = _mul(_dot(u[0], u[3]), _dot(u[1], u[2]))
    kp = _add(_mul("2-2*eps", pair), _mul("5", cross), _mul("-4", nested))
    km = _add(_mul("-2+2*eps", pair), _mul("4", cross), _mul("-5", nested))
    t2 = _mul("-2", av[4], _add(
        _mul(_add(_mul(av[0], av[2]), _mul(av[1], av[3])), kp),
        _mul(_add(_mul(av[0], av[3]), _mul(av[1], av[2])), km)))
    return {"T1_hat": t1, "T2_hat": t2, "K_plus": kp, "K_minus": km,
            "pnum1": _mul("-1/8", "2-3*eps", t1),
            "pnum2": _mul("-1/4", t2)}

def compute_Pnum(legs, edges, alpha_vars=None):
    """Return the two final dimensionless numerators, including all eps factors."""
    if tuple(edges) != tuple(l.edge for l in legs):
        raise ValueError("double_tgv_edge_order_disagrees_with_legs")
    result = compute_numerators(legs, alpha_vars)
    return result["pnum1"], result["pnum2"]

def compute_j_double_tgv(legs, edges, jacobian, x_vars):
    if tuple(edges) != tuple(l.edge for l in legs):
        raise ValueError("double_tgv_edge_order_disagrees_with_legs")
    return {"B5": compute_B5(legs), "U": compute_U(),
            **compute_numerators(legs), "jacobian": jacobian, "x_vars": list(x_vars),
            "alpha_vars": list(ALPHAS),
            "alpha_monomial": "(alpha1*alpha2*alpha3*alpha4*alpha5)^(-eps)",
            "denom_power_1": "3 - 3*eps", "denom_power_2": "2 - 3*eps"}