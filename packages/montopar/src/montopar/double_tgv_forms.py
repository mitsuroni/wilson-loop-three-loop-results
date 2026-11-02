"""Explicit compact and analytically reduced double-TGV input forms.

No CAS runtime: definitions use the existing expression-string convention.
The raw exporter remains the compatibility baseline. Reductions are analytic
boundary-value identities, not identities at a fixed auxiliary regulator.
"""
from copy import deepcopy
import re

from montopar.double_tgv import ALPHAS, _point, _validate_legs
from montopar.geometry import compute_Delta_prime, dot_product_factor
from montopar.j_factor import expand_polynomial

PREFIX = "dtg_"
NAME = re.compile(r"\b[A-Za-z_][A-Za-z0-9_]*\b")


def _add(*values):
    values = [v for v in values if v != "0"]
    return " + ".join(f"({v})" for v in values) if len(values) > 1 else values[0] if values else "0"


def _mul(*values):
    if "0" in values:
        return "0"
    values = [v for v in values if v != "1"]
    return "*".join(f"({v})" for v in values) if len(values) > 1 else values[0] if values else "1"


def _sub(a, b):
    return "0" if a == b else _add(a, _mul("-1", b))


class _Definitions:
    """Ordered scalar abbreviations; every definition has an explicit RHS."""
    def __init__(self):
        self.items = []

    def bind(self, name, expression):
        if expression in ("0", "1", "-1"):
            return expression
        symbol = PREFIX + name
        self.items.append({"symbol": symbol, "expression": expression})
        return symbol


def polynomial_fields(block):
    """Expand compact scalar fields using Montopar's existing polynomial parser.

    Useful to consumers and checks; no eval, implicit derivative or CAS needed.
    Reject duplicate/forward definitions rather than silently changing scope.
    """
    values = {}
    external = set(block.get("variables", []))
    for item in block["definitions"]:
        symbol, expression = item["symbol"], item["expression"]
        if not re.fullmatch(r"dtg_[A-Za-z0-9_]+", symbol) or symbol in values:
            raise ValueError("invalid_double_tgv_definition")
        unresolved = {x for x in NAME.findall(expression) if x.startswith(PREFIX)} - values.keys() - external
        if unresolved:
            raise ValueError("double_tgv_forward_or_missing_definition")
        expression = NAME.sub(lambda m: "(" + values[m[0]] + ")" if m[0] in values else m[0], expression)
        values[symbol] = expand_polynomial(expression)
    out = {}
    for key, expression in block["fields"].items():
        if {x for x in NAME.findall(expression) if x.startswith(PREFIX)} - values.keys() - external:
            raise ValueError("double_tgv_missing_field_definition")
        expression = NAME.sub(lambda m: "(" + values[m[0]] + ")" if m[0] in values else m[0], expression)
        out[key] = expand_polynomial(expression)
    return out


def _core(legs, alpha=None):
    """Build factored gradients with their exact degree-one alpha5 splitting."""
    _validate_legs(legs)
    if any(l.s_var.startswith(PREFIX) for l in legs):
        raise ValueError("double_tgv_reserved_form_variable")
    defs = _Definitions()
    aa = list(ALPHAS) if alpha is None else list(alpha)
    a1,a2,a3,a4,c = aa
    L, R = defs.bind("L", _add(a1,a2)), defs.bind("R", _add(a3,a4))
    U = defs.bind("U", _add(_mul(L,R),_mul(c,_add(L,R))))
    points = [_point(l) for l in legs]
    distances = {}
    Bterms = []
    for i in range(4):
        for j in range(i+1,4):
            dij = defs.bind(f"d{i+1}{j+1}", compute_Delta_prime(legs[i].edge,legs[j].edge,legs[i].s_var,legs[j].s_var))
            distances[i,j] = dij
            cij = _add(R,c) if j < 2 else _add(L,c) if i >= 2 else c
            Bterms.append(_mul(aa[i],aa[j],cij,dij))
    B = defs.bind("B",_add(*Bterms))
    # g_i = g_i0 + c*g_i1; g is the dimensionful gradient divided by vref.
    gradients = []
    for degree in (0,1):
        vectors = []
        for i in range(4):
            components = []
            for k in range(4):
                terms = []
                for j in range(4):
                    if i == j:
                        continue
                    same_vertex = (i < 2) == (j < 2)
                    coefficient = "1" if degree == 1 else (R if i < 2 else L) if same_vertex else "0"
                    terms.append(_mul(aa[j],coefficient,_sub(points[i][k],points[j][k])))
                components.append(defs.bind(f"g{degree}_{i+1}_{k+1}",_mul("-2",aa[i],_add(*terms))))
            vectors.append(components)
        gradients.append(vectors)

    def dot_vec(x,y):
        return _add(*(_mul(x[i],dot_product_factor(i+1,j+1),y[j]) for i in range(4) for j in range(4)))

    def dot_edge(edge,x):
        return _add(*(_mul(dot_product_factor(edge,k+1),x[k]) for k in range(4)))

    vertices = []
    for degree in (0,1):
        pair = []
        for side,(i,j) in enumerate(((0,1),(2,3))):
            gi,gj = gradients[degree][i],gradients[degree][j]
            ui,uj = legs[i].edge,legs[j].edge
            p = dot_product_factor(ui,uj)
            A = defs.bind(f"v{side}{degree}_a",_add(_mul("2",dot_edge(uj,gi)),dot_edge(uj,gj)))
            Bv = defs.bind(f"v{side}{degree}_b",_add(dot_edge(ui,gi),_mul("2",dot_edge(ui,gj))))
            pair.append([defs.bind(f"v{side}{degree}_{k+1}",_add(
                _mul(p,_sub(gj[k],gi[k])), A if ui == k+1 else "0", _mul("-1",Bv) if uj == k+1 else "0")) for k in range(4)])
        vertices.append(pair)
    vL0,vR0 = vertices[0]
    vL1,vR1 = vertices[1]
    N0 = defs.bind("N0",_mul("1/8","2-3*eps",dot_vec(vL0,vR0)))
    N1 = defs.bind("N1",_mul("1/8","2-3*eps",_add(dot_vec(vL0,vR1),dot_vec(vL1,vR0))))
    N2 = defs.bind("N2",_mul("1/8","2-3*eps",dot_vec(vL1,vR1)))
    # Also expose T1 without dividing by (2-3eps), so special eps has no fake pole.
    T1 = defs.bind("T1",_mul("-1",_add(dot_vec(vL0,vR0),
          _mul(c,_add(dot_vec(vL0,vR1),dot_vec(vL1,vR0))),_mul(f"({c})^2",dot_vec(vL1,vR1)))))
    p = lambda i,j: dot_product_factor(legs[i].edge,legs[j].edge)
    X,Y = _mul(p(0,2),p(1,3)),_mul(p(0,3),p(1,2))
    H0 = defs.bind("H0",_add(
        _mul(_sub(a2,a1),_sub(a4,a3),_add(_mul("2","2-2*eps",p(0,1),p(2,3)),X,Y)),
        _mul("9",L,R,_sub(X,Y))))
    M1 = defs.bind("M1",_mul("1/4",H0))
    fields = {"B5":B,"U":U,"T1_hat":T1,"T2_hat":_mul("-1",c,H0),
              "pnum1":_add(N0,_mul(c,N1),_mul(f"({c})^2",N2)),"pnum2":_mul(c,M1)}
    return defs.items, fields, (N0,N1,N2,M1), distances


def _normalization(raw):
    return {"relative_to":"C * J_raw_s", "color_projection":"none",
            "euler_compensation_count":0,"orbit_multiplicity_included":False,
            "crossing_assembled":False,"formal":"not_run",
            "physical_zero":raw["J_raw_s"]["physical_zero"]}


def build_compact(orbit,raw):
    definitions,fields,_,_ = _core(orbit.double_tgv_legs)
    zero = raw["J_raw_s"]["physical_zero"]
    B,U = fields["B5"],fields["U"]
    weight = f"({'*'.join(ALPHAS)})^(-eps)*({U})^(-1-2*eps)"
    def density(b):
        return f"({weight})*(({fields['pnum1']})*({b})^(-3+3*eps)+({fields['pnum2']})*({b})^(-2+3*eps))"
    out = {"schema_version":"double-tgv-compact-1.0","representation":"ordered_s_factored",
           "definitions":definitions,"definition_semantics":"ordered explicit scalar substitution",
           "fields":fields,"combined_C":raw["C"]["combined"],"relative_C":"1",
           "integrand":"0" if zero else density(B),
           "geometric_integrand_before_color_screening":density(B),
           "regulated_geometric_integrand":density(f"{B}+I*eta_reg*{U}"),
           "eta_reg":"dimensionless positive coordinate regulator, not alpha5; take common boundary value",
           "normalization":_normalization(raw),"equivalence":"pointwise after explicit substitution, including the declared finite regulator"}
    for key in ("variables","bounds","iterated_order_outer_to_inner","ordered_chains","alpha_vars","alpha_measure","delta_constraints","integration_order","jacobian"):
        out[key] = deepcopy(raw["J_raw_s"][key])
    return out


def _term(raw,term_id,ratio,kernel,variables,bounds,order,source):
    zero = raw["J_raw_s"]["physical_zero"]
    return {"term_id":term_id,"relative_C":ratio,
            "combined_C":f"({raw['C']['combined']}) * ({ratio})",
            "integrand":"0" if zero else kernel,"geometric_integrand_before_color_screening":kernel,
            "variables":list(variables),"bounds":deepcopy(bounds),
            "iterated_order_outer_to_inner":list(order),"jacobian":"included in integrand",
            "source":source}


def _paired_edges(orbit,raw):
    legs = orbit.double_tgv_legs
    if not (legs[0].edge == legs[1].edge and legs[2].edge == legs[3].edge and legs[0].edge != legs[2].edge):
        return False
    chains = raw["J_raw_s"]["ordered_chains"]
    return len(chains) == 2 and {frozenset(ch) for ch in chains} == {
        frozenset(l.s_var for l in legs[:2]),frozenset(l.s_var for l in legs[2:])}


def _same_edge_reduced(orbit,raw):
    legs = orbit.double_tgv_legs
    variables = ["dtg_rL","dtg_rR","dtg_x","dtg_y"]
    bounds = {v:["0","1"] for v in variables}
    chains = raw["J_raw_s"]["ordered_chains"]
    pairs = []
    orientation = 1
    for side,group in enumerate((legs[:2],legs[2:])):
        chain = next(ch for ch in chains if group[0].s_var in ch)
        orientation *= 1 if chain[0] == group[0].s_var else -1
        r = variables[side]
        pairs.append([(2,{chain[0]:r,chain[1]:r},"contact"),
                      (-1,{chain[0]:"0",chain[1]:r},"start"),
                      (-1,{chain[0]:r,chain[1]:"1"},"end")])
    d = dot_product_factor(legs[0].edge,legs[2].edge)
    terms = []
    for i,(cl,ml,fl) in enumerate(pairs[0]):
        for j,(cr,mr,fr) in enumerate(pairs[1]):
            m = {**ml,**mr}
            meanL = f"(dtg_x*({m[legs[0].s_var]})+(1-dtg_x)*({m[legs[1].s_var]}))"
            meanR = f"(dtg_y*({m[legs[2].s_var]})+(1-dtg_y)*({m[legs[3].s_var]}))"
            distance = compute_Delta_prime(legs[0].edge,legs[2].edge,meanL,meanR)
            ratio = f"({orientation*cl*cr})*({d})*Beta[eps,eps]/(4*eps*(1-3*eps))"
            kernel = f"(dtg_x*(1-dtg_x)*dtg_y*(1-dtg_y))^(-eps)*({distance})^(-1+3*eps)"
            term = _term(raw,f"pair-{i}{j}",ratio,kernel,variables,bounds,variables,
                         {"operation":"two ordered contour boundaries and two Beta parameter integrations",
                          "left_face":fl,"right_face":fr,"parent_to_face":m,"orientation":orientation,
                          "boundary_coefficient":cl*cr})
            term["distance"] = distance
            terms.append(term)
    return {"method":"paired_same_edge_beta","definitions":[],"terms":terms,
            "dimension_before":8,"dimension_after":4,
            "applicability":{"two_same_vertex_same_edge_pairs":True,"independent_ordered_pairs":True},
            "parameter_factor":"Beta[eps,eps]/(2*eps)",
            "operator_to_raw_ratio":f"({d})/(2*(1-3*eps))",
            "starting_region":"z1>0, v14<0, 0<Re(eps)<1/3; then continue the complete boundary combination"}


def _central_reduced(orbit,raw):
    aa = ["dtg_a*dtg_r","dtg_a*(1-dtg_r)","(1-dtg_a)*dtg_t","(1-dtg_a)*(1-dtg_t)","dtg_c"]
    definitions,fields,coefficients,distances = _core(orbit.double_tgv_legs,aa)
    # Discard c-dependent abbreviations: N0/N1/N2/M1 and distances are c-free.
    needed = set(NAME.findall(" ".join([*coefficients,*distances.values()])))
    retained = []
    for item in reversed(definitions):
        if item["symbol"] in needed:
            retained.append(item)
            needed.update(NAME.findall(item["expression"]))
    definitions = list(reversed(retained))
    d = distances
    definitions += [
        {"symbol":"dtg_h","expression":"dtg_a*(1-dtg_a)"},
        {"symbol":"dtg_massL","expression":f"dtg_r*(1-dtg_r)*({d[0,1]})"},
        {"symbol":"dtg_massR","expression":f"dtg_t*(1-dtg_t)*({d[2,3]})"},
        {"symbol":"dtg_massD","expression":f"dtg_r*dtg_t*({d[0,2]})+dtg_r*(1-dtg_t)*({d[0,3]})+(1-dtg_r)*dtg_t*({d[1,2]})+(1-dtg_r)*(1-dtg_t)*({d[1,3]})"},
        {"symbol":"dtg_Q0","expression":"dtg_a*dtg_massL+(1-dtg_a)*dtg_massR"},
        {"symbol":"dtg_Q1","expression":"dtg_a^2*dtg_massL+(1-dtg_a)^2*dtg_massR+dtg_h*dtg_massD"}]
    variables = raw["J_raw_s"]["variables"]+["dtg_a","dtg_r","dtg_t"]
    bounds = {**deepcopy(raw["J_raw_s"]["bounds"]),**{v:["0","1"] for v in variables[-3:]}}
    order = raw["J_raw_s"]["iterated_order_outer_to_inner"]+variables[-3:]
    terms = []
    for (n,m),coefficient in zip(((2,0),(2,1),(2,2),(1,1)),coefficients):
        if coefficient == "0":
            continue
        nu,A,D = f"{1+n}-3*eps",f"{m+1}-eps",str(n-m+1)
        ratio = f"Beta[{A},{D}]"
        special = f"(dtg_Q0)^(-({nu}))*Hypergeometric2F1[{nu},{A},{n+2}-eps,1-dtg_Q1/dtg_Q0]"
        kernel = f"(dtg_r*(1-dtg_r)*dtg_t*(1-dtg_t))^(-eps)*({coefficient})*dtg_h^({m-n}-2*eps)*({special})"
        term = _term(raw,f"central-{n}{m}",ratio,kernel,variables,bounds,order,
                     {"operation":"Cheng-Wu and exact central parameter Euler integral","denominator_shift":n,"central_power":m})
        term.update({"coefficient":coefficient,"hypergeometric":{"a":nu,"b":A,"c":f"{n+2}-eps","argument":"1-dtg_Q1/dtg_Q0","normalization":"ordinary Gauss 2F1"}})
        terms.append(term)
    return {"method":"central_parameter_gauss","definitions":definitions,"terms":terms,
            "dimension_before":8,"dimension_after":7,
            "ordered_chains":deepcopy(raw["J_raw_s"]["ordered_chains"]),
            "change_of_variables":{"simplex_alpha":[f"({v})/(1+dtg_c)" for v in aa],
                "simplex_jacobian":"dtg_a*(1-dtg_a)/(1+dtg_c)^5",
                "central_substitution":"dtg_c=dtg_h*x/(1-x)","central_jacobian":"dtg_h/(1-x)^2",
                "B_projective":"dtg_h*dtg_Q0+dtg_c*dtg_Q1","U_projective":"dtg_h+dtg_c",
                "finite_regulator_projective":"B_projective+I*eta_reg*U_projective*(1+dtg_c)"},
            "starting_region":"fixed interior remaining parameters with Q0,Q1>0 and Re(eps)<1; outer integration in a common convergence domain or by joint meromorphic continuation",
            "endpoint_rules":["Use the complete Euler/Gauss boundary value when Q0 tends to zero; do not plug zero into separate factors.",
                              "Q0=Q1 has hypergeometric value 1, with no difference-quotient pole.",
                              "Individual expanded terms need not have a common outer convergence strip; keep their common combination."]}


def build_reduced(orbit,raw):
    if all(raw["J_raw_s"][key] == "0" for key in ("pnum1","pnum2")):
        out = {"method":"zero_geometric_numerator","definitions":[],"terms":[],
               "dimension_before":8,"dimension_after":0,
               "starting_region":"the complete geometric numerator is identically zero"}
    else:
        out = _same_edge_reduced(orbit,raw) if _paired_edges(orbit,raw) else _central_reduced(orbit,raw)
    out.update({"schema_version":"double-tgv-reduced-1.0","representation":"sum_of_explicit_CJ_terms",
                "definition_semantics":"ordered explicit scalar substitution, shared by all terms",
                "normalization":_normalization(raw),"source_representation":"J_raw_s",
                "equivalence":{"kind":"analytic integral identity in stated starting region, then common meromorphic continuation",
                               "finite_auxiliary_regulator_identity":False,"complex_continuation_verified":False},
                "sum_rule":"C * J_raw_s = sum(term.combined_C * integral(term.integrand)); do not multiply raw C again",
                "evaluation_status":"algebraic_zero" if raw["J_raw_s"]["physical_zero"] else "not_evaluated"})
    return out
