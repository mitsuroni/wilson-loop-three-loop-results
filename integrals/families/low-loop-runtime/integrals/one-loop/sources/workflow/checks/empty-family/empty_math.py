"""Computational adapters for the supplied pure-exchange empty inputs.
No Lean, Kernel or proof-store implementation is changed.
"""
from pathlib import Path
import sys
import itertools
import math
import sympy as sp
from sympy.printing.mathematica import mathematica_code

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from preprocessing.ordered import arithmetic, integral_ast

eps = sp.Symbol("eps", positive=True)
z = sp.Symbol("z", positive=True)
ell = sp.Symbol("ell", real=True)
R2, R12, G = sp.symbols("R2 R12 G")
G_exact = sp.exp(-sp.EulerGamma * eps) * sp.gamma(1 - eps)
R2_exact = sp.gamma(1 + eps)**2 / sp.gamma(1 + 2*eps)
R12_exact = sp.gamma(1 + eps)*sp.gamma(1 + 2*eps) / sp.gamma(1 + 3*eps)
SUBS = {G: G_exact, R2: R2_exact, R12: R12_exact}


def to_sympy(text, variables):
    names = {"eps": eps, "z1": z, **{v: sp.Symbol(v, positive=True) for v in variables}}
    def visit(node):
        if node["kind"] == "rational":
            return sp.Rational(node["numerator"], node["denominator"])
        if node["kind"] == "var":
            return names[node["name"]]
        args = [visit(a) for a in node["args"]]
        return {"add": lambda: args[0]+args[1], "sub": lambda: args[0]-args[1],
                "mul": lambda: args[0]*args[1], "div": lambda: args[0]/args[1],
                "neg": lambda: -args[0], "rpow": lambda: args[0]**args[1]}[node["op"]]()
    return visit(arithmetic(text, set(names)))


def word_integral(word):
    m = len(word)
    if not word:
        return sp.Integer(1)
    if set(word) in ({"A"}, {"B"}):
        return 1/(sp.factorial(m)*eps**m)
    if word == "AB":
        return R2/eps**2
    if word == "BA":
        return (1-R2)/eps**2
    if word in ("AAB", "ABB"):
        return R12/(2*eps**3)
    if word in ("ABA", "BAB"):
        return (R2-R12)/eps**3
    if word in ("BAA", "BBA"):
        return (1-2*R2+R12)/(2*eps**3)
    raise ValueError("unimplemented ordered word: "+word)


def analyze(row):
    o, n = row["orbit"], row["loop_order"]
    j = o["J"]
    integral_ast(j)  # existing closed reader checks inventory, arithmetic and scope
    assert o["dispatch"] == "pure"
    assert o["integral_structure"]["direct_opposite_edge_links"]["channel_partition"] == []
    assert len(j["variables"]) == 2*n == 2*len(j["propagator_polynomials"])
    variables = {v: sp.Symbol(v, positive=True) for v in j["variables"]}
    letters, matches, scale = {}, [], sp.Integer(1)
    for name, raw in j["propagator_polynomials"].items():
        f = to_sympy(raw, variables)
        candidates = [(a,b,q) for a,b in itertools.permutations(variables,2)
                      for q in (sp.Integer(1),z)
                      if sp.expand(f-q*(1-variables[b])*variables[a]) == 0]
        assert len(candidates) == 1, (name, candidates)
        a,b,q = candidates[0]
        assert a not in letters and b not in letters, "repeated endpoint requires new mathematics"
        letters[a], letters[b] = "A","B"
        scale *= q
        matches.append({"propagator": name, "A": a, "B": b, "scale": str(q),
                        "polynomial_residual": "0"})
    assert sp.expand(to_sympy(j["numerator"], variables)-scale) == 0
    assert set(letters) == set(variables)
    child, tops = {}, []
    for v, (low,high) in j["bounds"].items():
        assert low == "0"
        if high == "1":
            tops.append(v)
        else:
            assert high in variables and high not in child, "domain must be disjoint chains"
            child[high] = v
    chains, seen = [], set()
    for top in sorted(tops):
        descending, v = [], top
        while True:
            assert v not in seen, "overlapping or cyclic domain"
            seen.add(v)
            descending.append(v)
            if v not in child: break
            v = child[v]
        ordered = list(reversed(descending))
        word = "".join(letters[v] for v in ordered)
        chains.append({"variables_low_to_high": ordered, "word": word})
    assert seen == set(variables)
    K = sp.factor(sp.prod(word_integral(c["word"]) for c in chains))
    r = int(sp.degree(scale,z))
    assert scale == z**r and 0 <= r <= n
    a1 = sp.Rational(o["color_preparation"]["a1"])
    assert str(a1) == o["color_preparation"]["coordinates"][0]
    ca, cf = sp.symbols("CA CF", positive=True)
    projected = to_sympy(o["color_preparation"]["projected_color"], ["CA","CF"])
    assert sp.expand(projected-a1*ca**(n-1)*cf) == 0
    assert o["C"]["components"]["color"] == o["color_preparation"]["projected_color"]
    sym = o["c4_symmetry"]["orbit_multiplicity"]
    expected_factor = "Exp[-EulerGamma*eps]" if n == 1 else f"Exp[-{n}*EulerGamma*eps]"
    expected_C = (f"({o['color_preparation']['projected_color']}) * ({expected_factor})"
                  f" * ((4*Pi*alphaS)^{n}) * (Gamma[1-eps]^{n}/(4*Pi^(2-eps))^{n})"
                  f" * ((-2*v14+I0)^({n}*eps)) * ((-2)^(-{n}))")
    assert o["C"]["combined"] == expected_C, "unexpected normalization factor"
    assert o["euler_preparation"]["gamma_power"] == n
    assert o["euler_preparation"]["scheme"] == "exp(-loop_order*EulerGamma*eps)"
    kappa = sp.Rational(sym)*a1/(2*(-2)**n)
    ratio = sp.cosh((sp.Rational(n,2)-r)*eps*ell)/sp.cosh(sp.Rational(n,2)*eps*ell)
    compact = kappa*G**n*K*ratio
    J_exact = z**(r*eps)*K.subs(SUBS)
    I_exact = compact.subs(SUBS).subs(ell,sp.log(z))
    return {"letters":letters,"propagators":matches,"chains":chains,"r":r,
            "a1":a1,"sym":sym,"kappa":kappa,"K":K,"ratio":ratio,
            "I_compact":compact,"J_exact":J_exact,"I_exact":I_exact}


def exp_series(log_coeff, order):
    coeff = [sp.Integer(1)]
    for m in range(1,order+1):
        coeff.append(sp.expand(sum(k*log_coeff.get(k,0)*coeff[m-k]
                                   for k in range(1,m+1))/m))
    return sum(coeff[k]*eps**k for k in range(order+1))


def truncate(expr, order):
    return sp.Poly(sp.expand(expr),eps).as_dict()


def laurent(a, n, through=2):
    order = 2*n+through
    g = exp_series({k:n*sp.zeta(k)/k for k in range(2,order+1)},order)
    r2 = exp_series({k:(-1)**k*sp.zeta(k)*(2-2**k)/k
                    for k in range(2,order+1)},order)
    r12 = exp_series({k:(-1)**k*sp.zeta(k)*(1+2**k-3**k)/k
                     for k in range(2,order+1)},order)
    def cut(expr):
        return sum(c*eps**powers[0] for powers,c in truncate(expr,order).items()
                   if powers[0] <= order)
    reg = cut(sp.cancel(eps**(2*n)*a["K"]).subs({R2:r2,R12:r12}))
    ratio = sp.series(a["ratio"],eps,0,order+1).removeO()
    reg = cut(cut(reg*g)*ratio*a["kappa"])
    return {k:sp.factor(reg.coeff(eps,k+2*n)) for k in range(-2*n,through+1)}


def latex_math(expr):
    return sp.latex(expr, symbol_names={eps:r"\varepsilon", ell:r"\ell"})


def expression(expr):
    return {"language":"WolframLanguage/InputForm",
            "expression":mathematica_code(expr), "latex":latex_math(expr)}
