"""J_n parametric integrand assembly.

J_n = ∫ d^{2n}x  ×  Jacobian  ×  ℓ_n  ×  Π_p f_p^{-1+ε}

Where:
  ℓ_n = Π_p (v_a·v_b / v_ref)       — dimensionless Lorentz factor
  f_p = X_p²/(2 v_ref)               — normalized squared distance
  Jacobian = hypercube transformation Jacobian
"""
import re
from collections import Counter
from fractions import Fraction

from montopar.models import (
    Orbit, HypercubeResult, JFactors, Pairing,
    TGVJFactors, FGVJFactors, MixedJFactors,
)
from montopar.tgv import compute_S, compute_Pnum
from montopar.fgv import compute_Q_prime, compute_lorentzZ, is_physical_zero

_ATOM_RE = re.compile(r"^[A-Za-z0-9_]+$")


# ---------------------------------------------------------------------------
# Polynomial expansion
# ---------------------------------------------------------------------------
#
# After substituting s -> x variables, each f_p is a product of sums (e.g.
# ``(1-x1*x3)*x2*x4*z1``).  The hepeflow-preprocessing reference stores the
# fully-expanded form (``x2*x4*z1 - x1*x2*x3*x4*z1``), produced by
# Mathematica's ``Expand``.  ``expand_polynomial`` reproduces that form so
# montopar output is byte-identical to the golden reference.
#
# Representation: a polynomial is a ``dict[monomial, coeff]`` where
# ``monomial`` is a ``frozenset``-like key derived from a ``Counter`` of
# variable-name → exponent, and ``coeff`` is an integer coefficient.

_TOKEN_RE = re.compile(r"\s*(\d+/\d+|\d+|[A-Za-z_][A-Za-z0-9_]*|[+\-*^()])")


def _tokenize(expr: str) -> list[str]:
    tokens = []
    pos = 0
    while pos < len(expr):
        m = _TOKEN_RE.match(expr, pos)
        if not m:
            if expr[pos].isspace():
                pos += 1
                continue
            raise ValueError(f"Unexpected character {expr[pos]!r} in {expr!r}")
        tokens.append(m.group(1))
        pos = m.end()
    return tokens


class _Parser:
    """Recursive-descent parser: expr -> polynomial dict."""

    def __init__(self, tokens: list[str]):
        self.tokens = tokens
        self.i = 0

    def peek(self) -> str | None:
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def next(self) -> str:
        tok = self.tokens[self.i]
        self.i += 1
        return tok

    def parse_expr(self) -> dict:
        poly = self.parse_term()
        while self.peek() in ("+", "-"):
            op = self.next()
            rhs = self.parse_term()
            poly = _poly_add(poly, rhs, sign=1 if op == "+" else -1)
        return poly

    def parse_term(self) -> dict:
        poly = self.parse_factor()
        while self.peek() == "*":
            self.next()
            rhs = self.parse_factor()
            poly = _poly_mul(poly, rhs)
        return poly

    def parse_factor(self) -> dict:
        tok = self.peek()
        if tok == "-":
            self.next()
            return _poly_scale(self.parse_factor(), -1)
        if tok == "+":
            self.next()
            return self.parse_factor()  # unary plus, no-op
        poly = self.parse_primary()
        if self.peek() == "^":
            self.next()
            exponent = self.next()
            if not exponent.isdigit():
                raise ValueError("Polynomial exponent must be a non-negative integer")
            poly = _poly_pow(poly, int(exponent))
        return poly

    def parse_primary(self) -> dict:
        tok = self.peek()
        if tok == "(":
            self.next()
            inner = self.parse_expr()
            if self.peek() != ")":
                raise ValueError("Missing ')'")
            self.next()
            return inner
        if tok is None:
            raise ValueError("Unexpected end of expression")
        # NUMBER or NAME
        self.next()
        return _poly_from_atom(tok)


def _poly_from_atom(tok: str) -> dict:
    if "/" in tok:
        num, den = tok.split("/", 1)
        return {frozenset(): Fraction(int(num), int(den))}
    if tok.isdigit():
        return {frozenset(): int(tok)}
    # variable name -> monomial with a single variable, exponent 1
    key = frozenset(Counter({tok: 1}).items())
    return {key: 1}


def _poly_scale(poly: dict, k: int) -> dict:
    return {m: c * k for m, c in poly.items() if c * k != 0}


def _poly_add(a: dict, b: dict, sign: int = 1) -> dict:
    out = dict(a)
    for m, c in b.items():
        out[m] = out.get(m, 0) + sign * c
    return {m: c for m, c in out.items() if c != 0}


def _poly_mul(a: dict, b: dict) -> dict:
    out: dict = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            merged = _merge_monomials(m1, m2)
            out[merged] = out.get(merged, 0) + c1 * c2
    return {m: c for m, c in out.items() if c != 0}


def _poly_pow(poly: dict, exponent: int) -> dict:
    result = {frozenset(): 1}
    base = poly
    power = exponent
    while power:
        if power & 1:
            result = _poly_mul(result, base)
        power //= 2
        if power:
            base = _poly_mul(base, base)
    return result


def _merge_monomials(m1: frozenset, m2: frozenset) -> frozenset:
    c1 = Counter(dict(m1))
    c2 = Counter(dict(m2))
    merged = c1 + c2
    return frozenset(merged.items())


def _monomial_key(m: frozenset) -> tuple:
    """Sort key for a monomial's variables: x-vars by index, then others."""
    counts = Counter(dict(m))
    x_items = []
    other_items = []
    for var, cnt in counts.items():
        if var.startswith("x") and var[1:].isdigit():
            x_items.append((int(var[1:]), cnt))
        else:
            other_items.append((var, cnt))
    x_items.sort(key=lambda t: t[0])
    other_items.sort(key=lambda t: t[0])
    return (tuple(x_items), tuple(other_items))


def _format_monomial(m: frozenset) -> str:
    x_key, other_key = _monomial_key(m)
    parts = []
    for idx, cnt in x_key:
        parts.append(f"x{idx}" if cnt == 1 else f"x{idx}^{cnt}")
    for var, cnt in other_key:
        parts.append(var if cnt == 1 else f"{var}^{cnt}")
    return "*".join(parts)


def _format_poly(poly: dict) -> str:
    if not poly:
        return "0"
    # Order terms by (total degree, monomial key).  Matches the hepeflow
    # convention of writing the lower-degree (positive) term first for
    # (1-A)*B-style f_p factors.
    ordered = sorted(
        poly.items(),
        key=lambda kv: (sum(dict(kv[0]).values()), _monomial_key(kv[0])),
    )
    out_parts = []
    for m, coeff in ordered:
        mono_str = _format_monomial(m)
        if mono_str == "":  # constant term
            out_parts.append(str(coeff))
        elif coeff == 1:
            out_parts.append(mono_str)
        elif coeff == -1:
            out_parts.append(f"-{mono_str}")
        else:
            out_parts.append(f"{coeff}*{mono_str}")
    # Join with explicit + / - signs.
    result = out_parts[0]
    for part in out_parts[1:]:
        if part.startswith("-"):
            result += f" - {part[1:]}"
        else:
            result += f" + {part}"
    return result


def _parse_polynomial(expr: str) -> dict:
    """Parse a polynomial expression into the internal exact representation."""
    tokens = _tokenize(expr)
    if not tokens:
        return {}
    parser = _Parser(tokens)
    poly = parser.parse_expr()
    if parser.i != len(tokens):
        raise ValueError(f"Trailing tokens in expression: {expr!r}")
    return poly


def expand_polynomial(expr: str) -> str:
    """Expand a product-of-sums polynomial string into canonical form."""
    return _format_poly(_parse_polynomial(expr))


def differentiate_polynomial(expr: str, var: str) -> str:
    """Differentiate a polynomial exactly with respect to one variable."""
    derivative: dict = {}
    for monomial, coeff in _parse_polynomial(expr).items():
        powers = Counter(dict(monomial))
        exponent = powers.get(var, 0)
        if exponent == 0:
            continue
        if exponent == 1:
            del powers[var]
        else:
            powers[var] = exponent - 1
        key = frozenset(powers.items())
        derivative[key] = derivative.get(key, 0) + coeff * exponent
    return _format_poly({m: c for m, c in derivative.items() if c != 0})


def linear_combination_polynomials(
    terms: list[tuple[str, str]],
) -> str:
    """Return the exact canonical sum of coefficient × polynomial terms."""
    result: dict = {}
    for coefficient, polynomial in terms:
        product = _poly_mul(
            _parse_polynomial(coefficient),
            _parse_polynomial(polynomial),
        )
        result = _poly_add(result, product)
    return _format_poly(result)



def dot_product_factor(edge_a: int, edge_b: int) -> str:
    """Compute v_a·v_b / v_ref for a given edge pair.

    N=4 kinematics:
      (1,2) → z, (3,4) → z
      (2,3) → 1, (4,1) → 1
      (1,3) → -(1+z), (2,4) → -(1+z)
    """
    # Canonical lookup with a < b
    a, b = (edge_a, edge_b) if edge_a < edge_b else (edge_b, edge_a)

    if (a, b) == (1, 2) or (a, b) == (3, 4):
        return "z1"
    elif (a, b) == (2, 3) or (a, b) == (1, 4):
        return "1"
    elif (a, b) == (1, 3) or (a, b) == (2, 4):
        return "-(1+z1)"
    else:
        raise ValueError(f"Unknown edge pair: ({edge_a}, {edge_b})")


def compute_f_p(edge_a: int, edge_b: int, s_a: str, s_b: str) -> str:
    """Compute f_p = X_p²/(2 v_ref) for a gluon connecting edge_a to edge_b.

    Adjacent edges (|a-b| = 1 mod 4):
      f_p = (1-s_a)*s_b * (v_a·v_b/v_ref)

    Opposite edges (|a-b| = 2 mod 4):
      {1,3}: f_p = z*(1-s_a)*(1-s_b) + s_a*s_b
      {2,4}: f_p = (1-s_a)*(1-s_b) + z*s_a*s_b
    """
    diff = abs(edge_a - edge_b) % 4

    if diff == 1 or diff == 3:
        # Adjacent edges: use contour ordering 1->2->3->4->1
        # diff=1 -> early=min, late=max; diff=3 -> early=max, late=min
        dot = dot_product_factor(edge_a, edge_b)
        if diff == 1:
            e_early, e_late = (edge_a, edge_b) if edge_a < edge_b else (edge_b, edge_a)
        else:
            e_early, e_late = (edge_a, edge_b) if edge_a > edge_b else (edge_b, edge_a)
        s_early = s_a if e_early == edge_a else s_b
        s_late = s_a if e_late == edge_a else s_b
        return f"(1-{s_early})*{s_late}*{dot}"

    elif diff == 2:
        # Opposite edges
        a, b = (edge_a, edge_b) if edge_a < edge_b else (edge_b, edge_a)
        if (a, b) == (1, 3):
            return f"z1*(1-{s_a})*(1-{s_b}) + {s_a}*{s_b}"
        else:  # (2, 4)
            return f"(1-{s_a})*(1-{s_b}) + z1*{s_a}*{s_b}"

    else:
        raise ValueError(f"Cannot compute f_p for edge pair ({edge_a}, {edge_b})")


def substitute_s_to_x(expr: str, s_to_x: dict[str, str]) -> str:
    """Replace all s variables in expr with their x expressions.

    Sorts s-variables by length (longest first) so that e.g. ``s311`` is
    replaced before ``s31`` is, preventing partial-match corruption.
    """
    result = expr
    for s_var in sorted(s_to_x.keys(), key=len, reverse=True):
        result = result.replace(s_var, s_to_x[s_var])
    return result


def compute_lorentz_factor(pairings: list[Pairing]) -> str:
    """Compute ℓ_n = product of (v_a·v_b/v_ref) for each pairing.

    Groups identical factors using power notation. Simple atomic factors
    (like ``z1``) are written without parentheses (``z1^2``); composite
    factors are wrapped in parentheses (``-(1+z1)^2``).
    """
    factors = [dot_product_factor(p.edge_a, p.edge_b) for p in pairings]
    counts = Counter(factors)

    parts = []
    for factor, count in sorted(counts.items()):
        if count == 1:
            parts.append(factor)
        else:
            if _ATOM_RE.fullmatch(factor):
                parts.append(f"{factor}^{count}")
            else:
                parts.append(f"({factor})^{count}")

    return "*".join(parts) if parts else "1"


def compute_j(orbit: Orbit, hc: HypercubeResult) -> JFactors:
    """Compute the J_n parametric integrand components.

    Args:
        orbit: Parsed orbit with edge distribution and pairings.
        hc: Hypercube transformation result (s→x mapping and Jacobian).

    Returns:
        JFactors with lorentz_z, jacobian, f_polynomials, f_product,
        lorentz_v14_power.
    """
    # Lorentz factor ℓ_n
    lorentz_z = compute_lorentz_factor(orbit.pairings)

    # f_p for each propagator, substituted s→x
    f_polynomials: dict[str, str] = {}
    f_factors: list[str] = []

    for i, pairing in enumerate(orbit.pairings, start=1):
        f_s = compute_f_p(pairing.edge_a, pairing.edge_b, pairing.s_a, pairing.s_b)
        f_x = substitute_s_to_x(f_s, hc.s_to_x)
        # Expand to canonical form to match the hepeflow-preprocessing
        # reference output (Mathematica ``Expand``).
        f_x = expand_polynomial(f_x)
        f_name = f"f{i}"
        f_polynomials[f_name] = f_x
        f_factors.append(f"{f_name}^(-1+eps)")

    f_product = " * ".join(f_factors)

    # Lorentz v14 power: one v_ref factor per pairing (denominator) → total
    # power equals the number of gluons / pairings.
    lorentz_v14_power = str(orbit.n_gluons)

    return JFactors(
        lorentz_z=lorentz_z,
        jacobian=hc.jacobian,
        f_product=f_product,
        f_polynomials=f_polynomials,
        lorentz_v14_power=lorentz_v14_power,
    )


def compute_j_expanded(j: JFactors) -> str:
    """将 f_p 多项式代入 J 被积函数，展开为完整表达式。

    输入: j.f_product = "f1^(-1+eps) * f2^(-1+eps)"
          j.f_polynomials = {"f1": "x2*z1 - x1*x2*z1", "f2": "..."}
    输出: "(x2*z1 - x1*x2*z1)^(-1+eps) * (...)^(-1+eps)"

    f_name 按长度降序替换，避免 "f1" 匹配到 "f11" 内部。
    """
    result = j.f_product
    for f_name in sorted(j.f_polynomials.keys(), key=len, reverse=True):
        f_poly = j.f_polynomials[f_name]
        result = result.replace(f_name, f"({f_poly})")
    return result


def compute_j_tgv(orbit: Orbit, hc: HypercubeResult) -> TGVJFactors:
    """Compute the J integrand for a TGV diagram.

    J = ∫ d³x ∫ d³α δ(1-Σα) Jacobian × Pnum / ((α₁α₂α₃)^ε · S^(2-2ε))

    Args:
        orbit: Parsed orbit with tgv_legs populated.
        hc: Hypercube transformation result (s→x mapping, Jacobian).

    Returns:
        TGVJFactors with pnum, S_polynomial, etc. in x-variables.
    """
    legs = orbit.tgv_legs

    # Compute S and Pnum in s-variables
    S_s = compute_S(legs)
    Pnum_s = compute_Pnum(legs)

    # s → x substitution + polynomial expansion
    S_x = expand_polynomial(substitute_s_to_x(S_s, hc.s_to_x))
    Pnum_x = expand_polynomial(substitute_s_to_x(Pnum_s, hc.s_to_x))

    return TGVJFactors(
        jacobian=hc.jacobian,
        pnum=Pnum_x,
        alpha_monomial="(alpha1*alpha2*alpha3)^(-eps)",
        S_polynomial=S_x,
        denom_power="2 - 2*eps",
        x_vars=hc.x_vars,
    )


def compute_j_fgv(orbit: Orbit, hc: HypercubeResult) -> FGVJFactors:
    """Compute the J integrand for an FGV diagram.

    J = ∫ d⁴x ∫ d⁴α δ(1-Σα) Jacobian × lorentzZ / ((α₁α₂α₃α₄)^ε · Q'^(2-3ε))

    Args:
        orbit: Parsed orbit with fgv_legs populated.
        hc: Hypercube transformation result.

    Returns:
        FGVJFactors with lorentzZ, Q_prime_polynomial, etc. in x-variables.
    """
    legs = orbit.fgv_legs

    # Q' in s-variables
    Qp_s = compute_Q_prime(legs)
    # lorentzZ is pure z — no s dependence
    lorentzZ_s = compute_lorentzZ(legs)

    # s → x substitution + polynomial expansion
    Qp_x = expand_polynomial(substitute_s_to_x(Qp_s, hc.s_to_x))
    # lorentzZ has no s-vars; expand to canonical form
    lorentzZ_x = expand_polynomial(lorentzZ_s)

    return FGVJFactors(
        jacobian=hc.jacobian,
        lorentzZ=lorentzZ_x,
        alpha_monomial="(alpha1*alpha2*alpha3*alpha4)^(-eps)",
        Q_prime_polynomial=Qp_x,
        denom_power="2 - 3*eps",
        x_vars=hc.x_vars,
        physical_zero=is_physical_zero(orbit.edge_distribution),
    )


def compute_j_mixed(
    orbit: "Orbit",
    hc: "HypercubeResult",
    pure_pairings: "list[Pairing]",
) -> "MixedJFactors":
    """Compute J for mixed (TGV + pure) diagrams.

    Args:
        orbit: Parsed orbit with tgv_legs populated.
        hc: Hypercube transformation result.
        pure_pairings: Pure gluon pairings (excludes TGV loops).

    Returns:
        MixedJFactors with pure and TGV parts in x-variables.
    """
    # Pure part
    lorentz_z = compute_lorentz_factor(pure_pairings)
    f_polynomials: dict[str, str] = {}
    for i, p in enumerate(pure_pairings, start=1):
        f_s = compute_f_p(p.edge_a, p.edge_b, p.s_a, p.s_b)
        f_x = expand_polynomial(substitute_s_to_x(f_s, hc.s_to_x))
        f_polynomials[f"f{i}"] = f_x

    f_product = " * ".join(
        f"{name}^(-1+eps)" for name in f_polynomials
    )

    # TGV part
    S_s = compute_S(orbit.tgv_legs)
    Pnum_s = compute_Pnum(orbit.tgv_legs)
    S_x = expand_polynomial(substitute_s_to_x(S_s, hc.s_to_x))
    Pnum_x = expand_polynomial(substitute_s_to_x(Pnum_s, hc.s_to_x))

    return MixedJFactors(
        lorentz_z=lorentz_z,
        jacobian=hc.jacobian,
        f_product=f_product,
        f_polynomials=f_polynomials,
        edge_pairs=[(p.edge_a, p.edge_b) for p in pure_pairings],
        pnum=Pnum_x,
        alpha_monomial="(alpha1*alpha2*alpha3)^(-eps)",
        S_polynomial=S_x,
        denom_power="2 - 2*eps",
        x_vars=hc.x_vars,
        n_pure=len(pure_pairings),
    )
