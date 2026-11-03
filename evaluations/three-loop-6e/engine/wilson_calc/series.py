"""Exact finite Laurent arithmetic over the existing scalar JSON AST.

``Series(c, lower_bound=a, through=k)`` knows every coefficient up to k;
all powers below a are zero.  Missing entries inside the known range are zero,
but entries above it are *unknown*.  This is an algebraic precision contract,
not a proof that a parameter integral admits a Laurent expansion.

The small scalar simplifier only uses rational arithmetic and ring identities.
It neither evaluates source text nor cancels a possibly zero symbolic divisor.
"""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
import json
import re


class SeriesError(ValueError):
    """Invalid input or an operation outside the supported coefficient ring."""


class InsufficientPrecision(SeriesError):
    """The requested coefficient depends on uncomputed input coefficients."""


_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_-]*\Z")
_INTEGER = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")
_NATURAL = re.compile(r"[1-9][0-9]*\Z")
_ARITY = {"add": 2, "sub": 2, "mul": 2, "div": 2, "neg": 1,
          "rpow": 2, "exp": 1, "log": 1, "gamma": 1, "zeta": 1,
          "pi": 0, "euler_gamma": 0}


def _integer(value, label):
    if type(value) is not int:
        raise SeriesError(label + " must be an integer")
    return value


def q(n, d=1):
    """Build a canonical rational in the existing AST (no floats or bools)."""
    if not isinstance(n, Fraction):
        _integer(n, "numerator")
    _integer(d, "denominator")
    if d == 0:
        raise SeriesError("zero denominator")
    value = Fraction(n, d)
    return {"kind": "rational", "numerator": str(value.numerator),
            "denominator": str(value.denominator)}


def _key(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def validate_scalar(term, variables=None):
    """Validate a supported subset of Wilson's scalar AST, without evaluation.

    If variables is provided, every free variable must occur in that inventory.
    Functions not already in Wilson's scalar grammar, binders, and raw strings
    are rejected.  Domain assumptions remain the caller's responsibility.
    """
    allowed = None if variables is None else set(variables)
    if allowed is not None and any(not isinstance(n, str) or not _NAME.fullmatch(n) for n in allowed):
        raise SeriesError("invalid variable inventory")
    count = 0

    def visit(node, depth):
        nonlocal count
        count += 1
        if count > 10000 or depth > 128:
            raise SeriesError("scalar AST exceeds validation budget")
        if not isinstance(node, dict):
            raise SeriesError("scalar must be an AST object")
        kind = node.get("kind")
        if kind == "rational":
            if set(node) != {"kind", "numerator", "denominator"}:
                raise SeriesError("invalid rational fields")
            n, d = node["numerator"], node["denominator"]
            if not isinstance(n, str) or not _INTEGER.fullmatch(n) or not isinstance(d, str) or not _NATURAL.fullmatch(d):
                raise SeriesError("invalid rational spelling")
            int(n); int(d)
        elif kind == "var":
            name = node.get("name")
            if set(node) != {"kind", "name"} or not isinstance(name, str) or not _NAME.fullmatch(name):
                raise SeriesError("invalid variable")
            if allowed is not None and name not in allowed:
                raise SeriesError("unbound variable: " + name)
        elif kind == "op":
            if set(node) != {"kind", "op", "args"} or node.get("op") not in _ARITY:
                raise SeriesError("unsupported scalar operation")
            args = node["args"]
            if not isinstance(args, list) or len(args) != _ARITY[node["op"]]:
                raise SeriesError("wrong scalar operation arity")
            for arg in args:
                visit(arg, depth + 1)
        else:
            raise SeriesError("unsupported scalar AST kind")
    visit(term, 0)


def free_variables(term):
    validate_scalar(term)
    if term["kind"] == "var":
        return {term["name"]}
    return set().union(*(free_variables(a) for a in term.get("args", [])))


def as_fraction(term):
    """Return an exact rational value only when syntactically rational."""
    if isinstance(term, dict) and term.get("kind") == "rational":
        validate_scalar(term)
        return Fraction(int(term["numerator"]), int(term["denominator"]))
    return None


def _binary(op, items):
    result = items[0]
    for item in items[1:]:
        result = {"kind": "op", "op": op, "args": [result, item]}
    return result


def _mul(*items):
    scalar, factors = Fraction(1), []
    def collect(item):
        nonlocal scalar
        value = as_fraction(item)
        if value is not None:
            scalar *= value
        elif item.get("kind") == "op" and item["op"] == "mul":
            for arg in item["args"]:
                collect(arg)
        else:
            factors.append(item)
    for item in items:
        collect(item)
    if not scalar:
        return q(0)
    factors.sort(key=_key)
    if scalar != 1 or not factors:
        factors.insert(0, q(scalar))
    return _binary("mul", factors)


def _add(*items):
    values, bases = defaultdict(Fraction), {}
    def collect(item, sign=Fraction(1)):
        value = as_fraction(item)
        if value is not None:
            values[""] += sign * value
        elif item.get("kind") == "op" and item["op"] == "add":
            for arg in item["args"]:
                collect(arg, sign)
        else:
            scalar, rest = Fraction(1), []
            def factors(x):
                nonlocal scalar
                val = as_fraction(x)
                if val is not None:
                    scalar *= val
                elif x.get("kind") == "op" and x["op"] == "mul":
                    for a in x["args"]:
                        factors(a)
                else:
                    rest.append(x)
            factors(item)
            base = _binary("mul", rest) if rest else q(1)
            key = _key(base)
            bases[key] = base
            values[key] += sign * scalar
    for item in items:
        collect(item)
    result = [q(values[""])] if values[""] else []
    result += [_mul(q(values[key]), bases[key]) for key in sorted(bases) if values[key]]
    return _binary("add", result) if result else q(0)


def normalize(term):
    """Canonical rational/ring simplification; no new AST operations."""
    validate_scalar(term)
    def visit(node):
        if node["kind"] == "rational":
            return q(Fraction(int(node["numerator"]), int(node["denominator"])))
        if node["kind"] == "var":
            return dict(node)
        op, args = node["op"], [visit(a) for a in node["args"]]
        if op == "add": return _add(*args)
        if op == "sub": return _add(args[0], _mul(q(-1), args[1]))
        if op == "neg": return _mul(q(-1), args[0])
        if op == "mul": return _mul(*args)
        if op == "div":
            denominator = as_fraction(args[1])
            if denominator == 0: raise SeriesError("division by zero")
            if denominator is not None: return _mul(args[0], q(1 / denominator))
        if op == "rpow":
            base, exponent = map(as_fraction, args)
            if exponent is not None and exponent.denominator == 1:
                if base == 0 and exponent < 0: raise SeriesError("negative power of zero")
                if exponent == 0: return q(1)
                if exponent == 1: return args[0]
                if base is not None: return q(base ** int(exponent))
        if op == "exp" and as_fraction(args[0]) == 0: return q(1)
        if op == "log" and as_fraction(args[0]) == 1: return q(0)
        return {"kind": "op", "op": op, "args": args}
    return visit(term)


def add_expr(*terms):
    return _add(*(normalize(t) for t in terms))


def mul_expr(*terms):
    return _mul(*(normalize(t) for t in terms))


def product_orders(through, left_lower, right_lower):
    """Sufficient input orders for a product, given established lower bounds."""
    for value in (through, left_lower, right_lower): _integer(value, "order")
    return through - right_lower, through - left_lower


@dataclass(frozen=True, init=False)
class Series:
    _terms: tuple
    lower_bound: int
    through: int
    variable: str

    def __init__(self, coefficients, lower_bound, through, *, variable="eps"):
        _integer(lower_bound, "lower_bound"); _integer(through, "through")
        if not isinstance(variable, str) or not _NAME.fullmatch(variable):
            raise SeriesError("invalid expansion variable")
        if not isinstance(coefficients, dict):
            raise SeriesError("coefficients must be a dictionary with integer powers")
        terms = []
        for power, term in coefficients.items():
            _integer(power, "coefficient power")
            if not lower_bound <= power <= through:
                raise SeriesError("coefficient outside declared precision range")
            term = normalize(term)
            if variable in free_variables(term):
                raise SeriesError("coefficient depends on the expansion variable")
            if as_fraction(term) != 0:
                terms.append((power, _key(term)))
        object.__setattr__(self, "_terms", tuple(sorted(terms)))
        object.__setattr__(self, "lower_bound", lower_bound)
        object.__setattr__(self, "through", through)
        object.__setattr__(self, "variable", variable)

    @property
    def coefficients(self):
        return {power: json.loads(term) for power, term in self._terms}

    def coefficient(self, power):
        _integer(power, "coefficient power")
        if power < self.lower_bound: return q(0)
        if power > self.through:
            raise InsufficientPrecision(f"coefficient {power} is unknown above {self.through}")
        return self.coefficients.get(power, q(0))

    def truncate(self, through):
        _integer(through, "through")
        if through > self.through and through >= self.lower_bound:
            raise InsufficientPrecision("truncation cannot increase precision")
        return Series({p: c for p, c in self.coefficients.items() if p <= through},
                      self.lower_bound, through, variable=self.variable)

    def _compatible(self, other):
        if not isinstance(other, Series) or self.variable != other.variable:
            raise SeriesError("series expansion variables differ")

    def add(self, other, *, through=None):
        self._compatible(other)
        known = min(max(self.through, self.lower_bound - 1),
                    max(other.through, other.lower_bound - 1))
        target = known if through is None else _integer(through, "through")
        if target > known: raise InsufficientPrecision("sum needs higher input coefficients")
        lower = min(self.lower_bound, other.lower_bound)
        return Series({p: add_expr(self.coefficient(p), other.coefficient(p))
                       for p in range(lower, target + 1)}, lower, target, variable=self.variable)

    def multiply(self, other, *, through=None):
        self._compatible(other)
        known = min(max(self.through, self.lower_bound - 1) + other.lower_bound,
                    max(other.through, other.lower_bound - 1) + self.lower_bound)
        target = known if through is None else _integer(through, "through")
        if target > known: raise InsufficientPrecision("product needs higher input coefficients")
        lower = self.lower_bound + other.lower_bound
        coefficients = {}
        for p in range(lower, target + 1):
            coefficients[p] = add_expr(*(mul_expr(self.coefficient(i), other.coefficient(p-i))
                                        for i in range(self.lower_bound, p-other.lower_bound+1)))
        return Series(coefficients, lower, target, variable=self.variable)

    def inverse(self, *, through=None):
        terms = self.coefficients
        if not terms: raise InsufficientPrecision("nonzero leading coefficient is not known")
        valuation = min(terms)
        leading = as_fraction(terms[valuation])
        if leading is None:
            raise SeriesError("inverse requires an established nonzero rational leading coefficient")
        known = self.through - 2 * valuation
        target = known if through is None else _integer(through, "through")
        if target > known: raise InsufficientPrecision("inverse needs higher input coefficients")
        degree = target + valuation
        inverse = {0: q(1 / leading)}
        for n in range(1, degree + 1):
            inverse[n] = mul_expr(q(-1 / leading), add_expr(*(
                mul_expr(self.coefficient(valuation+k), inverse[n-k]) for k in range(1, n+1))))
        return Series({n-valuation: c for n, c in inverse.items() if n <= degree},
                      -valuation, target, variable=self.variable)

    def exponential(self, *, through=None):
        if self.through < 0 and self.lower_bound <= 0:
            raise InsufficientPrecision("constant exponent coefficient is unknown")
        if any(p < 0 for p, _ in self._terms):
            raise SeriesError("exponential of a pole is not a Laurent series with finite lower bound")
        known = max(self.through, self.lower_bound - 1)
        target = known if through is None else _integer(through, "through")
        if target > known: raise InsufficientPrecision("exponential needs higher input coefficients")
        if target < 0: return Series({}, 0, target, variable=self.variable)
        constant = normalize({"kind": "op", "op": "exp", "args": [self.coefficient(0)]})
        out = {0: constant}
        for n in range(1, target + 1):
            out[n] = mul_expr(q(1, n), add_expr(*(mul_expr(q(k), self.coefficient(k), out[n-k])
                                                 for k in range(1, n+1))))
        return Series(out, 0, target, variable=self.variable)

    def to_json(self):
        return {"format": "wilson-laurent-series/0.1", "variable": self.variable,
                "lower_bound": self.lower_bound, "through": self.through,
                "coefficients": {str(p): c for p, c in self.coefficients.items()}}

    @classmethod
    def from_json(cls, value):
        fields = {"format", "variable", "lower_bound", "through", "coefficients"}
        if not isinstance(value, dict) or set(value) != fields or value["format"] != "wilson-laurent-series/0.1":
            raise SeriesError("unsupported Laurent series payload")
        coefficients = value["coefficients"]
        if not isinstance(coefficients, dict) or any(not isinstance(p, str) or not _INTEGER.fullmatch(p) or str(int(p)) != p for p in coefficients):
            raise SeriesError("invalid coefficient powers")
        return cls({int(p): c for p, c in coefficients.items()}, value["lower_bound"],
                   value["through"], variable=value["variable"])


def add(left, right, *, through=None): return left.add(right, through=through)
def multiply(left, right, *, through=None): return left.multiply(right, through=through)
def inverse(value, *, through=None): return value.inverse(through=through)
def exponential(value, *, through=None): return value.exponential(through=through)
