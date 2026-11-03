"""Exact sparse GPL expressions and shuffle operations, not an integrator.

GPL letters, endpoints, and coefficients reuse Wilson's existing scalar AST.
Validation is structural: it does not prove contour admissibility, endpoint
convergence, or connection to an original integral.  Tangential endpoint values
stay explicitly marked and are never silently promoted to ordinary limits.
"""
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import json
from math import comb

from .series import (q, normalize, validate_scalar, free_variables, as_fraction,
                     add_expr, mul_expr, SeriesError)


FORMAT = "wilson-hlog-sum/0.1"
KERNEL = "dt/(t-a)"


class HLogError(ValueError):
    pass


def _key(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _scalar(value):
    try:
        return normalize(value)
    except (SeriesError, TypeError, ValueError) as exc:
        raise HLogError(str(exc)) from exc


@dataclass(frozen=True, init=False)
class GPL:
    _data: str

    def __init__(self, word, endpoint, *, basepoint=None, branch="positive-real", endpoint_kind="ordinary"):
        if not isinstance(word, (list, tuple)):
            raise HLogError("GPL word must be a sequence of scalar AST letters")
        if not isinstance(branch, str) or not branch or len(branch) > 256:
            raise HLogError("GPL branch must be a nonempty recorded convention")
        if endpoint_kind not in ("ordinary", "tangential"):
            raise HLogError("unsupported GPL endpoint kind")
        if len(word) > 512:
            raise HLogError("GPL word exceeds payload budget")
        data = {"word": [_scalar(a) for a in word], "endpoint": _scalar(endpoint),
                "basepoint": _scalar(q(0) if basepoint is None else basepoint),
                "branch": branch, "endpoint_kind": endpoint_kind}
        # These elementary divergent limits are known without a general endpoint
        # algorithm. More complicated endpoint admissibility belongs to rules.
        if endpoint_kind == "ordinary" and data["word"]:
            at = _key(data["endpoint"])
            if len(data["word"]) == 1 and at == _key(data["word"][0]):
                raise HLogError("ordinary GPL endpoint has a known logarithmic divergence")
            if all(as_fraction(a) == 0 for a in data["word"]) and as_fraction(data["endpoint"]) == 0:
                raise HLogError("all-zero GPL at zero needs an explicit tangential endpoint")
        object.__setattr__(self, "_data", _key(data))

    def to_json(self):
        return json.loads(self._data)

    @classmethod
    def from_json(cls, value):
        expected = {"word", "endpoint", "basepoint", "branch", "endpoint_kind"}
        if not isinstance(value, dict) or set(value) != expected:
            raise HLogError("invalid GPL fields")
        return cls(value["word"], value["endpoint"], basepoint=value["basepoint"],
                   branch=value["branch"], endpoint_kind=value["endpoint_kind"])

    @property
    def word(self): return self.to_json()["word"]

    @property
    def variables(self):
        data = self.to_json()
        return set().union(*(free_variables(a) for a in [*data["word"], data["endpoint"], data["basepoint"]]))

    @property
    def is_unit(self): return not self.word

    def path(self):
        data = self.to_json()
        return _key({k: v for k, v in data.items() if k != "word"})


def unit_gpl():
    return GPL([], q(0))


@lru_cache(maxsize=128)
def _shuffle(left, right):
    if not left: return ((right, 1),)
    if not right: return ((left, 1),)
    result = Counter()
    for word, count in _shuffle(left[1:], right): result[(left[0],) + word] += count
    for word, count in _shuffle(left, right[1:]): result[(right[0],) + word] += count
    return tuple(sorted(result.items()))


def shuffle_words(left, right, *, max_terms=100000):
    """Return [(word_AST_list, multiplicity)], retaining repeated-letter counts."""
    if not isinstance(left, (list, tuple)) or not isinstance(right, (list, tuple)):
        raise HLogError("shuffle operands must be words")
    if len(left) + len(right) > 64:
        raise HLogError("shuffle exceeds operation word-length budget")
    if type(max_terms) is not int or max_terms < 1:
        raise HLogError("positive shuffle term budget required")
    encoded = [tuple(_key(_scalar(x)) for x in word) for word in (left, right)]
    letters = set(encoded[0] + encoded[1])
    if len(letters) <= 1:
        word = encoded[0] + encoded[1]
        return [([json.loads(a) for a in word], comb(len(word), len(left)))]
    # Conservative bound before expanding or populating the memoization cache.
    if comb(len(left) + len(right), len(left)) > max_terms:
        raise HLogError("shuffle exceeds term budget")
    return [([json.loads(a) for a in word], count) for word, count in _shuffle(*encoded)]


@dataclass(frozen=True, init=False)
class HLogTerm:
    _coefficient: str
    gpl: GPL
    periods: tuple

    def __init__(self, coefficient, gpl, periods=()):
        coefficient = _scalar(coefficient)
        if not isinstance(gpl, GPL): raise HLogError("term requires a GPL object")
        if not isinstance(periods, (list, tuple)) or any(not isinstance(p, GPL) for p in periods):
            raise HLogError("periods must be GPL objects")
        if any(p.variables for p in periods):
            raise HLogError("GPL periods must not contain free parameters")
        object.__setattr__(self, "_coefficient", _key(coefficient))
        object.__setattr__(self, "gpl", unit_gpl() if gpl.is_unit else gpl)
        object.__setattr__(self, "periods", tuple(sorted((p for p in periods if not p.is_unit), key=lambda p: p._data)))

    @property
    def coefficient(self): return json.loads(self._coefficient)

    def to_json(self):
        return {"coefficient": self.coefficient, "gpl": self.gpl.to_json(),
                "periods": [p.to_json() for p in self.periods]}

    @classmethod
    def from_json(cls, value):
        if not isinstance(value, dict) or set(value) != {"coefficient", "gpl", "periods"} or not isinstance(value["periods"], list):
            raise HLogError("invalid HLOG term fields")
        return cls(value["coefficient"], GPL.from_json(value["gpl"]),
                   [GPL.from_json(p) for p in value["periods"]])


@dataclass(frozen=True, init=False)
class HLogSum:
    terms: tuple
    parameters: tuple

    def __init__(self, terms=(), *, parameters=("z",)):
        if not isinstance(parameters, (tuple, list)) or any(not isinstance(p, str) for p in parameters) or len(set(parameters)) != len(parameters):
            raise HLogError("invalid or duplicate parameter inventory")
        try: validate_scalar(q(0), variables=parameters)
        except (SeriesError, TypeError) as exc: raise HLogError(str(exc)) from exc
        if not isinstance(terms, (tuple, list)) or any(not isinstance(t, HLogTerm) for t in terms):
            raise HLogError("terms must be HLogTerm objects")
        grouped, originals = {}, {}
        for term in terms:
            try: validate_scalar(term.coefficient, variables=parameters)
            except SeriesError as exc: raise HLogError(str(exc)) from exc
            if not term.gpl.variables <= set(parameters): raise HLogError("unbound GPL parameter")
            key = (term.gpl._data, tuple(p._data for p in term.periods))
            grouped[key] = add_expr(grouped.get(key, q(0)), term.coefficient)
            originals[key] = term
        result = tuple(HLogTerm(grouped[k], originals[k].gpl, originals[k].periods)
                       for k in sorted(grouped) if as_fraction(grouped[k]) != 0)
        object.__setattr__(self, "terms", result)
        object.__setattr__(self, "parameters", tuple(sorted(parameters)))

    def add(self, other):
        if not isinstance(other, HLogSum): raise HLogError("HLOG sum required")
        return HLogSum([*self.terms, *other.terms], parameters=sorted(set(self.parameters) | set(other.parameters)))

    def scale(self, coefficient):
        """Multiply by a scalar AST, retaining the declared parameter inventory."""
        coefficient = _scalar(coefficient)
        try: validate_scalar(coefficient, variables=self.parameters)
        except SeriesError as exc: raise HLogError(str(exc)) from exc
        return HLogSum([HLogTerm(mul_expr(coefficient, term.coefficient), term.gpl, term.periods)
                        for term in self.terms], parameters=self.parameters)

    def multiply(self, other, *, max_terms=100000):
        if not isinstance(other, HLogSum): raise HLogError("HLOG sum required")
        if type(max_terms) is not int or max_terms < 1: raise HLogError("positive term budget required")
        terms = []
        for left in self.terms:
            for right in other.terms:
                coefficient = mul_expr(left.coefficient, right.coefficient)
                periods = left.periods + right.periods
                if left.gpl.is_unit or right.gpl.is_unit:
                    terms.append(HLogTerm(coefficient, right.gpl if left.gpl.is_unit else left.gpl, periods))
                else:
                    if left.gpl.path() != right.gpl.path():
                        raise HLogError("shuffle requires identical endpoint, basepoint, branch and endpoint kind")
                    path = left.gpl.to_json()
                    for word, count in shuffle_words(left.gpl.word, right.gpl.word, max_terms=max_terms-len(terms)):
                        terms.append(HLogTerm(mul_expr(coefficient, q(count)),
                                             GPL(word, path["endpoint"], basepoint=path["basepoint"],
                                                 branch=path["branch"], endpoint_kind=path["endpoint_kind"]), periods))
                        if len(terms) > max_terms: raise HLogError("HLOG product exceeds term budget")
                if len(terms) > max_terms: raise HLogError("HLOG product exceeds term budget")
        return HLogSum(terms, parameters=sorted(set(self.parameters) | set(other.parameters)))

    def to_json(self):
        return {"format": FORMAT, "kernel": KERNEL, "parameters": list(self.parameters),
                "terms": [t.to_json() for t in self.terms]}

    @classmethod
    def from_json(cls, value):
        expected = {"format", "kernel", "parameters", "terms"}
        if not isinstance(value, dict) or set(value) != expected or value["format"] != FORMAT or value["kernel"] != KERNEL:
            raise HLogError("unsupported HLOG payload or kernel convention")
        if not isinstance(value["terms"], list): raise HLogError("HLOG terms must be a list")
        return cls([HLogTerm.from_json(t) for t in value["terms"]], parameters=value["parameters"])


def shuffle_product(left, right, *, max_terms=100000):
    return left.multiply(right, max_terms=max_terms)
