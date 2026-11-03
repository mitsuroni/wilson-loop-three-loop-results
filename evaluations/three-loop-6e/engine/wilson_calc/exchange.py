"""Structural reduction of ordered pure-exchange integrals.

The input is the existing ``ordered-parameters/1`` J mapping.  Graph IDs,
result tables, C, colour projections and Euler factors are deliberately absent
from this module.  Reducing a domain is not the same as evaluating its remaining
opposite-edge kernels, and the returned statuses preserve that distinction.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from functools import lru_cache
import hashlib
import itertools
import json
import math
import re
from typing import Mapping

import sympy as sp

from preprocessing.ordered import arithmetic, integral_ast

EPS = sp.Symbol("eps")
Z1 = sp.Symbol("z1", positive=True)
AT, BT, AS, BS, BAT, BAS = sp.symbols("AT BT AS BS BAT BAS")
WEIGHT_SYMBOLS = (AT, BT, AS, BS, BAT, BAS)

# Source pointers describe mathematical/algorithm provenance, not runtime
# imports or proof receipts.  In particular, no historical result is loaded.
SOURCES = (
    {"path": "workflow/checks/2ne-generation/audit.py", "sha256":
     "45f28ccff94785c8f790678becd4c1a14bb61f031d0730c7543522f681039979"},
    {"path": "workflow/checks/2ne-generation/ordered_words.py", "sha256":
     "c0429825338175a5080b4c3f4b8d1d1d7ab20540e9a408f5bdb2c21a56a2eafd"},
    {"path": "workflow/checks/single-opposite/reduce_inputs.py", "sha256":
     "ee53e5320bbbd6d5b0142adc4ae3fa046198052ce4ba86a1a23de4ed13a9f6fb"},
    {"path": "workflow/checks/empty-family/empty_math.py", "sha256":
     "c844a76f2eb7a1591d7cccd732821be488890d1196704ca397dce6eb0cd49373"},
)


class ExchangeResourceLimit(RuntimeError):
    """A finite word computation exceeded its configured term budget."""


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def arithmetic_to_sympy(text, variables=()):
    """Translate the existing closed AST; never parse source with sympify/eval."""
    symbols = {name: sp.Symbol(name, positive=True) for name in variables}
    symbols.update(eps=EPS, z1=Z1)

    def visit(node):
        if node["kind"] == "rational":
            return sp.Rational(node["numerator"], node["denominator"])
        if node["kind"] == "var":
            return symbols[node["name"]]
        args = [visit(arg) for arg in node["args"]]
        operations = {"add": lambda: args[0] + args[1],
                      "sub": lambda: args[0] - args[1],
                      "mul": lambda: args[0] * args[1],
                      "div": lambda: args[0] / args[1],
                      "neg": lambda: -args[0],
                      "rpow": lambda: args[0] ** args[1]}
        return operations[node["op"]]()

    return visit(arithmetic(text, set(symbols)))


@lru_cache(maxsize=4096)
def _shuffle_tuple(left, right):
    if not left:
        return ((right, 1),)
    if not right:
        return ((left, 1),)
    counts = Counter()
    for word, count in _shuffle_tuple(left[1:], right):
        counts[left[0] + word] += count
    for word, count in _shuffle_tuple(left, right[1:]):
        counts[right[0] + word] += count
    return tuple(sorted(counts.items()))


def shuffle_words(left, right):
    """Shuffle A/B words, including multiplicities of repeated letters."""
    _require(isinstance(left, str) and isinstance(right, str)
             and set(left + right) <= {"A", "B"}, "unknown endpoint word")
    return dict(_shuffle_tuple(left, right))


def reflect_word(word):
    _require(isinstance(word, str) and set(word) <= {"A", "B"},
             "unknown endpoint word")
    return word[::-1].translate(str.maketrans("AB", "BA"))


def compress_gaps(gaps, *, max_terms=10000):
    """Convert interval words between marked endpoints into lower-end words.

    For k marked points a key contains k lower-end words and one whole-edge
    word.  The formula is path composition followed by shuffle; this operation
    keeps the marked points and does not integrate an opposite-edge kernel.
    """
    gaps = tuple(gaps)
    _require(bool(gaps) and all(isinstance(w, str) and set(w) <= {"A", "B"}
                               for w in gaps), "invalid gap words")
    _require(type(max_terms) is int and max_terms > 0, "positive word term budget required")
    count_marks = len(gaps) - 1
    result = {tuple("" for _ in range(count_marks + 1)): 1}
    for gap_index, word in enumerate(gaps):
        next_result = Counter()
        cuts = (0,) if gap_index == 0 else range(len(word) + 1)
        for cut in cuts:
            left, right = word[:cut][::-1], word[cut:]
            for words, coefficient in result.items():
                left_terms = [(words, coefficient * (-1) ** cut)]
                if gap_index:
                    left_terms = []
                    for new_word, multiple in shuffle_words(words[gap_index - 1], left).items():
                        key = list(words)
                        key[gap_index - 1] = new_word
                        left_terms.append((tuple(key), coefficient * (-1) ** cut * multiple))
                for key, coefficient_left in left_terms:
                    for new_word, multiple in shuffle_words(key[gap_index], right).items():
                        new_key = list(key)
                        new_key[gap_index] = new_word
                        next_result[tuple(new_key)] += coefficient_left * multiple
                        if len(next_result) > max_terms:
                            raise ExchangeResourceLimit("endpoint word term budget exceeded")
        result = {key: value for key, value in next_result.items() if value}
    return result


def _unsupported(j, code, detail, descriptor=None):
    return {"status": "unsupported", "descriptor": descriptor,
            "residuals": [{"kind": code, "detail": detail, "J": deepcopy(j)}],
            "sources": deepcopy(SOURCES), "formal_state": "not_run"}


def _edge_map(names, endpoint_edges, variable_bindings):
    if endpoint_edges is not None:
        _require(set(endpoint_edges) == set(names), "endpoint edge inventory mismatch")
        _require(all(type(edge) is int and 1 <= edge <= 4 for edge in endpoint_edges.values()),
                 "endpoint edges must be integers from 1 to 4")
        return dict(endpoint_edges)
    reverse = {value: key for key, value in (variable_bindings or {}).items()}
    _require(len(reverse) == len(variable_bindings or {}), "variable bindings are not injective")
    result = {}
    for name in names:
        match = re.fullmatch(r"s([1-4])([1-9][0-9]*)", reverse.get(name, name))
        if match is None:
            return None
        result[name] = int(match[1])
    return result


def _check_topology(topology, propagators, names, edges, variable_bindings):
    if topology is None:
        return
    _require(not topology.get("internal_vertices"), "pure exchange cannot have internal vertices")
    ordinals = Counter()
    colour_names = {}
    for vertex in topology.get("vertices", []):
        edge = int(vertex["line"])
        ordinals[edge] += 1
        name = f"s{edge}{ordinals[edge]}"
        name = (variable_bindings or {}).get(name, name)
        _require(name in names and edges[name] == edge, "topology endpoint binding mismatch")
        _require(vertex["color"] not in colour_names, "repeated topology colour endpoint")
        colour_names[vertex["color"]] = name
    _require(set(colour_names.values()) == set(names), "topology endpoint inventory mismatch")
    links = topology.get("propagators", [])
    _require(len(links) == len(propagators), "topology propagator inventory mismatch")
    for observed, link in zip(propagators, links):
        try:
            expected = [colour_names[colour] for colour in link["color"]]
        except KeyError as exc:
            raise ValueError("unbound topology colour endpoint") from exc
        _require(len(expected) == 2 and set(expected) == set(observed["variables"]),
                 "topology pairing mismatch")


def analyze_exchange(j, *, loop_order=None, endpoint_edges=None, topology=None,
                     variable_bindings=None, max_word_terms=10000):
    """Derive a domain/kernel descriptor from a prepared J, without orbit lookup.

    Invalid binding/representation data raises ValueError.  A well-formed
    integral outside the known geometric rules returns ``unsupported`` with
    its original J intact.  ``reduced`` certifies only the stated research-level
    domain bridge and endpoint-word representation, not an HLOG evaluation.
    """
    integral_ast(j)
    names = j["variables"]
    factors = j["propagator_polynomials"]
    n = len(factors)
    _require(n > 0 and len(names) == 2 * n, "pure exchange needs two endpoints per propagator")
    _require(loop_order is None or (type(loop_order) is int and loop_order == n),
             "loop order differs from propagator count")
    edges = _edge_map(names, endpoint_edges, variable_bindings)
    if edges is None:
        return _unsupported(j, "missing_endpoint_geometry", "Supply endpoint_edges for non-geometric variable names.")
    symbols = {name: sp.Symbol(name, positive=True) for name in names}
    letters, propagators = {}, []
    adjacent_z_count = 0
    for factor_id, raw in factors.items():
        expression = arithmetic_to_sympy(raw, names)
        used = [name for name in names if symbols[name] in expression.free_symbols]
        if len(used) != 2 or expression.free_symbols - {Z1, *(symbols[name] for name in used)}:
            return _unsupported(j, "unsupported_propagator", f"{factor_id} is not a two-endpoint kernel.")
        used.sort(key=lambda name: (edges[name], names.index(name)))
        left, right = (symbols[name] for name in used)
        kernel = (1 - left) * (1 - right) + Z1 * left * right
        kernel_tilde = Z1 * (1 - left) * (1 - right) + left * right
        record = {"id": factor_id, "variables": used, "original_polynomial": raw}
        if sp.expand(expression - kernel) == 0 or sp.expand(expression - kernel_tilde) == 0:
            channel = "".join(str(edges[name]) for name in used)
            if channel not in ("13", "24"):
                return _unsupported(j, "geometry_mismatch", f"{factor_id} has an opposite kernel on non-opposite edges.")
            record.update(kind="opposite", channel=channel,
                          kernel="D" if sp.expand(expression - kernel) == 0 else "D_tilde")
            assignment = {name: "O" for name in used}
        else:
            matches = [(a, b, q) for a, b in itertools.permutations(used)
                       for q in (sp.Integer(1), Z1)
                       if sp.expand(expression - q * symbols[a] * (1 - symbols[b])) == 0]
            if len(matches) != 1:
                return _unsupported(j, "unsupported_propagator", f"{factor_id} has no supported adjacent/opposite kernel.")
            a, b, q = matches[0]
            if (edges[a] - edges[b]) % 4 not in (1, 3):
                return _unsupported(j, "geometry_mismatch", f"{factor_id} has an adjacent kernel on non-adjacent edges.")
            assignment = {a: "A", b: "B"}
            record.update(kind="adjacent", A=a, B=b, q=str(q))
            adjacent_z_count += int(q == Z1)
        if set(letters) & set(assignment):
            return _unsupported(j, "repeated_endpoint", "An endpoint appears in more than one propagator.")
        letters.update(assignment)
        propagators.append(record)
    _require(set(letters) == set(names), "uncovered endpoint")
    opposite = [item for item in propagators if item["kind"] == "opposite"]
    r = len(opposite)
    expected_numerator = (-(1 + Z1)) ** r * Z1 ** adjacent_z_count
    if sp.expand(arithmetic_to_sympy(j["numerator"], names) - expected_numerator) != 0:
        return _unsupported(j, "unsupported_numerator", "The numerator is outside the pure-exchange contraction rule.")
    chains = []
    for edge in range(1, 5):
        members = {name for name in names if edges[name] == edge}
        child, tops = {}, []
        for name in members:
            low, high = j["bounds"][name]
            lower = arithmetic_to_sympy(low, names)
            upper = arithmetic_to_sympy(high, names)
            if lower != 0:
                return _unsupported(j, "unsupported_domain", "Endpoint-word rules require zero lower bounds.")
            if upper == 1:
                tops.append(name)
            elif str(upper) in members and upper == symbols[str(upper)]:
                parent = str(upper)
                if parent in child:
                    return _unsupported(j, "unsupported_domain", "An edge is not one ordered endpoint chain.")
                child[parent] = name
            else:
                return _unsupported(j, "unsupported_domain", "An upper bound is not a parent on the same edge.")
        if members and len(tops) != 1:
            return _unsupported(j, "unsupported_domain", "Each occupied edge needs exactly one chain.")
        descending = []
        if members:
            current = tops[0]
            while True:
                _require(current not in descending, "cyclic endpoint chain")
                descending.append(current)
                if current not in child:
                    break
                current = child[current]
        _require(set(descending) == members, "incomplete endpoint chain")
        chain = descending[::-1]
        word = "".join(letters[name] for name in chain)
        marks = [name for name in chain if letters[name] == "O"]
        gaps = word.split("O")
        compression = compress_gaps(gaps, max_terms=max_word_terms)
        chains.append({"edge": edge, "variables_low_to_high": chain, "word": word,
                       "marked_endpoints": marks, "gap_words": gaps,
                       "compression": [{"words_at_marked_endpoints": list(key[:-1]),
                                        "whole_edge_word": key[-1], "coefficient": value}
                                       for key, value in sorted(compression.items())]})
    _check_topology(topology, propagators, names, edges, variable_bindings)
    channels = {}
    for channel in ("13", "24"):
        left, right = [chains[int(edge) - 1]["marked_endpoints"] for edge in channel]
        partner = {item["variables"][0]: item["variables"][1] for item in opposite
                   if item["channel"] == channel}
        _require(set(partner) == set(left) and set(partner.values()) == set(right),
                 "opposite channel endpoint inventory mismatch")
        channels[channel] = {"count": len(left), "left": left, "right": right,
                             "permutation": [right.index(partner[name]) + 1 for name in left]}
    remaining_bounds = {}
    for chain in chains:
        marks = chain["marked_endpoints"]
        for index, name in enumerate(marks):
            remaining_bounds[name] = ["0", marks[index + 1] if index + 1 < len(marks) else "1"]
    m = n - r
    prefactor = (-(1 + Z1)) ** r * Z1 ** (adjacent_z_count * EPS) / EPS ** (2 * m)
    source_hash = hashlib.sha256(json.dumps(j, sort_keys=True, ensure_ascii=False,
                                          separators=(",", ":")).encode("utf-8")).hexdigest()
    descriptor = {"loop_order": n, "adjacent_count": m, "opposite_count": r,
                  "adjacent_z_count": adjacent_z_count, "endpoint_edges": edges,
                  "letters": letters, "propagators": propagators, "edges": chains,
                  "channels": channels,
                  "channel_partition": sorted(value["count"] for value in channels.values() if value["count"]),
                  "integration_order_outer_to_inner": list(j["integration_order_outer_to_inner"]),
                  "remaining_order_outer_to_inner": [name for name in j["integration_order_outer_to_inner"]
                                                       if name in remaining_bounds],
                  "remaining_bounds": remaining_bounds, "prefactor": prefactor,
                  "input_J_sha256": source_hash,
                  "word_normalization": "Each A/B density includes eps; interval words are ordered from low to high."}
    bridge = {"kind": "pure_exchange_endpoint_words", "valid_domain": "z1 > 0, Re(eps) > 0",
              "reason": "Distinct endpoints: the absolute-value bound on the enclosing cube factors into integrable adjacent densities and two-dimensional D kernels. Grouping preserves all ordered O endpoints.",
              "epsilon_expansion": "Not established by this domain reduction; each evaluator must supply its own meromorphic/holomorphic continuation.",
              "source": "benchmarks/2ne/theory/generation/ordered-words.md",
              "formal_state": "not_run"}
    return {"status": "reduced", "descriptor": descriptor, "prefactor": prefactor,
            "domain_bridge": bridge, "sources": deepcopy(SOURCES), "residuals": [],
            "formal_state": "not_run"}


def _gamma_ratio(a, b):
    return sp.gamma(1 + a * EPS) * sp.gamma(1 + b * EPS) / sp.gamma(1 + (a + b) * EPS)


def whole_word_value(word):
    """Exact *normalised* whole-edge word, or None for an unsupported word.

    Pure powers work at any length.  A^p B^q follows from a Beta integral.
    A^p B A^q follows by binomially expanding its upper A interval; reflection
    also covers words with one A.  All remaining length-four words follow
    from a Beta integral and shuffle products.  The identities hold first
    for Re(eps)>0 on the original ordered finite domain; the Gamma expressions
    give their meromorphic continuation.  Longer unsupported words remain.
    """
    _require(isinstance(word, str) and set(word) <= {"A", "B"}, "unknown endpoint word")
    if not word:
        return sp.Integer(1)
    if len(set(word)) == 1:
        return sp.Rational(1, math.factorial(len(word)))
    if word.count("A") == 1 and word.count("B") != 1:
        return whole_word_value(reflect_word(word))
    if word.count("B") == 1:
        p, q = word.index("B"), len(word) - word.index("B") - 1
        return sp.Add(*[(-1) ** k * sp.binomial(q, k) * _gamma_ratio(p + k, 1)
                        / (sp.factorial(p) * sp.factorial(q)) for k in range(q + 1)])
    match = re.fullmatch(r"(A+)(B+)", word)
    if match:
        p, q = len(match[1]), len(match[2])
        return _gamma_ratio(p, q) / (sp.factorial(p) * sp.factorial(q))
    if len(word) == 4 and word.count("A") == 2:
        u = _gamma_ratio(1, 1)  # W_AB; W_BA=1-W_AB.
        aabb = _gamma_ratio(2, 2) / 4
        # Integrate the two lower B and two upper A variables first:
        # W_BBAA = eps/2 * integral_0^1 x^(eps-1)(1-x^eps)
        #                               * (1-(1-x)^eps)^2 dx.
        # Expanding those two factors and using Beta gives the expression
        # below.  Shuffle then gives AB^2, BA^2 and AB*BA; reflection of the
        # whole unit interval identifies W_ABBA=W_BAAB.
        bbaa = sp.Rational(1, 4)-u+_gamma_ratio(1, 2)-aabb
        abab = u**2/2-2*aabb
        baba = (1-u)**2/2-2*bbaa
        abba = (u*(1-u)-abab-baba)/4
        return {"ABAB": abab, "ABBA": abba, "BAAB": abba,
                "BABA": baba, "BBAA": bbaa}[word]
    return None


def _descriptor(analysis):
    descriptor = analysis.get("descriptor", analysis)
    _require(isinstance(descriptor, Mapping) and "opposite_count" in descriptor,
             "a reduced exchange descriptor is required")
    return descriptor


def evaluate_empty(analysis):
    """Evaluate an all-adjacent J as an exact Gamma/power expression."""
    descriptor = _descriptor(analysis)
    if descriptor["opposite_count"]:
        return {"status": "unsupported", "residuals": [{"kind": "opposite_kernels_remaining",
                "descriptor": descriptor}], "formal_state": "not_run"}
    values, pending = [], []
    for edge in descriptor["edges"]:
        value = whole_word_value(edge["word"])
        if value is None:
            pending.append({"kind": "unsupported_whole_edge_word", "edge": edge["edge"],
                            "word": edge["word"], "bounds": [0, 1],
                            "normalization": "eps^length times the original iterated density"})
        else:
            values.append(value)
    if pending:
        return {"status": "partial", "known_factor": descriptor["prefactor"] * sp.prod(values),
                "residuals": pending, "descriptor": descriptor, "formal_state": "not_run"}
    expression = descriptor["prefactor"] * sp.prod(values)
    return {"status": "computed", "J": expression, "expression": expression,
            "epsilon_lower_bound": -2 * descriptor["loop_order"],
            "continuation": "Finite Gamma ratios regular at eps=0; explicit eps powers carry a sufficient pole bound.",
            "domain": "z1 > 0", "sources": deepcopy(SOURCES), "residuals": [],
            "formal_state": "not_run"}


def _interval_weight(word, side, a, b, ba):
    if not word:
        return sp.Integer(1)
    if set(word) == {"A"}:
        return (a if side == "lower" else 1 - a) ** len(word) / sp.factorial(len(word))
    if set(word) == {"B"}:
        return (1 - b if side == "lower" else b) ** len(word) / sp.factorial(len(word))
    g11 = _gamma_ratio(1, 1)
    rules = {("AB", "lower"): a * (1 - b) - ba, ("BA", "lower"): ba,
             ("AB", "upper"): g11 - a + ba,
             ("BA", "upper"): a + b - a * b - g11 - ba}
    return rules.get((word, side))


def reflect_single_weight(expression):
    """Reflect the complete single-opposite weight, retaining BA differences."""
    g11 = _gamma_ratio(1, 1)
    return sp.expand(expression.subs({AT: BT, BT: AT, AS: BS, BS: AS,
                                     BAT: AT + BT - AT * BT - g11 - BAT,
                                     BAS: AS + BS - AS * BS - g11 - BAS}, simultaneous=True))


def single_weight_terms(expression):
    """Read powers t^(a eps)(1-t)^(b eps)s^(c eps)(1-s)^(d eps).

    Raises ValueError while an ordered BA insertion remains; it must not be
    silently discarded before physical crossing or independent integration.
    """
    _require(not expression.has(BAT, BAS), "ordered BA insertion still requires integration or an established cyclic identity")
    try:
        polynomial = sp.Poly(sp.expand(expression), AT, BT, AS, BS)
    except sp.PolynomialError as exc:
        raise ValueError("single-opposite weight is not a polynomial in endpoint powers") from exc
    return [{"powers": tuple(powers), "coefficient": coefficient}
            for powers, coefficient in polynomial.terms() if coefficient != 0]


def single_opposite_weights(analysis):
    """Extract actual single-kernel endpoint weights; no graph-specific atlas."""
    descriptor = _descriptor(analysis)
    if descriptor["opposite_count"] != 1:
        return {"status": "unsupported", "residuals": [{"kind": "not_single_opposite",
                "descriptor": descriptor}], "formal_state": "not_run"}
    opposite = next(item for item in descriptor["propagators"] if item["kind"] == "opposite")
    if opposite["kernel"] not in ("D", "D_tilde"):
        return {"status": "unsupported", "residuals": [{"kind": "unknown_single_kernel",
                "descriptor": descriptor}], "formal_state": "not_run"}
    transport = "z-inversion" if opposite["kernel"] == "D_tilde" else "identity"
    t_name, s_name = opposite["variables"]
    weight, pending = sp.Integer(1), []
    for edge in descriptor["edges"]:
        if not edge["marked_endpoints"]:
            value = whole_word_value(edge["word"])
        else:
            mark = edge["marked_endpoints"][0]
            a, b, ba = (AT, BT, BAT) if mark == t_name else (AS, BS, BAS)
            lower, upper = edge["gap_words"]
            first = _interval_weight(lower, "lower", a, b, ba)
            last = _interval_weight(upper, "upper", a, b, ba)
            value = None if first is None or last is None else first * last
        if value is None:
            pending.append({"kind": "unsupported_interval_word", "edge": edge["edge"],
                            "word": edge["word"], "gap_words": edge["gap_words"]})
        else:
            weight *= value
    weight = sp.expand(weight)
    result = {"status": "reduced", "normalized_weight": weight,
              "prefactor": descriptor["prefactor"], "variables": {"t": t_name, "s": s_name},
              "kernel": opposite["kernel"]+"_z(t,s)^(eps-1)",
              "transport": transport, "terms": [], "residuals": pending,
              "sources": deepcopy(SOURCES), "formal_state": "not_run"}
    if pending:
        result["status"] = "partial"
        result["known_weight_factor"] = result.pop("normalized_weight")
    elif weight.has(BAT, BAS):
        result["status"] = "partial"
        result["residuals"] = [{"kind": "ordered_BA_insertion",
            "definition": "B_A(x)=eps^2 integral_(0<u<v<x) (1-u)^(eps-1) v^(eps-1) du dv",
            "detail": "BAT=B_A(t), BAS=B_A(s). Keep the raw J insertion; do not delete it using an I-only identity."}]
    else:
        result["terms"] = single_weight_terms(weight)
    return result
