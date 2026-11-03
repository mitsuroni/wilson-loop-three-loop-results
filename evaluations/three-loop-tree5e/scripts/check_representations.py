"""Independent Fraction/AST back-substitution of published coefficient views.

No producer module, SymPy algebra, HyperInt parser or shuffle helper is
imported.  Admitted generator dictionaries are explicit reused inputs.
Final reference records may be supplied only after candidate freezing.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re

PACK = Path(__file__).resolve().parents[1]
PERIODS = ("zeta_2", "zeta_3", "zeta_5", "log_2", "A4", "A6")
ZERO = (0,) * 6
LETTERS = {F(0), F(-1), F(-1, 2), F(-2)}
ONE = {((), ZERO): F(1)}
LOG = {((F(0),), ZERO): F(1)}


def clean(v):
    return {key: F(value) for key, value in v.items() if value}


def add(*vectors):
    answer = defaultdict(F)
    for vector in vectors:
        for key, value in vector.items():
            answer[key] += value
    return clean(answer)


def scale(vector, factor):
    return clean({key: value * F(factor) for key, value in vector.items()})


@lru_cache(maxsize=None)
def shuffle(left, right):
    """Count order-preserving interleavings, including repeated letters."""
    if not left:
        return {right: 1}
    if not right:
        return {left: 1}
    answer = Counter()
    for tail, multiplicity in shuffle(left[1:], right).items():
        answer[(left[0],) + tail] += multiplicity
    for tail, multiplicity in shuffle(left, right[1:]).items():
        answer[(right[0],) + tail] += multiplicity
    return dict(answer)


def multiply(left, right):
    answer = defaultdict(F)
    for (wl, pl), cl in left.items():
        for (wr, pr), cr in right.items():
            periods = tuple(a + b for a, b in zip(pl, pr))
            for word, multiplicity in shuffle(wl, wr).items():
                answer[word, periods] += cl * cr * multiplicity
    return clean(answer)


def power(vector, exponent):
    if type(exponent) is not int or exponent < 0:
        raise ValueError("nonnegative integer power required")
    answer = ONE
    for _ in range(exponent):
        answer = multiply(answer, vector)
    return answer


def scalar(vector):
    if not vector:
        return F(0)
    if set(vector) == {((), ZERO)}:
        return vector[(), ZERO]
    raise ValueError("expected an exact rational scalar")


def zeta(index):
    if index in (2, 3, 5):
        exponents = list(ZERO)
        exponents[(2, 3, 5).index(index)] = 1
        return {((), tuple(exponents)): F(1)}
    if index == 4:
        return scale(power(zeta(2), 2), F(2, 5))
    if index == 6:
        return scale(power(zeta(2), 3), F(8, 35))
    raise ValueError(f"unsupported zeta index: {index}")


def rows_vector(rows):
    """Admit either three-zeta or six-period JSON coordinate rows."""
    answer = defaultdict(F)
    for row in rows:
        if "Hlog_word" in row:
            # Original T1/T2 output encodes pi powers separately.  The stored
            # coefficient already contains the log^n -> n!*Hlog conversion.
            word = tuple(F(a) for a in row["Hlog_word"])
            pip = row["pi_power"]
            if type(pip) is not int or pip < 0 or pip % 2:
                raise ValueError("legacy pi power must be nonnegative and even")
            periods = (pip // 2, row["zeta3_power"], row["zeta5_power"], 0, 0, 0)
            if not set(word) <= LETTERS or any(type(n) is not int or n < 0 for n in periods):
                raise ValueError("invalid legacy Hlog/period row")
            answer[word, periods] += F(row["coefficient"]) * 6 ** (pip // 2)
            continue
        word = tuple(F(a) for a in row["word"])
        if not set(word) <= LETTERS:
            raise ValueError(f"unsupported Hlog alphabet: {word}")
        if "period_powers" in row:
            periods = tuple(row["period_powers"])
            if len(periods) != 6:
                raise ValueError("six period powers required")
        else:
            periods = tuple(row["zeta_powers"])
            if len(periods) != 3:
                raise ValueError("three zeta powers required")
            periods += (0, 0, 0)
        if any(type(n) is not int or n < 0 for n in periods):
            raise ValueError("nonnegative period powers required")
        answer[word, periods] += F(row["coefficient"])
    return clean(answer)


def default_atoms():
    atoms = {name: {((), tuple(int(i == j) for i in range(6))): F(1)}
             for j, name in enumerate(PERIODS)}
    atoms.update({"ell": LOG, "rho": LOG, "log_z": LOG, "logz": LOG})
    for n in (2, 3, 4, 5, 6):
        atoms[f"Z{n}"] = zeta(n)
        atoms[f"zeta_{n}"] = zeta(n)
    atoms["L2"] = atoms["log_2"]
    return atoms


class ExpressionReader:
    def __init__(self, atoms=None):
        self.atoms = default_atoms()
        self.atoms.update(atoms or {})
        self.used_names = set()

    @staticmethod
    def _name(node):
        return node.id if isinstance(node, ast.Name) else None

    def rational_node(self, node):
        return scalar(self.node(node))

    def integer_node(self, node):
        value = self.rational_node(node)
        if value.denominator != 1:
            raise ValueError("integer expected")
        return value.numerator

    def word(self, node):
        if isinstance(node, ast.Call) and self._name(node.func) == "Tuple":
            nodes = node.args
        elif isinstance(node, (ast.Tuple, ast.List)):
            nodes = node.elts
        else:
            raise ValueError("Hlog word must be an explicit tuple/list")
        word = tuple(self.rational_node(item) for item in nodes)
        if not set(word) <= LETTERS:
            raise ValueError(f"unsupported Hlog letters: {word}")
        return word

    def node(self, node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return scale(ONE, node.value)
        if isinstance(node, ast.Name):
            if node.id not in self.atoms:
                raise ValueError(f"unbound coefficient generator: {node.id}")
            self.used_names.add(node.id)
            return self.atoms[node.id]
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return scale(self.node(node.operand), -1)
            if isinstance(node.op, ast.UAdd):
                return self.node(node.operand)
        if isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Pow):
                exponent = self.integer_node(node.right)
                if self._name(node.left) in ("pi", "Pi"):
                    if exponent < 0 or exponent % 2:
                        raise ValueError("only nonnegative even pi powers are admitted")
                    return scale(power(zeta(2), exponent // 2), 6 ** (exponent // 2))
                base = self.node(node.left)
                if exponent < 0:
                    return scale(ONE, scalar(base) ** exponent)
                return power(base, exponent)
            left, right = self.node(node.left), self.node(node.right)
            if isinstance(node.op, ast.Add):
                return add(left, right)
            if isinstance(node.op, ast.Sub):
                return add(left, scale(right, -1))
            if isinstance(node.op, ast.Mult):
                return multiply(left, right)
            if isinstance(node.op, ast.Div):
                return scale(left, 1 / scalar(right))
        if isinstance(node, ast.Call) and not node.keywords:
            name = self._name(node.func)
            if name in ("zeta", "Zeta") and len(node.args) == 1:
                return zeta(self.integer_node(node.args[0]))
            if name in ("Rational", "Fraction") and len(node.args) == 2:
                return scale(ONE, self.rational_node(node.args[0]) / self.rational_node(node.args[1]))
            if name in ("log", "ln") and len(node.args) == 1:
                if self._name(node.args[0]) in ("z", "z1"):
                    return LOG
                if self.rational_node(node.args[0]) == 2:
                    return self.atoms["log_2"]
                raise ValueError("only log(z) and log(2) are admitted")
            if name in ("Hlog", "G") and len(node.args) == 2:
                if self._name(node.args[0]) not in ("z", "z1"):
                    raise ValueError("Hlog endpoint must be z")
                return {(self.word(node.args[1]), ZERO): F(1)}
        raise ValueError("unsupported coefficient expression: " + ast.dump(node)[:180])

    def read(self, expression):
        if isinstance(expression, int):
            return scale(ONE, expression)
        if not isinstance(expression, str):
            raise ValueError("coefficient must be an exact expression string")
        expression = expression.strip()
        expression = expression.replace("zeta[1,-3]", "A4").replace("zeta[1, -3]", "A4")
        expression = expression.replace("zeta[1,-5]", "A6").replace("zeta[1, -5]", "A6")
        expression = re.sub(r"zeta\[(\d+)\]", r"zeta(\1)", expression)
        expression = expression.replace("^", "**")
        return self.node(ast.parse(expression, mode="eval").body)


def coefficient_vector(row, atoms=None):
    if isinstance(row, list):
        return rows_vector(row)
    if isinstance(row, dict):
        if "hlog_terms" in row:
            return rows_vector(row["hlog_terms"])
        if "expression" in row:
            return ExpressionReader(atoms).read(row["expression"])
        raise ValueError("coefficient object has neither Hlog terms nor expression")
    return ExpressionReader(atoms).read(row)


def raw_m_dictionary():
    run = PACK / "snapshots/workflow/checks/single-opposite/runs/20260916T090538777022Z"
    atoms, sources = {}, []
    for path in sorted(run.glob("M_*_*.txt")):
        text = path.read_text(encoding="utf-8").split("=", 1)[-1].strip()
        atoms[path.stem] = ExpressionReader().read(text)
        sources.append({"path": path.relative_to(PACK).as_posix(),
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "generator": path.stem})
    return atoms, sources


def record_atoms(record, raw_m):
    definitions = record.get("basis_definitions", {})
    if isinstance(definitions, list):
        definitions = {item["name"]: item for item in definitions}
    atoms = dict(raw_m)
    reused, raw_matches = [], []
    for name, definition in definitions.items():
        terms = definition.get("hlog_terms") if isinstance(definition, dict) else None
        if terms is None:
            raise ValueError(f"generator {name} has no explicit Hlog coordinates")
        vector = rows_vector(terms)
        if name in raw_m:
            if vector != raw_m[name]:
                raise ValueError(f"published generator disagrees with raw M dictionary: {name}")
            raw_matches.append(name)
        atoms[name] = vector
        reused.append(name)
    return atoms, reused, raw_matches


def _record_map(data):
    if isinstance(data, list):
        return {record["orbit_id"]: record for record in data}
    if "orbits" in data:
        return _record_map(data["orbits"])
    if "results" in data and isinstance(data["results"], (dict, list)):
        return _record_map(data["results"])
    return data


def check_all(candidate, references):
    """Check named views and old Hlogs against the already frozen candidate."""
    candidate, references = _record_map(candidate), _record_map(references)
    if set(candidate) != set(references):
        raise ValueError("candidate/reference orbit sets differ")
    expected = {f"5e-3g-{i:02}" for i in range(1, 45)}
    if set(candidate) != expected:
        raise ValueError("complete 44-orbit P6 candidate and reference sets required")
    raw_m, raw_sources = raw_m_dictionary()
    report = {"format": "wilson-P6-representation-backsubstitution/1",
              "arithmetic": "independent Fraction AST and shuffle implementation; no producer or SymPy imports",
              "period_order": list(PERIODS),
              "admitted_dictionary_reuse": {"status": "reused evaluated dictionaries, not reintegrated",
                                           "raw_M_sources": raw_sources,
                                           "per_reference_basis_definitions": True},
              "orbits": {}, "coefficient_count": 0, "passed": True,
              "formal": "not_run"}
    for oid in sorted(candidate):
        cand, ref = candidate[oid], references[oid]
        atoms, admitted, raw_matches = record_atoms(ref, raw_m)
        candidate_coeffs = cand["laurent"]["coefficients"]
        expressions = ref["laurent"]["coefficients"]
        if set(candidate_coeffs) != {str(k) for k in range(-6, 1)} or set(expressions) != set(candidate_coeffs):
            raise ValueError(f"complete seven-position Laurent sets required: {oid}")
        reference_hlogs = ref.get("hlog_coefficients", {})
        rows, used = {}, set()
        for k in range(-6, 1):
            key = str(k)
            target = coefficient_vector(candidate_coeffs[key])
            reader = ExpressionReader(atoms)
            view = expressions[key]
            named = coefficient_vector(view, atoms) if isinstance(view, (dict, list)) else reader.read(view)
            used.update(reader.used_names)
            named_ok = named == target
            hlog_ok = None if key not in reference_hlogs else rows_vector(reference_hlogs[key]) == target
            rows[key] = {"named_backsubstitution_equals_frozen_candidate": named_ok,
                         "reference_hlog_equals_frozen_candidate": hlog_ok,
                         "canonical_term_count": len(target)}
            report["coefficient_count"] += 1
            if not named_ok or hlog_ok is False:
                report["passed"] = False
        report["orbits"][oid] = {"coefficients": rows,
                                 "admitted_basis_generators": sorted(admitted),
                                 "basis_M_matches_raw_dictionary": sorted(raw_matches),
                                 "used_named_atoms": sorted(used)}
    report["orbit_count"] = len(report["orbits"])
    report["negative_controls"] = self_check()
    report["passed"] = report["passed"] and all(report["negative_controls"].values())
    return report


def self_check():
    reader = ExpressionReader()
    checks = {
        "shuffle_multiplicity": reader.read("log(z)^2") == reader.read("2*Hlog(z,[0,0])"),
        "mixed_shuffle": reader.read("Hlog(z,[0])*Hlog(z,[-1])") == reader.read("Hlog(z,[0,-1])+Hlog(z,[-1,0])"),
        "even_zeta_reduction": reader.read("zeta(4)+pi^4/90") == reader.read("4*zeta(2)^2/5"),
        "wide_letters": reader.read("Hlog(z,[-1/2,-2])") == {((F(-1, 2), F(-2)), ZERO): F(1)},
        "mutation_detected": reader.read("Hlog(z,[0])+1") != reader.read("Hlog(z,[0])"),
    }
    for label, expression in (("unbound_generator_rejected", "UnadmittedM+1"),
                              ("unexpected_letter_rejected", "Hlog(z,[1])"),
                              ("nonrational_denominator_rejected", "1/Hlog(z,[0])")):
        try:
            reader.read(expression)
            checks[label] = False
        except ValueError:
            checks[label] = True
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--references", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        result = self_check()
        print(json.dumps(result, indent=2))
        return 0 if all(result.values()) else 1
    if args.candidate is None or args.references is None:
        parser.error("--candidate and --references are required after candidate freezing")
    result = check_all(json.loads(args.candidate.read_text(encoding="utf-8")),
                       json.loads(args.references.read_text(encoding="utf-8")))
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
