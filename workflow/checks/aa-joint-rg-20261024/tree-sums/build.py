"""Fresh current 6e and 5e sums for the AA double-pole assembly.

Run with D:/py/python.exe -X utf8 workflow/checks/aa-joint-rg-20261024/tree-sums/build.py.
Only this directory is written. No old writer or external CAS is run.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from fractions import Fraction
import sys

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
POWERS = tuple(range(-6, -1))
SOURCES = {}
CHECKS = []


def binding(path):
    p = path if isinstance(path, Path) else ROOT / path
    key = p.relative_to(ROOT).as_posix()
    digest = hashlib.sha256(p.read_bytes()).hexdigest()
    if key in SOURCES:
        assert SOURCES[key] == digest, ("input changed during run", key)
    SOURCES[key] = digest
    return p


def read(path):
    return json.loads(binding(path).read_text(encoding="utf-8"))


def check(name, passed, **extra):
    CHECKS.append(dict(name=name, passed=bool(passed), **extra))
    assert passed, name


lib = binding("workflow/checks/single-opposite/hlog_algebra.py")
spec = importlib.util.spec_from_file_location("aa_tree_hlog", lib)
ha = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ha)


def vec(terms):
    out = {}
    for term in terms:
        periods = term.get("zeta_powers", term.get("period_powers"))
        assert periods is not None and not any(periods[3:]), term
        word = tuple(Fraction(x) for x in term["word"])
        assert all(x in (0, -1) for x in word), term
        key = (tuple(int(x) for x in word), tuple(periods[:3]))
        out[key] = out.get(key, Fraction(0)) + Fraction(term["coefficient"])
    return ha.clean(out)


def parse_coeff(d, k):
    c = d["laurent"]["coefficients"][str(k)]
    base = f"/laurent/coefficients/{k}"
    if isinstance(c, dict) and "hlog_terms" in c:
        return vec(c["hlog_terms"]), base + "/hlog_terms"
    if "hlog_coefficients" in d:
        terms = d["hlog_coefficients"][str(k)]
        if all("word" in t for t in terms):
            return vec(terms), f"/hlog_coefficients/{k}"
    if isinstance(c, dict):
        assert "expression" in c
        expr = (parse_mathematica(c["expression"])
                if c["language"].startswith("Wolfram")
                else sp.sympify(c["expression"]))
        expr = expr.replace(lambda t: t.is_Function and t.func.__name__ == "Zeta",
                            lambda t: sp.zeta(*t.args))
        expr = expr.subs(sp.Symbol("Pi"), sp.pi)
        expr = expr.replace(lambda t: t.func == sp.log and str(t.args[0]) == "z",
                            lambda t: sp.Symbol("rho"))
        pointer = base + "/expression"
    else:
        expr = sp.sympify(c, locals={"E": sp.Symbol("E")})
        pointer = base
    atoms = {sp.Symbol("ell"): ha.LOG, sp.Symbol("rho"): ha.LOG}
    atoms.update(BASIS_ATOMS)
    return ha.from_sympy(expr, atoms), pointer


M01, M11, M21 = sp.symbols("M01 M11 M21")
z2, z3 = sp.symbols("zeta_2 zeta_3")
moment_source = "workflow/checks/single-opposite/runs/20260916T090538777022Z"
moments = {name: ha.parse(binding(f"{moment_source}/{old}.txt"))
           for name, old in [(M01, "M_0_1"), (M11, "M_1_1"), (M21, "M_2_1")]}
BASIS_ATOMS = {**moments, z2: ha.period(2), z3: ha.period(3)}
for old, new in [("M_0_1", M01), ("M_1_1", M11), ("M_2_1", M21)]:
    BASIS_ATOMS[sp.Symbol(old)] = moments[new]
compact_basis = [sp.Integer(1), z2, M01, z3, M11, z2**2, z2*M01, M01**2, M21]
compact_vectors = [ha.from_sympy(expr, BASIS_ATOMS) for expr in compact_basis]


def compact(v):
    keys = sorted(set(v).union(*(set(x) for x in compact_vectors)))
    mat = sp.Matrix([[x.get(key, 0) for x in compact_vectors] for key in keys])
    rhs = sp.Matrix([v.get(key, 0) for key in keys])
    sol, free = mat.gauss_jordan_solve(rhs)
    assert not free
    expr = sp.expand(sum(c*x for c, x in zip(sol, compact_basis)))
    assert ha.from_sympy(expr, BASIS_ATOMS) == v
    return expr


def selected_6e():
    rows = []
    for p in sorted((ROOT / "benchmarks/2ne/families").glob("*/catalog.json")):
        d = read(p)
        for i, r in enumerate(d["items"]):
            if not r["orbit_id"].startswith("6e-"):
                continue
            key = "latest_result" if "latest_result" in r else "result"
            rows.append(dict(id=r["orbit_id"], path=r[key],
                             selector=p.relative_to(ROOT).as_posix(),
                             selector_pointer=f"/items/{i}/{key}",
                             declared_sha256=r.get("result_sha256")))
    check("6e current catalog coverage", sorted(r["id"] for r in rows)
          == [f"6e-{i:02}" for i in range(1, 40)])
    return sorted(rows, key=lambda r: r["id"])


def selected_5e():
    p = "benchmarks/tgv/5e-3g/current.json"
    d = read(p)
    rows = [dict(id=oid, path=f"benchmarks/tgv/5e-3g/{r['result']}",
                 selector=p, selector_pointer=f"/orbits/{oid}/result",
                 declared_sha256=r["sha256"])
            for oid, r in sorted(d["orbits"].items())]
    check("5e current selector coverage", [r["id"] for r in rows]
          == [f"5e-3g-{i:02}" for i in range(1, 45)])
    return rows


def assemble(name, rows):
    total = {k: {} for k in POWERS}
    ledger = []
    selected_vectors = {}
    for r in rows:
        d = read(r["path"])
        digest = SOURCES[r["path"]]
        if r["declared_sha256"]:
            check("selected result hash", digest == r["declared_sha256"], id=r["id"])
        coeffs = {}
        selected_vectors[r["id"]] = {}
        for k in POWERS:
            v, pointer = parse_coeff(d, k)
            total[k] = ha.add(total[k], v)
            selected_vectors[r["id"]][k] = v
            coeffs[str(k)] = dict(pointer=pointer, hlog_terms=ha.serialize(v))
        ledger.append({**r, "sha256": digest, "assembly_weight": "1",
                       "normalization": d.get("normalization"),
                       "I_definition": d.get("I_definition"),
                       "coefficients": coeffs})
    # Derive fresh compact answers from moment Hlog tables, before comparison.
    coeffs = {str(k): str(compact(total[k])) for k in POWERS}
    for k in POWERS:
        check("compact coefficient back-substitution", True, family=name, power=k)
    previous_path = f"benchmarks/recursion-audit/results/{name}-effective-sum/r001.json"
    previous = read(previous_path)
    comparisons = []
    for k in POWERS:
        reference = vec(previous["hlog_coefficients"][str(k)])
        delta = ha.add(total[k], ha.scale(reference, -1))
        check("existing independent recursion sum", not delta, family=name, power=k)
        comparisons.append(dict(power=k, residual=ha.serialize(delta), passed=not delta))
    # Retain symmetric sources: omitting a nonzero source must change the sum.
    control_id = "6e-02" if name == "6e" else "5e-3g-35"
    control_power = -6 if name == "6e" else -4
    omitted = ha.add(total[control_power], ha.scale(selected_vectors[control_id][control_power], -1))
    residual = ha.add(omitted, ha.scale(total[control_power], -1))
    check("omitted symmetric source detected", bool(residual), family=name,
          omitted_id=control_id, power=control_power, residual=ha.serialize(residual))
    return dict(integral_count=len(rows), powers=list(POWERS), coefficients=coeffs,
                hlog_coefficients={str(k): ha.serialize(total[k]) for k in POWERS},
                source_ids=[r["id"] for r in rows], sources=ledger,
                assembly_weight="1 per selected physical I; no symmetry deduplication",
                comparison=dict(path=previous_path, sha256=SOURCES[previous_path],
                                pointer="/hlog_coefficients", rows=comparisons))


def main():
    binding(Path(__file__))
    for p in ["benchmarks/recursion-audit/sum-6e.md", "benchmarks/recursion-audit/sum-5e.md",
              "benchmarks/2ne/theory/unified-weight-six.md"]:
        binding(p)
    families = {"6e": assemble("6e", selected_6e()), "5e": assemble("5e", selected_5e())}
    check("source files stable throughout run", all(
        hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h for p, h in SOURCES.items()))
    result = dict(format="wilson-aa-tree-sums/1", date="2026-10-24", task="WILSON-2.13",
                  passed=True, formal="not_run", expression_language="SymPy exact expressions",
                  convention=dict(stage="scale-stripped bare connected exponent w3, CF CA^2 coefficient",
                      rho="log(z), never the scale logarithm L", M01="M_{0,1}(z)",
                      M11="M_{1,1}(z)", M21="M_{2,1}(z)",
                      zeta_2="Riemann zeta(2)", zeta_3="Riemann zeta(3)",
                      hlog="G(word;z), first letter outermost, G(0^n;z)=log(z)^n/n!",
                      normalization="a1/fixed color, C4 orbit multiplicity, crossing and exactly one Euler compensation already included in selected physical I; no additional web or log(W) subtraction",
                      outside="Only common CF CA^2, coupling power and S3 scale sum are outside these family sums"),
                  families=families, checks=CHECKS, source_sha256=SOURCES,
                  independence="Fresh algebraic sum of selected per-diagram coefficients. Existing independent recursion family sums are comparison targets only, loaded after fresh compact answers. Original integrations were reused, not rerun. No RG/cusp target was used.",
                  limitations=["Only powers -6 through -2 are assembled here.",
                               "This does not separately prove individual cancellation directions.",
                               "No all-epsilon, complex-branch or formal integral proof is claimed."])
    (HERE / "tree-sums.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(passed=True, source_files=len(SOURCES), checks=len(CHECKS),
                         families={k: dict(count=v['integral_count'], coefficients=v['coefficients'])
                                   for k,v in families.items()}), indent=2))


if __name__ == "__main__":
    main()
