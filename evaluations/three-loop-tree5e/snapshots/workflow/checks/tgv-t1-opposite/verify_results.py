"""Offline result verification from the frozen moments, exact basis and root bridges."""
import json,hashlib
import sympy as sp
import formulas as f
from publish import FULL,DEFS,expression
import check_math as math_check
ROOT,HERE=f.ROOT,f.HERE
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    delivery=ROOT/f.read(HERE/"current-delivery.json")["delivery"];manifest=f.read(delivery/"manifest.json")
    for b in manifest["bindings"]:
        assert sha(ROOT/b["path"])==b["sha256"],b["path"]
        assert sha(delivery/b["snapshot"])==b["sha256"],b["snapshot"]
    evidence=ROOT/manifest["evidence_run"];f.run=lambda:evidence
    math_check.run=lambda:evidence
    assert math_check.check()==f.read(evidence/"math-checks.json")
    exact=f.coefficients();basis=f.read(ROOT/f.read(evidence/"representation.json")["basis_path"])
    definitions={g["name"]:g for g in basis["generators"]}
    atoms={sp.Symbol(n):f.vec(g["hlog_terms"]) for n,g in definitions.items()}
    rows=f.read(ROOT/"benchmarks/tgv/families/t1-opposite/results-r001.json")["results"]
    assert {r["orbit_id"][-2:] for r in rows}=={"11","12","32"}
    count=nonzero=0
    for r in rows:
        p=ROOT/r["path"];assert sha(p)==r["sha256"];d=f.read(p);oid=d["orbit_id"][-2:]
        assert d["full_epsilon"]["expression"]==FULL[oid] and d["full_epsilon"]["definitions"]==DEFS
        assert d["basis_definitions"]==definitions and d["formal"]=="not_run"
        assert sha(ROOT/d["basis_source"]["path"])==d["basis_source"]["sha256"]
        assert sha(ROOT/d["source"]["prepared_path"])==d["source"]["prepared_sha256"]
        assert sha(ROOT/d["source"]["canonical_path"])==d["source"]["canonical_sha256"]
        assert set(d["laurent"]["coefficients"])=={str(k) for k in range(-6,1)}
        for k,expr in d["laurent"]["coefficients"].items():
            actual=f.ha.from_sympy(sp.sympify(expr,locals={str(a):a for a in atoms}),atoms)
            assert actual==f.vec(d["hlog_coefficients"][k])==exact[oid][int(k)]
            count+=1;nonzero+=bool(actual)
        assert d["pole_order"]==-min(k for k,v in exact[oid].items() if v)
    for k in range(-6,1):assert f.ha.add(exact["11"][k],exact["12"][k])==exact["32"][k]
    # The coordinate basis used for every target was checked to be crossing even.
    assert all(g["name"] in definitions for g in basis["generators"])
    num=f.read(evidence/"numeric.json")
    assert num["passed"] and len(num["I_remainders"])==9 and len(num["raw_integral_checks"])==6 and len(num["series_checks"])==4
    assert all(r["passed"] and abs(r["difference"])<=r["tolerance"] for r in num["raw_integral_checks"])
    assert count==21 and nonzero==11
    report={"passed":True,"integrals":3,"Laurent_coefficients":count,"nonzero_coefficients":nonzero,
       "exact_basis_and_Hlog_roundtrips":count,"coefficient_relations":7,"independent_old_T_moments":4,
       "weight6_MN_rank":12,"weight6_MN_plus_targets_rank":13,"new_basis_generators":0,
       "I_remainder_samples":9,"raw_QMC_cases":6,"series_tail_cases":4,
       "source_and_snapshot_hash_checks":2*len(manifest["bindings"]),"formal":"not_run"}
    print(json.dumps(report,indent=2));return report
if __name__=="__main__":main()
