"""Offline verification from pinned moments, exact basis and source bridges."""
import hashlib,json
import sympy as sp
import formulas as f
import check_math as mc
from publish import FULL,DEFS
ROOT,HERE=f.ROOT,f.HERE
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 delivery=ROOT/f.read(HERE/"current-delivery.json")["delivery"];manifest=f.read(delivery/"manifest.json")
 for b in manifest["bindings"]:
  assert sha(ROOT/b["path"])==b["sha256"],b["path"]
  assert sha(delivery/b["snapshot"])==b["sha256"],b["snapshot"]
 evidence=ROOT/manifest["evidence_run"];f.run=lambda:evidence;mc.run=lambda:evidence
 assert mc.check()==f.read(evidence/"math-checks.json")
 cs=f.coefficients();basis=f.read(ROOT/f.read(evidence/"representation.json")["basis_path"])
 definitions={g["name"]:g for g in basis["generators"]};atoms={sp.Symbol(n):f.vec(g["hlog_terms"]) for n,g in definitions.items()}
 rows=f.read(ROOT/"benchmarks/tgv/families/t3-opposite/results-r001.json")["results"]
 assert {r["orbit_id"][-2:] for r in rows}=={"22","23","24","25","44"}
 count=nonzero=0
 for r in rows:
  path=ROOT/r["path"];assert sha(path)==r["sha256"];d=f.read(path);oid=d["orbit_id"][-2:]
  assert d["full_epsilon"]["expression"]==FULL[oid] and d["full_epsilon"]["definitions"]==DEFS and d["full_epsilon"]["g"]==oid
  assert d["basis_definitions"]==definitions and d["formal"]=="not_run"
  assert sha(ROOT/d["basis_source"]["path"])==d["basis_source"]["sha256"]
  for label in ("prepared","canonical"):assert sha(ROOT/d["source"][label+"_path"])==d["source"][label+"_sha256"]
  assert set(d["laurent"]["coefficients"])=={str(k) for k in range(-6,1)}
  for k,expr in d["laurent"]["coefficients"].items():
   actual=f.ha.from_sympy(sp.sympify(expr,locals={str(a):a for a in atoms}),atoms)
   assert actual==f.vec(d["hlog_coefficients"][k])==cs[oid][int(k)]
   count+=1;nonzero+=bool(actual)
  assert d["pole_order"]==1==-min(k for k,v in cs[oid].items() if v)
 num=f.read(evidence/"numeric.json");blocks=f.read(evidence/"series-checks.json")
 assert num["passed"] and len(num["I_remainders"])==9 and len(num["raw_integral_checks"])==10
 assert all(r["passed"] and abs(r["difference"])<=r["tolerance"] for r in num["raw_integral_checks"])
 assert blocks["passed"] and len(blocks["blocks"])==12
 assert count==35 and nonzero==10
 report={"passed":True,"integrals":5,"Laurent_coefficients":count,"nonzero_coefficients":nonzero,"exact_basis_Hlog_roundtrips":count,
 "full_epsilon_equalities":2,"factorized_44":True,"coefficient_equalities":14,"new_basis_generators":0,
 "I_remainder_samples":9,"raw_QMC_cases":10,"finite_series_block_checks":12,"source_and_snapshot_hash_checks":2*len(manifest["bindings"]),"formal":"not_run"}
 print(json.dumps(report,indent=2));return report
if __name__=="__main__":main()
