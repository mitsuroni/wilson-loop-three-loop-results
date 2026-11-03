"""Offline result and frozen dependency verification; no CAS process or external source needed."""
from pathlib import Path
import json,hashlib,math
import sympy as sp
import formulas as f
from publish import DEFINITIONS,FULL
ROOT,HERE=f.ROOT,f.HERE
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    delivery=ROOT/json.loads((HERE/"current-delivery.json").read_text())["delivery"]
    manifest=json.loads((delivery/"manifest.json").read_text())
    for b in manifest["bindings"]:
        assert sha(ROOT/b["path"])==b["sha256"],b["path"]
        assert sha(delivery/b["snapshot"])==b["sha256"],b["snapshot"]
    evidence=ROOT/manifest["evidence_run"];f.run=lambda:evidence
    exact=f.coefficients()
    index=json.loads((ROOT/"benchmarks/tgv/families/t1-empty/results-r001.json").read_text())
    count=0;nonzero=0
    for rec in index["results"]:
        p=ROOT/rec["path"];assert sha(p)==rec["sha256"]
        d=json.loads(p.read_text(encoding="utf-8"));oid=d["orbit_id"][-2:]
        assert d["full_epsilon"]["expression"]==FULL[oid] and d["full_epsilon"]["definitions"]==DEFINITIONS
        assert sha(ROOT/d["source"]["prepared_path"])==d["source"]["prepared_sha256"]
        assert d["formal"]=="not_run"
        for power,text in d["laurent"]["coefficients"].items():
            val=sp.sympify(text,locals={"ell":f.ell})
            assert sp.expand(val-exact[oid][int(power)])==0
            h=0
            for term in d["hlog_coefficients"][power]:
                n=len(term["Hlog_word"]);assert term["Hlog_word"]==[0]*n
                h+=sp.Rational(term["coefficient"])*sp.pi**term["pi_power"]*sp.zeta(3)**term["zeta3_power"]*sp.zeta(5)**term["zeta5_power"]*f.ell**n/math.factorial(n)
            assert sp.expand(h-val)==0
            count+=1;nonzero+=int(val!=0)
    assert count==49
    print(json.dumps({"passed":True,"integrals":7,"Laurent_coefficients":count,"nonzero_coefficients":nonzero,
        "Hlog_roundtrips":count,"dependency_hash_checks":2*len(manifest["bindings"]),
        "external_repository_required":False,"formal":"not_run"},indent=2))
if __name__=="__main__":main()
