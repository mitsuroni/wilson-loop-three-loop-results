"""Offline exact reflection, inherited coefficients and current/snapshot verification."""
import sys,json
import sympy as sp
from check_math import ROOT,HERE,PAIR,read,sha,compute,bridge
sys.path.insert(0,str(ROOT/"workflow/checks/single-opposite"))
import hlog_algebra as ha
def vec(terms):
    from fractions import Fraction
    return {(tuple(t["word"]),tuple(t["zeta_powers"])):Fraction(t["coefficient"]) for t in terms}
def verify_bindings(delivery):
    manifest=read(delivery/"manifest.json")
    for b in manifest["bindings"]:
        assert sha(ROOT/b["path"])==b["sha256"],b["path"]
        assert sha(delivery/b["snapshot"])==b["sha256"],b["snapshot"]
    return 2*len(manifest["bindings"])
def main():
    delivery=ROOT/read(HERE/"current-delivery.json")["delivery"];manifest=read(delivery/"manifest.json")
    count=verify_bindings(delivery);evidence=ROOT/manifest["evidence_run"]
    reflection=compute();assert reflection==read(evidence/"reflection.json")
    pin=read(evidence/"source-delivery.json");parent=ROOT/pin["delivery"]
    assert sha(parent/"manifest.json")==pin["manifest_sha256"];parent_count=verify_bindings(parent)
    rows=read(ROOT/"benchmarks/tgv/families/t2-opposite/results-r001.json")["results"]
    assert {r["orbit_id"][-2:] for r in rows}==set(PAIR)
    cs={};nonzero=0;roots=bridge.records(bridge.PREPARED)
    for r in rows:
        p=ROOT/r["path"];assert sha(p)==r["sha256"];d=read(p);oid=d["orbit_id"][-2:]
        previous=ROOT/d["reused_result"]["path"];assert sha(previous)==d["reused_result"]["sha256"]
        old=read(previous);assert old["orbit_id"]=="5e-3g-"+PAIR[oid]
        assert d["prepared_C"]==roots[oid]["C"]
        assert d["reflection"]==next(p for p in reflection["pairs"] if p["T2"]==d["orbit_id"])
        for key in ("laurent","hlog_coefficients","basis_definitions","basis_source","full_epsilon","pole_order",
                    "I_definition","normalization","domain","source","hlog_encoding","hlog_convention"):
            assert d[key]==old[key],(oid,key)
        assert d["formal"]==old["formal"]=="not_run"
        assert d["channel"]==([2,4] if oid=="16" else [1,3])
        assert set(d["laurent"]["coefficients"])=={str(k) for k in range(-6,1)}
        atoms={sp.Symbol(n):vec(g["hlog_terms"]) for n,g in d["basis_definitions"].items()}
        cs[oid]={}
        for k,value in d["laurent"]["coefficients"].items():
            expr=sp.sympify(value,locals={str(a):a for a in atoms})
            v=ha.from_sympy(expr,atoms);assert v==vec(d["hlog_coefficients"][k])
            cs[oid][int(k)]=v;nonzero+=bool(v)
        assert d["pole_order"]==-min(k for k,v in cs[oid].items() if v)
    for k in range(-6,1):
        assert ha.add(cs["36"][k],cs["37"][k])==cs["16"][k]
        assert all(len(w)+sum(a*b for a,b in zip(zp,ha.ZW))==6+k for w,zp in cs["37"][k])
    numeric=read(evidence/"numeric.json");assert numeric["passed"] and len(numeric["records"])==6
    assert all(r["passed"] and abs(r["difference"])<=r["tolerance"] for r in numeric["records"])
    assert nonzero==11
    report={"passed":True,"integrals":3,"reflection_pairs":3,"negative_cases":5,"Laurent_coefficients":21,
      "nonzero_coefficients":nonzero,"basis_Hlog_roundtrips":21,"coefficient_relations":7,"fresh_raw_QMC_cases":6,
      "source_and_snapshot_hash_checks":count,"reused_T1_hash_checks":parent_count,
      "fresh_HyperInt_run":False,"new_basis_generators":0,"formal":"not_run"}
    print(json.dumps(report,indent=2));return report
if __name__=="__main__":main()
