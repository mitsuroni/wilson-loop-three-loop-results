"""Offline exact bridge, Laurent/Hlog, and current/snapshot provenance verification."""
import math
import sympy as sp
from check_math import ROOT,HERE,PAIR,read,sha,compute,dump
def bindings(manifest,delivery):
    count=0
    for b in manifest["bindings"]:
        assert sha(ROOT/b["path"])==b["sha256"],b["path"]
        assert sha(delivery/b["snapshot"])==b["sha256"],b["snapshot"]
        count+=2
    return count
def main():
    delivery=ROOT/read(HERE/"current-delivery.json")["delivery"]
    manifest=read(delivery/"manifest.json");count=bindings(manifest,delivery)
    evidence=ROOT/manifest["evidence_run"]
    bridge=compute();assert bridge==read(evidence/"reflection.json")
    reused_counts={}
    for name,r in read(evidence/"reuse.json").items():
        d=ROOT/r["delivery"];assert sha(d/"manifest.json")==r["manifest_sha256"]
        reused_counts[name]=bindings(read(d/"manifest.json"),d)
    index=read(ROOT/"benchmarks/tgv/families/t2-empty/results-r001.json")
    assert {r["orbit_id"][-2:] for r in index["results"]}==set(PAIR)
    ell=sp.Symbol("ell");coeff={};nonzero=0
    for row in index["results"]:
        p=ROOT/row["path"];assert sha(p)==row["sha256"];d=read(p);t=d["orbit_id"][-2:];s=PAIR[t]
        old=read(ROOT/d["reused_result"]["path"])
        assert old["orbit_id"]=="5e-3g-"+s
        assert sha(ROOT/d["reused_result"]["path"])==d["reused_result"]["sha256"]
        for key in ("full_epsilon","laurent","hlog_coefficients","I_definition","source","normalization","domain","pole_order"):
            assert d[key]==old[key],(t,key)
        assert d["formal"]==old["formal"]=="not_run"
        assert d["reflection"]==next(r for r in bridge["pairs"] if r["T2"]==d["orbit_id"])
        coeff[t]={}
        assert set(d["laurent"]["coefficients"])=={str(i) for i in range(-6,1)}
        for k,v in d["laurent"]["coefficients"].items():
            val=sp.sympify(v,locals={"ell":ell});hv=sp.Integer(0)
            for term in d["hlog_coefficients"][k]:
                n=len(term["Hlog_word"]);assert term["Hlog_word"]==[0]*n
                hv+=sp.Rational(term["coefficient"])*sp.pi**term["pi_power"]*sp.zeta(3)**term["zeta3_power"]*sp.zeta(5)**term["zeta5_power"]*ell**n/math.factorial(n)
            assert sp.expand(hv-val)==0
            assert sp.expand(val.subs(ell,-ell)-val)==0
            coeff[t][int(k)]=val;nonzero+=int(val!=0)
        assert -min(k for k,v in coeff[t].items() if v!=0)==d["pole_order"]
    for k in range(-6,1):
        assert sp.expand(coeff["07"][k]-2*coeff["06"][k])==0
        assert sp.expand(coeff["08"][k]-coeff["05"][k]-coeff["06"][k])==0
        assert sp.expand(coeff["31"][k]-coeff["20"][k]-coeff["21"][k])==0
    num=read(evidence/"numeric.json")
    assert num["passed"] and len(num["records"])==14
    assert all(r["passed"] and abs(r["difference"])<=r["tolerance"] for r in num["records"])
    report={"passed":True,"integrals":7,"exact_reflection_pairs":7,"negative_cases":4,
      "Laurent_coefficients":49,"nonzero_coefficients":nonzero,"Hlog_roundtrips":49,
      "crossing_coefficients":49,"linear_coefficient_relations":21,"fresh_raw_QMC_cases":14,
      "dependency_hash_checks":count,"reused_dependency_hash_checks":reused_counts,
      "fresh_HyperInt_run":False,"formal":"not_run"}
    print(__import__("json").dumps(report,indent=2))
    return report
if __name__=="__main__":main()
