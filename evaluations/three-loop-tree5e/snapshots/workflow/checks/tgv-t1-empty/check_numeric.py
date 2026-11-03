"""High precision full-epsilon/remainder checks plus direct original seven-dimensional QMC."""
import json,sys,time
from pathlib import Path
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.stats import qmc
from formulas import ROOT,HERE,IDS,run,coefficients,numerical_full,moment,e,ell
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-mixed"))
from numeric import make_function,evaluate
def main():
    records=[];start=time.time()
    # Finite analytic moments are checked by full functions at small nonzero e.
    for ep in ("0.02","0.01","0.005"):
        full=numerical_full(ep,"1.3",65);eps=mp.mpf(ep)
        fc=sum(mp.mpf(str(sp.N(moment("F",k),70)))*eps**k for k in range(6))
        zc=sum(mp.mpf(str(sp.N(moment("Z",k),70)))*eps**k for k in range(1,5))
        records.append({"kind":"moment_remainder","eps":ep,"F_residual":str(full["F"]-fc),"F_residual_over_eps6":str((full["F"]-fc)/eps**6),
            "Z_residual":str(full["Z"]-zc),"Z_residual_over_eps5":str((full["Z"]-zc)/eps**5)})
        print("MOMENTS",ep,"F/eps6",records[-1]["F_residual_over_eps6"],"Z/eps5",records[-1]["Z_residual_over_eps5"],flush=True)
    for key in ("F","Z"):
        residuals=[mp.mpf(r[key+"_residual"]) for r in records]
        power=6 if key=="F" else 5
        ratios=[abs(residuals[i]/residuals[i+1]) for i in range(2)]
        assert all(mp.mpf(2)**(power-1)<x<mp.mpf(2)**(power+1) for x in ratios),(key,ratios)
    cs=coefficients();numeric=[]
    for ep in ("0.02","0.01","0.005"):
        f=numerical_full(ep,"1.3",75);eps=mp.mpf(ep)
        for oid in IDS:
            ev=sp.lambdify(ell,sum(v*e**k for k,v in cs[oid].items()).subs(e,sp.Rational(ep)),"mpmath")
            approx=ev(mp.log(mp.mpf("1.3")))
            residual=f["I"][oid]-approx
            numeric.append({"orbit_id":"5e-3g-"+oid,"eps":ep,"z":"1.3","residual":str(residual),"residual_over_eps":str(residual/eps)})
        print("I_REMAINDERS",ep,flush=True)
    for oid in IDS:
        vals=[mp.mpf(r["residual"]) for r in numeric if r["orbit_id"].endswith("-"+oid)]
        assert all(mp.mpf("1.4")<abs(vals[i]/vals[i+1])<mp.mpf("2.8") for i in range(2)),(oid,vals)
    # Full original J_raw, alpha simplex as Dirichlet measure; no IBP moments used here.
    raw=json.loads((ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json").read_text(encoding="utf-8"))["orbits"]
    raw={r["orbit_id"][-2:]:r for r in raw if r["orbit_id"][-2:] in IDS}
    direct=[];ep,zval=0.65,0.4;full=numerical_full("0.65","0.4",55)
    for oid in IDS:
        fn=make_function(raw[oid]["J_raw"],ep,zval)
        vals=[]
        for seed in range(8):
            pts=qmc.Sobol(len(fn[0])+2,scramble=True,seed=9200+seed).random_base2(17)
            vals.append(evaluate(fn,pts,ep))
        mean=float(np.mean(vals));se=float(np.std(vals,ddof=1)/np.sqrt(len(vals)));exact=float(full["J"][oid])
        diff=mean-exact;tolerance=8*se+1e-4*max(abs(exact),1e-10)
        passed=abs(diff)<=tolerance
        direct.append({"orbit_id":"5e-3g-"+oid,"eps":ep,"z":zval,"raw_J_mean":mean,
            "closed_J":exact,"difference":diff,"empirical_standard_error":se,"tolerance":tolerance,
            "relative_difference":diff/exact,"samples":vals,"points_per_scramble":2**17,"scrambles":8,"passed":passed})
        print("RAW_QMC",oid,"difference",diff,"SE",se,"passed",passed,flush=True)
    assert all(r["passed"] for r in direct)
    report={"passed":True,"moment_checks":records,"I_remainders":numeric,"raw_integral_checks":direct,
      "error_kind":"QMC errors empirical, not rigorous; small-e checks use 65/75 decimal precision",
      "seconds":time.time()-start,"formal":"not_run"}
    (run()/"numeric.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("NUMERIC_COMPLETE",time.time()-start,flush=True)
if __name__=="__main__":main()
