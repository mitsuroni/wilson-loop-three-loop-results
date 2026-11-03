"""Fresh raw T2[1] QMC versus the frozen full-function T1[1] evaluations."""
import sys,time,json
import numpy as np
from scipy.stats import qmc
from check_math import ROOT,HERE,PAIR,run,read,dump,sha,bridge
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-mixed"))
from numeric import make_function
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-t2-empty"))
from stable_qmc import evaluate
def main():
    started=time.time();pin=read(run()/"source-delivery.json")
    delivery=ROOT/pin["delivery"];manifest=read(delivery/"manifest.json")
    assert sha(delivery/"manifest.json")==pin["manifest_sha256"]
    path=ROOT/pin["evidence_run"]/"numeric.json"
    binding=next(b for b in manifest["bindings"] if b["path"]==path.relative_to(ROOT).as_posix())
    assert sha(path)==binding["sha256"]
    prior=read(path);assert prior["passed"]
    records=prior["raw_integral_checks"]
    raw=bridge.records(bridge.SOURCE);rows=[]
    for case,(ep,z) in enumerate((("0.65","0.4"),("0.7","1.3"))):
        for target,source in PAIR.items():
            reference=next(r for r in records if r["orbit_id"]==source and r["eps"]==ep and r["z"]==z)
            fn=make_function(raw[target]["J_raw"],float(ep),float(z));samples=[]
            for seed in range(8):
                points=qmc.Sobol(7,scramble=True,seed=31100+100*case+13*int(target)+seed,bits=30).random_base2(17)+2.0**-31
                samples.append(evaluate(fn,points,float(ep)))
            exact=float(reference["analytic_J"]);mean=float(np.mean(samples));se=float(np.std(samples,ddof=1)/np.sqrt(8))
            diff=mean-exact;tol=8*se+1e-4*max(abs(exact),1e-10)
            rows.append({"orbit_id":"5e-3g-"+target,"source_T1":"5e-3g-"+source,"eps":ep,"z":z,
              "raw_J_mean":mean,"analytic_J":reference["analytic_J"],"difference":diff,"empirical_standard_error":se,
              "relative_difference":diff/exact,"tolerance":tol,"passed":abs(diff)<=tol,"samples":samples,
              "scrambles":8,"points_per_scramble":2**17,
              "sampler":"30-bit Sobol plus uniform 2^-31 half-cell shift, complementary beta quantiles; no clipping"})
            print(target,ep,z,"relative",diff/exact,"SE",se,"passed",rows[-1]["passed"],flush=True)
    report={"passed":all(r["passed"] for r in rows),"records":rows,
      "analytic_reference":{"path":path.relative_to(ROOT).as_posix(),"sha256":sha(path),"kind":"previously evaluated 55-digit full T1 functions, not old QMC means"},
      "max_absolute_relative_difference":max(abs(r["relative_difference"]) for r in rows),
      "error_kind":"empirical QMC error, not a rigorous enclosure","seconds":time.time()-started,"formal":"not_run"}
    dump(run()/"numeric.json",report);assert report["passed"]
    print("NUMERIC_COMPLETE",report["seconds"],flush=True)
if __name__=="__main__":main()
