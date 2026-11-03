"""Fresh seven-dimensional raw T2 evaluations against the mapped analytic T1 functions."""
import json,sys,time
import numpy as np
from scipy.stats import qmc
from check_math import ROOT,HERE,PAIR,records,SOURCE,run,dump,load_module
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-mixed"))
from numeric import make_function
from stable_qmc import evaluate
formulas=load_module("t2_t1formulas",ROOT/"workflow/checks/tgv-t1-empty/formulas.py")
def main():
    started=time.time();raw=records(SOURCE);reports=[]
    for index,(ep,zval) in enumerate((("0.65","0.4"),("0.7","2.3"))):
        full=formulas.numerical_full(ep,zval,55)
        for target,source in PAIR.items():
            fn=make_function(raw[target]["J_raw"],float(ep),float(zval));values=[]
            for seed in range(8):
                points=qmc.Sobol(len(fn[0])+2,scramble=True,seed=15000+100*index+11*int(target)+seed,bits=30).random_base2(17)+2.0**-31
                values.append(evaluate(fn,points,float(ep)))
            mean=float(np.mean(values));se=float(np.std(values,ddof=1)/np.sqrt(len(values)))
            exact=float(full["J"][source]);diff=mean-exact
            tol=8*se+1e-4*max(abs(exact),1e-10)
            rec={"orbit_id":"5e-3g-"+target,"reference_T1":"5e-3g-"+source,"eps":ep,"z":zval,
                 "raw_mean":mean,"analytic_J":str(full["J"][source]),"difference":diff,
                 "empirical_standard_error":se,"relative_difference":diff/exact,"tolerance":tol,
                 "sampling":"30-bit Sobol plus a uniform half-cell shift 2-31; open cube, no deleted samples","points_per_scramble":2**17,"scrambles":8,"samples":values,"passed":abs(diff)<=tol}
            reports.append(rec);print(target,ep,zval,"relative",diff/exact,"SE",se,"passed",rec["passed"],flush=True)
    report={"passed":all(x["passed"] for x in reports),"records":reports,
      "error_kind":"empirical Sobol scramble spread; not a rigorous bound",
      "max_absolute_relative_difference":max(abs(x["relative_difference"]) for x in reports),
      "seconds":time.time()-started,"formal":"not_run"}
    dump(run()/"numeric.json",report);assert report["passed"]
    print("NUMERIC_COMPLETE",report["seconds"],flush=True)
if __name__=="__main__":main()
