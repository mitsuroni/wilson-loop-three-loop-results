"""Full functions versus coefficients, independent raw J, and integrated-series tails."""
import sys,time,json
import numpy as np
import mpmath as mp
from scipy.stats import qmc
from formulas import ROOT,HERE,run,read,dump,coefficients
from full_functions import all_values,integrated_series
from numeric_hlog import HlogEvaluator
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-mixed"))
from numeric import make_function
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-t2-empty"))
from stable_qmc import evaluate
def main():
    start=time.time();cs=coefficients();mp.mp.dps=65;z="0.4";ev=HlogEvaluator(z,256)
    coeff={oid:{k:ev.vector(v) for k,v in rows.items()} for oid,rows in cs.items()}
    remainders=[]
    for ep in ("0.03","0.015","0.0075"):
        full=all_values(ep,z,65);e=mp.mpf(ep)
        for oid in cs:
            approx=mp.fsum(v*e**k for k,v in coeff[oid].items())
            residual=full["I"][oid]-approx
            remainders.append({"orbit_id":oid,"eps":ep,"z":z,"full_I":str(full["I"][oid]),"residual":str(residual),"residual_over_eps":str(residual/e)})
        print("REMAINDERS",ep,[r["residual_over_eps"] for r in remainders[-3:]],flush=True)
    for oid in cs:
        vals=[mp.mpf(r["residual"]) for r in remainders if r["orbit_id"]==oid]
        assert all(mp.mpf("1.4")<abs(vals[i]/vals[i+1])<mp.mpf("2.8") for i in range(2)),(oid,vals)
    source={r["orbit_id"][-2:]:r for r in read(ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json")["orbits"]}
    direct=[];series=[]
    for case,(ep,zval) in enumerate((("0.65","0.4"),("0.7","1.3"))):
        full=all_values(ep,zval,55)
        for oid in cs:
            fn=make_function(source[oid]["J_raw"],float(ep),float(zval));samples=[]
            for seed in range(8):
                points=qmc.Sobol(7,scramble=True,seed=24100+100*case+13*int(oid)+seed,bits=30).random_base2(17)+2.0**-31
                samples.append(evaluate(fn,points,float(ep)))
            mean=float(np.mean(samples));se=float(np.std(samples,ddof=1)/np.sqrt(8));exact=float(full["J"][oid]);delta=mean-exact
            tolerance=8*se+1e-4*max(abs(exact),1e-10)
            direct.append({"orbit_id":oid,"eps":ep,"z":zval,"raw_J_mean":mean,"analytic_J":str(full["J"][oid]),
              "difference":delta,"empirical_standard_error":se,"relative_difference":delta/exact,"tolerance":tolerance,
              "passed":abs(delta)<=tolerance,"scrambles":8,"points_per_scramble":2**17,"samples":samples,
              "sampler":"30-bit Sobol plus uniform 2^-31 half-cell shift; stable complementary beta quantiles; no clipping"})
            print("RAW",oid,ep,zval,delta/exact,"SE",se,"passed",direct[-1]["passed"],flush=True)
        mp.mp.dps=65;partial=integrated_series(ep,zval,32)
        for name in ("T","K"):
            gap=full[name]-partial[name];bound=partial[name+"_tail_bound"]
            assert abs(gap)<=bound and gap!=0
            gaps=[abs(full[name]-mp.mpf(partial["partials"][n][name])) for n in (8,16,32)]
            assert gaps[0]>gaps[1]>gaps[2]
            series.append({"moment":name,"eps":ep,"z":zval,"N":32,"full":str(full[name]),"partial":str(partial[name]),
              "remainder":str(gap),"analytic_tail_bound":str(bound),"partial_errors":list(map(str,gaps)),"passed":True})
        print("SERIES",ep,zval,"passed",flush=True)
    assert all(r["passed"] for r in direct)
    report={"passed":True,"I_remainders":remainders,"raw_integral_checks":direct,"series_checks":series,
      "max_raw_relative_difference":max(abs(r["relative_difference"]) for r in direct),
      "seconds":time.time()-start,"error_kind":"raw QMC has empirical error; high-precision reduced integrals are not interval arithmetic","formal":"not_run"}
    dump(run()/"numeric.json",report);print("NUMERIC_COMPLETE",report["seconds"],flush=True)
if __name__=="__main__":main()
