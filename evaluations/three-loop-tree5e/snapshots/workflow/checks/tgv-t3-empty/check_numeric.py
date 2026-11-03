"""Fresh raw seven-dimensional evaluations and independent full-epsilon remainders."""
import sys,time,json
import numpy as np
import mpmath as mp
from scipy.stats import qmc
from formulas import ROOT,HERE,run,read,dump,coefficients
from full_functions import values,Ls,Z,H
from numeric_hlog import HlogEvaluator
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-mixed"))
from numeric import make_function
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-t2-empty"))
from stable_qmc import evaluate
def main():
 start=time.time();cs=coefficients();mp.mp.dps=65;z="0.4";ev=HlogEvaluator(z,256)
 cv={oid:{k:ev.vector(v) for k,v in c.items()} for oid,c in cs.items()}
 rem=[]
 for ep in ("0.03","0.015","0.0075"):
  full=values(ep,z,65);e=mp.mpf(ep)
  for oid in ("09","10","26"):
   residual=full["I"][oid]-mp.fsum(v*e**k for k,v in cv[oid].items())
   rem.append({"orbit_id":oid,"eps":ep,"z":z,"full_I":str(full["I"][oid]),"residual":str(residual),"residual_over_eps":str(residual/e)})
  print("REMAINDERS",ep,[r["residual_over_eps"] for r in rem[-3:]],flush=True)
  dump(run()/"numeric-progress.json",{"remainders":rem})
 for oid in ("09","10","26"):
  vals=[mp.mpf(r["residual"]) for r in rem if r["orbit_id"]==oid]
  assert all(mp.mpf("1.4")<abs(vals[i]/vals[i+1])<mp.mpf("2.8") for i in range(2)),(oid,vals)
 source={r["orbit_id"][-2:]:r for r in read(ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json")["orbits"]}
 direct=[]
 for case,(ep,zval) in enumerate((("0.65","0.4"),("0.7","1.3"))):
  full=values(ep,zval,55)
  for oid in cs:
   fn=make_function(source[oid]["J_raw"],float(ep),float(zval));samples=[]
   for seed in range(8):
    pts=qmc.Sobol(7,scramble=True,seed=42100+100*case+13*int(oid)+seed,bits=30).random_base2(17)+2.0**-31
    samples.append(evaluate(fn,pts,float(ep)))
   mean=float(np.mean(samples));se=float(np.std(samples,ddof=1)/np.sqrt(8));exact=float(full["J"][oid]);delta=mean-exact
   tol=8*se+1e-4*max(abs(exact),1e-10)
   rec={"orbit_id":oid,"eps":ep,"z":zval,"raw_J_mean":mean,"analytic_J":str(full["J"][oid]),
    "difference":delta,"empirical_standard_error":se,"relative_difference":delta/exact,"tolerance":tol,
    "passed":abs(delta)<=tol,"scrambles":8,"points_per_scramble":2**17,"samples":samples}
   direct.append(rec);print("RAW",oid,ep,zval,delta/exact,"SE",se,"passed",rec["passed"],flush=True)
   dump(run()/"numeric-progress.json",{"remainders":rem,"raw_integral_checks":direct})
 report={"passed":all(r["passed"] for r in direct),"I_remainders":rem,"raw_integral_checks":direct,
  "max_raw_relative_difference":max(abs(r["relative_difference"]) for r in direct),
  "seconds":time.time()-start,"error_kind":"empirical QMC errors; reduced quadrature at 55/65 digits, not rigorous enclosures",
  "sampler":"30-bit Sobol + 2^-31 half-cell shift; stable beta complements; no clipping or discarded samples","formal":"not_run"}
 dump(run()/"numeric.json",report);assert report["passed"]
 print("NUMERIC_COMPLETE",report["seconds"],flush=True)
if __name__=="__main__":main()
