"""Independent full-function and original seven-dimensional numerical checks."""
import sys,time,json
import numpy as np
import mpmath as mp
from scipy.stats import qmc
from scipy.special import gamma
from formulas import ROOT,HERE,run,read,dump,coefficients
from full_functions import values,full44,middleH
from numeric_hlog import HlogEvaluator
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-mixed"))
from numeric import make_function
sys.path.insert(0,str(ROOT/"workflow/checks/tgv-t2-empty"))
from stable_qmc import evaluate
def main():
 start=time.time();mp.mp.dps=65;ev=HlogEvaluator(".4",256)
 cv={i:{k:ev.vector(v) for k,v in c.items()} for i,c in coefficients().items()}
 remainders=[]
 for ep in ("0.03","0.015","0.0075"):
  e=float(ep);z=.4
  lo=values(e,z,256);hi=values(e,z,512)
  g=np.exp(-3*float(mp.euler)*e)*gamma(1-e)*gamma(1-2*e)
  factor=g/(32*e**2*(1+z**(3*e)))
  lo["I"]["23"]=factor*middleH(e,z,1024);hi["I"]["23"]=factor*middleH(e,z,2048)
  for oid in ("22","23"):
   approx=float(mp.fsum(v*mp.mpf(ep)**k for k,v in cv[oid].items()))
   residual=float(hi["I"][oid]-approx);refinement=float(hi["I"][oid]-lo["I"][oid])
   assert abs(refinement)<.04*abs(residual),(oid,ep,residual,refinement)
   remainders.append({"orbit_id":oid,"eps":ep,"z":"0.4","full_I":float(hi["I"][oid]),"residual":residual,"residual_over_eps":residual/e,
    "refinement_difference":refinement,"quadrature_orders":[256,512] if oid=="22" else [1024,2048],"precision":"float64 quadrature; coefficients by 65-digit Hlog ODE series"})
  full=full44(ep,".4",65);em=mp.mpf(ep);res=full["I"]-mp.fsum(v*em**k for k,v in cv["44"].items())
  remainders.append({"orbit_id":"44","eps":ep,"z":"0.4","full_I":str(full["I"]),"residual":str(res),"residual_over_eps":str(res/em),"precision":"65-digit reduced quadrature"})
  print("REMAINDER",ep,[(r["orbit_id"],r["residual_over_eps"]) for r in remainders[-3:]],flush=True)
 for oid in ("22","23","44"):
  vals=[mp.mpf(r["residual"]) for r in remainders if r["orbit_id"]==oid]
  assert all(mp.mpf("1.4")<abs(vals[i]/vals[i+1])<mp.mpf("2.8") for i in range(2)),(oid,vals)
 sources={d["orbit_id"][-2:]:d for d in read(ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json")["orbits"]}
 raw=[];refinements=[];middle_bridges=[]
 for case,(ep,zval) in enumerate((("0.65","0.4"),("0.7","1.3"))):
  e,z=float(ep),float(zval);lo=values(e,z,128,True);hi=values(e,z,256,True);f44=full44(ep,zval,55)
  lo["J"]["44"]=hi["J"]["44"]=float(f44["J"])
  delta_mid=hi["Hmiddle"]-4*(1-2*e)*e**2*(hi["J"]["23"]+hi["J"]["24"])
  assert abs(delta_mid)<5e-6*max(1,abs(hi["Hmiddle"]))
  middle_bridges.append({"eps":ep,"z":zval,"Hmiddle_minus_raw_sum":delta_mid,"passed":True})
  for oid in ("22","23","24","25","44"):
   exact=hi["J"][oid];refine=abs(exact-lo["J"][oid]);assert refine<.001*abs(exact),(oid,refine,exact)
   refinements.append({"orbit_id":oid,"eps":ep,"z":zval,"low":lo["J"][oid],"high":exact,"absolute_difference":refine})
   fn=make_function(sources[oid]["J_raw"],e,z);samples=[]
   for seed in range(8):
    pts=qmc.Sobol(7,scramble=True,seed=52100+100*case+13*int(oid)+seed,bits=30).random_base2(17)+2.0**-31
    samples.append(evaluate(fn,pts,e))
   mean=float(np.mean(samples));se=float(np.std(samples,ddof=1)/np.sqrt(8));delta=mean-exact
   tol=8*se+1e-4*abs(exact)+4*refine
   rec={"orbit_id":oid,"eps":ep,"z":zval,"raw_J_mean":mean,"reduced_full_J":exact,
    "difference":delta,"empirical_standard_error":se,"relative_difference":delta/exact,"tolerance":tol,"passed":abs(delta)<=tol,
    "quadrature_refinement":refine,"scrambles":8,"points_per_scramble":2**17,"samples":samples}
   raw.append(rec);print("RAW",oid,ep,zval,delta/exact,"passed",rec["passed"],flush=True)
   dump(run()/"numeric-progress.json",{"I_remainders":remainders,"raw_integral_checks":raw})
 report={"passed":all(r["passed"] for r in raw),"I_remainders":remainders,"raw_integral_checks":raw,"full_quadrature_refinements":refinements,
  "middle_partition_checks":middle_bridges,"max_raw_relative_difference":max(abs(r["relative_difference"]) for r in raw),
  "seconds":time.time()-start,"error_kind":"empirical QMC and quadrature refinement differences; not rigorous enclosures",
  "sampler":"30-bit Sobol plus 2^-31 uniform half-cell shift; stable complementary beta quantiles; no deleted samples","formal":"not_run"}
 dump(run()/"numeric.json",report);assert report["passed"];print("NUMERIC_COMPLETE",report["seconds"],flush=True)
if __name__=="__main__":main()
