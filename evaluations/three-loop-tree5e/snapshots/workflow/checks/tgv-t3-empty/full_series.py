"""Explicit, fully integrated full-epsilon blocks; positive z."""
import sympy as sp
import numpy as np
from scipy.special import roots_jacobi
from formulas import ROOT,HERE,run,dump
t,b,x,y=sp.symbols("t b x y")
def block(n,e,z,oid):
 F=z*(1-t)*(1-b*x-(1-b)*y)+t*(b*x+(1-b)*y)
 poly=sp.Poly(sp.expand(((1+z)*t-z)*(2*b-1)*(1-F/(1+z))**n),t,b,x,y)
 terms=[]
 for (a,i,c,d),coeff in poly.terms():
  beta=sp.gamma(i+1-e)*sp.gamma(1-e)/sp.gamma(i+2-2*e)
  if oid=="09":m=(1/sp.Integer(c+d+2)-1/(c+d+2+e))/((a+1)*(c+1))
  elif oid=="10":m=(1/sp.Integer(c+1)-1/(c+1+e))/((a+1)*(c+d+2+e))
  elif oid=="26":m=(1/sp.Integer(a+1)-1/(a+1+e))/((c+1)*(c+d+2))
  else:raise ValueError(oid)
  terms.append(coeff*beta*m)
 return sum(terms)
def jacobi_grid(left,right,n=10):
 nodes,weights=roots_jacobi(n,right,left)
 return (nodes+1)/2,weights/2**(left+right+1)
def direct_polynomial(n,e,z,oid):
 # x=r*y; Jacobi weights integrate powers before evaluating the remaining polynomial.
 if oid=="09":regions=[(1,0,0,1),(-1,0,0,1+e)]
 elif oid=="10":regions=[(1,0,0,1+e),(-1,0,e,1+e)]
 else:regions=[(1,0,0,1),(-1,e,0,1)]
 ans=0.
 for sign,at,ar,ay in regions:
  specs=[jacobi_grid(at,0),jacobi_grid(-e,-e),jacobi_grid(ar,0),jacobi_grid(ay,0)]
  T,B,R,Y=np.meshgrid(*(s[0] for s in specs),indexing="ij")
  weights=np.prod(np.array(np.meshgrid(*(s[1] for s in specs),indexing="ij")),axis=0)
  U=B*R*Y+(1-B)*Y;F=z*(1-T)*(1-U)+T*U
  ans+=sign*np.sum(weights*((1+z)*T-z)*(2*B-1)*(1-F/(1+z))**n)
 return float(ans)
def main():
 reports=[];e=sp.Rational(13,20);z=sp.Rational(2,5)
 for n in range(5):
  for oid in ("09","10","26"):
   exact=float(sp.N(block(n,e,z,oid),40));direct=direct_polynomial(n,float(e),float(z),oid);err=direct-exact
   assert abs(err)<2e-11*max(1,abs(exact)),(n,oid,exact,direct)
   reports.append({"n":n,"orbit_id":oid,"eps":str(e),"z":str(z),"integrated_block":exact,"independent_Jacobi_quadrature":direct,"difference":err,"passed":True})
 report={"passed":True,"blocks":reports,"test":"tensor Gauss-Jacobi on x=r*y versus exact finite monomial sum","finite_block_checks":len(reports),"full_series_convergence":"proved by dominated binomial expansion in method.md; not a rate or interval bound","formal":"not_run"}
 dump(run()/"series-checks.json",report);print("SERIES_BLOCKS_PASSED",len(reports))
if __name__=="__main__":main()
