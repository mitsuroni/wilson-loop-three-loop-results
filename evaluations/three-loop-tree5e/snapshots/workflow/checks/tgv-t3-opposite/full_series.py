"""Fully integrated two-denominator series blocks and independent finite-block checks."""
import sympy as sp
import numpy as np
from scipy.special import roots_jacobi
from formulas import ROOT,HERE,run,dump
t,p,x,y,q,b=sp.symbols("t p x y q b")
def block(k,l,e,z,oid):
 F=lambda s,u:z*(1-s)*(1-u)+s*u
 poly=sp.Poly(sp.expand(((1+z)*t-z)*(2*b-1)*(1-F(t,b*x+(1-b)*y)/(1+z))**k*(1-F(p,q)/(1+z))**l),t,p,x,y,q,b)
 result=0
 for (a,h,c,d,j,i),coeff in poly.terms():
  beta=sp.gamma(i+1-e)*sp.gamma(1-e)/sp.gamma(i+2-2*e)
  tp=1/sp.Integer((a+1)*(a+h+2)) if oid in ("22","23") else 1/sp.Integer((h+1)*(a+h+2))
  if oid=="22":edge=1/sp.Integer((c+1)*(c+d+2)*(c+d+j+3))
  elif oid in ("23","24"):edge=1/sp.Integer((c+1)*(c+j+2)*(c+d+j+3))
  elif oid=="25":edge=1/sp.Integer((j+1)*(j+c+2)*(j+c+d+3))
  else:raise ValueError(oid)
  result+=coeff*beta*tp*edge
 return result
def grid(a,d,n=6):
 nodes,w=roots_jacobi(n,d,a);return (nodes+1)/2,w/2**(a+d+1)
def direct(k,l,e,z,oid):
 specs=[grid(1,0),grid(0,0),grid(2,0),grid(1,0),grid(0,0),grid(-e,-e)]
 P,R,O,S,V,B=np.meshgrid(*(s[0] for s in specs),indexing="ij")
 W=np.prod(np.array(np.meshgrid(*(s[1] for s in specs),indexing="ij")),axis=0)
 T=P*R
 if oid=="22":Q=O;Y=O*S;X=Y*V
 elif oid=="23":Y=O;Q=O*S;X=Q*V
 else:raise ValueError(oid)
 U=B*X+(1-B)*Y
 FT=z*(1-T)*(1-U)+T*U;FP=z*(1-P)*(1-Q)+P*Q
 return float(np.sum(W*((1+z)*T-z)*(2*B-1)*(1-FT/(1+z))**k*(1-FP/(1+z))**l))
def main():
 reports=[];e=sp.Rational(13,20);z=sp.Rational(2,5)
 for total in range(3):
  for k in range(total+1):
   l=total-k
   for oid in ("22","23"):
    exact=float(sp.N(block(k,l,e,z,oid),40));num=direct(k,l,float(e),float(z),oid);delta=num-exact
    assert abs(delta)<3e-12*max(1,abs(exact)),(k,l,oid,num,exact)
    reports.append({"orbit_id":oid,"k":k,"l":l,"eps":str(e),"z":str(z),"exact_block":exact,"direct_Gauss_Jacobi":num,"difference":delta,"passed":True})
 dump(run()/"series-checks.json",{"passed":True,"blocks":reports,"finite_block_checks":len(reports),"method":"independent six-dimensional tensor Gauss-Jacobi versus fully integrated monomial sum","formal":"not_run"})
 print("SERIES_BLOCKS_PASSED",len(reports))
if __name__=="__main__":main()
