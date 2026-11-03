"""Full-epsilon polynomial blocks for weighted T4; no residual integrals in the sum."""
import sympy as sp
import numpy as np
from scipy.special import roots_jacobi
from check_math import S,P,x,y,w,a1,a2,a3,z,ROOT,HERE,dump
from functools import lru_cache

def block(n,ep,zv,oid):
 poly=sp.Poly(sp.expand((P*(1-S/(1+z))**n).subs(z,zv)),x,y,w,a1,a2,a3)
 @lru_cache(None)
 def dm(i,j,k):return sp.N(sp.gamma(i+1-ep)*sp.gamma(j+1-ep)*sp.gamma(k+1-ep)/sp.gamma(i+j+k+3-3*ep),45)
 @lru_cache(None)
 def wm(a,b,c):
  A=sp.Rational(1,a+1);B=sp.Rational(1,b+1);C=sp.Rational(1,c+1)
  Be=sp.gamma(b+1)*sp.gamma(1+ep)/sp.gamma(b+2+ep)
  if oid==13:return A*Be*(C-1/(c+1+ep))
  if oid==14:return A*(B-Be)*(C-1/(c+1+ep))
  if oid==15:return A*(B-Be)/(c+1+ep)
  return (A-1/(a+1+ep))*B*C
 return sum(cf*dm(i,j,k)*sp.N(wm(a,b,c),45) for (a,b,c,i,j,k),cf in poly.terms())

def grid(left,right,n=9):
 xx,ww=roots_jacobi(n,right,left);return (xx+1)/2,ww/2**(left+right+1)

def direct(n,ep,zv,oid):
 cases={13:[(1,0,1,0),(-1,0,1,1)],14:[(1,0,0,0),(-1,0,1,0),(-1,0,0,1),(1,0,1,1)],15:[(1,0,0,1),(-1,0,1,1)],40:[(1,0,0,0),(-1,1,0,0)]}
 f=sp.lambdify((x,y,w,a1,a2,a3,z),P*(1-S/(1+z))**n,'numpy',cse=True)
 ans=0.
 for sign,C,A,B in cases[oid]:
  specs=[grid(C*ep,0),grid(0,A*ep),grid(B*ep,0),grid(-ep,1-2*ep),grid(-ep,-ep)]
  xx,yy,ww,al,bt=np.meshgrid(*(q[0] for q in specs),indexing='ij')
  weights=np.prod(np.array(np.meshgrid(*(q[1] for q in specs),indexing='ij')),axis=0)
  ans+=sign*np.sum(weights*f(xx,yy,ww,al,(1-al)*bt,(1-al)*(1-bt),zv))
 return float(ans)

def main():
 rows=[];ep=sp.Rational(13,20);zv=sp.Rational(2,5)
 for n in range(3):
  for oid in (13,14,15,40):
   v=float(block(n,ep,zv,oid));q=direct(n,float(ep),float(zv),oid);err=q-v
   assert abs(err)<3e-10*max(1,abs(v)),(n,oid,v,q)
   rows.append({'n':n,'representative':oid,'eps':str(ep),'z_D':str(zv),'monomial_sum':v,'Gauss_Jacobi':q,'difference':err,'passed':True})
   print('BLOCK',n,oid,err,flush=True)
 dump(HERE/'series-checks.json',{'passed':True,'blocks':rows,'formal':'not_run','scope':'finite blocks only; convergence argument is in method.md, no numerical infinite-tail bound'})
if __name__=='__main__':main()