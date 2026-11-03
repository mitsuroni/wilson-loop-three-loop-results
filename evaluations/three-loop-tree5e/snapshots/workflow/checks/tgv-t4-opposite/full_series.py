"""Fully integrated two-kernel full-epsilon series blocks and independent cubature."""
from functools import lru_cache
import sympy as sp,numpy as np
from scipy.special import roots_jacobi
from check_math import ROOT,HERE,dump,S,P,x,y,w,z,a1,a2,a3
p,q=sp.symbols('p q');Dp=(1-p)*(1-q)+z*p*q

def block(k,l,e,zv,oid):
 pol=sp.Poly(sp.expand((P*(1-S/(1+z))**k*(1-Dp/(1+z))**l).subs(z,zv)),x,y,w,p,q,a1,a2,a3)
 @lru_cache(None)
 def am(i,j,n):return sp.N(sp.gamma(i+1-e)*sp.gamma(j+1-e)*sp.gamma(n+1-e)/sp.gamma(i+j+n+3-3*e),45)
 ans=0
 for (a,b,c,h,j,i1,i2,i3),cf in pol.terms():
  mm=sp.Rational(1,(a+1)*(a+h+2)*(b+1)*(c+j+2)*(c+1 if oid==27 else j+1))
  ans+=cf*mm*am(i1,i2,i3)
 return float(ans)

def grid(left,right,n=5):
 xx,ww=roots_jacobi(n,right,left);return (xx+1)/2,ww/2**(left+right+1)

def direct(k,l,e,zv,oid):
 specs=[grid(0,0),grid(0,0),grid(0,0),grid(1,0),grid(1,0),grid(-e,1-2*e),grid(-e,-e)]
 A,Y,B,Pv,Qv,al,bt=np.meshgrid(*(s[0] for s in specs),indexing='ij')
 weights=np.prod(np.array(np.meshgrid(*(s[1] for s in specs),indexing='ij')),axis=0)
 X=A*Pv
 if oid==27:W=B*Qv;pureq=Qv
 else:W=Qv;pureq=B*Qv
 f=sp.lambdify((x,y,w,p,q,a1,a2,a3,z),P*(1-S/(1+z))**k*(1-Dp/(1+z))**l,'numpy',cse=True)
 return float(np.sum(weights*f(X,Y,W,Pv,pureq,al,(1-al)*bt,(1-al)*(1-bt),zv)))

def main():
 rows=[];e=sp.Rational(13,20);zv=sp.Rational(2,5)
 for total in range(3):
  for k in range(total+1):
   l=total-k
   for oid in (27,28):
    a=block(k,l,e,zv,oid);b=direct(k,l,float(e),float(zv),oid);err=b-a;assert abs(err)<3e-10*max(1,abs(a)),(oid,k,l,a,b)
    rows.append({'id':oid,'k':k,'l':l,'eps':str(e),'z_D':str(zv),'monomial_sum':a,'independent_7D_Gauss_Jacobi':b,'difference':err,'passed':True});print('BLOCK',oid,k,l,err,flush=True)
 dump(HERE/'series-checks.json',{'passed':True,'finite_blocks':rows,'convergence_scope':'initial absolute convergence strip; finite signed blocks first; no numerical infinite-tail bound','formal':'not_run'})
if __name__=='__main__':main()