"""Independent Gauss-Jacobi quadrature of the reduced full-epsilon kernels."""
from functools import lru_cache
import numpy as np
import mpmath as mp
from scipy.special import roots_jacobi,beta,gamma
@lru_cache(None)
def grid(n,left=0.,right=0.):
 x,w=roots_jacobi(n,right,left)
 return (x+1)/2,w/2**(left+right+1)
def F(z,t,u):return z*(1-t)*(1-u)+t*u
def partial(e,z,t,q):
 f=F(z,t,q);w=z/(1+z);d=q-w
 if abs(d)<1e-12:return (1+z)*(1-t)*q**(e-1)
 return f**e*np.expm1(e*np.log(q/f))/(e*d)
def Z(e,z,u):
 w=z/(1+z);a=z*(1-u);d=u-w
 with np.errstate(divide="ignore",invalid="ignore"):
  val=a**e*np.expm1(e*np.log(u/a))/(e*d)
 return np.where(abs(d)<1e-12,(1+z)*u**(e-1),val)
def rawH(e,z,oid,n):
 T,wt=grid(n);T=T[:,None];wt=wt[:,None]
 B=beta(1-e,1-e);value=0.
 if oid=="22":pieces=[(1,1.,0.,-e,-e,"lo",0),( -B,1.,0.,0.,0.,"lo",0)]
 else:
  pieces=[(1,1-e,0.,-e,0.,"lo",1),(-1,1.,0.,-e,-e,"lo",0),
          (1,0.,1-e,0.,-e,"hi",2),(-1,0.,1.,-e,-e,"hi",0)]
 for sign,ql,qr,rl,rr,side,extra in pieces:
  qs,wq=grid(n,ql,qr);rs,wr=grid(n,rl,rr);R=rs[None,:]
  subtotal=0.
  for q,qw in zip(qs,wq):
   U=q*R if side=="lo" else q+(1-q)*R
   factor=(1-U)**(-e) if extra==1 else U**(-e) if extra==2 else 1.
   vals=(1+z)*F(z,T,U)**(2*e-1)*partial(e,z,T,q)*factor
   subtotal+=qw*np.sum(vals*wt*wr[None,:])
  value+=sign*subtotal
 return float(value)
def middleH(e,z,n):
 value=0.
 pieces=[(1,1-e,0.,-e,0.,"lo",1),(-1,1.,0.,-e,-e,"lo",0),
         (1,0.,1-e,0.,-e,"hi",2),(-1,0.,1.,-e,-e,"hi",0)]
 for sign,ql,qr,rl,rr,side,extra in pieces:
  qs,wq=grid(n,ql,qr);rs,wr=grid(n,rl,rr);R=rs[None,:];Q=qs[:,None]
  U=Q*R if side=="lo" else Q+(1-Q)*R
  factor=(1-U)**(-e) if extra==1 else U**(-e) if extra==2 else 1.
  vals=Z(2*e,z,U)*Z(e,z,Q)*factor
  value+=sign*np.sum(vals*wq[:,None]*wr[None,:])
 return float(value)
def full44(ep,zval,dps=60):
 mp.mp.dps=dps;e=mp.mpf(ep);z=mp.mpf(zval);w=z/(1+z)
 def h0(delta):
  a=delta/2
  def f(b):
   W=mp.log(b*(1-b));L=mp.log(z)+mp.log1p(-b)-mp.log(b)
   if abs(b-w)<mp.mpf(10)**(-dps/2):return -(w*(1-w))**(a-1)
   return mp.exp(a*W)*mp.sinh(a*L)/(a*(b-w))
  return mp.quad(f,[0,w/2,w,(1+w)/2,1])
 C=mp.pi*(mp.cos(mp.pi*e)-mp.cosh(e*mp.log(z)))/(e*mp.sin(mp.pi*e))
 diff=C-mp.beta(1-e,1-e)*h0(2*e)
 j=z**(3*e/2)*h0(e)*diff/(4*(1-2*e)*e**2)
 G=mp.exp(-3*mp.euler*e)*mp.gamma(1-e)*mp.gamma(1-2*e)
 return {"J":j,"I":G*h0(e)*diff/(32*e**2*mp.cosh(3*e*mp.log(z)/2))}
def values(e,z,n=128,raw_middle=False):
 e=float(e);z=float(z);h=rawH(e,z,"22",n);hi=rawH(e,1/z,"22",n);km=middleH(e,z,n)
 g=float(mp.exp(-3*mp.euler*e))*gamma(1-e)*gamma(1-2*e)
 j={"22":h/(4*(1-2*e)*e**2),"25":z**(3*e)*hi/(4*(1-2*e)*e**2)}
 if raw_middle:
  h23=rawH(e,z,"23",n);h24=rawH(e,1/z,"23",n)
  j.update({"23":h23/(4*(1-2*e)*e**2),"24":z**(3*e)*h24/(4*(1-2*e)*e**2)})
 ii={"22":g*(h+z**(3*e)*hi)/(32*e**2*(1+z**(3*e))),"23":g*km/(32*e**2*(1+z**(3*e)))}
 ii["25"]=ii["22"];ii["24"]=ii["23"]
 return {"J":j,"I":ii,"H22":h,"H22_inverse":hi,"Hmiddle":km}
