"""Independent numerical evaluation of reduced full-epsilon functions, z>0."""
import mpmath as mp
def quad(f,z):
 w=z/(1+z)
 return mp.quad(f,[0,w/2,w,(1+w)/2,1])
def Ls(e,u):
 B=mp.beta(1-e,1-e)
 A=u**e*mp.betainc(1-e,-2*e,0,1-u)
 T=(u*(1-u))**(-e)
 aa=mp.gamma(1-e)*mp.gamma(1+e)
 Q=aa*(mp.hyp2f1(1+2*e,-e,1,u)+2*e*u*mp.hyp2f1(1+2*e,1-e,2,u))
 return B*mp.expm1(e*mp.log(u))+e*A,T-e*A-Q
def Z(e,z,u):
 w=z/(1+z);a=z*(1-u)
 if abs(u-w)<mp.mpf(10)**(-mp.mp.dps/2):
  return (1+z)*u**(2*e-1)
 return (u**(2*e)-a**(2*e))/(2*e*(u-w))
def Zweight(e,z,u):
 a=z*(1-u);v=((1+z)*u-z)/a;lam=1-2*e
 return (1+z)*a**(-lam)*(mp.hyp2f1(lam,1,2,-v)-mp.hyp2f1(lam,1+e,2+e,-v)/(1+e))
def H(e,z,oid):
 B=mp.beta(1-e,1-e)
 if oid in ("09","10"):
  ix=0 if oid=="09" else 1
  return quad(lambda u:Z(e,z,u)*Ls(e,u)[ix],z)
 return quad(lambda u:Zweight(e,z,u)*((u*(1-u))**(-e)-B),z)
def values(ep,zval,dps=55):
 mp.mp.dps=dps;e=mp.mpf(ep);z=mp.mpf(zval)
 hs={oid:H(e,z,oid) for oid in ("09","10","26")}
 his={oid:H(e,1/z,oid) for oid in hs}
 j={oid:-v/(4*(1-2*e)*e**4) for oid,v in hs.items()}
 for a,b in (("39","09"),("38","10"),("30","26")):j[a]=-z**(3*e)*his[b]/(4*(1-2*e)*e**4)
 G=mp.exp(-3*mp.euler*e)*mp.gamma(1-e)*mp.gamma(1-2*e)
 ii={oid:-G*(hs[oid]+z**(3*e)*his[oid])/(32*e**4*(1+z**(3*e))) for oid in hs}
 for a,b in (("39","09"),("38","10"),("30","26")):ii[a]=ii[b]
 return {"J":j,"I":ii,"H":hs,"H_inverse_z":his}
