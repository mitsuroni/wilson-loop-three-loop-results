"""Independent ODE power/log series and Mobius continuation for real Hlogs."""
from functools import lru_cache
from fractions import Fraction as Q
from math import factorial
import mpmath as mp
import wide_hlog as wh
from check_math import ROOT,HERE,read,dump

def mpq(q):return mp.mpf(q.numerator)/q.denominator

class SeriesEvaluator:
 def __init__(self,N):self.N=N
 @lru_cache(None)
 def series(self,w):
  N=self.N
  if not w:return [{0:mp.mpf(1)}]+[{} for _ in range(N)]
  tail=self.series(w[1:]);ans=[{} for _ in range(N+1)]
  if w[0]==0:
   for k,c in tail[0].items():ans[0][k+1]=c/(k+1)
   numer=tail
  else:
   a=mpq(w[0]);numer=[{} for _ in range(N+1)];running={}
   for n in range(1,N+1):
    running={k:(running.get(k,0)-tail[n-1].get(k,0))/a for k in set(running)|set(tail[n-1])};numer[n]=running.copy()
  for n in range(1,N+1):
   for k,c in numer[n].items():
    for j in range(k+1):ans[n][k-j]=ans[n].get(k-j,0)+c*(-1)**j*factorial(k)/factorial(k-j)/mp.mpf(n)**(j+1)
  return ans
 def word(self,w,z):
  z=mp.mpf(z);lz=mp.log(z);val=mp.mpf(0)
  for row in reversed(self.series(w)):val=val*z+sum(c*lz**k for k,c in row.items())
  return val
 def vector(self,vec,z,periods):
  return mp.fsum(mpq(c)*mp.fprod(v**p for v,p in zip(periods,ps))*self.word(w,z) for (w,ps),c in vec.items())

@lru_cache(None)
def mobius_word(w):
 if not w:return {():Q(1)}
 a=w[0];choices=[(Q(1),Q(-1))] if a==-1 else [(a/(1+a),Q(1)),(Q(1),Q(-1))]
 out={}
 for h,c in choices:
  for tail,d in mobius_word(w[1:]).items():out[(h,)+tail]=out.get((h,)+tail,Q(0))+c*d
 return {k:v for k,v in out.items() if v}

def mobius_vec(v):
 out={}
 for (w,ps),c in v.items():
  for ww,m in mobius_word(w).items():out[ww,ps]=out.get((ww,ps),Q(0))+c*m
 return {k:v for k,v in out.items() if v}

def main():
 mp.mp.dps=70
 def A(n):return mp.quad(lambda t:mp.log1p(t)**2*(-mp.log(t))**(n-2)/(2*mp.factorial(n-2)*t) if t else 0,[0,mp.mpf('.5'),1])
 # A4 is sum (-1)^n H_{n-1}/n^3; A6 has denominator n^5.
 a4=A(3);a6=A(5);periods=[mp.zeta(2),mp.zeta(3),mp.zeta(5),mp.log(2),a4,a6]
 r=ROOT/read(HERE/'current-omega_basis.json')['run'];omega=wh.parse(r/'Omega6.txt');prim=wh.parse(r/'Tprimitive.txt');pm=mobius_vec(prim)
 ev=SeriesEvaluator(640);small=SeriesEvaluator(384)
 elementary=[]
 for letter in (Q(-1),Q(-1,2),Q(-2)):
  zz=mp.mpf('.2');got=small.word((letter,),zz);expected=mp.log1p(-zz/mpq(letter));err=abs(got-expected);assert err<mp.mpf('1e-60');elementary.append({'letter':str(letter),'error':str(err)})
 periodchecks=[]
 for value,w in [(a4,(Q(0),Q(0),Q(-1),Q(-1))),(a6,(Q(0),)*4+(Q(-1),Q(-1)))]:
  vv=mobius_vec({(w,wh.ZERO):Q(1)});got=small.vector(vv,mp.mpf('.5'),periods);err=abs(got-value);assert err<mp.mpf('1e-60');periodchecks.append({'word':[str(a) for a in w],'quadrature':str(value),'ODE':str(got),'error':str(err)})
 checks=[]
 for ztext in ('.2','.4'):
  z=mp.mpf(ztext);u=z/(1+z)
  canonical=ev.vector(omega,z,periods)
  primitive=ev.vector(pm,u,periods)+ev.vector(pm,1-u,periods)
  err=abs(canonical-primitive);assert err<mp.mpf('1e-42'),(ztext,err)
  checks.append({'z':ztext,'canonical_Hlog':str(canonical),'Tz_plus_Tinv':str(primitive),'absolute_difference':str(err),'passed':True});print('OMEGA',ztext,mp.nstr(err,6),flush=True)
 dump(HERE/'omega-numeric.json',{'passed':True,'working_precision':70,'series_terms':640,'elementary_checks':elementary,'period_checks':periodchecks,'Omega_definition_checks':checks,'period_values':{k:str(v) for k,v in zip(wh.PERIODS,periods)},'error_kind':'independent convergent ODE series versus integral-defined periods and Mobius-transformed primitive; not interval-certified','formal':'not_run'})
if __name__=='__main__':main()