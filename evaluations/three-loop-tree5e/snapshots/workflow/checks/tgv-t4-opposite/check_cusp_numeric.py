"""Check the newly retained w=0 cusp, including its small-epsilon expansion."""
import sys,numpy as np,mpmath as mp
from check_math import ROOT,HERE,read,dump
from check_numeric import points,alpha,polys,compare
from formulas import ha
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
# Use a separately named module because this family also has numeric_hlog.py.
from check_math import module
oldnum=module('cusp_old_hlog',ROOT/'workflow/checks/single-opposite/numeric_hlog.py')

def full(e,z):
 e,z=mp.mpf(e),mp.mpf(z)
 fa=(mp.gamma(1-e)*mp.gamma(1+e)-mp.gamma(1-e)**2/mp.gamma(1-2*e))/(2*e**2)
 def f(p):
  if p==0:return 0
  if p==1:return z**(e-1)/e
  ell=mp.log(z*p/(1-p));qp=(1-p)**(e-1)*(mp.expm1(e*ell)/(e*mp.expm1(ell)) if ell else 1)
  return qp*(-mp.expm1(2*e*mp.log1p(-p)))
 b=mp.quad(f,[0,1/(1+z),1])
 return fa*b/(2*e**2)

def main():
 mp.mp.dps=70;r=ROOT/read(HERE/'current-cusp_export.json')['run'];ev=oldnum.HlogEvaluator('.4',384)
 cm=ev.vector(ha.parse(r/'CuspMnorm.txt'));cz=ev.vector(ha.parse(r/'CuspZnorm.txt'))
 rem=[]
 for ep in ['.03','.015','.0075']:
  e=mp.mpf(ep);z=mp.mpf('.4');got=-(1+z)*full(e,z);res=got-cm/e-cz
  rem.append({'eps':ep,'z':'.4','normalized_full_cusp':str(got),'residual':str(res),'residual_over_eps':str(res/e)})
 assert all(mp.mpf('1.7')<abs(mp.mpf(rem[i]['residual'])/mp.mpf(rem[i+1]['residual']))<mp.mpf('2.3') for i in range(2))
 rows=[]
 for case,(e,z) in enumerate(((.65,.4),(.7,1.3))):
  exact=float(full(str(e),str(z)));samples=[]
  for seed in range(8):
   pt=points(6,82100+100*case+seed,17);al,norm=alpha(pt,e);xr,y,p,q=[pt[:,i] for i in range(4)];x=p*xr
   ss=polys(x,y,0,al,z)[0];dp=(1-p)*(1-q)+z*p*q
   samples.append(float(norm*np.mean(p*ss**(-1+2*e)*dp**(-1+e))))
  mean=np.mean(samples);se=np.std(samples,ddof=1)/np.sqrt(8);delta=mean-exact;tol=8*se+1e-4*abs(exact)
  row={'eps':e,'z_D':z,'raw_boundary':float(mean),'Gamma_reduced':exact,'difference':float(delta),'empirical_SE':float(se),'tolerance':float(tol),'passed':bool(abs(delta)<=tol),'samples':samples};rows.append(row);print('CUSP',e,z,delta/exact,row['passed'],flush=True)
 d={'passed':all(r['passed'] for r in rows),'raw_boundary_checks':rows,'cusp_remainder_samples':rem,'working_precision':70,'scope':'only the explicit cusp, not the full I remainder','formal':'not_run'}
 dump(HERE/'cusp-numeric.json',d);assert d['passed']
if __name__=='__main__':main()