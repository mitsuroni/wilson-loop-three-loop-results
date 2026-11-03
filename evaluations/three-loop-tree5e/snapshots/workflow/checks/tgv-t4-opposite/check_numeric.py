"""Independent full-epsilon ordered-domain, boundary/recoil and factorization checks."""
import sys,time
import numpy as np
import mpmath as mp
from scipy.special import betaincinv,gammaln
from scipy.stats import qmc
from check_math import ROOT,HERE,read,dump
sys.path.insert(0,str(ROOT/'workflow/checks/tgv-mixed'))
from numeric import make_function
sys.path.insert(0,str(ROOT/'workflow/checks/tgv-t2-empty'))
from stable_qmc import evaluate

def points(dim,seed,power=17):return qmc.Sobol(dim,scramble=True,seed=seed,bits=30).random_base2(power)+2**-31

def alpha(pt,e):
 a=1-e;uu,vv=pt[:,-2],pt[:,-1]
 aa=betaincinv(a,2*a,uu);ac=betaincinv(2*a,a,1-uu);b=betaincinv(a,a,vv);bc=betaincinv(a,a,1-vv)
 al=(aa,ac*b,ac*bc);assert all(np.all(v>0) for v in al)
 return al,np.exp(3*gammaln(a)-gammaln(3*a))

def polys(x,y,w,al,z):
 a,b,c=al
 S=a*b*y*(1-x)+a*c*((1-x)*(1-w)+z*x*w)+z*b*c*w*(1-y)
 P=(a*b*(-(1+z)*(1-x)+z*y)+a*c*(1+3*z-(1+z)*x+z*(1+z)*w)+z*b*c*(1-y-(1+z)*w))/4
 R=(1+z)*c*(a*(1-x+z*w)-z*b*w)/2
 return S,P,R

def ordered(e,z,oid,pt,ibp=False):
 al,norm=alpha(pt,e);xi,y,wi,p,q=[pt[:,i] for i in range(5)]
 x=p*xi
 if oid==27:w=q*wi;span=q;wlo=0;wup=q
 else:w=q+(1-q)*wi;span=1-q;wlo=q;wup=1
 jac=p*span;pure=((1-p)*(1-q)+z*p*q)**(-1+e)
 S,P,R=polys(x,y,w,al,z)
 if not ibp:val=P*S**(-2+2*e)*jac
 else:
  bnd=z*span*(polys(p,y,w,al,z)[0]**(-1+2*e)-polys(0,y,w,al,z)[0]**(-1+2*e))
  bnd+=(1+z)*jac*(polys(x,1,w,al,z)[0]**(-1+2*e)-polys(x,0,w,al,z)[0]**(-1+2*e))
  bnd-=p*(polys(x,y,wup,al,z)[0]**(-1+2*e)-polys(x,y,wlo,al,z)[0]**(-1+2*e))
  val=bnd/(4*(1-2*e))+R*S**(-2+2*e)*jac
 vals=-(1+z)*norm*val*pure;assert np.all(np.isfinite(vals));return float(np.mean(vals))

def core(e,z,pt):
 al,norm=alpha(pt,e);S,P,_=polys(pt[:,0],pt[:,1],pt[:,2],al,z)
 return float(norm*np.mean(P*S**(-2+2*e)))

def pure(e,z):
 with mp.workdps(55):
  e,z=mp.mpf(str(e)),mp.mpf(str(z))
  def ff(p):
   if p==0:return 1/e
   if p==1:return z**(e-1)/e
   L=mp.log(z*p/(1-p))
   return (1-p)**(e-1)*(mp.expm1(e*L)/(e*mp.expm1(L)) if L else 1)
  return -(1+z)*mp.quad(ff,[0,1/(1+z),1])

def stats(xs):return float(np.mean(xs)),float(np.std(xs,ddof=1)/np.sqrt(len(xs)))
def compare(label,a,b):
 am,ae=stats(a);bm,be=stats(b);diff=am-bm;tol=8*np.hypot(ae,be)+1e-4*max(abs(am),abs(bm))
 return {'label':label,'left':am,'right':bm,'difference':diff,'empirical_SE_left':ae,'empirical_SE_right':be,'relative_difference':diff/bm,'tolerance':float(tol),'passed':bool(abs(diff)<=tol),'left_samples':list(a),'right_samples':list(b)}

def main():
 start=time.time();src={g['orbit_id'][-2:]:g for g in read(ROOT/'benchmarks/tgv/5e-3g/sources/canonical-r001.json')['orbits']}
 ibps=[];raws=[];factors=[];refs=[];computed={}
 for case,(e,z) in enumerate(((.65,.4),(.7,1.3))):
  jd={}
  for zi,zz in enumerate((z,1/z)):
   for oid in (27,28):
    aa=[ordered(e,zz,oid,points(7,71000+1000*case+100*zi+10*oid+k)) for k in range(8)]
    bb=[ordered(e,zz,oid,points(7,73000+1000*case+100*zi+10*oid+k),True) for k in range(8)]
    jd[oid,zi]=aa;r=compare(f'IBP {oid} eps={e} zD={zz}',aa,bb);ibps.append(r);print(r['label'],r['relative_difference'],r['passed'],flush=True)
  input_samples={}
  for oid in ('27','28','29','41','42'):
   fn=make_function(src[oid]['J_raw'],e,z);aa=[evaluate(fn,points(7,75600+1000*case+10*int(oid)+k),e) for k in range(8)];input_samples[oid]=aa
   if oid in ('27','28'):bb=[z**(3*e)*v for v in jd[int(oid),1]]
   elif oid=='29':bb=jd[27,0]
   else:continue
   r=compare(f'RAW {oid} eps={e} z={z}',aa,bb);raws.append(r);print(r['label'],r['relative_difference'],r['passed'],flush=True)
  jp=float(pure(e,z));cc=[core(e,z,points(5,78400+1000*case+k)) for k in range(8)]
  fac=compare(f'PARTITION 41+42 eps={e} z={z}',[a+b for a,b in zip(input_samples['41'],input_samples['42'])],[jp*c for c in cc]);factors.append(fac);print(fac['label'],fac['relative_difference'],fac['passed'],flush=True)
  fn=make_function(src['41']['J_raw'],e,1/z);aa=[evaluate(fn,points(7,79500+1000*case+k),e) for k in range(8)]
  ref=compare(f'REFLECTION 42/41 eps={e} z={z}',input_samples['42'],[z**(3*e)*a for a in aa]);refs.append(ref);print(ref['label'],ref['relative_difference'],ref['passed'],flush=True)
  computed[str(case)]={'eps':e,'z':z,'root_samples':input_samples,'extra_41_inverse_z':aa,'pure_J':jp,'bare_T4_samples':cc}
  dump(HERE/'numeric-progress.json',{'ordered_IBP':ibps,'original_root_comparisons':raws,'factorization':factors,'reflection':refs})
 report={'passed':all(r['passed'] for r in ibps+raws+factors+refs),'ordered_IBP':ibps,'original_root_comparisons':raws,'factorization':factors,'reflection':refs,'root_sampling':computed,'raw_root_sample_cases':12,'scrambles':8,'points_per_scramble':2**17,'seconds':time.time()-start,'error_kind':'empirical QMC estimates, no rigorous enclosure','small_e_remainder_check':'not_run','formal':'not_run'}
 dump(HERE/'numeric.json',report);assert report['passed'];print('NUMERIC_COMPLETE',report['seconds'],flush=True)
if __name__=='__main__':main()