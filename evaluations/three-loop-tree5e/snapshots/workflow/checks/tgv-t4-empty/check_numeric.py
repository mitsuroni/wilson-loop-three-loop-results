"""Original seven-dimensional roots, endpoint reduction and full weighted IBP.
All errors are empirical Sobol scramble estimates, never certified intervals.
"""
import sys,time
import numpy as np
from scipy.special import betaincinv,gammaln
from scipy.stats import qmc
from check_math import ROOT,HERE,read,dump,IDS,PAIRS
sys.path.insert(0,str(ROOT/'workflow/checks/tgv-mixed'))
from numeric import make_function
sys.path.insert(0,str(ROOT/'workflow/checks/tgv-t2-empty'))
from stable_qmc import evaluate
BASE=(13,14,15,40)

def points(dim,seed,power=17):return qmc.Sobol(dim,scramble=True,seed=seed,bits=30).random_base2(power)+2**-31

def alpha(pts,e):
 a=1-e;u,v=pts[:,-2],pts[:,-1]
 t=betaincinv(a,2*a,u);tc=betaincinv(2*a,a,1-u)
 b=betaincinv(a,a,v);bc=betaincinv(a,a,1-v)
 als=t,tc*b,tc*bc
 assert all(np.all(x>0) for x in als)
 return als,np.exp(3*gammaln(a)-gammaln(3*a))

def polys(x,y,w,als,z):
 a,b,c=als
 S=a*b*y*(1-x)+a*c*((1-x)*(1-w)+z*x*w)+z*b*c*w*(1-y)
 P=(a*b*(-(1+z)*(1-x)+z*y)+a*c*(1+3*z-(1+z)*x+z*(1+z)*w)+z*b*c*(1-y-(1+z)*w))/4
 R=(1+z)*c*(a*(1-x+z*w)-z*b*w)/2
 return S,P,R

def weights(x,y,w,e):
 A=(1-y)**e;B=w**e
 return (A*(1-B),(1-A)*(1-B),(1-A)*B,1-x**e)

def reduced(e,z,pts):
 al,norm=alpha(pts,e);x,y,w=pts[:,0],pts[:,1],pts[:,2]
 SS,PP,RR=polys(x,y,w,al,z);W=weights(x,y,w,e)
 assert np.all(SS>0)
 return norm*np.array([np.mean(PP*SS**(-2+2*e)*ww) for ww in W])

def ibp(e,z,pts):
 al,norm=alpha(pts,e);x,y,w=pts[:,0],pts[:,1],pts[:,2]
 SS,PP,RR=polys(x,y,w,al,z)
 reco=np.array([np.mean(RR*SS**(-2+2*e)*ww) for ww in weights(x,y,w,e)])
 boundary=np.zeros(4)
 for axis,coeff in [(0,z),(1,1+z),(2,-1)]:
  for end,sgn in [(0,-1),(1,1)]:
   vals=[x,y,w];vals[axis]=end
   ss,_,_=polys(*vals,al,z)
   boundary+=sgn*coeff*np.array([np.mean(ss**(-1+2*e)*ww) for ww in weights(*vals,e)])
 # The e factor in each derivative cancels the inverse-CDF Jacobian 1/e.
 yy=1-y**(1/e);ss,_,_=polys(x,yy,w,al,z)
 dy0=np.mean(ss**(-1+2*e));dy1=np.mean(ss**(-1+2*e)*w**e)
 ww=w**(1/e);ss,_,_=polys(x,y,ww,al,z)
 dw0=np.mean(ss**(-1+2*e));dw1=np.mean(ss**(-1+2*e)*(1-y)**e)
 xx=x**(1/e);ss,_,_=polys(xx,y,w,al,z);dx=np.mean(ss**(-1+2*e))
 deriv=np.array([(1+z)*(dy0-dy1)-dw1,-(1+z)*(dy0-dy1)-dw0+dw1,-(1+z)*dy1+dw0-dw1,z*dx])
 return norm*((boundary+deriv)/(4*(1-2*e))+reco)

def main():
 start=time.time();sources={d['orbit_id'][-2:]:d for d in read(ROOT/'benchmarks/tgv/5e-3g/sources/canonical-r001.json')['orbits']}
 rawchecks=[];ibpchecks=[]
 for case,(e,zv) in enumerate(((.65,.4),(.7,1.3))):
  kvals={};kse={}
  for zz in (zv,1/zv):
   rr=[];ii=[]
   for k in range(8):
    rr.append(reduced(e,zz,points(5,61100+1000*case+(0 if zz==zv else 100)+k)))
    ii.append(ibp(e,zz,points(5,62200+1000*case+(0 if zz==zv else 100)+k)))
   rm=np.mean(rr,axis=0);im=np.mean(ii,axis=0);rs=np.std(rr,axis=0,ddof=1)/np.sqrt(8);is_=np.std(ii,axis=0,ddof=1)/np.sqrt(8)
   kvals[zz]=rm;kse[zz]=rs
   for j,oid in enumerate(BASE):
    diff=im[j]-rm[j];tol=8*np.hypot(rs[j],is_[j])+1e-4*abs(rm[j])
    row={'id':oid,'eps':e,'z_D':zz,'raw_K':float(rm[j]),'IBP_K':float(im[j]),'difference':float(diff),'relative_difference':float(diff/rm[j]),'tolerance':float(tol),'passed':bool(abs(diff)<=tol),'raw_scrambles':[float(v[j]) for v in rr],'IBP_scrambles':[float(v[j]) for v in ii]}
    ibpchecks.append(row);print('IBP',oid,e,zz,row['relative_difference'],row['passed'],flush=True)
  for oid in IDS:
   if oid in PAIRS:base=int(PAIRS[oid]);ref=zv**e/e**2*kvals[zv][BASE.index(base)];se=zv**e/e**2*kse[zv][BASE.index(base)]
   else:base=int(oid);ref=zv**(2*e)/e**2*kvals[1/zv][BASE.index(base)];se=zv**(2*e)/e**2*kse[1/zv][BASE.index(base)]
   fn=make_function(sources[oid]['J_raw'],e,zv)
   samples=[evaluate(fn,points(7,63400+1000*case+10*int(oid)+k),e) for k in range(8)]
   mean=np.mean(samples);err=np.std(samples,ddof=1)/np.sqrt(8);diff=mean-ref;tol=8*np.hypot(err,se)+1e-4*abs(ref)
   row={'id':oid,'eps':e,'z':zv,'original_J':float(mean),'reduced_J':float(ref),'difference':float(diff),'relative_difference':float(diff/ref),'tolerance':float(tol),'passed':bool(abs(diff)<=tol),'scrambles':samples,'empirical_SE_raw':float(err),'empirical_SE_reduced':float(se)}
   rawchecks.append(row);print('RAW',oid,e,zv,row['relative_difference'],row['passed'],flush=True)
   dump(HERE/'numeric-progress.json',{'raw':rawchecks,'weighted_IBP':ibpchecks})
 report={'passed':all(r['passed'] for r in rawchecks+ibpchecks),'raw':rawchecks,'weighted_IBP':ibpchecks,'points_per_scramble':2**17,'scrambles':8,'seconds':time.time()-start,'error_kind':'empirical QMC estimates, no rigorous enclosure','coefficient_remainder_check':'not_run; these positive-epsilon samples verify full domain reduction and weighted IBP, not individual Laurent coefficients','formal':'not_run'}
 dump(HERE/'numeric.json',report);assert report['passed'];print('NUMERIC_COMPLETE',report['seconds'],flush=True)
if __name__=='__main__':main()