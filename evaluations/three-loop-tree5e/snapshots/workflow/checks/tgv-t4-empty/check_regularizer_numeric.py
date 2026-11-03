"""Independent numerical comparison of raw and Euler-regularized T_w at nonzero epsilon."""
import sympy as sp,numpy as np,time,sys
sys.setrecursionlimit(30000)
from check_math import ROOT,HERE,read,dump
from check_numeric import alpha,points,polys

def main():
 start=time.time()
 rr,ss_,tt,uu,vv,zz,ee,aa,bb=sp.symbols('r s t u v z eps aa bb');ss=(rr,ss_,tt,uu,vv,zz,ee,aa,bb)
 FF=rr*tt*uu*vv*zz+rr*ss_*uu*vv+rr*tt*vv*zz+ss_*tt*vv*zz+rr*ss_*uu+ss_*vv*zz+rr*uu+rr
 NN=(1+vv)*(rr*(zz*tt-1)*(1+uu)+zz*ss_*(1+tt))
 E1=lambda f:rr*sp.diff(f,rr)+vv*sp.diff(f,vv)
 E2=lambda f:vv*sp.diff(f,vv)-tt*sp.diff(f,tt)
 l1=-ee+(-1-ee)*rr/(1+rr+ss_)+(-1-2*ee)*vv/(1+vv)+(-2+2*ee)*E1(FF)/FF+bb*ee/(1+vv)
 l2=(1+2*ee)*tt/(1+tt)-(1+2*ee)*vv/(1+vv)+(-2+2*ee)*E2(FF)/FF+bb*ee/(1+vv)
 rg1=(((1+bb)*ee-2)*NN-E1(NN)-NN*l1)/((1+bb)*ee)
 rg=((2+bb)*ee*rg1-E2(rg1)-rg1*l2)/((2+bb)*ee)
 pp=ROOT/read(HERE/'current-regularizer_points.json')['run']/'regularizer-points.txt'
 cert=[]
 point_lines=pp.read_text().splitlines();assert len(point_lines)==6, "expected six exported operator points"
 for line in point_lines:
  k,B,val=line.split();k=int(k);B=int(B)
  pt={rr:k,ss_:k+1,tt:k+2,uu:k+3,vv:k+4,zz:sp.Rational(2,5),ee:sp.Rational(13,20),aa:0,bb:B}
  got=rg.subs(pt);assert got==sp.Rational(val),(k,B)
  cert.append({'k':k,'B':B,'exact_difference':'0'})
 fn=sp.lambdify(ss,rg,'numpy',cse=True)
 rows=[]
 for case,(e,z) in enumerate(((.65,.4),(.7,1.3))):
  samples={key:[] for key in [(a,b) for a in range(2) for b in range(2)]}
  for seed in range(8):
   pt=points(5,65200+case*100+seed,16);al,norm=alpha(pt,e)
   x,y,w=[pt[:,k].astype(np.longdouble) for k in range(3)];a1,a2,a3=[a.astype(np.longdouble) for a in al]
   rr=a1/a3;s=a2/a3;t=x/(1-x);u=y/(1-y);v=w/(1-w)
   S,_,_=polys(x,y,w,(a1,a2,a3),np.longdouble(z));sw=-a1*a3*(1-x)+z*a1*a3*x+z*a2*a3*(1-y)
   scale=a3*a3*(1-x)*(1-y)*(1-w)
   for A,B in samples:
    W=(1-y)**(A*e)*w**(B*e)
    original=sw*S**(-2+2*e)*W
    transformed=fn(rr,s,t,u,v,np.longdouble(z),np.longdouble(e),A,B)*scale*S**(-2+2*e)*W
    assert np.all(np.isfinite(transformed))
    samples[A,B].append((float(norm*np.mean(original)),float(norm*np.mean(transformed))))
  for (A,B),vals in samples.items():
   vv=np.asarray(vals);means=vv.mean(axis=0);diff=vv[:,1]-vv[:,0];delta=diff.mean();se=diff.std(ddof=1)/np.sqrt(8);tol=8*se+1e-4*max(abs(means))
   row={'A':A,'B':B,'eps':e,'z_D':z,'original_Tw':float(means[0]),'regularized_Tw':float(means[1]),'difference':float(delta),'empirical_SE_paired':float(se),'tolerance':float(tol),'passed':bool(abs(delta)<=tol),'samples':vals}
   rows.append(row);print('REG',A,B,e,z,delta,abs(delta)/max(abs(means)),row['passed'],flush=True)
 report={'passed':all(r['passed'] for r in rows),'checks':rows,'independent_operator_point_certificates':cert,'points_per_scramble':2**16,'scrambles':8,'seconds':time.time()-start,'arithmetic':'float64 beta quantiles; numpy longdouble rational kernel evaluation','scope':'full nonzero-epsilon regularizer, not Laurent remainder verification','formal':'not_run'}
 dump(HERE/'regularizer-numeric.json',report);assert report['passed']
if __name__=='__main__':main()