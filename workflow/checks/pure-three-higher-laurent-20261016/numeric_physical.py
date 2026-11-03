"""New physical epsilon^3: independent original ordered-kernel quadrature.

The route integrates one linear kernel coordinate analytically and retains its
log-interval variance. It does not use the new GPL-word endpoint integration.
"""
from functools import lru_cache
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import mpmath as mp
import numpy as np
from scipy.special import roots_legendre
import hlog9 as h

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
from numeric_hlog import HlogEvaluator

@lru_cache(None)
def grid(n):
    x,w=roots_legendre(n);x=(x.astype(np.longdouble)+1)/2;b=1-x;d=x**3+b**3
    return x**3/d,w.astype(np.longdouble)/2*3*x*x*b*b/d**2

def primitive(lo,hi,s,z):
    slope=(1+z)*s-1;dl=(1-lo)*(1-s)+z*lo*s;dh=(1-hi)*(1-s)+z*hi*s
    diff=np.log(dh)-np.log(dl);ratio=(hi-lo)*slope/dl
    diff=np.where(abs(ratio)<.2,np.log1p(np.maximum(ratio,np.longdouble('-.5'))),diff)
    shape=np.broadcast_shapes(np.shape(diff),np.shape(slope))
    value=np.divide(diff,slope,out=np.broadcast_to((hi-lo)/dl,shape).copy(),where=np.broadcast_to(abs(slope)>1e-12,shape))
    return value,(np.log(dh)+np.log(dl))/2,diff

def interval_polys(center,width,omit_variance=False):
    variance=0 if omit_variance else width*width/24
    return [np.ones_like(center),center,center*center/2+variance,center**3/6+center*variance]

def raw(z,n,omit_variance=False):
    z=np.longdouble(str(z));x,wt=grid(n)
    u=x[:,None,None];w=x[None,:,None];v=w*x[None,None,:]
    a,ac,ad=primitive(0,u,v,z)
    b,bc,bd=primitive(u,1,w,z);c,cc,cd=primitive(w,1,u,z)
    center=ac+bc+cc-1.5*np.log(z)
    variance=0 if omit_variance else (ad*ad+bd*bd+cd*cd)/24
    poly=[np.ones_like(center),center,center*center/2+variance,center**3/6+center*variance]
    measure=wt[:,None,None]*wt[None,:,None]*wt[None,None,:]*w
    U=[float(-(1+z)**3*np.sum(measure*a*b*c*pp)) for pp in poly]
    r=w+(1-w)*x[None,None,:];b,bc,bd=primitive(u,1,r,z)
    am=[np.sum(a*pp*wt[None,None,:]*w,axis=2) for pp in interval_polys(ac,ad,omit_variance)]
    bm=[np.sum(b*pp*wt[None,None,:]*(1-w),axis=2) for pp in interval_polys(bc,bd,omit_variance)]
    mid=(1-x[:,None])*(1-x[None,:])+z*x[:,None]*x[None,:]
    phi=np.log(mid)-1.5*np.log(z)
    kernels=[sum(am[i]*bm[j]*phi**(q-i-j)/math.factorial(q-i-j)
                 for i in range(q+1) for j in range(q-i+1)) for q in range(4)]
    A=[float(-(1+z)**3*np.sum(pp/mid*wt[:,None]*wt[None,:])) for pp in kernels]
    return np.array(A),np.array(U)

def physical(z,n,omit_variance=False):
    A,U=raw(z,n,omit_variance);Ai,Ui=raw(1/mp.mpf(str(z)),n,omit_variance)
    B=(U+Ui)/2
    p2=float(3*mp.zeta(2)/2-9*mp.log(mp.mpf(str(z)))**2/8);p3=float(mp.zeta(3))
    a=A.copy();b=B.copy()
    a[2]+=p2*A[0];b[2]+=p2*B[0]
    a[3]+=p2*A[1]+p3*A[0];b[3]+=p2*B[1]+p3*B[0]
    return {'A':A,'B':B,'A_inverse':Ai,'I20':-a/16,'I21':-b/32}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run');args=ap.parse_args()
    run=ROOT/args.run;path=run/'coefficients.json'
    data=json.loads(path.read_text(encoding='utf-8'))
    vectors={name:h.lift(data['bare'][name][3]) for name in ('A','B')}
    vectors.update({name:h.lift(data['physical'][oid][3]) for name,oid in (('I20','6e-20'),('I21','6e-21'))})
    records=[];negatives=[]
    with mp.workdps(55):
        z53=mp.nsum(lambda n:(mp.zeta(3)-mp.zeta(3,n))/n**5,[1,mp.inf])
        periods=[mp.zeta(2),mp.zeta(3),mp.zeta(5),mp.zeta(7),z53,mp.zeta(9)]
        saved_periods=[str(x) for x in periods]
        for zz in ('0.4','0.7'):
            ev=HlogEvaluator(zz,384)
            print('original quadratures z='+zz,flush=True)
            low=physical(zz,64);high=physical(zz,108)
            expected={}
            for name,v in vectors.items():
                value=float(mp.fsum(mp.mpf(c.numerator)/c.denominator*mp.fprod(a**p for a,p in zip(periods,ps))*ev.word(word)
                                    for (word,ps),c in v.items()))
                expected[name]=value
                error=abs(high[name][3]-value)
                refinement=abs(high[name][3]-low[name][3])
                tolerance=max(10*refinement,3e-6*max(1,abs(value)))
                assert error<tolerance,(zz,name,value,high[name][3],error,tolerance)
                records.append({'z':zz,'target':name,'epsilon_power':3,'Hlog':value,
                                'quadrature64':float(low[name][3]),'quadrature108':float(high[name][3]),
                                'absolute_error':float(error),'refinement_change':float(refinement),'tolerance':float(tolerance)})
            assert abs(high['A'][3]-high['A_inverse'][3])<1e-5*(1+abs(high['A'][3]))
            wrong=physical(zz,108,True)
            for name in ('A','B'):
                residual=abs(wrong[name][3]-expected[name])
                assert residual>1e-4*max(1,abs(expected[name]))
                negatives.append({'z':zz,'target':name,'mutation':'omit interval log variance at epsilon cubed',
                                  'nonzero_difference':float(residual),'rejected':True})
            print('new physical coefficient checks passed at z='+zz,flush=True)
    report={'passed':True,'formal':'not_run','cases':records,'negative_controls':negatives,
            'period_values':saved_periods,'Hlog_digits':55,'Hlog_series_terms':384,
            'coefficient_file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'route':'Original three/four-dimensional ordered-kernel endpoint integrals; analytic linear primitive moments then graded tensor quadrature.',
            'error_scope':'Observed quadrature refinement, not rigorous intervals or convergence-rate bounds.'}
    dest=run/'numeric-verification.json';assert not dest.exists()
    dest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'cases':len(records),'negative_controls':len(negatives),
                      'max_absolute_error':max(r['absolute_error'] for r in records)}))

if __name__=='__main__':main()
