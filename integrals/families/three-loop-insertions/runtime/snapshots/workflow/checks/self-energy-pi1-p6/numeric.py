"""P6 independent kernel-product quadrature and finite-part remainder checks."""
from pathlib import Path
import sys,json
import mpmath as mp
import numpy as np
from scipy.special import roots_jacobi
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'workflow/checks/self-energy-pi1-baselines'))
import numeric as baseline
mp.mp.dps=55
H=baseline.H0
records=[]
def bare(e,z,n):
    r,wr=roots_jacobi(n,0,e);r=(r+1)/2;wr=wr/2**(e+1)
    t,wt=leggauss(n);t=(t+1)/2;wt=wt/2
    h=1-t[None,:]+z*t[None,:]-(1+z)*r[:,None]*t[None,:]*(1-t[None,:])
    return -2*(1+z)*float(np.einsum('i,ij,j->',wr,h**(e-1),wt))
def full(e,z):
    G=mp.exp(-mp.euler*e)*mp.gamma(1-e);L=mp.log(z)
    return G**3*H(str(e),str(z))*(H(str(2*e),str(z))+mp.cosh(e*L)/(2*e*(1+e)))/(64*e*(1-2*e)*mp.cosh(3*e*L/2))
def moment(q,z):
    # Integrate one bare endpoint exactly, then Taylor-expand its divided difference.
    b=1/(1+z);L=mp.log(z)
    def f(t):
        if t==0 or t==1:return mp.mpf(0)
        A=mp.log(z*t)-L/2;B=mp.log1p(-t)-L/2
        if abs(t-b)<mp.mpf('1e-45'):return -A**q/(mp.factorial(q)*b*(1-b))
        return -(1+z)*(A**(q+1)-B**(q+1))/(mp.factorial(q+1)*((1+z)*t-1))
    return mp.quad(f,[0,b,1])
def main():
    for es in ['.17','.23','.31']:
        for zs in ['.5','1','2']:
            e,z=mp.mpf(es),mp.mpf(zs)
            vals=[float(z**(2*e))*baseline.quadrature(float(e),float(1/z),n)*bare(float(e),float(z),n) for n in [128,256]]
            J=(z**e*H(str(2*e),str(z))+(1+z**(2*e))/(4*e*(1+e)))*z**(e/2)*H(str(e),str(z))
            err=abs(mp.mpf(vals[-1])-J)
            assert err<mp.mpf('2e-7'),(es,zs,err)
            records.append({'type':'two independent source-kernel quadratures after exact factorization','eps':es,'z':zs,'orders':[128,256],'values':vals,'absolute_error':str(err),'passed':True})
        print('source products eps='+es+' passed',flush=True)
    for zs in ['.5','1','2']:
        z=mp.mpf(zs);L=mp.log(z);t0,t1,t2=[moment(q,z) for q in range(3)]
        assert abs(t0-baseline.m01(z))<mp.mpf('1e-40')
        assert abs(t1-baseline.m11(z)/2)<mp.mpf('1e-40')
        for q,v in enumerate([t0,t1,t2]):records.append({'type':'existing H0 coefficient by defining moment','q':q,'z':zs,'value':str(v),'passed':True})
        coeff={-2:t0/128,-1:(t0+t1+2*t0*t0)/128,0:(t2+t1+(3+mp.mpf('1.5')*mp.zeta(2)-5*L*L/8)*t0+4*t0*t0+6*t0*t1)/128}
        errors=[];herrors=[]
        for den in [64,128,256,512]:
            e=mp.mpf(1)/den
            errors.append(abs(full(e,z)-mp.fsum(c*e**p for p,c in coeff.items())))
            herrors.append(abs(H(str(e),str(z))-t0-e*t1-e*e*t2))
        ratios=[errors[i+1]/errors[i] for i in range(3)]
        hratios=[herrors[i+1]/herrors[i] for i in range(3)]
        assert abs(ratios[-1]-.5)<.03,(zs,ratios)
        assert abs(hratios[-1]-.125)<.01,(zs,hratios)
        records.append({'type':'I O(eps) and H0 O(eps^3) remainder','z':zs,'I_errors':list(map(str,errors)),'I_ratios':list(map(str,ratios)),'H0_ratios':list(map(str,hratios)),'passed':True})
        print('moments and remainders z='+zs+' passed',flush=True)
    (HERE/'numeric-checks.json').write_text(json.dumps({'count':len(records),'checks':records,'mpmath_dps':55,'quadrature':'float64; not certified interval','formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls recorded')
if __name__=='__main__':main()
