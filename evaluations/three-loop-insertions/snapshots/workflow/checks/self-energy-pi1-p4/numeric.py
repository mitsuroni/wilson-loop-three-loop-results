from pathlib import Path
from functools import lru_cache
import json
import mpmath as mp
import numpy as np
from scipy.special import roots_jacobi
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
mp.mp.dps=50
records=[]
@lru_cache(None)
def H(e,z,c=0,d=0):
    e,z=mp.mpf(e),mp.mpf(z)
    if z>1:return H(str(e),str(1/z),d,c)
    return -(1+z)*z**(-e/2)*mp.meijerg([[e,0,-c*e],[]],[[0,e,(d+1)*e],[]],z)/(mp.gamma(1-e)*mp.gamma(1+e)*mp.gamma(1+(c+d+1)*e))
def full(e,z):
    L=mp.log(z);G=mp.exp(-mp.euler*e)*mp.gamma(1-e)
    h0=H(str(e),str(z));a=H(str(e),str(z),2,0);b=H(str(e),str(z),0,2)
    return G**3*(mp.cosh(e*L)*(h0-(a+b)/2)-mp.sinh(e*L)*(a-b)/2)/(128*e**3*(1-2*e)*(1+e)*mp.cosh(3*e*L/2))
def divided(e,z,y):
    b=1/(1+z)
    if abs(y-b)<mp.mpf('1e-40'):return (z/(1+z))**(e-1)
    return ((z*y)**e-(1-y)**e)/(e*((1+z)*y-1))
def moment(q,z,c=0,d=0):
    b=1/(1+z);L=mp.log(z)
    def f(t):
        if t==0 or t==1:return mp.mpf(0)
        shift=c*mp.log(t)+d*mp.log1p(-t)
        A=mp.log(z*t)-L/2+shift;B=mp.log1p(-t)-L/2+shift
        if abs(t-b)<mp.mpf('1e-40'):return -A**q/(mp.factorial(q)*b*(1-b))
        return -(1+z)*(A**(q+1)-B**(q+1))/(mp.factorial(q+1)*((1+z)*t-1))
    return mp.quad(f,[0,b,1])
def triangle(e,z,n):
    t,wt=roots_jacobi(n,0,2*e);t=(t+1)/2;wt/=2**(2*e+1)
    r,wr=leggauss(n);r=(r+1)/2;wr/=2
    y=t[:,None]*r[None,:];den=(1+z)*y-1
    f=((z*y)**e-(1-y)**e)/(e*den)
    return -(1+z)*z**(2*e)/(2*e*(1+e))*float(np.einsum('i,ij,j->',wt,f,wr))
def main():
    for es,zs in [('.19','.5'),('.27','1'),('.33','2')]:
        e,z=mp.mpf(es),mp.mpf(zs);b=1/(1+z)
        j4=z**(mp.mpf('2.5')*e)*(H(str(e),str(z))-H(str(e),str(z),2,0))/(4*e**2*(1+e))
        j9=z**(e/2)*(H(str(e),str(z))-H(str(e),str(z),0,2))/(4*e**2*(1+e))
        for label,power,weight,expected in [('04',z**(2*e),lambda t:1-t**(2*e),j4),('09',mp.mpf(1),lambda t:1-(1-t)**(2*e),j9)]:
            integral=-(1+z)*power/(4*e**2*(1+e))*mp.quad(lambda t:weight(t)*divided(e,z,t),[0,b,1])
            error=abs(integral-expected);assert error<mp.mpf('1e-35'),(label,error)
            records.append({'type':'reduced source versus MeijerG','orbit':label,'eps':es,'z':zs,'absolute_error':str(error),'passed':True})
        values=[triangle(float(e),float(z),n) for n in [128,256]];error=abs(mp.mpf(values[-1])/j4-1)
        assert error<mp.mpf('2e-6'),error
        records.append({'type':'04 original remaining ordered triangle','eps':es,'z':zs,'values':values,'relative_error':str(error),'passed':True})
        print('source checks '+es+' '+zs+' passed',flush=True)
    for zs in ['.5','1','2']:
        z=mp.mpf(zs);L=mp.log(z);D={};E={}
        for q in [1,2,3]:
            a=moment(q,z,2,0);b=moment(q,z,0,2);D[q]=moment(q,z)-(a+b)/2;E[q]=L*(a-b)/2
        coeff={-2:D[1]/128,-1:(D[1]+D[2]-E[1])/128,0:(D[3]+D[2]-E[1]-E[2]+(3+3*mp.zeta(2)/2-5*L**2/8)*D[1])/128}
        errors=[]
        for den in [64,128,256,512]:
            e=mp.mpf(1)/den;errors.append(abs(full(e,z)-sum(c*e**p for p,c in coeff.items())))
        ratios=[errors[k+1]/errors[k] for k in range(3)];assert abs(ratios[-1]-.5)<.03,(zs,ratios)
        records.append({'type':'independent H0/H2 moments and O(eps) remainder','z':zs,'D':{str(k):str(v) for k,v in D.items()},'eta':{str(k):str(v) for k,v in E.items()},'errors':list(map(str,errors)),'ratios':list(map(str,ratios)),'passed':True})
        print('moments/remainder '+zs+' passed',flush=True)
    (HERE/'numeric-checks.json').write_text(json.dumps({'count':len(records),'checks':records,'mpmath_dps':50,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls passed')
if __name__=='__main__':main()
