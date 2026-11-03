"""Three-dimensional source quadratures, independent Hlog rho1 and remainder tests."""
from pathlib import Path
from functools import lru_cache
import json,math
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
mp.mp.dps=55;records=[]
@lru_cache(None)
def H(e,z,c=0,d=0):
    e,z=mp.mpf(e),mp.mpf(z)
    if z>1:return H(str(e),str(1/z),d,c)
    return -(1+z)*z**(-e/2)*mp.meijerg([[e,0,-c*e],[]],[[0,e,(d+1)*e],[]],z)/(mp.gamma(1-e)*mp.gamma(1+e)*mp.gamma(1+(c+d+1)*e))
def boundary(e,z):
    l=mp.log(z);p=H(str(e),str(z),2,0);m=H(str(e),str(z),0,2)
    A=H(str(3*e),str(z))+mp.cosh(e*l)*(H(str(e),str(z))-p-m)-mp.sinh(e*l)*(p-m)
    return A/(2*e*(1+e))
@lru_cache(None)
def quad(e,z,n):
    nodes,w=leggauss(n);nodes=(nodes+1)/2;w=w/2
    v=nodes[None,:,None];r=nodes[None,None,:];y=v*r;B=1-y;h=(1+z)*y-1
    values=np.zeros(6)
    for start in range(0,n,8):
        u=nodes[start:start+8,None,None];weights=w[start:start+8,None,None]*w[None,:,None]*w[None,None,:]*v
        logB=np.log(B);q=np.log1p(u*h/B);logD=logB+q
        outer=(1-u)*(1-v)+z*u*v;logouter=np.log(outer)
        def divided(a):
            if a==0:return q/h
            return np.exp(a*logB)*np.expm1(a*q)/(a*h)
        common=(1+z)**2
        t21=common*z**(-1.5*e)*divided(2*e)*outer**(e-1)
        t12=common*z**(-1.5*e)*divided(e)*outer**(2*e-1)
        se_prefix=-(1+z)/(1+e)*divided(2*e)-z*(2*e-1)/(2*(1+e))*divided(2*e-1)
        j7=se_prefix*(-1-z)*outer**(e-1)
        se_outer=(-2*(1+z)*outer+(1-2*e)*z)/(2*(1+e))*outer**(-2+2*e)
        j8=se_outer*(-1-z)*divided(e)
        k0=common*divided(0)/outer
        rho1=k0*(logouter+(logD+logB)/2-np.log(z))
        for i,arr in enumerate([t21,t12,j7,j8,k0,rho1]):values[i]+=np.sum(weights*arr)
    return tuple(map(float,values))
def hlog(word,z):
    if not word:return mp.mpf(1)
    if all(x==0 for x in word):return mp.log(z)**len(word)/mp.factorial(len(word))
    assert word.count(-1)==1 and set(word)<=set([0,-1])
    a=word.index(-1);b=len(word)-a-1
    return mp.quad(lambda t:mp.log(z/t)**a*mp.log(t)**b/(mp.factorial(a)*mp.factorial(b)*(1+t)),[0,z])
@lru_cache(None)
def rho1(z):
    z=mp.mpf(z);basis=json.loads((ROOT/'benchmarks/2ne/families/same-channel-two/structural-basis-functions.json').read_text())['functions']['rho_id_0_0_1']
    return mp.fsum(mp.mpf(t['coefficient'])*hlog(t['word'],z)*mp.fprod(mp.zeta(k)**p for k,p in zip([2,3,5],t['zeta_powers'])) for t in basis['hlog_terms'])
def moment(k,n,z):
    b=z/(1+z);l=mp.log(z)
    def f(t):
        if t==0 or t==1:return mp.mpf(0)
        W=mp.log(t)+mp.log1p(-t);T=l+mp.log1p(-t)-mp.log(t)
        if abs(t-b)<mp.mpf('1e-45'):return -W**k/(b*(1-b)) if n==1 else mp.mpf(0)
        return W**k*T**n/(t-b)
    return mp.quad(f,[0,b,1])
def coeff(z):
    l=mp.log(z);m01=-(l*l+mp.pi**2)/2;m11=moment(1,1,z)
    return (m01*m01/384-m11/256,(3*rho1(str(z))-moment(2,1,z)/2-moment(0,3,z)/3+l*l*m01/4+2*m01*m01/3-m11/2)/128)
def main():
    for es,zs in [('.17','.5'),('.23','1'),('.31','2')]:
        e,z=mp.mpf(es),mp.mpf(zs);lo=quad(float(e),float(z),64);hi=quad(float(e),float(z),128)
        reduced=mp.mpf(hi[0]+hi[1])+boundary(e,z)
        original=mp.mpf(hi[2]+hi[3])/z**(1.5*e);err=abs(original-reduced)
        assert err<mp.mpf('3e-4'),(es,zs,err)
        records.append({'type':'two full source kernels versus complete metric and eight-corner reduction','eps':es,'z':zs,'low':lo[:4],'high':hi[:4],'absolute_residual':str(err),'passed':True})
        print('source/reduction '+es+' '+zs+' '+str(err),flush=True)
    for zs in ['.5','1','2']:
        z=mp.mpf(zs);m01=-(mp.log(z)**2+mp.pi**2)/2
        lo=quad(0,float(z),64);hi=quad(0,float(z),128)
        errs=[abs(mp.mpf(hi[4])-m01*m01/6),abs(mp.mpf(hi[5])-rho1(str(z)))]
        assert max(errs)<mp.mpf('5e-4'),(zs,errs)
        records.append({'type':'reused rho0/rho1 Hlog versus independent ordered integral','z':zs,'errors':list(map(str,errs)),'low':lo[4:],'high':hi[4:],'passed':True})
        pole,finite=coeff(z);errors=[];uncertainties=[]
        for den in [32,64,128]:
            e=mp.mpf(1)/den;lo=quad(float(e),float(z),96);hi=quad(float(e),float(z),192)
            G=mp.exp(-mp.euler*e)*mp.gamma(1-e);pref=G**3/(128*e*(1-2*e)*mp.cosh(1.5*e*mp.log(z)))
            value=pref*(mp.mpf(hi[0]+hi[1])+boundary(e,z))
            errors.append(value-pole/e-finite);uncertainties.append(abs(pref*((hi[0]+hi[1])-(lo[0]+lo[1]))))
        ratios=[abs(errors[i+1]/errors[i]) for i in range(2)]
        assert mp.mpf('.35')<ratios[-1]<mp.mpf('.65'),(zs,ratios,errors,uncertainties)
        records.append({'type':'O(eps) remainder with quadrature resolution comparison','z':zs,'eps_denominators':[32,64,128],'remainders':list(map(str,errors)),'halving_ratios':list(map(str,ratios)),'resolution_differences':list(map(str,uncertainties)),'passed':True})
        print('rho/finite remainder '+zs+' passed',flush=True)
    (HERE/'numeric-checks.json').write_text(json.dumps({'count':len(records),'checks':records,'mpmath_dps':55,'limitation':'ordered three-dimensional quadratures use float64; resolution differences are not certified errors; mixed series not summed to high precision','formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls passed')
if __name__=='__main__':main()
