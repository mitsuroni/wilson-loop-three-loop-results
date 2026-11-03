from pathlib import Path
from functools import lru_cache
import json
import sympy as s
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
mp.mp.dps=55;records=[]
@lru_cache(None)
def H(e,z,c=0,d=0):
    e,z=mp.mpf(e),mp.mpf(z);c,d=mp.mpf(c),mp.mpf(d)
    if z>1:return H(str(e),str(1/z),str(d),str(c))
    delta=2*e
    return -(1+z)*z**(-e)*mp.meijerg([[delta,0,-c*delta],[]],[[0,delta,(d+1)*delta],[]],z)/(mp.gamma(1-delta)*mp.gamma(1+delta)*mp.gamma(1+(c+d+1)*delta))
def J(e,z):
    Q=mp.gamma(1+e)*mp.gamma(1+2*e)/mp.gamma(1+3*e)
    return z**(2*e)*(H(str(e),str(z))-H(str(e),str(z),'.5',0))/e**2+z**e*(z**(2*e)+3*(1-Q))/(12*e**3*(1+e))
def full(e,z):
    G=mp.exp(-mp.euler*e)*mp.gamma(1-e)
    return G**3*(J(e,z)+z**(3*e)*J(e,1/z))/(32*e*(1-2*e)*(1+z**(3*e)))
def source(e,z,n,reverse=False):
    nodes,w=leggauss(n);nodes=(nodes+1)/2;w=w/2
    rho=nodes[:,None]**(1/(2*e));t=nodes[None,:]
    x=rho*t;y=1-rho*(1-t)
    h=1-t+z*t-(1+z)*rho*t*(1-t)
    K=(-2*(1+z)*rho*h+(1-2*e)*z)/(2*e*(1+e))*h**(-2+2*e)
    weight=(2-(1-x)**e-(1-y)**e)/2 if reverse else (2-x**e-y**e)/2
    return (1 if reverse else z**e)/e**2*float(np.einsum('i,ij,j->',w,K*weight,w))
def moment(k,n,z):
    b=z/(1+z);l=mp.log(z)
    def f(t):
        if t==0 or t==1:return mp.mpf(0)
        W=mp.log(t)+mp.log1p(-t);T=l+mp.log1p(-t)-mp.log(t)
        if abs(t-b)<mp.mpf('1e-45'):return -W**k/(b*(1-b)) if n==1 else mp.mpf(0)
        return W**k*T**n/(t-b)
    return mp.quad(f,[0,b,1])
def main():
    for es,zs in [('.19','.5'),('.27','1'),('.33','2')]:
        e,z=mp.mpf(es),mp.mpf(zs)
        for reverse in [False,True]:
            vals=[source(float(e),float(z),n,reverse) for n in [128,256]]
            expected=z**(3*e)*J(e,1/z) if reverse else J(e,z)
            err=abs(mp.mpf(vals[-1])-expected)
            assert err<mp.mpf('2e-5'),(es,zs,reverse,err)
            records.append({'type':'complete weighted source K quadrature','orbit':'10' if reverse else '03','eps':es,'z':zs,'values':vals,'absolute_error':str(err),'passed':True})
        print('weighted source '+es+' '+zs+' passed',flush=True)
    result=json.loads((BASE/'results/4e-Pi1-03/r001.json').read_text(encoding='utf-8'))
    names=sorted(set().union(*(s.sympify(v['expression']).free_symbols for v in result['laurent']['coefficients'].values())),key=str)
    fs={int(p):s.lambdify(names,s.sympify(v['expression']),'mpmath') for p,v in result['laurent']['coefficients'].items()}
    for zs in ['.5','1','2']:
        z=mp.mpf(zs);values=[]
        for name in names:
            n=str(name);values.append(mp.log(z) if n=='L' else moment(int(n[1]),int(n[2]),z))
        coeff={p:f(*values) for p,f in fs.items()};errors=[]
        for den in [64,128,256,512]:
            e=mp.mpf(1)/den;errors.append(abs(full(e,z)-sum(c*e**p for p,c in coeff.items())))
        ratios=[errors[i+1]/errors[i] for i in range(3)];assert abs(ratios[-1]-.5)<.03,(zs,ratios)
        records.append({'type':'independent M defining moments and finite-part remainder','z':zs,'moments':{str(k):str(v) for k,v in zip(names,values)},'errors':list(map(str,errors)),'ratios':list(map(str,ratios)),'passed':True})
        print('M moments/remainder '+zs+' passed',flush=True)
    (HERE/'numeric-checks.json').write_text(json.dumps({'count':len(records),'checks':records,'mpmath_dps':55,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls passed')
if __name__=='__main__':main()
