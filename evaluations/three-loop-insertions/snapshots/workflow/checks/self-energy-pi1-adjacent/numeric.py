"""Independent quadrature of actual ordered source kernels plus epsilon remainders."""
from pathlib import Path
import json, sys
import numpy as np
from scipy.special import roots_jacobi
import mpmath as mp
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parent))
import compute as co

mp.mp.dps=65
records=[]
orbits={o['orbit_id'][-2:]:o for o in json.loads((co.BASE/'prepared/4e-Pi1.json').read_text(encoding='utf-8'))['orbits']}
funcs={key:s.lambdify((co.e,co.z,co.x,co.a,co.b,co.c,co.d),co.parse(orbits[key]['J']['integrand']),'numpy',cse=True) for key in ['01','02','05','06']}
def grid(powers,n):
    nodes=[];weights=[]
    for i,power in enumerate(powers):
        x,w=roots_jacobi(n,0,power)
        shape=[1]*len(powers);shape[i]=n
        nodes.append(((x+1)/2).reshape(shape));weights.append((w/2**(power+1)).reshape(shape))
    return nodes, np.prod(np.broadcast_arrays(*weights),axis=0)
def ordered(key,e,z,n):
    f=funcs[key]
    if key in ['01','02']:
        p,q=(2*e,e) if key=='01' else (e,2*e)
        powers=[p+q-1,q-1,p+q-1,p-1]
        (v,t,w,r),weights=grid(powers,n)
        # s31=1-v, s32=1-v*t; s41=w*r, s42=w. |Jac|=v*w.
        full=f(e,z,0,1-v,1-v*t,w*r,w)*v*w
        density=v**powers[0]*t**powers[1]*w**powers[2]*r**powers[3]
        return float(np.sum(weights*full/density))
    p,q=(e,2*e) if key=='05' else (2*e,e)
    alpha,beta=q,p
    value=0.
    # r=1-s31, t=s32: r+t>1. Partition into r<1/2, t<1/2,
    # and [1/2,1]^2; there is no overlap except boundaries of measure zero.
    for region in ['small-r','small-t','center']:
        hpower=p if region=='small-r' else q if region=='small-t' else 0
        (u,v,h,k),weights=grid([alpha-1,beta-1,hpower,0],n)
        if region=='small-r':r=h/2;t=1-r*k;jac=r/2
        elif region=='small-t':t=h/2;r=1-t*k;jac=t/2
        else:r=(1+h)/2;t=(1+k)/2;jac=1/4
        full=f(e,z,1-u,1-r,t,v,0)*jac
        density=u**(alpha-1)*v**(beta-1)*h**hpower
        value+=float(np.sum(weights*full/density))
    return value
def jclosed(key,e,z):
    if key in ['01','02']:return z**(3*e)/(18*e**4*(1+e))
    Q=mp.gamma(1+e)*mp.gamma(1+2*e)/mp.gamma(1+3*e)
    return z**((1 if key=='05' else 2)*e)*(1-Q)/(4*e**4*(1+e))
def iclosed(key,e,z):
    g=mp.exp(-mp.euler*e)*mp.gamma(1-e)
    if key in ['01','02']:return g**3/(576*e**5*(1-2*e)*(1+e))
    Q=mp.gamma(1+e)*mp.gamma(1+2*e)/mp.gamma(1+3*e)
    return g**3*(1-Q)*mp.cosh(e*mp.log(z)/2)/(128*e**5*(1-2*e)*(1+e)*mp.cosh(3*e*mp.log(z)/2))
def main():
    for key in ['01','02','05','06']:
        for e0 in ['0.17','0.29','0.41']:
            for z0 in ['0.4','1','2.5']:
                ee,zz=mp.mpf(e0),mp.mpf(z0)
                low=ordered(key,float(ee),float(zz),8);high=ordered(key,float(ee),float(zz),16)
                exact=jclosed(key,ee,zz)
                error=float(abs((high-exact)/exact));drift=abs((high-low)/high)
                assert error<2e-11 and drift<2e-9,(key,e0,z0,error,drift)
                records.append({'type':'actual-four-dimensional-source-quadrature','orbit':key,'eps':e0,'z':z0,'orders':[8,16],'relative_error':error,'order_drift':drift,'passed':True})
        print('original ordered source '+key+' passed',flush=True)
    # Build an actual dimensionful prepared-C value from saved components.
    for key in ['01','05']:
        o=orbits[key];comp=o['C']['components']
        check_eps=mp.mpf('.23');zz=mp.mpf('1.7');v14=mp.mpf('-1.3');al=mp.mpf('.2');ca=mp.mpf(3);cf=mp.mpf(4)/3;tf=mp.mpf('.5')
        for nf in [0,4]:
            names={'eps':check_eps,'alpha_s':al,'CA':ca,'CF':cf,'TF':tf,'nf':nf,'v14':v14,'I0':0,'Pi':mp.pi,'EulerGamma':mp.euler,'Gamma':mp.gamma,'exp':mp.exp}
            vals=[]
            for component in comp.values():
                vals.append(eval(component.replace('^','**').replace('Gamma[1-eps]','Gamma(1-eps)'),{'__builtins__':{}},names))
            C=mp.fprod(vals)
            B1=((10-6*check_eps)*ca-8*(1-check_eps)*tf*nf)/(3-2*check_eps)
            external=(al/mp.pi)**3*ca*cf*B1*(-2*mp.pi*v14)**(3*check_eps)*(1+zz**(3*check_eps))
            assembled=(mp.mpf(o['I_target']['sym'])/2)*C*(jclosed(key,check_eps,zz)+zz**(3*check_eps)*jclosed(key,check_eps,1/zz))/external
            error=abs(assembled/iclosed(key,check_eps,zz)-1)
            assert error<mp.mpf('1e-55'),error
            records.append({'type':'prepared-C-reconstruction','orbit':key,'nf':nf,'relative_error':str(error),'passed':True})
    for key in ['01','05']:
        result=json.loads((co.BASE/f'results/4e-Pi1-{key}/r001.json').read_text(encoding='utf-8'))
        coeff={int(k):s.lambdify(co.z,s.sympify(v['expression'],locals={'z':co.z}),'mpmath') for k,v in result['laurent']['coefficients'].items()}
        for z0 in ['0.4','1','2.5']:
            zz=mp.mpf(z0);errors=[]
            for den in [64,128,256,512]:
                ee=mp.mpf(1)/den
                poly=mp.fsum(fun(zz)*ee**k for k,fun in coeff.items())
                errors.append(abs(iclosed(key,ee,zz)-poly))
            ratios=[errors[i+1]/errors[i] for i in range(3)]
            assert abs(ratios[-1]-mp.mpf('.5'))<mp.mpf('.02'),(key,z0,ratios)
            records.append({'type':'O(eps)-remainder-halving','orbit':key,'z':z0,'eps_denominators':[64,128,256,512],'absolute_errors':list(map(str,errors)),'ratios':list(map(str,ratios)),'passed':True})
    out={'count':len(records),'checks':records,'source_quadrature':'double precision Gauss-Jacobi after exact ordered-domain maps; noncertified','closed_and_remainders':'65 decimal digit mpmath','formal':'not_run'}
    (co.HERE/'numeric-checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls passed')

if __name__=='__main__':main()
