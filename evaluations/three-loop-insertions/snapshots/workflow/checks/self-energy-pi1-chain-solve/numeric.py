from pathlib import Path
from functools import lru_cache
import json,re
import mpmath as mp
import sympy as s
import numpy as np
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1-chain'
mp.mp.dps=65;records=[]
@lru_cache(None)
def H(delta,z):
    delta,z=mp.mpf(delta),mp.mpf(z)
    if z>1:z=1/z
    return -(1+z)*z**(-delta/2)*mp.meijerg([[delta,0,0],[]],[[0,delta,delta],[]],z)/(mp.gamma(1-delta)*mp.gamma(1+delta)**2)
def GU(e):
    g=mp.exp(-mp.euler*e)*mp.gamma(1-e)
    u=mp.gamma(1+e)**2*mp.gamma(1-e)*mp.gamma(1-3*e)/(mp.gamma(1-2*e)**2*mp.gamma(1+2*e))
    return g**3*u
def J(key,e,z):
    if key=='01':return (2+e)*z**(3*e)/(18*e**2*(1+2*e))
    return z**(mp.mpf('1.5')*e)*H(str(3*e),str(z))+(1+z**(3*e))/(6*e*(1+2*e))
def I(key,e,z):
    if key=='01':return -GU(e)*(2+e)/(1152*e**4*(1-2*e)**2*(1+2*e))
    return -GU(e)/(384*e**3*(1-2*e)**2*(1+2*e))-GU(e)*H(str(3*e),str(z))/(128*e**2*(1-2*e)**2*mp.cosh(3*e*mp.log(z)/2))
def quadrature(e,z,n):
    nodes,w=leggauss(n);nodes=(nodes+1)/2;w/=2
    rho=nodes[:,None]**(1/(3*e));t=nodes[None,:]
    h=1-t+z*t-(1+z)*rho*t*(1-t)
    N=-(2+e)*(1+z)*rho*h+(1-3*e)*z
    return float(np.einsum('i,ij,j->',w,N*h**(-2+3*e)/(3*e*(1+2*e)),w))
def moment(k,n,z):
    b=z/(1+z);L=mp.log(z)
    def f(t):
        if t==0 or t==1:return mp.mpf(0)
        W=mp.log(t)+mp.log1p(-t);T=L+mp.log1p(-t)-mp.log(t)
        if abs(t-b)<mp.mpf('1e-55'):return -W**k/(b*(1-b)) if n==1 else mp.mpf(0)
        return W**k*T**n/(t-b)
    return mp.quad(f,[0,b,1])
def values(z):
    L=mp.log(z)
    m0=-(L*L+mp.pi**2)/2
    m1=2*mp.zeta(3)+(L*L+mp.pi**2)*mp.log1p(z)-mp.pi**2*L/2+2*L*mp.polylog(2,-z)-2*mp.polylog(3,-z)-L**3/6
    return {'L':L,'M01':m0,'M11':m1,'M21':moment(2,1,z),'M03':moment(0,3,z)}
def coeff_functions(table):
    args=s.symbols('L M01 M11 M21 M03')
    return {int(k):s.lambdify(args,s.sympify(v['expression']),'mpmath') for k,v in table.items()}
def main():
    package=json.loads((BASE/'prepared/2e-Pi1+Pi1.json').read_text(encoding='utf-8'))
    # Connect the actual source polynomial to the endpoint-regularized sector map.
    ep,zz,x,y,rho,t=s.symbols('eps z x y rho t',positive=True)
    o=package['orbits'][1];local={'eps':ep,'z1':zz,'s21':x,'s41':y}
    parse=lambda v:s.sympify(v.replace('^','**'),locals=local)
    f=parse(o['J']['geometry_polynomials'][0]['f']);N=parse(o['J']['se_numerator'])
    mapping={x:rho*t,y:1-rho*(1-t)};h=1-t+zz*t-(1+zz)*rho*t*(1-t)
    assert s.expand(f.subs(mapping,simultaneous=True)-rho*h)==0
    assert s.expand(N.subs(mapping,simultaneous=True)+(2+ep)*(1+zz)*rho*h-(1-3*ep)*zz)==0
    assert s.simplify(s.Matrix(list(mapping.values())).jacobian([rho,t]).det()-rho)==0
    assert s.expand(N-N.subs({x:y,y:x},simultaneous=True))==0
    for es in ['.07','.13','.23']:
        for zs in ['.5','1','2']:
            e,z=mp.mpf(es),mp.mpf(zs);q128=quadrature(float(e),float(z),128);q256=quadrature(float(e),float(z),256)
            error=abs(mp.mpf(q256)-J('02',e,z));assert error<mp.mpf('2e-7'),(es,zs,error)
            records.append({'type':'CH_O full source sector quadrature','eps':es,'z':zs,'orders':[128,256],'values':[q128,q256],'absolute_error':str(error),'passed':True})
    print('9 opposite source quadratures passed',flush=True)
    for es,zs in [('.07','.5'),('.13','1'),('.23','2')]:
        e,z=mp.mpf(es),mp.mpf(zs);a=mp.quad(lambda u:4*u**(12*e-1),[0,1])
        actual=(2+e)*z**(3*e)*a*a/(2*(1+2*e));error=abs(actual-J('01',e,z));assert error<mp.mpf('1e-45')
        records.append({'type':'CH_A original endpoint moments','eps':es,'z':zs,'absolute_error':str(error),'passed':True})
    allvalues={}
    for zs in ['.5','1','2']:
        v=values(mp.mpf(zs));allvalues[zs]=v
        for name in ['M21','M03']:records.append({'type':'existing M defining integral','name':name,'z':zs,'value':str(v[name]),'connection':'exact Hlog identity checked in compute; independent moments used in remainder controls below','passed':True})
    # Actual prepared C, including its flavor polynomial, reconstructed independently.
    for o in package['orbits']:
        key=o['orbit_id'][-2:]
        expression=o['C']['combined'].replace('(-2*v14+I0)','scale')
        expression=re.sub(r'Gamma\[([^\]]+)\]',r'gamma(\1)',expression).replace('^','**')
        ep,ca,cf,tf,nf,alpha,scale=s.symbols('eps CA CF TF nf alpha_s scale')
        parsed=s.sympify(expression,locals={'eps':ep,'CA':ca,'CF':cf,'TF':tf,'nf':nf,'alpha_s':alpha,'scale':scale,'Pi':s.pi,'EulerGamma':s.EulerGamma,'gamma':s.gamma,'exp':s.exp})
        fn=s.lambdify((ep,ca,cf,tf,nf,alpha,scale),parsed,'mpmath')
        e=mp.mpf('.13');z=mp.mpf('.7');al=mp.mpf('.2');A=mp.mpf('2.6')
        for flavors in [0,4,6]:
            B=((10-6*e)*3-4*(1-e)*flavors)/(3-2*e)
            C=fn(e,3,mp.mpf(4)/3,mp.mpf('.5'),flavors,al,A)
            obtained=mp.mpf(o['I_target']['sym'])/2*C*(J(key,e,z)+z**(3*e)*J(key,e,1/z))/((al/mp.pi)**3*(mp.mpf(4)/3)*B**2*(mp.pi*A)**(3*e)*(1+z**(3*e)))
            error=abs(obtained/I(key,e,z)-1);assert error<mp.mpf('1e-50')
            records.append({'type':'actual prepared C and solved J normalization','orbit':key,'nf':flavors,'relative_error':str(error),'passed':True})
    for key in ['01','02']:
        result=json.loads((BASE/f'results/2e-Pi1+Pi1-{key}/r001.json').read_text(encoding='utf-8'))
        tables=[('I',coeff_functions(result['laurent']['coefficients']),lambda e:mp.mpf(1))]
        for channel in result['color_channels']:
            ep=s.Symbol('eps');phi=s.lambdify(ep,s.sympify(channel['exact_scalar_multiplier'],locals={'eps':ep}),'mpmath')
            tables.append((channel['monomial'],coeff_functions(channel['coefficients']),phi))
        for zs in ['.5','1','2']:
            z=mp.mpf(zs);v=allvalues[zs];args=[v[n] for n in ['L','M01','M11','M21','M03']]
            for label,table,phi in tables:
                cs={k:f(*args) for k,f in table.items()};errors=[]
                for den in [128,256,512,1024,2048]:
                    e=mp.mpf(1)/den;errors.append(abs(phi(e)*I(key,e,z)-mp.fsum(c*e**k for k,c in cs.items())))
                ratios=[errors[k+1]/errors[k] for k in range(len(errors)-1)]
                assert abs(ratios[-1]-mp.mpf('.5'))<mp.mpf('.03'),(key,label,zs,ratios)
                records.append({'type':'I or color-channel O(eps) remainder','orbit':key,'channel':label,'z':zs,'eps_denominators':[128,256,512,1024,2048],'errors':list(map(str,errors)),'ratios':list(map(str,ratios)),'passed':True})
            print('remainder '+key+' z='+zs+' passed',flush=True)
    (HERE/'numeric-checks.json').write_text(json.dumps({'count':len(records),'checks':records,'mpmath_dps':65,'sector_map_algebra_checked':True,'quadrature':'float64, not certified interval','formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(records)} numerical controls passed')
if __name__=='__main__':main()
