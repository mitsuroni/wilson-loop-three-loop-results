"""Independent original ordered-kernel and real-moment controls for the packet."""
from pathlib import Path
import argparse,json,math,time,sys
from functools import lru_cache

sys.dont_write_bytecode=True
import mpmath as mp
import numpy as np
import sympy as s
from family_engines import parse,M,N,eps,ell

H=Path(__file__).resolve().parent
CHECKS=[]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def check(name,value,target,tol=0,**meta):
    if tol:
        err=abs(mp.mpf(value)-mp.mpf(target));good=err<mp.mpf(tol)
        row={'name':name,'passed':bool(good),'value':mp.nstr(value,35),'target':mp.nstr(target,35),'absolute_error':mp.nstr(err,8),'tolerance':str(tol),**meta}
    else:
        residual=s.simplify(value-target);good=residual==0
        row={'name':name,'passed':bool(good),'residual':str(residual),**meta}
    CHECKS.append(row);assert good,row

def exact_polynomial_J(oid,zv,ev):
    p=read(H/'sources/benchmarks/2ne/prepared/4e.json')
    row=next(x for x in p['orbits'] if x['orbit_id']==oid)['J']
    from sympy.parsing.mathematica import parse_mathematica
    parse_m=lambda t:parse_mathematica(t).subs(s.Symbol('z1'),zv)
    f=s.expand(parse_m(row['numerator'])*parse_m(row['propagator_polynomials']['f1'])**(ev-1)*parse_m(row['propagator_polynomials']['f2'])**(ev-1))
    for var in reversed(row['integration_order_outer_to_inner']):
        lo,hi=row['bounds'][var]
        f=s.integrate(f,(s.Symbol(var),s.sympify(lo),s.sympify(hi)))
    return s.simplify(f)

def vm_poly(n,zv):
    u,v=s.symbols('u v');D=(1-u)*(1-v)+zv*u*v
    return s.simplify(-(1+zv)*zv**(-s.Rational(n,2))*s.integrate(s.expand(D**(n-1)),(u,0,1),(v,0,1)))

def vn_integral(e,z):
    r=z/(1+z);lz=mp.log(z)
    def f(b):
        if b==0 or b==1:return mp.mpf(0)
        L=lz+mp.log1p(-b)-mp.log(b);R=L-lz;W=mp.log(b)+mp.log1p(-b)
        if abs(b-r)<mp.mpf('1e-35'):return mp.mpf(0)
        return mp.exp(e*W)*R*(mp.cosh(e*L)-1)/(e**2*(b-r))
    return mp.quad(f,[0,r/2,r,(1+r)/2,1])

def original_kernel_checks():
    for z in [s.Rational(2,3),s.Rational(3,2)]:
        for e in [1,2]:
            for oid in [f'4e-{i:02}' for i in range(1,7)]:
                J=exact_polynomial_J(oid,z,e);Ji=exact_polynomial_J(oid,1/z,e)
                mult=2 if oid=='4e-04' else 1 if oid=='4e-06' else 4
                physical=-s.Rational(mult,16)*(J+z**(2*e)*Ji)/(1+z**(2*e))
                c2=(z**e+z**(-e))/2;c1=(z**s.Rational(e,2)+z**(-s.Rational(e,2)))/2
                V,V2=vm_poly(e,z),vm_poly(2*e,z)
                if oid=='4e-01':target=-s.Rational(1,16*e**4)
                elif oid in {'4e-02','4e-05'}:target=(V2-c1*V)/(4*e**2*c2)
                elif oid=='4e-03':target=(s.gamma(1+e)**2/s.gamma(1+2*e)-1)/(4*e**4*c2)
                elif oid=='4e-06':target=-V**2/(16*c2)
                else:
                    zz=mp.mpf(str(z.p))/int(z.q);VN=vn_integral(mp.mpf(e),zz)
                    Q=mp.mpf(str(s.N(V**2/2,55)))-2*VN+2*mp.harmonic(2*e)/e*mp.mpf(str(s.N(V2,55)))
                    target=-Q/(8*mp.mpf(str(s.N(c2,55))))
                    check(f'{oid} original ordered polynomial kernel to closed Q eps={e},z={z}',mp.mpf(str(s.N(physical,55))),target,'2e-30',Gamma_factor_removed=True)
                    continue
                check(f'{oid} original ordered polynomial kernel to evaluated mother eps={e},z={z}',physical,target,Gamma_factor_removed=True)
    u,v,z=s.symbols('u v z')
    for name,f,d,expected in [('A',z*(1-u)*v,z,2*z*z*(1-u)*v),('O',(1-u)*(1-v)+z*u*v,-(1+z),None)]:
        if expected is None:expected=-2*(1+z)*f+(1-2*eps)*z
        full=2*(1+eps)*d*f+f*s.diff(f,u,v)+(2*eps-1)*s.diff(f,u)*s.diff(f,v)
        check('Pi1 '+name+' complete transverse numerator from mixed derivative',full,expected)
        for ev in [1,2]:
            zv=s.Rational(2,3)
            J=s.integrate(s.expand((expected*f**(2*eps-2)/(2*(1+eps))).subs({eps:ev,z:zv})),(u,0,1),(v,0,1))
            target=zv**(2*ev)/(4*ev**2*(1+ev)) if name=='A' else zv**ev*vm_poly(2*ev,zv)+(1+zv**(2*ev))/(4*ev*(1+ev))
            check(f'Pi1 {name} original kernel/corners eps={ev}',J,target)
    for ev in [mp.mpf(1)/4,mp.mpf(1)/3]:
        def f(b):
            if b==0 or b==1:return mp.mpf(0)
            return b**(-ev-1)*(1-b)**(-ev)*(2*b-1)*(-mp.expm1((-1+2*ev)*mp.log1p(-b)))/(-1+2*ev)
        j=mp.quad(f,[0,mp.mpf('.5'),mp.mpf('.9'),1])/(16*ev**3)
        target=(2*mp.gamma(1+ev)*mp.gamma(1-ev)-mp.beta(1-ev,1-ev))/(16*(1-2*ev)*ev**4)
        check('TGV01 original surviving Beta integral eps='+str(ev),j,target,'2e-10',endpoint_difference_kept_whole=True)

@lru_cache(None)
def moment(k,j,ztext):
    z=mp.mpf(ztext);r=z/(1+z)
    def f(b):
        if b==0 or b==1:return mp.mpf(0)
        W=mp.log(b)+mp.log1p(-b);L=mp.log(z)+mp.log1p(-b)-mp.log(b)
        if abs(b-r)<mp.mpf('1e-35'):return -W**k/(r*(1-r)) if j==1 else mp.mpf(0)
        return W**k*L**j/(b-r)
    return mp.quad(f,[0,r/2,r,(1+r)/2,1])

@lru_cache(None)
def grid(n):
    x,w=np.polynomial.legendre.leggauss(n);x=(x+1)/2;w=w/2
    den=x**3+(1-x)**3;t=x**3/den;tc=(1-x)**3/den
    lt=np.log(t);lc=np.log(tc)
    weights=w*3*x**2*(1-x)**2/den**2
    return t,tc,lt,lc,weights

def tau(z,n):
    t,tc,lt,lc,w=grid(n);D=tc[:,None]*tc[None,:]+z*t[:,None]*t[None,:]
    xp=np.log(D)+lt[:,None]+lt[None,:]-math.log(z)/2
    xm=np.log(D)+lc[:,None]+lc[None,:]-math.log(z)/2
    return -(1+z)/48*np.sum(w[:,None]*w[None,:]*(xp**4+xm**4)/D)

def q04(z,n):
    u,uc,lu,luc,w=grid(n);v,vc,lv,lvc,_=grid(n)
    u=u[:,None];uc=uc[:,None];v=v[None,:];vc=vc[None,:]
    D=uc*vc+z*u*v;dl=np.log(D)-lvc[None,:];du=math.log(z)+lu[:,None]-np.log(D)
    a=(1+z)*v-1;b=(1+z)*u-1
    K=(1+z)**2*dl*du/(a*b)
    phi=np.log(D)+(lvc[None,:]+lu[:,None]-math.log(z))/2
    ww=w[:,None]*w[None,:]
    return [np.sum(ww*K),np.sum(ww*K*phi),np.sum(ww*K*(phi**2/2+(dl**2+du**2)/24))]

def evaluate(expr,ztext,tau_value=0):
    val={ell:mp.log(mp.mpf(ztext)),s.Symbol('tau_D4'):mp.mpf(float(tau_value))}
    for a in expr.free_symbols:
        name=str(a)
        if name.startswith('M_'):
            _,k,j=name.split('_');val[a]=moment(int(k),int(j),ztext)
        elif name.startswith('N_'):
            _,k,j=name.split('_');k,j=int(k),int(j);val[a]=moment(k,j+1,ztext)-val[ell]*moment(k,j,ztext)
    return s.lambdify(tuple(val),expr,'mpmath')(*val.values())

def moment_checks():
    targets=read(H/'sources/workflow/checks/tgv04-high-order/runs/positive-epsilon/numeric-targets.json')
    for text in ['0.5','1.0','2.0']:
        z=mp.mpf(text);m01=moment(0,1,text);l=mp.log(z)
        check('original M01 definition z='+text,m01,-(l*l+mp.pi**2)/2,'1e-30')
        check('original M03 definition z='+text,moment(0,3,text),-m01*m01,'1e-28')
        m11=(l*l+mp.pi**2)*mp.log1p(z)-l**3/6+2*l*mp.polylog(2,-z)-mp.pi**2*l/2-2*mp.polylog(3,-z)+2*mp.zeta(3)
        check('original M11 definition z='+text,moment(1,1,text),m11,'1e-29')
        tq1,tq2=tau(float(z),320),tau(float(z),640)
        check('independent tauD4 mesh convergence z='+text,tq2,tq1,'3e-7',empirical_not_interval=True)
        target=next(t for t in targets if float(t['z'])==float(z))
        for q in [1,2]:
            expr=parse(read(H/'3e-3g-04/result.json')['coefficients'][str(q)]['sympy'])
            val=evaluate(expr,text,tq2)
            check(f'TGV04 eps{q} original moments/tau versus frozen integrated HLOG z={text}',val,mp.mpf(str(target[f'I{q}'])),'3e-7' if q==2 else '3e-11',old_HLOG_reference_reused=True)
        qs1,qs2=q04(float(z),160),q04(float(z),320)
        iv=[-qs2[0]/8,-qs2[1]/8,-(qs2[2]+float(mp.zeta(2)-l*l/2)*qs2[0])/8]
        for q in range(3):
            check(f'4e04 Q{q} original two-variable mesh convergence z={text}',qs2[q],qs1[q],'5e-5',empirical_not_interval=True)
            expr=parse(read(H/'4e-04/result.json')['coefficients'][str(q)]['sympy'])
            check(f'4e04 physical eps{q} original Q insertion versus independent MN moments z={text}',iv[q],evaluate(expr,text),'8e-6',empirical_not_interval=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();out=a.output_dir.resolve()
    public_root=H.parents[1]
    if out.is_relative_to(public_root) and not out.is_relative_to(public_root/'replay-output'):
        raise SystemExit('Inside this repository choose replay-output/<new-name>; immutable science/docs/package trees cannot be output directories.')
    if out.exists() and any(out.iterdir()):raise SystemExit('Choose a fresh output directory')
    out.mkdir(parents=True,exist_ok=True);tick=time.monotonic();mp.mp.dps=50
    original_kernel_checks();moment_checks()
    record={'status':'passed','checks':CHECKS,'check_count':len(CHECKS),'seconds':time.monotonic()-tick,'mpmath_dps':mp.mp.dps,'new_integer_epsilon_original_pure_integrations':True,'new_Pi1_mixed_derivative_and_original_integrations':True,'new_TGV01_reduced_Beta_integral_controls':True,'new_4e04_subtracted_2D_coefficient_integrations':True,'new_original_MN_tau_controls':True,'new_full_TGV04_5D_CAS_integrations':False,'empirical_mesh_error_not_interval':True,'formal':'not_run'}
    (out/'integral-checks.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:record[k] for k in ['status','check_count','seconds','formal']}))

if __name__=='__main__':main()
