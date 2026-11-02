"""Independent moment quadrature and direct evaluation of the original WRen.m."""
import json
import mpmath as mp
import sympy as s
from calculate import ROOT,HERE,DEST,data,encode,sha
from reference import fetch,raw_W

def moment11(z):
    split=z/(1+z)
    def f(b):
        if b==0 or b==1:return mp.mpf(0)
        R=mp.log(b)+mp.log1p(-b)
        if b==split:return -R/(split*(1-split))
        return R*(mp.log(z)+mp.log1p(-b)-mp.log(b))/(b-split)
    return mp.quadts(f,[0,split,1])

def main():
    fetch();result=data(DEST/'result-r001.json');mp.mp.dps=70;rows=[]
    def close(name,a,b,tol='1e-55'):
        error=abs(a-b)/max(1,abs(a),abs(b));ok=error<mp.mpf(tol)
        rows.append({'name':name,'passed':bool(ok),'relative_error':mp.nstr(error,10),'tolerance':tol})
        assert ok,rows[-1]
    a,ss,tt,mu,CA,CF,nf,TR=s.symbols('as4pi s t mu CA CF nf TR')
    # No symbolic polylog identities from reference.py are applied on this route.
    raw=s.expand(raw_W());p1=raw.coeff(a,1)/4;p2=raw.coeff(a,2)/16
    paper=s.lambdify((ss,tt,mu,CA,CF,nf,TR),p2-p1**2/2,'mpmath')
    L,rho,M11,T=s.symbols('L rho M11 T')
    F1=s.lambdify((L,rho),s.sympify(result['F1_over_CF']['sympy']),'mpmath')
    F2=s.lambdify((L,rho,M11,CA,T),s.sympify(result['F2_over_CF']['sympy']),'mpmath')
    z=s.Symbol('z');m11=s.lambdify(z,s.sympify(result['definitions']['M11']['sympy']),'mpmath')
    k2=s.lambdify((CA,CF,T),s.sympify(result['k2']['sympy']),'mpmath')
    h2=s.lambdify((CA,CF,T),s.sympify(result['h2']['sympy']),'mpmath')
    for zs in ['.2','.7','1','1.4','5']:
        zz=mp.mpf(zs);r=mp.log(zz);m=moment11(zz)
        close('M11_integral_vs_classical_z'+zs,m,m11(zz))
        close('M11_crossing_z'+zs,m,m11(1/zz))
        for ls in ['-.4','.3']:
            ell=mp.mpf(ls);mup=mp.exp((ell-r/2)/2)
            for label,ca,cf,nn in [('CA',mp.mpf(1),mp.mpf(1),mp.mpf(0)),
                                   ('T',mp.mpf(0),mp.mpf(1),mp.mpf(2)),
                                   ('SU3_nf5',mp.mpf(3),mp.mpf(4)/3,mp.mpf(5))]:
                flavor=nn/2
                ours=cf*F2(ell,r,m,ca,flavor)
                theirs=paper(mp.sqrt(zz),mp.mpf(1),mup,ca,cf,nn,mp.mpf(1)/2)
                close('original_WRen_'+label+'_z'+zs+'_L'+ls,theirs,ours)
            ca,cf,flavor=mp.mpf(3),mp.mpf(4)/3,mp.mpf(5)/2
            b0=(11*ca-4*flavor)/12
            step=mp.mpf('1e-20')
            derivative=cf*(F2(ell+step,r,m,ca,flavor)-F2(ell-step,r,m,ca,flavor))/(2*step)
            value=derivative-b0*cf*F1(ell,r)
            close('finite_RG_scale_difference_z'+zs+'_L'+ls,value,-2*k2(ca,cf,flavor)*ell-h2(ca,cf,flavor),'1e-37')
    report={'working_dps':70,'checks':rows,'passed':len(rows),'code_sha256':sha(__file__),
        'result_sha256':sha(DEST/'result-r001.json'),'reference_source_sha256':sha(HERE/'reference-source.json'),
        'limits':'Empirical precision controls; no certified intervals or full original 12-integral re-integration.'}
    (HERE/'numeric-checks.json').write_bytes(encode(report))
    print(json.dumps({'numeric_checks':len(rows),'maximum_relative_error':max(float(r['relative_error']) for r in rows)}))

if __name__=='__main__':main()
