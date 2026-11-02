"""Read low-loop finite functions using original M11 quadrature and polylogs."""
from pathlib import Path
import argparse, hashlib, json, sys
sys.dont_write_bytecode=True
import mpmath as mp
import sympy as s

HERE=Path(__file__).resolve().parent
L,rho,M11,CA,T,CF,B,Li2,Li3=s.symbols('L rho M_1_1 CA T CF B Li2 Li3')

def raw_m11(z):
    r=z/(1+z)
    def f(b):
        if b==r:return mp.log(r*(1-r))*(-1/(r*(1-r)))
        lb=mp.log1p((r-b)/(b*(1-r)))
        return mp.log(b*(1-b))*lb/(b-r)
    return mp.quad(f,[0,r/2,r,(1+r)/2,1])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--coefficients',type=Path,default=HERE/'coefficients.json');ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    out=a.output_dir.resolve()
    root=HERE.parents[1]
    if any(out.is_relative_to((root/name).resolve()) for name in ['data','packages','results','integrals']):
        raise SystemExit('Output must be outside the frozen source/data/result/article trees.')
    if out.exists():raise SystemExit('Use a new output directory.')
    fixture=HERE/'reference/fixture.json'
    hashes={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in [a.coefficients,fixture]}
    data=json.loads(a.coefficients.read_text(encoding='utf-8-sig'));ref=json.loads(fixture.read_text(encoding='utf-8-sig'))
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==digest for p,digest in hashes.items()),'Input changed during reading'
    f1=s.lambdify((L,rho),s.sympify(data['one_loop']['F']),'mpmath')
    f2=s.lambdify((L,rho,M11,CA,T),s.sympify(data['two_loop']['F']),'mpmath')
    r2=s.lambdify((L,rho,B,Li2,Li3,CA,T),s.sympify(ref['F2_over_CF']),'mpmath')
    rows=[];errors=[];values={}
    for dps in [40,60]:
        with mp.workdps(dps):
            for frac in [(2,3),(1,1),(3,2)]:
                z=mp.mpf(frac[0])/frac[1];ell=mp.log(z);m=raw_m11(z)
                closed=2*mp.zeta(3)+(ell**2+mp.pi**2)*mp.log1p(z)-mp.pi**2*ell/2+2*ell*mp.polylog(2,-z)-2*mp.polylog(3,-z)-ell**3/6
                err=abs(m-closed)/(1+abs(closed));assert err<mp.mpf(10)**(-dps+7),('M11',z,err)
                errors.append(mp.nstr(err,12))
                for ltext in ['0','0.37']:
                    lv=mp.mpf(ltext);one=f1(lv,ell)
                    result={'dps':dps,'z':f'{frac[0]}/{frac[1]}','L':ltext,'M11_original_quadrature':mp.nstr(m,dps),'F1_over_CF':mp.nstr(one,dps)}
                    for label,ca,t in [('CA',1,0),('T',0,1)]:
                        own=f2(lv,ell,m,ca,t);external=r2(lv,ell,mp.log1p(z),mp.polylog(2,-z),mp.polylog(3,-z),ca,t)
                        error=abs(own-external)/(1+abs(external));assert error<mp.mpf(10)**(-dps+6),('finite reference',label,error)
                        errors.append(mp.nstr(error,12));result['F2_'+label+'_over_CF']=mp.nstr(own,dps)
                        values[(dps,frac,ltext,label)]=own
                    rows.append(result)
            for ltext in ['0','0.37']:
                for label in ['CA','T']:
                    x=values[(dps,(2,3),ltext,label)];y=values[(dps,(3,2),ltext,label)]
                    assert abs(x-y)/(1+abs(x))<mp.mpf(10)**(-dps+6),('crossing',dps,ltext,label)
    changes=[]
    with mp.workdps(70):
        for key,low in values.items():
            dps,frac,ltext,label=key
            if dps==40:
                difference=abs(low-values[(60,frac,ltext,label)])/(1+abs(values[(60,frac,ltext,label)]))
                assert difference<mp.mpf('1e-33')
                changes.append(difference)
        largest=mp.nstr(max(changes),15)
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==digest for p,digest in hashes.items()),'Input changed during calculation'
    report={'status':'passed','samples':rows,'precisions':[40,60],'maximum_scaled_40_to_60_change':largest,'M11_and_external_finite_scaled_errors':errors,'source_coefficients_sha256':hashes[a.coefficients],'reference_fixture_sha256':hashes[fixture],'scope':'New original one-dimensional M11 definition, stable function reading, finite-reference and inversion samples. Not new original two-loop multivariable integration or strict interval proof.','formal':'not_run'}
    out.mkdir(parents=True);(out/'numeric-samples.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','samples':len(rows),'maximum_scaled_change':largest,'formal':'not_run'}))

if __name__=='__main__':main()
