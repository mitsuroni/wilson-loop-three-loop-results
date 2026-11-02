"""Current-only low-loop assembly from fourteen canonical individual I tables."""
from pathlib import Path
import argparse,hashlib,json,sys,time
sys.dont_write_bytecode=True
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[2]
eps,L,rho,CA,T=s.symbols('eps L rho CA T');M=lambda k,j:s.Symbol(f'M_{k}_{j}')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return s.sympify(v['sympy'] if isinstance(v,dict) else v)
def clean(v):
    v=s.sympify(v).subs({M(0,3):-M(0,1)**2,s.Symbol('N_0_2'):M(0,1)**2/3+4*s.zeta(2)*M(0,1)},simultaneous=True)
    return s.expand(v.subs(M(0,1),-(rho**2+s.pi**2)/2))
def encode(v):return {str(p):str(clean(x)) for p,x in v.items()}
def scaled(values,n,lo,hi):
    return {p:clean(sum(v*n**(p-q)*((L+rho/2)**(p-q)+(L-rho/2)**(p-q))/s.factorial(p-q) for q,v in values.items() if q<=p)) for p in range(lo,hi+1)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();out=a.output_dir.resolve()
    if out.exists() or (out.is_relative_to(R) and not any(out.is_relative_to(x) for x in [H,R/'.wilson',R/'replay-output'])):raise ValueError('Use a fresh output outside the scientific source trees.')
    tick=time.monotonic();inputs={};checks=[]
    def read(rel):
        p=R/rel;inputs[rel]=sha(p);return json.loads(p.read_text(encoding='utf-8-sig'))
    def check(name,v,negative=False):
        v=clean(v);assert v!=0 if negative else v==0,(name,v)
        checks.append({'name':name,'negative_control':negative,'residual':str(v)})
    catalog=read('benchmarks/current/catalog.json');one={};two={};rows=[]
    for row in catalog['entries']:
        n=row['loop_order']
        if n>2:continue
        ref=row['current_coefficients'];data=read(ref['path']);assert inputs[ref['path']]==ref['sha256']
        assert data['physical_convention']=='wilson-qcd-action-path-fourier/1.0' and data['normalization']['assembly_sign']==1
        assert not any(k in data for k in ['legacy_to_canonical_multiplier','family_sign','legacy_sign'])
        cs={int(p):parse(v) for p,v in data['coefficients'].items()}
        if n==1:one[row['orbit_id']]=cs;weight=s.S.One
        else:
            weight=CA*parse(data['color_weights_over_CF']['CA'])+T*parse(data['color_weights_over_CF']['T'])
            series=s.series(weight,eps,0,7).removeO().expand()
            two[row['orbit_id']]={p:clean(sum(v*series.coeff(eps,p-q) for q,v in cs.items() if q<=p)) for p in range(-4,3)}
        rows.append({'orbit_id':row['orbit_id'],'loop_order':n,'assembly_sign':1,'physical_weight':str(weight),'result':ref['path'],'sha256':ref['sha256']})
    assert len(one)==2 and len(two)==12
    U1={p:clean(sum(x[p] for x in one.values())) for p in range(-2,4)}
    U2={p:clean(sum(x[p] for x in two.values())) for p in range(-4,3)}
    B1=scaled(U1,1,-2,3);B2=scaled(U2,2,-4,2);b0=(11*CA-4*T)/12
    CT={p:clean(-b0*B1.get(p+1,0)) for p in range(-4,3)};P2={p:clean(B2[p]+CT[p]) for p in range(-4,3)}
    F1,F2=B1[0],P2[0];k2=-s.diff(P2[-1],L);h2=-2*P2[-1].subs(L,0)
    check('finite one-loop RG',s.diff(F1,L)+2*L)
    check('finite two-loop RG',s.diff(F2,L)-b0*F1+2*k2*L+h2)
    check('cusp',k2-(CA*(s.Rational(67,36)-s.pi**2/12)-5*T/9))
    check('noncusp',h2-(CA*(s.Rational(101,27)-7*s.zeta(3)/2-11*s.pi**2/72)+T*(-s.Rational(28,27)+s.pi**2/18)))
    for p,v in {-4:0,-3:3*b0/2,-2:b0*L-k2/2,-1:-k2*L-h2/2}.items():check('pole'+str(p),P2[p]-v)
    result={'format':'wilson-current-low-loop-assembly/1','physical_convention':'wilson-qcd-action-path-fourier/1.0','one_loop':{'U':encode(U1),'bare':encode(B1),'post_counterterm':encode(B1),'F':str(F1),'logZ':encode({p:B1[p] for p in [-2,-1]})},'two_loop':{'U':encode(U2),'bare':encode(B2),'counterterm':encode(CT),'post_counterterm':encode(P2),'F':str(F2),'logZ':encode({p:P2[p] for p in range(-4,0)}),'k_over_CF':str(k2),'h_over_CF':str(h2)},'source_rows':rows,'formal':'not_run'}
    out.mkdir(parents=True);(out/'constructed-before-holdout.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    holdout=read('benchmarks/current/low-loops/coefficients.json')
    for loop in ['one_loop','two_loop']:
        for stage in ['U','bare','post_counterterm']:
            for power,value in result[loop][stage].items():check(loop+' '+stage+' '+power,parse(value)-parse(holdout[loop][stage][power]))
        check(loop+' finite',parse(result[loop]['F'])-parse(holdout[loop]['F']))
    check('negative duplicate current self-energy sign',two['2e-Pi1-01'][-3]*2,True)
    check('negative missing finite CT',CT[0],True)
    for rel,digest in inputs.items():assert sha(R/rel)==digest,rel
    result.update(status='passed',checks=checks,inputs=inputs,seconds=time.monotonic()-tick,legacy_family_sign_applied=False)
    (out/'coefficients.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','integrals':14,'checks':len(checks),'negative_controls':2,'seconds':result['seconds']}))
if __name__=='__main__':main()
