"""Complete current four-color scales/CT/RG from directly assembled native sums.

Old totals are read only as exact function-representation targets and holdouts,
after the independently constructed native expression has been saved.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[2]
spec=importlib.util.spec_from_file_location('current_combine_maths',H/'consumer-agent/assemble_matter.py')
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
k=m.k;L,rho,CA,CF,T=m.L,m.rho,m.CA,m.CF,m.T
COLORS={'AA':CA**2,'Af':CA*T,'Ff':CF*T,'ff':T**2};POWERS=range(-6,1)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def enc(v):return {str(p):str(s.expand(x)) for p,x in v.items()}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True)
    for name in ['low','matter','aa']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args();out=a.output_dir.resolve()
    if out.exists() or (out.is_relative_to(R) and not any(out.is_relative_to(x) for x in [H,R/'.wilson',R/'replay-output'])):raise ValueError('Use a new output outside scientific source trees')
    out.mkdir(parents=True);tick=time.monotonic();inputs={};checks=[]
    def read(p):
        p=Path(p).resolve();inputs[str(p)]=sha(p);return json.loads(p.read_text(encoding='utf-8-sig'))
    dictionary=read(R/'workflow/checks/four-leg-bootstrap-20261101/basis/dictionary.json')
    k.MOMENTS={s.Symbol(x['name']):k.vec_expr(x['hlog_terms']) for x in dictionary['atoms']};k.MOMENTS[s.Symbol('ell')]=rho
    def parse(v):return k.math.parse(v)
    def check(name,value,negative=False):
        poly=s.Poly(s.expand(value),CA,CF,T)
        v=s.expand(sum(m.normL(c)*CA**a*CF**b*T**t for (a,b,t),c in poly.terms()))
        assert v!=0 if negative else v==0,(name,str(v)[:800])
        checks.append({'name':name,'residual':str(v),'negative_control':negative})
    low,matter,aa=read(a.low),read(a.matter),read(a.aa)
    assert low['status']==matter['status']=='passed' and aa['status']=='passed_fresh_current_AA_reassembly'
    assert matter['legacy_family_sign_applied'] is False and aa['construction_did_not_consume_old_aggregate_or_AA_total']
    assert aa['source_orbit_coverage']=={'tree':83,'insertions':31,'fourleg':35,'total':149}
    w1={int(p):parse(v) for p,v in low['one_loop']['bare'].items()};w2={int(p):parse(v) for p,v in low['two_loop']['bare'].items()}
    b0,b1=(11*CA-4*T)/12,(34*CA**2-20*CA*T-12*CF*T)/48
    CT={p:s.expand(-2*b0*w2.get(p+1,0)+b0*b0*w1.get(p+2,0)-b1*w1.get(p+1,0)/2) for p in POWERS}
    native={}
    for c,color in COLORS.items():
        U={int(p):parse(v) for p,v in (aa['U_constant_MN'] if c=='AA' else matter['channels'][c]['U']).items()}
        bare={p:s.expand(sum(v*3**(p-q)*((L+rho/2)**(p-q)+(L-rho/2)**(p-q))/s.factorial(p-q) for q,v in U.items() if q<=p)) for p in POWERS}
        ct={p:s.expand(CT[p]).coeff(CA,2 if c=='AA' else 1 if c=='Af' else 0).coeff(CF,1 if c=='Ff' else 0).coeff(T,2 if c=='ff' else 1 if c in ['Af','Ff'] else 0) for p in POWERS}
        post={p:s.expand(bare[p]+ct[p]) for p in POWERS};F=post[0]
        cusp=-s.Rational(3,2)*s.diff(post[-1],L);noncusp=-3*post[-1].subs(L,0)
        native[c]={'U':enc(U),'bare':enc(bare),'counterterm':enc(ct),'post_counterterm':enc(post),'F':str(F),'k3':str(cusp),'h3':str(noncusp)}
        check(c+' finite RG',s.diff(F,L)-s.expand(2*b0*parse(low['two_loop']['F'])+b1*parse(low['one_loop']['F'])).coeff(CA,2 if c=='AA' else 1 if c=='Af' else 0).coeff(CF,1 if c=='Ff' else 0).coeff(T,2 if c=='ff' else 1 if c in ['Af','Ff'] else 0)+2*cusp*L+noncusp)
        if c!='AA':
            for stage in ['bare','counterterm','post_counterterm']:
                for p in POWERS:check(c+' independent '+stage+str(p),parse(native[c][stage][str(p)])-parse(matter['channels'][c][stage][str(p)]))
    dump(out/'constructed-before-holdout.json',{'native_channels':native,'old_total_used_to_construct':False,'inputs_before_holdout':inputs})
    target=read(R/'benchmarks/renormalization/three-loop/first-complete-20261102/coefficients.json');functions={}
    for c,rec in native.items():
        target_c=target['channels'][c];fn={}
        for stage in ['U','bare','counterterm','post_counterterm']:
            fn[stage]={}
            for p,v in rec[stage].items():
                candidate=parse(target_c[stage][p]);check(c+' exact function representative '+stage+' '+p,parse(v)-candidate)
                fn[stage][p]=str(s.expand(candidate))
        for field in ['k3','h3']:
            candidate=parse(target_c[field]);check(c+' extracted '+field+' reference',parse(rec[field])-candidate);fn[field]=str(s.expand(candidate))
        fn['F']=fn['post_counterterm']['0'];fn['f_at_L0']=str(s.expand(parse(fn['F']).subs(L,0)))
        fn['logZ']={p:v for p,v in fn['post_counterterm'].items() if int(p)<0};functions[c]=fn
    for c in ['AA','Af','Ff','ff']:
        for p in [-6,-5]:check(c+' highest pole cancellation '+str(p),parse(functions[c]['post_counterterm'][str(p)]))
    check('negative duplicated native sign',2*parse(functions['Af']['U']['-4']),True)
    check('negative missing lower-order finite term',2*b0*parse(low['two_loop']['F'])+b1*parse(low['one_loop']['F']),True)
    for path,digest in inputs.items():assert sha(Path(path))==digest,path
    common={'physical_convention':'wilson-qcd-action-path-fourier/1.0','normalization':{'quantity':'coefficient in connected log W; external CF removed','lambda':'alpha_s_MSbar/pi','T':'TF*nf','L':'(Ls+Lt)/2','rho':'log(z)','family_sign_conversion_allowed':False,'Euler_and_orbit_restoration_count_at_assembly':0},'formal':'not_run'}
    for n,name in [(1,'one_loop'),(2,'two_loop')]:
        dump(out/f'{n}loop.json',{'format':'wilson-current-loop-result/1','loop_order':n,**common,'coefficients':low[name]})
    dump(out/'3loop.json',{'format':'wilson-current-loop-result/1','loop_order':3,**common,'color_basis_over_CF':{c:str(v) for c,v in COLORS.items()},'channels':functions,'native_common_z_HLOG_channels':native,'source_orbit_count':149,'fourleg_joint_objects':9,'function_representation_scope':'AA U independently solved in a predefined MN basis; matter/native expressions converted to accepted function representatives only after full exact HLOG equality. Old totals never construct the native sum.','validation_scope':'Direct current individual/joint coefficient assembly, scales/CT/RG and reference equality; original integrations and exact dictionaries reused at their prior research scope.'})
    report={'status':'passed','checks':len(checks),'negative_controls':sum(x['negative_control'] for x in checks),'checks_detail':checks,'inputs':inputs,'seconds':time.monotonic()-tick,'source_count':149,'extra_family_signs':False,'formal':'not_run'}
    dump(out/'checks.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['inputs','checks_detail']}))
if __name__=='__main__':main()
