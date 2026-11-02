"""Assemble existing two-loop I, ordinary counterterms, strict poles and finite RG."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sys
import sympy as s
from sympy.parsing.mathematica import parse_mathematica

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DEST=ROOT/'benchmarks/renormalization/two-loop'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def data(p): return json.loads(Path(p).read_text(encoding='utf8'))
def encode(x): return (json.dumps(x,ensure_ascii=False,indent=2)+'\n').encode('utf8')
def obj(x): return {'sympy':str(x),'latex':s.latex(x)}
def immutable(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    blob=encode(x)
    if p.exists(): assert p.read_bytes()==blob, 'New revision required: '+str(p)
    else: p.write_bytes(blob)

e,L,rho,CA,CF,T,sigma=s.symbols('eps L rho CA CF T sigma')
M01,M11,M21,M03=s.symbols('M01 M11 M21 M03')
z=s.Symbol('z')
moments={(0,1):M01,(1,1):M11,(2,1):M21,(0,3):M03}
elementary={M01:-(rho**2+s.pi**2)/2,M03:-(rho**2+s.pi**2)**2/4}

def read_expression(rec):
    x=rec['expression']
    out=parse_mathematica(x) if rec['language'].startswith('Wolfram') else s.sympify(x)
    out=out.replace(s.Function('Zeta'),s.zeta).replace(s.Function('Gamma'),s.gamma)
    out=out.replace(lambda a:a.is_Function and str(a.func) in ['WilsonM','M'],
                    lambda a:moments[tuple(map(int,a.args[:2]))])
    return out.subs(s.log(z),rho)

def main():
    rows=[];dependencies={}
    def bind(rel):
        p=ROOT/rel;dependencies[rel]=sha(p);return data(p)
    def check(name,ok,detail=None):
        rows.append(dict(name=name,passed=bool(ok),**({'detail':detail} if detail is not None else {})))
        assert ok,name
    def zero(name,x): check(name,s.expand(x)==0)
    bridge=bind('benchmarks/renormalization/two-loop/assembly-bridge-r001.json')
    audit=bind('workflow/checks/rg-two-loop-audit/checks.json')
    for rel,hashval in audit['dependencies'].items():
        dependencies[rel]=sha(ROOT/rel);check('audit_dependency:'+rel,dependencies[rel]==hashval)
    check('audited_signs',all(row['physical_multiplier_of_existing_I']==(-1 if row['family']=='2e-Pi1' else 1)
                            for row in bridge['rows']))
    atoms={rho:ha.LOG}
    base='workflow/checks/single-opposite/runs/20260916T090538777022Z'
    for (k,q),symbol in moments.items():
        rel=f'{base}/M_{k}_{q}.txt';dependencies[rel]=sha(ROOT/rel);atoms[symbol]=ha.parse(ROOT/rel)
    check('M01_elementary',atoms[M01]==ha.from_sympy(elementary[M01],atoms))
    check('M03_elementary',atoms[M03]==ha.from_sympy(elementary[M03],atoms))
    cs={};coeffrows={}
    for row in bridge['rows']:
        oid=row['orbit_id'];result=bind(row['result']);coeff={}
        for field in ['package','manifest']:
            inp=result['input'];bind(inp[field]);check(oid+':'+field,dependencies[inp[field]]==inp[field+'_sha256'])
        check(oid+':no_extra_C4',result['normalization']['sym']==row['m'])
        for power,record in result['laurent']['coefficients'].items():
            chosen=record.get('MN_alternative',record.get('compact',record.get('M',record)))
            value=read_expression(chosen);coeff[int(power)]=value
            if 'hlog_terms' in record:
                vector={(tuple(x['word']),tuple(x['zeta_powers'])):Fraction(x['coefficient']) for x in record['hlog_terms']}
                check(oid+':Hlog_connection:'+power,ha.clean(vector)==ha.from_sympy(value,atoms))
            check(oid+':no_unknown_symbols:'+power,value.free_symbols<={rho,M01,M11,M21,M03})
        cs[oid]=sum(v*e**p for p,v in coeff.items())
        coeffrows[oid]={str(p):obj(v) for p,v in coeff.items()}
    pure=s.expand(sum(cs[k] for k in cs if k.startswith('4e-')))
    tgv=s.expand(sum(cs[k] for k in cs if k.startswith('3e-')))
    se=s.expand(sum(cs[k] for k in cs if k.startswith('2e-Pi1')))
    B=((10-6*e)*CA-8*(1-e)*T)/(3-2*e)
    Bseries=s.series(B,e,0,4).removeO()
    b0=(11*CA-4*T)/12
    combined=s.expand(CA*(pure+tgv)+sigma*Bseries*se)
    c={p:s.expand(combined.coeff(e,p).subs(elementary)) for p in range(-4,1)}
    zero('fourth_pole_cancels_for_either_SE_sign',c[-4])
    for p,x in c.items():
        check('M21_M03_cancel_in_combined:'+str(p),not(x.has(M21) or x.has(M03)))
    scale={n:s.expand(2**n*((L+rho/2)**n+(L-rho/2)**n)/s.factorial(n)) for n in range(5)}
    bare={p:s.expand(sum(c[j]*scale[p-j] for j in range(-4,p+1))) for p in range(-4,1)}
    one=bind('benchmarks/renormalization/one-loop/result-r001.json')
    onecs={int(p):s.expand(s.sympify(x['sympy']).subs({s.Symbol('L_s'):L+rho/2,s.Symbol('L_t'):L-rho/2}).subs(elementary))
           for p,x in one['laurent_unit_coefficients'].items()}
    check('one_loop_through_eps1',1 in onecs)
    ct={p:-b0*onecs.get(p+1,0) for p in range(-4,1)}
    post={p:s.expand(bare[p]+ct[p]) for p in range(-4,1)}
    # Exact finite convolution of the two complete scale powers, independently at graph level.
    for p in range(-4,1):
        direct=0
        for row in bridge['rows']:
            poly=cs[row['orbit_id']]
            color=CA if row['family']!='2e-Pi1' else -Bseries
            direct+=s.expand(sum(scale[j]*e**j for j in range(5))*color*poly).coeff(e,p)
        zero('per_graph_vs_group_convolution:'+str(p),direct.subs(elementary)-bare[p].subs(sigma,-1))
    Z1=sum(onecs[p]*e**p for p in [-2,-1]);F1=onecs[0]
    cases={}
    referenceK=CA*(s.Rational(67,36)-s.pi**2/12)-s.Rational(5,9)*T
    for sign in [-1,1]:
        case={p:s.expand(x.subs(sigma,sign)) for p,x in post.items()}
        Z2=sum(case[p]*e**p for p in range(-4,0));F2=case[0]
        Dz=s.expand(s.diff(Z2,L)-2*e*Z2-b0*Z1)
        Df=s.expand(s.diff(F2,L)-b0*F1)
        cuspZ=-s.diff(case[-1],L);hZ=-2*case[-1].subs(L,0)
        cuspF=-s.diff(Df,L)/2;hF=-Df.subs(L,0)
        failure={str(p):s.expand(Dz.coeff(e,p)) for p in range(-4,0)}
        rgshape=s.expand(Df+2*referenceK*L+hF)
        if sign==-1:
            for p,residual in failure.items():zero('Ddim_pole_cancellation:'+p,residual)
            zero('two_cusp_extractions_agree',cuspZ-cuspF)
            zero('two_h_extractions_agree',hZ-hF)
            zero('Ddim_finite_RG_connection',Dz-2*cuspZ*L-hZ)
            zero('finite_RG_connection',Df+2*cuspZ*L+hZ)
            zero('CA_cusp_reference',s.expand(cuspZ-referenceK).coeff(CA))
            zero('T_cusp_reference',s.expand(cuspZ-referenceK).coeff(T))
            check('noncusp_geometry_independent',hZ.free_symbols<={CA,T})
            zero('Z2_cubic_template',case[-3]-3*b0/2)
            zero('Z2_double_template',case[-2]-b0*L+cuspZ/2)
            zero('Z2_simple_template',case[-1]+cuspZ*L+hZ/2)
            for p,x in case.items():zero('crossing_post:'+str(p),x-x.subs(rho,-rho))
            f2=F2.subs(L,0)
            zero('finite_L_polynomial',F2-(-b0*L**3/3-cuspZ*L**2+(b0*(rho**2/4+s.pi**2/3)-hZ)*L+f2))
        else:
            check('negative_control_wrong_SE_locality',s.diff(case[-1],L,2)!=0)
            check('negative_control_wrong_SE_finite_RG',rgshape!=0)
            check('negative_control_wrong_SE_nf_alone',s.expand(s.diff(case[-1],L,2)).coeff(T)!=0)
        cases[str(sign)]={'post_counterterm_coefficients_over_CF':{str(p):obj(x) for p,x in case.items()},
            'Ddim_pole_residuals_over_CF':{p:obj(x) for p,x in failure.items()},
            'finite_RG_over_CF':obj(Df),'cusp_from_simple_pole_over_CF':obj(cuspZ),
            'cusp_from_finite_RG_over_CF':obj(cuspF),'h_from_simple_pole_over_CF':obj(hZ),
            'h_from_finite_RG_over_CF':obj(hF),'cusp_reference_residual_over_CF':obj(s.expand(cuspZ-referenceK)),
            'local_cusp_RG_consistent':sign==-1}
    best=cases['-1'];F2=post[0].subs(sigma,-1).expand()
    k2=s.sympify(best['cusp_from_simple_pole_over_CF']['sympy'])
    h2=s.sympify(best['h_from_simple_pole_over_CF']['sympy'])
    f2=F2.subs(L,0)
    d2=CA*(-s.Rational(607,162)-s.pi**4/48+11*s.zeta(3)/36+335*s.pi**2/432)+T*(s.Rational(82,81)-25*s.pi**2/108-s.zeta(3)/9)
    zero('compact_f2',f2-(k2*rho**2/4-b0*M11/2+d2))
    compact_bare={-4:0,-3:-b0/4,-2:-k2/4,-1:b0*(rho**2/4+s.pi**2/6)-h2/4,
                  0:k2*rho**2/4-b0*M11/2+d2/2-b0*s.zeta(3)/3}
    for power,value in compact_bare.items():zero('compact_bare_full_power:'+str(power),c[power].subs(sigma,-1)-value)
    zero('one_loop_eps1_used_in_CT',onecs[1]-(-L**3/3+L*(rho**2/4+s.pi**2/3)-M11/2-2*s.zeta(3)/3))
    # Destructive alternatives must fail without modifying any integral.
    missing=F2-ct[0]
    check('negative_control_missing_counterterm',s.expand(s.diff(missing,L)-b0*F1+2*k2*L+h2)!=0)
    check('negative_control_missing_beta_derivative',s.expand(s.diff(F2,L)+2*k2*L+h2)!=0)
    Bzero=s.expand(CA*(pure+tgv)-B.subs(e,0)*se)
    check('negative_control_four_dimensional_B1',s.expand(Bzero.coeff(e,-2)-c[-2].subs(sigma,-1))!=0)
    # Retain M11 as an explicit classical polylogarithm, not an unknown integral.
    M11exact=2*s.zeta(3)+(s.log(z)**2+s.pi**2)*s.log(1+z)-s.pi**2*s.log(z)/2+2*s.log(z)*s.polylog(2,-z)-2*s.polylog(3,-z)-s.log(z)**3/6
    body={'format':'wilson-two-loop-rg-result/1','revision':'r001',
        'scope':'fundamental massless QCD, Feynman gauge, non-self-intersecting null quadrilateral',
        'domain':'s,t<0 as analytic invariants; fixed z=s/t>0; Feynman branches retained before continuation',
        'definitions':{'lambda':'alpha_s_MSbar/pi','T':'TF*nf','TF':'1/2',
            'L':'(L_s+L_t)/2','rho':'L_s-L_t=log(z)','mu_x_squared':'exp(2*EulerGamma)*mu_MSbar_squared/4',
            'X_s':'mu_x_squared*(-s+i0)','X_t':'mu_x_squared*(-t+i0)',
            'beta0':str(4*b0),'b0':str(b0),'M11':obj(M11exact)},
        'old_I_unchanged':True,'physical_multiplier_of_old_Pi1_I':-1,
        'group_coefficients':{name:{str(p):obj(poly.coeff(e,p)) for p in range(-4,1)} for name,poly in [('pure',pure),('TGV',tgv),('Pi1_without_B1',se)]},
        'B1_expanded_through_eps3':obj(Bseries),
        'bare_w2_full_scale_form':{'prefactor':'CF*(X_s^(2eps)+X_t^(2eps))',
            'coefficient_series':{str(p):obj(x.subs(sigma,-1)) for p,x in c.items()},'remainder':'O(eps), fixed geometry'},
        'bare_w2_coefficients_over_CF':{str(p):obj(x.subs(sigma,-1)) for p,x in bare.items()},
        'ordinary_counterterm_over_CF':{str(p):obj(s.expand(x)) for p,x in ct.items()},
        'sign_comparison':cases,
        'F1_over_CF':obj(F1),'F2_over_CF':obj(F2),'f2_over_CF':obj(f2),
        'log_Z4_two_loop_coefficient':obj(CF*sum(post[p].subs(sigma,-1)*e**p for p in range(-4,0))),
        'k2':obj(CF*k2),'h2':obj(CF*h2),'gamma4_two_loop':obj(2*CF*h2),
        'W_R_two_loop_coefficient':obj(CF*F2+CF**2*F1**2/2),
        'status':'computed through finite term; two routes of RG/cusp agree',
        'formal':'not_run','MCP':'no Wilson proof MCP tool available in this session',
        'limits':['general complex z continuation not certified','old upstream Pi1 convention not changed',
                  'auxiliary symbolic checks are not a QFT formal proof'],
        'dependencies':dict(sorted(dependencies.items()))}
    for rel,h in dependencies.items():check('unchanged:'+rel,sha(ROOT/rel)==h)
    immutable(DEST/'result-r001.json',body)
    (HERE/'exact-checks.json').write_bytes(encode({'checks':rows,'passed':len(rows),
        'negative_controls':sum(x['name'].startswith('negative_control') for x in rows),
        'result_sha256':sha(DEST/'result-r001.json'),'code_sha256':sha(__file__),'formal':'not_run'}))
    print(json.dumps({'checks':len(rows),'k2':str(CF*k2),'gamma4_two_loop':str(2*CF*h2)},ensure_ascii=False))

if __name__=='__main__':main()
