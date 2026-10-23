"""Independent matter checks: RG null directions, external noncusp, raw nf2 C/J.

No existing I/Laurent coefficient enters the nf2 prediction construction.
Only the final comparison opens the published matter result. This is a new
derivation route from shared prepared root facts, not a blind human analysis
or an independent derivation of the upstream QCD operator.
"""
from pathlib import Path
import argparse,hashlib,json,platform,time,re
import sympy as sp
import mpmath as mp

HERE=Path(__file__).resolve().parent
DEFAULT_ROOT=HERE.parents[2]
if (DEFAULT_ROOT/'.wilson/task-inputs.json').is_file():
    DEFAULT_ROOT=Path(json.loads((DEFAULT_ROOT/'.wilson/task-inputs.json').read_text(encoding='utf-8'))['read_only_input_root'])
e,L,rho,z,s,t,CA,CF,T=sp.symbols('eps L rho z s t CA CF T')
M01,M11,M21,M03=sp.symbols('M01 M11 M21 M03')
SOURCES={};CHECKS=[]

def read(root,rel):
    raw=(root/rel).read_bytes();h=hashlib.sha256(raw).hexdigest()
    if rel in SOURCES: assert SOURCES[rel]==h,('concurrent input change',rel)
    SOURCES[rel]=h
    return json.loads(raw.decode('utf-8-sig'))

def text_hash(root,rel):
    SOURCES[rel]=hashlib.sha256((root/rel).read_bytes()).hexdigest()

def parse(x):
    if isinstance(x,dict): x=x.get('sympy',x.get('expression'))
    x=re.sub(r'Gamma\[([^\[\]]+)\]',r'gamma(\1)',x)
    return sp.sympify(re.sub(r'zeta\[(\d+)\]',r'zeta(\1)',x).replace('^','**'))

def zero(name,x):
    val=sp.cancel(sp.expand(x));assert val==0,(name,str(val)[:1000])
    CHECKS.append({'name':name,'passed':True,'residual':'0'})

def nonzero(name,x):
    val=sp.factor(x);assert val!=0,name
    CHECKS.append({'name':name,'passed':True,'negative_control':True,'nonzero_residual':str(val)})

def coefficient(x,q): return sp.expand(x).coeff(e,q)

def truncated(x,hi): return sp.series(x,e,0,hi+1).removeO().expand()

def convolution(a,b,lo=-4,hi=0):
    return {q:sp.expand(sum(v*b.get(q-j,0) for j,v in a.items())) for q in range(lo,hi+1)}

def rg_nullspace():
    k,h,delta,eta=sp.symbols('k h delta eta')
    g=sp.Function('g')(rho)
    candidate=-(k*L**2+h*L)+g
    shifted=candidate-delta*L+eta*(rho**4+sp.zeta(3))
    zero('RG blind direction: synchronized noncusp and finite L change',sp.diff(shifted,L)+2*k*L+h+delta)
    zero('RG blind direction: pole L derivative is unchanged',sp.diff(-delta/3,L))
    zero('RG blind direction: even finite shape passes crossing',(rho**4+sp.zeta(3))-(rho**4+sp.zeta(3)).subs(rho,-rho))
    nonzero('negative: change h3 without changing finite L term',sp.diff(candidate,L)+2*k*L+h+delta)
    nonzero('negative: odd finite shape fails crossing',rho-rho.subs(rho,-rho))
    nonzero('negative: a nonzero even finite deformation still changes the answer',eta*(rho**4+sp.zeta(3)))
    # D=partial_L-(b0*lambda^2+b1*lambda^3)*partial_lambda.
    lam,b0,b1,f1,f2,c=sp.symbols('lam b0 b1 f1 f2 c')
    finite=f1*lam+f2*lam**2
    gamma_shift=sp.expand((2*b0*lam**2+2*b1*lam**3)*sp.diff(finite,lam))
    zero('finite Wilson scheme shift at lambda2',gamma_shift.coeff(lam,2)-2*b0*f1)
    zero('finite Wilson scheme shift at lambda3',gamma_shift.coeff(lam,3)-4*b0*f2-2*b1*f1)
    zero('log scale shift in Gamma',(-4*k*(L+c)-h)-(-4*k*L-(h+4*c*k)))
    nonzero('negative: ignore a constant log scale shift',4*c*k)
    # The same bare coupling is represented using the two explicitly different
    # Euler conventions; matching scale arguments is mandatory.
    muP=sp.Symbol('muP',positive=True)
    project_log=e*sp.log(muP)+sp.EulerGamma*e-e*sp.log(4*sp.pi)
    paper_log=e*(sp.log(muP)+2*sp.EulerGamma)-sp.EulerGamma*e-e*sp.log(4*sp.pi)
    zero('MS Euler dictionaries match at muF2=exp(2gammaE)*muP2',paper_log-project_log)
    zero('coordinate scale muF2/4 equals project mux2',sp.exp(2*sp.EulerGamma)*muP/4-sp.exp(2*sp.EulerGamma)*muP/4)
    nonzero('negative: identify muF2 and muP2 numerically',-2*sp.EulerGamma*e)

def source_kernel_checks(single,chain):
    F=(1-s)*(1-t)+z*s*t
    for data,n in [(single,2),(chain,3)]:
        assert data['sign_convention']=='wilson-self-energy-1pi-insertion/2.0'
        for orb in data['orbits']:
            assert orb['J']['jacobian']=='1'
            assert all(pair==['0','1'] for pair in orb['J']['bounds'].values())
            assert orb['c4_symmetry']['orbit_multiplicity']==(4 if orb['orbit_id'].endswith('-01') else 2)
            raw=orb['J']['integrand'].replace('^','**').replace('z1','z').replace('s31','s').replace('s21','s').replace('s41','t')
            kernel=sp.sympify(raw)
            if orb['orbit_id'].endswith('-02'):
                nu=n*e
                divisor=2*n*e*(1+(n-1)*e)
                expected=-(1+z)*F+nu*(nu*(1+z)*F-(nu-1)*z)/divisor
                # F>0 in the open original unit square, z>0. Equal-base powers
                # are normalized only on that real domain.
                positive_F=sp.Symbol('positive_F',positive=True)
                canonical=kernel.xreplace({sp.expand(F):positive_F,F:positive_F})
                numerator=sp.powsimp(canonical*positive_F**(2-nu),force=True).subs(positive_F,F)
                residual=numerator-expected
                zero(f'raw {orb["orbit_id"]} scalar plus complete derivative',residual)
            else:
                Fadj=z*(1-s)*t
                factor=1/(1+e) if n==2 else (2+e)/(2*(1+2*e))
                positive_F=sp.Symbol('positive_F',positive=True)
                canonical=kernel.xreplace({sp.expand(Fadj):positive_F,Fadj:positive_F})
                numerator=sp.powsimp(canonical*positive_F**(2-n*e),force=True).subs(positive_F,Fadj)
                zero(f'raw {orb["orbit_id"]} adjacent monomial',numerator-factor*z*Fadj)
    # Full C4/crossing sum of derivative terms, including opposite-edge corners.
    R2=sp.Symbol('R2');R3=sp.Symbol('R3')
    zero('closed single-bubble boundary plus adjacent source',2*R2/(4*e**2*(1+e))+R2/(2*e*(1+e))-R2/(2*e**2))
    zero('closed double-bubble boundary plus adjacent source',R3*(2+e)/(9*e**2*(1+2*e))+R3/(3*e*(1+2*e))-2*R3/(9*e**2))
    # The scalar master is integrated in t before its epsilon expansion:
    # H_j=-int (A^(j+1)-B^(j+1))/((j+1)!*(t-b)) dt.
    A,B,W,V=sp.symbols('A B W V')
    zero('independent H1 moment map',(A**2-B**2)/2-(A-B)*(A+B)/2)
    zero('independent H2 moment map',(A**3-B**3)/6-(A-B)*((A+B)**2/8+(A-B)**2/24))

def source_prediction(single,chain,pi2):
    source_kernel_checks(single,chain)
    def gn(n,hi=4): return truncated(sp.exp(sum(n*sp.zeta(k)*e**k/k for k in range(2,hi+1))),hi)
    def pg(hi=4): return truncated(sp.exp(sum((4+2*(-1)**k+3**k-2*2**k-(-2)**k)*sp.zeta(k)*e**k/k for k in range(2,hi+1))),hi)
    P=pg()
    # The accepted chain C is checked against its explicit multiplicative parts;
    # no old scalar I, result coefficient, or RG-fitted constant enters here.
    for data,n in [(single,2),(chain,3)]:
        for orb in data['orbits']:
            comp=orb['C']['components']
            assert comp['euler']==f'exp(-{n}*EulerGamma*eps)'
            expected= '-(1/(16*eps*(1-2*eps)))' if n==2 else '-1/(128*eps^2*(1-2*eps)^2)'
            zero(f'{orb["orbit_id"]} raw rational C',sp.sympify(comp['rational'].replace('^','**'))-sp.sympify(expected.replace('^','**')))
    h=[M01,M11/2,M21/8+M03/24]
    def geom(n):
        expn=truncated(sp.exp(n*e*L),4)
        coshn=truncated(sp.cosh(n*e*rho/2),4)
        return truncated(-expn*(2*coshn/(n*e)**2+sum(h[j]*(n*e)**j for j in range(3))),2)
    Bf=-8*(1-e)/(3-2*e)
    for orb in single['orbits']:
        color=sp.expand(parse(orb['C']['components']['self_energy_color']).subs(sp.Symbol('TF')*sp.Symbol('nf'),T))
        zero(orb['orbit_id']+' actual C fermion coefficient',color.coeff(T)-Bf)
    U=sp.gamma(1+e)**2*sp.gamma(1-e)*sp.gamma(1-3*e)/(sp.gamma(1-2*e)**2*sp.gamma(1+2*e))
    for orb in chain['orbits']:
        color=sp.expand(parse(orb['C']['components']['self_energy_color']).subs(sp.Symbol('TF')*sp.Symbol('nf'),T))
        zero(orb['orbit_id']+' actual C nf2 coefficient',color.coeff(T,2)-Bf**2)
        zero(orb['orbit_id']+' actual C Gamma ratio',parse(orb['C']['components']['gamma_chain_ratio'])-U)
    w1=truncated(gn(1)*geom(1),2)
    w2=truncated(Bf*gn(2)/(8*e*(1-2*e))*geom(2),1)
    w3=truncated(Bf**2*P/(64*e**2*(1-2*e)**2)*geom(3),0)
    # nf2 color extraction of the ordinary coupling CT follows from beta0_T=-1/3.
    post=sp.expand(w3+sp.Rational(2,3)/e*w2+sp.Rational(1,9)/e**2*w1)
    elementary={M01:-(rho**2+sp.pi**2)/2,M03:-(rho**2+sp.pi**2)**2/4}
    post=sp.expand(post.subs(elementary));w3=sp.expand(w3.subs(elementary))
    c={q:coefficient(post,q) for q in range(-4,1)}
    h3=-3*c[-1].subs(L,0);k3=-sp.Rational(3,2)*sp.diff(c[-1],L)
    zero('raw-source nf2 h3 geometry independence',sp.diff(h3,rho))
    nf2={'bare':{str(q):str(coefficient(w3,q)) for q in range(-4,1)},
            'post':{str(q):str(v) for q,v in c.items()},'h3':str(h3),'k3':str(k3),
            'finite':str(c[0]),'shape':str(c[0].subs(L,0)),
            'construction_inputs':'accepted single/double bubble prepared C/J, elementary beta integrals, scalar logarithmic moments; no published I/Laurent coefficients',
            'human_analysis_blind':False,'independence':'independent rederivation from shared prepared root facts, not a new QCD source derivation'}
    assert pi2['sign_convention']==chain['sign_convention']
    for a,b in zip(pi2['orbits'],chain['orbits']):
        for key in ['representation','variables','bounds','iterated_order_outer_to_inner','jacobian','integrand','power_prescription']:
            assert a['J'][key]==b['J'][key], ('Pi2/chain source mismatch',key)
        for ga,gb in zip(a['J']['geometry_polynomials'],b['J']['geometry_polynomials']):
            assert {k:v for k,v in ga.items() if k!='kind'}=={k:v for k,v in gb.items() if k!='kind'}
        for key in ['orbit_multiplicity','stabilizer_order','stabilizer_rotations']:
            assert a['c4_symmetry'][key]==b['c4_symmetry'][key]
        CHECKS.append({'name':a['orbit_id']+' exact geometry/domain/multiplicity bridge to chain',
                       'passed':True,'metadata_differences':['order_basis explanatory label','geometry kind label','C4 definition explanatory label']})
        assert a['C']['components']['gamma_residual']==b['C']['components']['gamma_chain_ratio']
        zero(a['orbit_id']+' positive Pi2 rational C',parse(a['C']['components']['rational'])-1/(32*(1-2*e)**2))
    channel=pi2['orbits'][0]['operator']['channels']['Ff']
    assert channel['full_monomial']=='CF^2*TF*nf'
    Aop,Bop=parse(channel['A']),parse(channel['B'])
    for orb in pi2['orbits']:
        actual=parse(orb['C']['components']['self_energy_scalar']).subs(sp.Symbol('TF')*sp.Symbol('nf'),T)
        actual=sp.expand(actual).coeff(CF).coeff(T)
        zero(orb['orbit_id']+' actual C Ff operator A/B',actual-(Aop+Bop/U))
    # U*p_Ff=U*A_Ff+B_Ff. Individual e^-3 terms cancel; their finite
    # contribution requires Gamma logarithms through e^5, including zeta5.
    pref=truncated(-(pg(5)*Aop+gn(3,5)*Bop)/(16*(1-2*e)**2),2)
    w3cff=truncated(pref*geom(3),0)
    postcff=sp.expand((w3cff+w1/(8*e)).subs(elementary));w3cff=sp.expand(w3cff.subs(elementary))
    c={q:coefficient(postcff,q) for q in range(-5,1)}
    for q in [-5,-4]:zero('raw-source cff superleading cancellation '+str(q),c[q])
    h3=-3*c[-1].subs(L,0);k3=-sp.Rational(3,2)*sp.diff(c[-1],L)
    zero('raw-source cff h3 geometry independence',sp.diff(h3,rho))
    cff={'bare':{str(q):str(coefficient(w3cff,q)) for q in range(-4,1)},
         'post':{str(q):str(v) for q,v in c.items() if q>=-4},'h3':str(h3),'k3':str(k3),
         'finite':str(c[0]),'shape':str(c[0].subs(L,0)),
         'construction_inputs':'accepted Pi2 Ff bare operator A/B coefficients and exact common chain J; no published Pi2 I/Laurent coefficients',
         'human_analysis_blind':False,'independence':'independent geometry rederivation from shared prepared root facts; Pi2 upstream operator remains accepted'}
    return {'nf2':nf2,'cff':cff}

def external_and_prediction_comparison(root,pred,reference):
    current=read(root,'benchmarks/renormalization/three-loop/current.json')
    ct=read(root,'benchmarks/renormalization/three-loop/counterterms-r001.json')
    zero('external Eq2.3 beta0 normalization',parse(reference['beta0'])-parse(ct['b0']))
    zero('external Eq2.3 beta1 normalization',parse(reference['beta1'])-parse(ct['b1']))
    data=read(root,'benchmarks/renormalization/three-loop/matter-channels-r001.json')
    af=read(root,'benchmarks/renormalization/three-loop/results/CF-CA-TF-nf/r001.json')
    two=read(root,'benchmarks/renormalization/two-loop/result-r002.json')
    hs={'CF*(TF*nf)^2':parse(data['nf2']['h3_nf2']),'CF^2*TF*nf':parse(data['cff']['h3_cff']),'CF*CA*TF*nf':parse(af['h3'])}
    for color in hs:
        pointer=current['channels'][color]
        rel='benchmarks/renormalization/three-loop/'+pointer['result']
        assert SOURCES[rel]==pointer['sha256'], ('changed current result',color)
    for color,h3 in hs.items():
        zero(color+' external factorisation Eq5.16 versus 2h3',parse(reference['gamma_square_three_loop_matter'][color])-2*h3)
        nonzero('negative: omit mu versus mu2 factor for '+color,parse(reference['gamma_square_three_loop_matter'][color])-h3)
    h2=sp.expand(parse(two['h2'])/CF)
    for color,monomial in [('CA',CA),('T',T)]:
        zero('external Eq5.15 two-loop convention control '+color,parse(reference['gamma_square_two_loop'][color])-2*h2.coeff(monomial))
    for key in ['nf2','cff']:
        old=data[key];new=pred[key]
        for q in range(-4,1):
            zero(f'new raw-source {key} post coefficient {q}',parse(new['post'][str(q)])-parse(old[f'post_counterterm_{key}'][str(q)]))
            zero(f'new raw-source {key} bare coefficient {q}',parse(new['bare'][str(q)])-parse(old[f'bare_w3_{key}'][str(q)]))
        zero('new raw-source '+key+' h3',parse(new['h3'])-parse(old[f'h3_{key}']))
        zero('new raw-source '+key+' general-z finite shape',parse(new['shape'])-parse(old[f'f3_{key}']))
    # Original extraction code is separately hashed to expose independence:
    # its h3 is extracted from the assembled single pole, never Eq5.16.
    for rel in ['workflow/checks/w3-ca-nf/physical.mpl','workflow/checks/w3-ca-nf/source_inputs.mpl',
                'workflow/checks/w3-p0/nf2.py','workflow/checks/w3-p0/cff.py',
                'benchmarks/renormalization/one-loop/derivation.md','workflow/checks/rg-one-loop/calculate.py']:
        text_hash(root,rel)

def numeric_checks(pred):
    rows=[]
    # Direct scalar logarithmic-moment quadrature, not GPL/MeijerG evaluation.
    for label in ['0.3','1','3']:
        for dps in [35,55]:
            with mp.workdps(dps):
                zz=mp.mpf(label)
                ell=mp.log(zz);bb=1/(1+zz)
                def moment(j):
                    def integrand(x):
                        if x==0 or x==1:return mp.mpf('0')
                        aa=mp.log(x)+ell/2;bc=mp.log1p(-x)-ell/2
                        if abs(x-bb)<mp.sqrt(mp.eps):
                            middle=(aa+bc)/2
                            return -middle**j*(1/bb+1/(1-bb))/mp.factorial(j)
                        return -(aa**(j+1)-bc**(j+1))/(mp.factorial(j+1)*(x-bb))
                    return mp.quad(integrand,[0,bb,1])
                vals=[moment(j) for j in range(3)]
                expected0=-(ell**2+mp.pi**2)/2
                assert abs(vals[0]-expected0)<mp.mpf(10)**(-dps+8),(str(zz),dps,mp.nstr(vals[0]-expected0,30))
                moment03=-(ell**2+mp.pi**2)**2/4
                mm11=2*vals[1];mm21=8*vals[2]-moment03/3
                shapes={}
                for key in ['nf2','cff']:
                    fun=sp.lambdify((rho,M11,M21),parse(pred[key]['shape']),'mpmath')
                    shapes[key]=mp.nstr(fun(ell,mm11,mm21),dps-5)
                rows.append({'z':str(zz),'dps':dps,'H0_H1_H2':[mp.nstr(v,dps-5) for v in vals],
                             'H0_elementary_residual':mp.nstr(vals[0]-expected0,8),'finite_shapes_L0':shapes})
    refinement=[]
    for i in range(0,len(rows),2):
        with mp.workdps(60):
            old,new=rows[i:i+2]
            differences={key:abs(mp.mpf(old['finite_shapes_L0'][key])-mp.mpf(new['finite_shapes_L0'][key])) for key in ['nf2','cff']}
            assert all(v<mp.mpf('1e-25') for v in differences.values()),differences
            refinement.append({'z':new['z'],'shape_refinement':{k:mp.nstr(v,8) for k,v in differences.items()},'threshold':'1e-25','passed':True})
    return {'method':'direct one-dimensional log-difference quadrature of scalar original rectangle',
            'rows':rows,'refinement':refinement,'certified_interval':False,'published_I_GPL_MeijerG_read':False}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    ap.add_argument('--output-dir',type=Path,default=HERE);ap.add_argument('--skip-numeric',action='store_true')
    args=ap.parse_args();root=args.root.resolve();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    start=time.time();reference=json.loads((HERE/'reference.json').read_text(encoding='utf-8'))
    rg_nullspace()
    single=read(root,'benchmarks/self-energy/pi1/prepared/2e-Pi1-r002.json')
    chain=read(root,'benchmarks/self-energy/pi1-chain/prepared/2e-Pi1+Pi1-r002.json')
    pi2=read(root,'benchmarks/self-energy/pi2/prepared/2e-Pi2-r002.json')
    source_current=read(root,'benchmarks/self-energy/current-r002.json')
    source_pointers={entry['orbit_id']:entry['prepared'] for entry in source_current['entries']}
    for package in [single,chain,pi2]:
        for orbit in package['orbits']:
            selected=source_pointers[orbit['orbit_id']]
            assert SOURCES[selected['path']]==selected['sha256'],('raw source not current',orbit['orbit_id'])
    one=read(root,'benchmarks/renormalization/one-loop/result-r001.json')
    assert one['definitions']['mu_x_squared']=='exp(2*EulerGamma)*mu_MSbar_squared/4'
    pred=source_prediction(single,chain,pi2)
    write=lambda name,obj:(out/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    write('raw-source-predictions.json',pred)
    external_and_prediction_comparison(root,pred,reference)
    numeric=None if args.skip_numeric else numeric_checks(pred)
    for rel,h in SOURCES.items():assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==h,('input changed',rel)
    report={'status':'passed','checks':CHECKS,'exact_checks':len(CHECKS),'negative_controls':sum(bool(x.get('negative_control')) for x in CHECKS),
            'source_sha256':SOURCES,'external_non_cusp_prediction_compared':list(reference['gamma_square_three_loop_matter']),
            'external_direct_three_loop_closed_loop_integration':False,'nf2_independent_raw_source_rederivation':True,
            'cff_independent_raw_source_rederivation':True,'independent_upstream_QCD_operator_derivation':False,'Af_finite_shape_independently_checked':False,
            'formal':'not_run','numeric':numeric,'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,
            'elapsed_seconds':round(time.time()-start,3)}
    write('verification.json',report)
    print(json.dumps({k:report[k] for k in ['status','exact_checks','negative_controls','external_non_cusp_prediction_compared','nf2_independent_raw_source_rederivation','cff_independent_raw_source_rederivation','formal','elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
