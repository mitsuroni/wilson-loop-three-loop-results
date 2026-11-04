"""Original two-endpoint chain/Pi2 evaluation and full colour convolution.

Uses native prepared C/J and existing evaluated M words. No result table,
operator-expansion table, selector, CAS process or original writer is opened.
"""
from __future__ import annotations
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import sympy as s

PACK=Path(__file__).resolve().parents[1]
E,Z=s.symbols('eps z',positive=True)
X,Y=s.symbols('x y',positive=True)
RHO=s.Symbol('rho',real=True)
M01,M11,M21,M03=s.symbols('M01 M11 M21 M03')
MOMENT_RUN='workflow/checks/single-opposite/runs/20260916T090538777022Z'

def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def parse(text):
    text=re.sub(r'Gamma\[([^\]]+)\]',r'gamma(\1)',text).replace('^','**')
    return s.sympify(text,locals={'eps':E,'z':Z,'z1':Z,'s31':X,'s21':X,'s41':Y,
        'gamma':s.gamma,'Pi':s.pi,'CA':s.Symbol('CA'),'CF':s.Symbol('CF'),
        'TF':s.Symbol('TF'),'nf':s.Symbol('nf')})

def coeffs(expr,lo=-6,hi=0):
    value=s.series(expr,E,0,hi+1).removeO().expand()
    return {p:s.cancel(value.coeff(E,p)) for p in range(lo,hi+1)}

def exp_jet(logarithm,order):
    rows={0:s.Integer(1)}
    for n in range(1,order+1):rows[n]=s.expand(sum(k*logarithm.get(k,0)*rows[n-k] for k in range(1,n+1))/n)
    return rows

def gamma_jet(slopes,order):
    logs={k:sum(mult*(-slope)**k for slope,mult in slopes)*s.zeta(k)/k for k in range(2,order+1)}
    return exp_jet(logs,order)

def product(a,b,lo=-6,hi=0):
    return {p:s.expand(sum(v*b.get(p-k,0) for k,v in a.items())) for p in range(lo,hi+1)}

def rational_coeffs(expr,lo,hi):
    """Exact finite Laurent division; avoids expanding a large rational sum."""
    n,d=s.fraction(s.cancel(expr))
    if n==0:return {p:s.Integer(0) for p in range(lo,hi+1)}
    pn,pd=s.Poly(n,E),s.Poly(d,E)
    n0=min(t[0][0] for t in pn.terms());d0=min(t[0][0] for t in pd.terms())
    offset=n0-d0;regular={}
    for j in range(max(0,hi-offset+1)):
        regular[j]=s.cancel((pn.nth(n0+j)-sum(pd.nth(d0+k)*regular[j-k] for k in range(1,j+1)))/pd.nth(d0))
    return {p:regular.get(p-offset,s.Integer(0)) for p in range(lo,hi+1)}

def generate(snapshot_root=None):
    root=Path(snapshot_root) if snapshot_root is not None else PACK/'snapshots'
    path=root/'workflow/checks/single-opposite/hlog_algebra.py'
    spec=importlib.util.spec_from_file_location('chains_original_hlog_algebra',path)
    ha=importlib.util.module_from_spec(spec);spec.loader.exec_module(ha)
    moments=[root/MOMENT_RUN/f'M_{k}_{j}.txt' for k,j in [(0,1),(1,1),(2,1),(0,3)]]
    vectors={atom:ha.parse(p) for atom,p in zip([M01,M11,M21,M03],moments)};vectors[RHO]=ha.LOG
    checks=[]
    def check(label,condition):
        if not condition:raise ValueError('chain/Pi2 check: '+label)
        checks.append({'name':label,'passed':True})
    def eq(label,a,b=0):
        residual=s.cancel(s.expand(a-b))
        check(label,residual==0 or s.simplify(residual)==0)
    def rows(expressions):
        return {str(p):{'expression':s.sstr(v),'language':'SymPy/str','latex':s.latex(v),
                       'hlog_terms':ha.serialize(ha.from_sympy(v,vectors))} for p,v in expressions.items()}
    # Native C contains G^3 U. U has no epsilon^1 or epsilon^2 term.
    uslopes=[(1,2),(-1,1),(-3,1),(-2,-2),(2,-1)]
    Uj=gamma_jet(uslopes,5);invUj=gamma_jet([(k,-m) for k,m in uslopes],5)
    Pj=gamma_jet([(-1,4),(1,2),(-3,1),(-2,-2),(2,-1)],4)
    Gj=gamma_jet([(-1,3)],4)
    eq('U no quadratic',Uj[2]);eq('U cubic',Uj[3],6*s.zeta(3))
    for k,v in product(Uj,invUj,0,5).items():eq('U inverse coefficient '+str(k),v,int(k==0))
    regP=sum(v*E**p for p,v in Pj.items());regG=sum(v*E**p for p,v in Gj.items())
    Hjet=M01+3*E*M11/2+9*E**2*(M21/8+M03/24)
    bare=[(2+E)/(1152*E**4*(1-2*E)**2*(1+2*E)),
          1/(384*E**3*(1-2*E)**2*(1+2*E))+Hjet/(128*E**2*(1-2*E)**2*s.cosh(3*E*RHO/2))]
    geom={i:coeffs(regP*b) for i,b in enumerate(bare)}
    plain={i:coeffs(regG*b) for i,b in enumerate(bare)}
    G=s.exp(-s.EulerGamma*E)*s.gamma(1-E)
    U=s.gamma(1+E)**2*s.gamma(1-E)*s.gamma(1-3*E)/(s.gamma(1-2*E)**2*s.gamma(1+2*E))
    delta=s.Symbol('delta')
    H0=-(1+Z)*Z**(-delta/2)*s.meijerg([delta,0,0],[],[0,delta,delta],[],Z)/(s.gamma(1-delta)*s.gamma(1+delta)**2)
    full_base=[G**3*U*(2+E)/(1152*E**4*(1-2*E)**2*(1+2*E)),
               G**3*U/(384*E**3*(1-2*E)**2*(1+2*E))+G**3*U*H0.subs(delta,3*E)/(128*E**2*(1-2*E)**2*s.cosh(3*E*s.log(Z)/2))]
    results={};operator_jets={}
    for family in ['2e-Pi1+Pi1','2e-Pi2']:
        source=root/'benchmarks/current/prepared/3loop'/f'{family}.json'
        package=load(source)
        for index,orbit in enumerate(package['orbits']):
            oid=orbit['orbit_id'];j=orbit['J'];c=orbit['C']['components']
            check(oid+' one Euler compensation',orbit['normalization']['euler_compensation_count']==1)
            eq(oid+' explicit Euler factor',parse(c['euler']).subs(s.Symbol('EulerGamma'),s.EulerGamma),s.exp(-3*s.EulerGamma*E))
            check(oid+' native physical sign',orbit['normalization']['assembly_sign']==1)
            check(oid+' original unit domain',all(b==['0','1'] for b in j['bounds'].values()) and j['jacobian']=='1')
            check(oid+' prescribed iterated order',j['iterated_order_outer_to_inner']==(['s31','s41'] if index==0 else ['s21','s41']))
            check(oid+' C4 orbit factor',orbit['c4_symmetry']['orbit_multiplicity']==(4 if index==0 else 2))
            D=Z*(1-X)*Y if index==0 else (1-X)*(1-Y)+Z*X*Y
            N=(2+E)*Z*D if index==0 else -(2+E)*(1+Z)*D+(1-3*E)*Z
            eq(oid+' actual source distance',parse(j['geometry_polynomials'][0]['f']),D)
            eq(oid+' complete numerator',parse(j['se_numerator']),N)
            eq(oid+' actual source kernel',parse(j['integrand']),N*D**(-2+3*E)/(2*(1+2*E)))
            if index==1:
                mixed=3*E*(s.diff(D,X,Y)*D+(3*E-1)*s.diff(D,X)*s.diff(D,Y))
                eq(oid+' full mixed-derivative reduction',N/(2*(1+2*E)),-(1+Z)*D+mixed/(6*E*(1+2*E)))
                for xx,yy,value in [(0,0,1),(1,1,Z),(0,1,0),(1,0,0)]:eq(oid+' ordinary boundary '+str((xx,yy)),D.subs({X:xx,Y:yy}),value)
            else:
                # The two actual endpoint powers integrate to 1/(3 eps).
                eq(oid+' endpoint primitive derivative',s.diff(X**(3*E)/(3*E),X),X**(3*E-1))
            Rnative=parse(c['rational'])
            if family=='2e-Pi2':Rnative=Rnative/(4*E**2)
            sign=s.cancel(Rnative*128*E**2*(1-2*E)**2)
            check(oid+' C-derived constant normalization',sign in (-1,1))
            eq(oid+' native residual Gamma ratio',parse(c.get('gamma_chain_ratio',c.get('gamma_residual'))),U)
            geoc={p:s.expand(sign*v) for p,v in geom[index].items()}
            channels=[];operator_coeffs={}
            if family=='2e-Pi1+Pi1':
                CA,TF,nf=s.symbols('CA TF nf');B=((10-6*E)*CA-8*(1-E)*TF*nf)/(3-2*E)
                eq(oid+' native full self-energy colour',parse(c['self_energy_color']),B**2)
                a=(10-6*E)/(3-2*E);b=8*(1-E)/(3-2*E)
                mults=[('AA','CA^2*CF',a*a),('Af','CA*CF*TF*nf',-2*a*b),('ff','CF*(TF*nf)^2',b*b)]
                for name,monomial,mult in mults:
                    phi=rational_coeffs(mult,0,4);final=product(phi,geoc)
                    operator_coeffs[name]={str(p):s.sstr(v) for p,v in phi.items()}
                    encoded=rows(final)
                    channels.append({'channel':name,'monomial':monomial,'color_monomial':monomial,
                        'exact_multiplier':s.sstr(mult),'multiplier_coefficients':operator_coeffs[name],
                        'coefficients':encoded,'laurent':{'from_power':-6,'through_power':0,'coefficients':encoded,'remainder':'O(eps)'}})
                external='CF*B1(eps)^2'
            else:
                operator=orbit['operator'];check(oid+' adopted two-loop 1PI source',operator['insertion_count']==1 and operator['subgraph_loop_order']==2 and operator['not_Pi1_squared'])
                colors=[]
                for name,monomial in [('AA','CA^2*CF'),('Af','CA*CF*TF*nf'),('Ff','CF^2*TF*nf')]:
                    A,B=[parse(operator['channels'][name][k]) for k in ['A','B']]
                    p_exact=parse(operator['channels'][name]['p']);eq(oid+' source p '+name,p_exact,A+B/U)
                    colors.append(parse(operator['channels'][name]['operator_monomial'])*p_exact)
                    key=(s.sstr(A),s.sstr(B))
                    if key not in operator_jets:
                        aj,bj=rational_coeffs(A,-3,2),rational_coeffs(B,-3,2)
                        phi={k+2:4*s.expand(aj[k]+sum(v*invUj.get(k-j,0) for j,v in bj.items())) for k in range(-3,3)}
                        psi={k+2:4*s.expand(bj[k]+sum(v*Uj.get(k-j,0) for j,v in aj.items())) for k in range(-3,3)}
                        operator_jets[key]=(phi,psi)
                    phi,psi=[dict(v) for v in operator_jets[key]]
                    eq(oid+' no artificial phi pole '+name,phi.pop(-1));eq(oid+' no artificial Uphi pole '+name,psi.pop(-1))
                    final=product(phi,geoc);direct=product(psi,{p:sign*v for p,v in plain[index].items()})
                    for power in range(-6,1):eq(oid+' independent direct C route '+name+str(power),final[power],direct[power])
                    operator_coeffs[name]={str(p):s.sstr(v) for p,v in phi.items()}
                    encoded=rows(final)
                    channels.append({'channel':name,'monomial':monomial,'color_monomial':monomial,
                        'exact_multiplier':s.sstr(4*E**2*p_exact),'multiplier_coefficients':operator_coeffs[name],
                        'coefficients':encoded,'laurent':{'from_power':-6,'through_power':0,'coefficients':encoded,'remainder':'O(eps)'},
                        'operator_source':'Exact adopted native full-D Pi2 A+B/U; phi regenerated, not read from an expansion table.'})
                eq(oid+' full native Pi2 scalar',parse(c['self_energy_scalar']),sum(colors))
                external='CF*4*eps^2*P2(eps)'
            encoded=rows(geoc)
            result={'format':'wilson-p5-original-evaluation/1','orbit_id':oid,'family':family+('-adjacent' if index==0 else '-opposite'),
                'quantity':'physical I (shared Pi2 geometry before external phi_c)' if family=='2e-Pi2' else 'physical I',
                'normalization':{'scheme':'native-physical-CJ','physical_convention':'wilson-qcd-action-path-fourier/1.0','assembly_sign':1,
                    'loop_order':3,'Euler_count':1,'sym':orbit['c4_symmetry']['orbit_multiplicity'],'a1':'1',
                    'external_color':external,'C_derived_residual_sign':str(sign),'includes_U':True,
                    'formula':'(sym/2)*R_native*(J(z)+z^(3eps)*J(1/z))/(1+z^(3eps))'},
                'source':{'path':source.relative_to(root).as_posix(),'sha256':sha(source),'json_pointer':'/orbits/'+str(index)},
                'J':j,'I':{'sympy':s.sstr(sign*full_base[index]),'language':'SymPy/str'},
                'laurent':{'from_power':-6,'through_power':0,'coefficients':encoded,'remainder':'O(eps)'},
                'color_channels':channels,'operator_coefficients':operator_coeffs,
                'subintegrals':{'adjacent':'Two original endpoint integrals 1/(3eps), multiplied in prescribed order.',
                    'opposite':'H0(3eps,z) plus complete mixed-derivative boundary (1+z^(3eps))/(6eps(1+2eps)).',
                    'H0_series':'M01+(3eps/2)M11+9eps^2(M21/8+M03/24)+O(eps^3)',
                    'regular_geometry_budget':4 if index==0 else 3,'operator_multiplier_budget':4 if index==0 else 3,
                    'Pi2_Gamma_ratio_budget_before_cancellation':5,'H0_jet_order':2,
                    'evaluated_M_reuse':[{'path':p.relative_to(root).as_posix(),'sha256':sha(p)} for p in moments]},
                'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/3; local meromorphic continuation at eps=0','complex_kinematics':'not newly verified'},
                'generation_scope':'Fresh source checks, endpoint/IBP boundary reduction, Gamma jets and colour convolutions; existing evaluated M words reused. Native C determines sign. No final table or new CAS read.',
                'formal':'not_run'}
            if family=='2e-Pi2':result['I_geom']=result['I'];result['I_geom_laurent']=result['laurent']
            results[oid]=result
    # Identity is a consequence of C/J and does not identify the operators.
    for index in [1,2]:
        a=results[f'2e-Pi1+Pi1-{index:02d}'];b=results[f'2e-Pi2-{index:02d}']
        for p in range(-6,1):eq('current shared geometric sign relation '+str((index,p)),
            s.sympify(a['laurent']['coefficients'][str(p)]['expression'])+s.sympify(b['laurent']['coefficients'][str(p)]['expression']))
    for r in results.values():r['generation_checks']={'count':len(checks),'checks':checks}
    return results
