"""Complete ordered Pi1 P5, explicit mixed-exponent series and reused coefficients."""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import sys,json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha
e,z,L=s.symbols('eps z L',positive=True);x,y,u,v=s.symbols('x y u v',positive=True)
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def eq(name,value):
    assert s.simplify(s.expand(value))==0,(name,value)
    checks.append({'name':name,'passed':True})
@lru_cache(None)
def qpoly(k,zz):
    a,b=s.symbols('a b')
    return s.Poly(s.expand(((zz+a+b)/(1+zz)-a*b)**k),a,b).terms()
def block(k,l,zz):
    return s.factor(sum(c1*c2/((a+1)*(a+c+2)*(b+1)*(b+d+2)) for (a,b),c1 in qpoly(k,zz) for (c,d),c2 in qpoly(l,zz)))
def polynomial_integral(expr):
    # Original outer-to-inner u,x,v,y, with 0<x<u, 0<y<v.
    return s.integrate(s.integrate(s.integrate(s.integrate(expr,(y,0,v)),(v,0,1)),(x,0,u)),(u,0,1))
def main():
    package=BASE/'prepared/4e-Pi1.json';manifest=BASE/'prepared/manifest.json';meta=json.loads(manifest.read_text(encoding='utf-8'))
    for d in meta['inputs']+meta['outputs']+[meta['normalization'],meta['preparer'],meta['color_backend']]:assert sha(ROOT/d['path'])==d['sha256']
    selected=[o for o in json.loads(package.read_text(encoding='utf-8'))['orbits'] if o['orbit_id'] in ['4e-Pi1-07','4e-Pi1-08']]
    D=lambda a,b:(1-a)*(1-b)+z*a*b
    N=lambda a,b:-2*(1+z)*D(a,b)+(1-2*e)*z
    for o in selected:
        j=o['J'];local={'eps':e,'z1':z,'s21':x,'s41':y,'s22':u,'s42':v}
        parse=lambda t:s.sympify(t.replace('^','**'),locals=local)
        a,b,c,d=(x,y,u,v) if o['orbit_id'].endswith('07') else (u,v,x,y)
        eq(o['orbit_id']+' SE distance',parse(j['geometry_polynomials'][0]['f'])-D(a,b))
        eq(o['orbit_id']+' bare distance',parse(j['geometry_polynomials'][1]['f'])-D(c,d))
        eq(o['orbit_id']+' complete numerator',parse(j['se_numerator'])-N(a,b))
        eq(o['orbit_id']+' complete source',parse(j['integrand'])-N(a,b)*D(a,b)**(-2+2*e)/(2*(1+e))*(-1-z)*D(c,d)**(-1+e))
        assert o['I_target']['sym']==2 and o['I_target']['a1']=='-1/2'
        assert j['bounds']=={'s21':['0','s22'],'s22':['0','1'],'s41':['0','s42'],'s42':['0','1']}
        assert j['iterated_order_outer_to_inner']==['s22','s21','s42','s41']
    eq('reflection distance',D(1-x,1-y)-z*D(x,y).subs(z,1/z))
    eq('reflection full numerator',N(1-x,1-y)-z**2*N(x,y).subs(z,1/z))
    mixed=2*e*(s.diff(D(x,y),x,y)*D(x,y)+(2*e-1)*s.diff(D(x,y),x)*s.diff(D(x,y),y))
    eq('complete SE derivative identity',N(x,y)/(2*(1+e))+(1+z)*D(x,y)-mixed/(4*e*(1+e)))
    # The two rectangles have four corners each. Keep all eight terms.
    H03,H01,H2p,H2m,T=s.symbols('H03 H01 H2p H2m Tmixed_sum')
    raw_boundary=(2*z**(3*e/2)*H03+(1+z**(2*e))*z**(e/2)*H01-2*z**(e/2)*(H2p-H2m)-2*z**(5*e/2)*(H2p+H2m))/(4*e*(1+e))
    A=H03+s.cosh(e*s.log(z))*(H01-2*H2p)-2*s.sinh(e*s.log(z))*H2m
    eq('eight boundary corners assembled',(raw_boundary-z**(3*e/2)*A/(2*e*(1+e))).rewrite(s.exp))
    G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
    full=G**3*(T+A/(2*e*(1+e)))/(128*e*(1-2*e)*s.cosh(3*e*s.log(z)/2))
    eq('sym=2 normalized assembly',(G**3*(z**(3*e/2)*T+raw_boundary)/(64*e*(1-2*e)*(1+z**(3*e)))-full).rewrite(s.exp))
    # Fully integrated polynomial blocks; no uncomputed coefficient is introduced.
    for aa,bb in [(1,1),(2,1),(1,2),(4,2)]:
        zz=s.Rational(2,3)
        actual=polynomial_integral(s.expand((1+zz)**2*D(x,y).subs(z,zz)**(aa-1)*D(u,v).subs(z,zz)**(bb-1)))
        summed=(1+zz)**(aa+bb)*sum(s.rf(1-aa,k)*s.rf(1-bb,l)/s.factorial(k)/s.factorial(l)*block(k,l,zz) for k in range(aa) for l in range(bb))
        eq(f'mixed terminating series {aa},{bb}',actual-summed)
    R0,R1=s.symbols('rho0 rho1');t00,t01,t02,t21,t22,eta21=s.symbols('tau0q0 tau0q1 tau0q2 tau2q1 tau2q2 eta2q1')
    # Mixed ordered kernels are symmetric under exchange of their regulator labels
    # after summing. Differentiate at the common zero: coefficient =3*rho1.
    h03=t00+3*e*t01+9*e**2*t02;h01=t00+e*t01+e**2*t02;h2p=t00+e*t21+e**2*t22
    AA=h03+s.cosh(e*L)*(h01-2*h2p)-2*s.sinh(e*L)*e*eta21/L
    expr=(2*R0+3*e*R1+AA/(2*e*(1+e)))/(128*e*(1-2*e)*s.cosh(3*e*L/2))
    series=s.series(expr,e,0,1).removeO().expand()
    single_path=ROOT/'benchmarks/2ne/families/single-opposite/basis-functions.json'
    double_path=ROOT/'benchmarks/2ne/families/same-channel-two/structural-basis-functions.json'
    single=json.loads(single_path.read_text());double=json.loads(double_path.read_text())['functions']
    vector=lambda v:{(tuple(a['word']),tuple(a['zeta_powers'])):Fraction(a['coefficient']) for a in v['hlog_terms']}
    vectors={q:vector(single[str(q)]) for q in [t00,t01,t02,t21,t22,eta21]};vectors.update({R0:vector(double['rho_id_0_0_0']),R1:vector(double['rho_id_0_0_1']),L:ha.LOG})
    run=ROOT/json.loads((ROOT/'workflow/checks/single-opposite/current-moments.json').read_text())['run']
    M01,M11,M21,M03=s.symbols('M01 M11 M21 M03')
    mfiles=[run/f'M_{k}_{n}.txt' for k,n in [(0,1),(1,1),(2,1),(0,3)]]
    vectors.update({m:ha.parse(p) for m,p in zip([M01,M11,M21,M03],mfiles)})
    expected={-1:M01**2/384-M11/256,0:(3*R1-M21/2-M03/3+L**2*M01/4+2*M01**2/3-M11/2)/128}
    coeff={}
    for power in range(-6,1):
        actual=ha.from_sympy(series.coeff(e,power),vectors);goal=expected.get(power,s.S.Zero)
        assert ha.add(actual,ha.scale(ha.from_sympy(goal,vectors),-1))=={},power
        checks.append({'name':'Laurent Hlog reduction '+str(power),'passed':True})
        coeff[str(power)]={'expression':str(goal),'language':'SymPy/str','latex':s.latex(goal),'hlog_terms':ha.serialize(actual)}
    dependencies=[single_path,double_path,ROOT/'benchmarks/2ne/families/same-channel-two/unified-structure.md',ROOT/'benchmarks/2ne/families/same-channel-two/method.md',*mfiles]
    series_definition={'T(a,b,z)':'((1+z)/sqrt(z))^(a+b)*Sum[(1-a)_k*(1-b)_l/(k!*l!)*B(k,l,z),{k,0,infinity},{l,0,infinity}]',
      'B(k,l,z)':'Sum[p_k(i,j;z)*p_l(m,n;z)/((i+1)*(i+m+2)*(j+1)*(j+n+2)),{i,j=0..k},{m,n=0..l}]',
      'p_k(i,j;z)':'Sum[k!*(-1)^h*z^(k-i-j+h)/((k-i-j+h)!*(i-h)!*(j-h)!*h!*(1+z)^(k-h)),{h,max(0,i+j-k),min(i,j)}]',
      'order':'each finite B block first, then the double sum; no arbitrary rearrangement of signed monomials',
      'domain':'z>0, initially 0<Re(eps)<1/2 with a=2eps,b=eps; locally normally convergent near eps=0'}
    for o in selected:
        result={'format':'wilson-computational-result/1','orbit_id':o['orbit_id'],'revision':'r001','family':'Pi1-P5','input':{'package':package.relative_to(ROOT).as_posix(),'package_sha256':sha(package),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},'normalization':{'scheme':'pi1-color-euler-normalization-1.0','loop_order':3,'sym':2,'a1':'-1/2','Euler_count':1,'external_color':'CA*CF*((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)'},
          'I':{'kind':'explicit convergent mixed-exponent series plus MeijerG','sympy_template':str(full),'bindings':{'Tmixed_sum':'T(2eps,eps,z)+T(eps,2eps,z)','H03':'H0(3eps,z)','H01':'H0(eps,z)','H2p':'(H2(eps,z)+H2(eps,1/z))/2','H2m':'(H2(eps,z)-H2(eps,1/z))/2'},'series':series_definition},
          'laurent':{'from_power':-6,'through_power':0,'leading_power':-1,'remainder':'O(eps)','coefficients':coeff},'functions':{'rho1':'existing rho_id_0_0_1 = [eps] Q_id;0,0(eps,z), explicit Hlog in dependency','L':'log(z)'},'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in dependencies],'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/2; meromorphic near zero','complex_kinematics':'not verified'},'formal':'not_run'}
        p=BASE/f"results/{o['orbit_id']}/r001.json";p.parent.mkdir(parents=True,exist_ok=True);data=json.dumps(result,indent=2,ensure_ascii=False)+'\n'
        if p.exists():assert p.read_text(encoding='utf-8')==data
        else:p.write_text(data,encoding='utf-8')
    (HERE/'exact-checks.json').write_text(json.dumps({'count':len(checks),'checks':checks,'positions':14,'nonzero':4,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)} exact checks; simple poles; 14 slots / 4 nonzero; two results frozen')
if __name__=='__main__':main()
