"""Solve the two frozen chain I targets, with full flavor-channel reconstruction."""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib,re
import sympy as s
from sympy.printing.mathematica import mathematica_code
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1-chain'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha
e,z,x,y,L=s.symbols('eps z x y L',positive=True)
M01,M11,M21,M03,T02=s.symbols('M01 M11 M21 M03 T02')
G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
U=s.gamma(1+e)**2*s.gamma(1-e)*s.gamma(1-3*e)/(s.gamma(1-2*e)**2*s.gamma(1+2*e))
R=-G**3*U/(128*e**2*(1-2*e)**2)
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def eq(name,value):
    assert s.simplify(s.expand_func(value))==0,(name,value)
    checks.append({'name':name,'passed':True})
def freeze(p,obj):
    data=(json.dumps(obj,ensure_ascii=False,indent=2)+'\n').encode('utf-8');p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==data,p
    else:p.write_bytes(data)
def parse(t):return s.sympify(t.replace('^','**'),locals={'eps':e,'z1':z,'s31':x,'s21':x,'s41':y})
def main():
    manifest=BASE/'prepared/manifest.json';meta=json.loads(manifest.read_text(encoding='utf-8'))
    for d in meta['inputs']+meta['outputs']+[meta['normalization'],meta['preparer'],meta['color_backend']]:assert sha(ROOT/d['path'])==d['sha256'],d['path']
    package=BASE/'prepared/2e-Pi1+Pi1.json';orbits=json.loads(package.read_text(encoding='utf-8'))['orbits']
    assert [o['I_target']['sym'] for o in orbits]==[4,2] and all(o['I_target']['a1']=='1' for o in orbits)
    assert all(o['normalization']['euler_compensation_count']==1 for o in orbits)
    for o in orbits:
        assert o['J']==o['J_raw_s']
        assert o['J']['iterated_order_outer_to_inner']==(['s31','s41'] if o['orbit_id'].endswith('01') else ['s21','s41'])
        assert all(b==['0','1'] for b in o['J']['bounds'].values())
    ja=orbits[0]['J'];f=z*(1-x)*y
    eq('CH_A actual source distance',parse(ja['geometry_polynomials'][0]['f'])-f)
    eq('CH_A complete source numerator',parse(ja['se_numerator'])-(2+e)*z*f)
    eq('CH_A actual integrand assembly',parse(ja['integrand'])-parse(ja['se_numerator'])*f**(-2+3*e)/(2*(1+2*e)))
    endpoint=s.integrate(x**(3*e-1),(x,0,1));eq('CH_A positive endpoint integral',endpoint-1/(3*e))
    JA=(2+e)*z**(3*e)*endpoint**2/(2*(1+2*e))
    IA=-G**3*U*(2+e)/(1152*e**4*(1-2*e)**2*(1+2*e))
    eq('CH_A source-to-I with m=4',2*R*(JA+z**(3*e)*JA.subs(z,1/z))/(1+z**(3*e))-IA)
    jb=orbits[1]['J'];D=(1-x)*(1-y)+z*x*y;d=-(1+z)
    NN=-(2+e)*(1+z)*D+(1-3*e)*z
    eq('CH_O actual source distance',parse(jb['geometry_polynomials'][0]['f'])-D)
    eq('CH_O complete source numerator',parse(jb['se_numerator'])-NN)
    eq('CH_O actual integrand assembly',parse(jb['integrand'])-NN*D**(-2+3*e)/(2*(1+2*e)))
    mixed=3*e*(s.diff(D,x,y)*D+(3*e-1)*s.diff(D,x)*s.diff(D,y))
    eq('CH_O full derivative identity',NN/(2*(1+2*e))-d*D-mixed/(6*e*(1+2*e)))
    for ss,tt,value in [(0,0,1),(1,1,z),(0,1,0),(1,0,0)]:eq(f'CH_O corner {ss},{tt}',D.subs({x:ss,y:tt})-value)
    H=s.Symbol('H0_3eps')
    JB=z**(3*e/2)*H+(1+z**(3*e))/(6*e*(1+2*e))
    eq('CH_O J crossing',s.powsimp(JB-z**(3*e)*JB.subs(z,1/z),force=True))
    IO_symbol=-G**3*U/(384*e**3*(1-2*e)**2*(1+2*e))+R*H/s.cosh(3*e*s.log(z)/2)
    eq('CH_O C4 source-to-I with m=2',(R*(JB+z**(3*e)*JB.subs(z,1/z))/(1+z**(3*e))-IO_symbol).rewrite(s.exp))
    delta=s.Symbol('delta')
    H0=-(1+z)*z**(-delta/2)*s.meijerg([delta,0,0],[],[0,delta,delta],[],z)/(s.gamma(1-delta)*s.gamma(1+delta)**2)
    IO=IO_symbol.subs(H,H0.subs(delta,3*e))
    # Source-bound evaluated moments, not fresh unevaluated coefficient definitions.
    basis_path=ROOT/'benchmarks/2ne/families/single-opposite/basis-functions.json'
    basis=json.loads(basis_path.read_text(encoding='utf-8'))
    vector=lambda key:{(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in basis[key]['hlog_terms']}
    run=ROOT/json.loads((ROOT/'workflow/checks/single-opposite/current-moments.json').read_text())['run']
    files=[run/f'M_{k}_{n}.txt' for k,n in [(0,1),(1,1),(2,1),(0,3)]]
    vectors={m:ha.parse(p) for m,p in zip([M01,M11,M21,M03],files)};vectors[L]=ha.LOG
    for name,residual in [
        ('H0 constant is M01',ha.add(vector('tau0q0'),ha.scale(vectors[M01],-1))),
        ('H0 first coefficient is M11/2',ha.add(vector('tau0q1'),ha.scale(vectors[M11],-s.Rational(1,2)))),
        ('H0 second coefficient in existing M',ha.add(vector('tau0q2'),ha.scale(vectors[M21],-s.Rational(1,8)),ha.scale(vectors[M03],-s.Rational(1,24))))]:
        assert residual=={},name;checks.append({'name':name,'passed':True})
    # log Gamma(1+t)=-gamma*t+sum_{k>=2}(-1)^k*zeta(k)*t^k/k.
    slopes=[(-1,4),(1,2),(-3,1),(-2,-2),(2,-1)]
    gamma_product=s.prod(s.gamma(1+k*e)**p for k,p in slopes)*s.exp(-3*s.EulerGamma*e)
    eq('complete Gamma product includes U',gamma_product/(G**3*U)-1)
    eq('Euler linear cancellation',-sum(k*p for k,p in slopes)-3)
    logreg=sum(s.zeta(k)*e**k/k*sum(p*(-slope)**k for slope,p in slopes) for k in range(2,5))
    regular=s.series(s.exp(logreg),e,0,5).removeO()
    eq('G3U quadratic coefficient',s.expand(regular).coeff(e,2)-3*s.zeta(2)/2)
    eq('G3U cubic coefficient',s.expand(regular).coeff(e,3)-7*s.zeta(3))
    eq('G3U quartic coefficient',s.expand(regular).coeff(e,4)-s.Rational(201,16)*s.zeta(4))
    seriesA=s.series(-regular*(2+e)/(1152*e**4*(1-2*e)**2*(1+2*e)),e,0,1).removeO().expand()
    hseries=M01+3*e*M11/2+9*e**2*(M21/8+M03/24)
    seriesO=s.series(-regular/(384*e**3*(1-2*e)**2*(1+2*e))-regular*hseries/(128*e**2*(1-2*e)**2*s.cosh(3*e*L/2)),e,0,1).removeO().expand()
    Z2,Z3,Z4=s.zeta(2),s.zeta(3),s.zeta(4)
    expectedA={-4:-s.Rational(1,576),-3:-s.Rational(5,1152),-2:-(18+3*Z2)/1152,
      -1:-(40+s.Rational(15,2)*Z2+14*Z3)/1152,
      0:-(112+27*Z2+35*Z3+s.Rational(201,8)*Z4)/1152}
    expectedO={-3:-s.Rational(1,384),-2:-(2+3*M01)/384,
      -1:-(8+3*Z2/2+12*M01+9*M11/2)/384,
      0:-(16+3*Z2+7*Z3+27*M21/8+9*M03/8+18*M11+(36+9*Z2/2-27*L**2/8)*M01)/384}
    contract_path=BASE/'prepared/normalization.json';contract=json.loads(contract_path.read_text(encoding='utf-8'))
    dependencies=[basis_path,ROOT/'benchmarks/2ne/families/single-opposite/basis.md',*files]
    for o,full,jfull,series,expected,pole in zip(orbits,[IA,IO],[JA,JB.subs(H,H0.subs(delta,3*e))],[seriesA,seriesO],[expectedA,expectedO],[4,3]):
        coeff={}
        for power in range(-6,1):
            value=s.simplify(series.coeff(e,power));eq(o['orbit_id']+' coefficient '+str(power),value-expected.get(power,0))
            coeff[str(power)]={'expression':str(value),'language':'SymPy/str','latex':s.latex(expected.get(power,s.S.Zero)),
                              'hlog_terms':ha.serialize(ha.from_sympy(value,vectors))}
        colored=[]
        for channel in contract['color_channels']:
            phi=s.sympify(channel['coefficient'],locals={'eps':e})
            phi_series=s.series(phi,e,0,pole+1).removeO()
            product=s.series(phi_series*series,e,0,1).removeO().expand()
            cc={}
            for power in range(-6,1):
                value=s.simplify(product.coeff(e,power))
                convolution=sum(s.expand(phi_series).coeff(e,k)*series.coeff(e,power-k) for k in range(max(0,power+pole)+1))
                eq(o['orbit_id']+' '+channel['monomial']+' convolution '+str(power),value-convolution)
                cc[str(power)]={'expression':str(value),'language':'SymPy/str','hlog_terms':ha.serialize(ha.from_sympy(value,vectors))}
            colored.append({'monomial':channel['monomial'],'exact_scalar_multiplier':str(phi),'coefficients':cc,'interpretation':'coefficients of phi(eps)*the same I; coupling and scale remain external'})
        result={'format':'wilson-computational-result/1','orbit_id':o['orbit_id'],'revision':'r001','family':'Pi1-chain-'+('adjacent' if pole==4 else 'opposite'),
          'input':{'package':package.relative_to(ROOT).as_posix(),'package_sha256':sha(package),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},
          'normalization':{'scheme':contract['schema'],'sym':o['I_target']['sym'],'a1':'1','loop_order':3,'Euler_count':1,'external_color':'CF*B1(eps)^2'},
          'J':{'expression':mathematica_code(jfull),'language':'WolframLanguage/InputForm'},
          'I':{'expression':mathematica_code(full),'language':'WolframLanguage/InputForm','sympy':str(full),'latex':s.latex(full)},
          'laurent':{'from_power':-6,'through_power':0,'leading_power':-pole,'remainder':'O(eps)','coefficients':coeff},
          'color_channels':colored,'functions':{'L':'log(z)','Mkn':'existing evaluated M_{k,n}; bound Hlog files in dependencies','H0':mathematica_code(H0)},
          'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in dependencies]+[{'path':contract_path.relative_to(ROOT).as_posix(),'sha256':sha(contract_path)}],
          'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/3; meromorphic near zero','complex_kinematics':'not verified'},
          'verification':{'exact':(HERE/'exact-checks.json').relative_to(ROOT).as_posix(),'numeric':(HERE/'numeric-checks.json').relative_to(ROOT).as_posix()},'formal':'not_run'}
        freeze(BASE/f"results/{o['orbit_id']}/r001.json",result)
    (HERE/'exact-checks.json').write_text(json.dumps({'count':len(checks),'checks':checks,'solved_orbits':2,'I_positions':14,'I_nonzero':9,'color_positions':42,'color_nonzero':27,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)} exact checks; 2 I results / 14 slots / 9 nonzero; 42 color-channel slots')
if __name__=='__main__':main()
