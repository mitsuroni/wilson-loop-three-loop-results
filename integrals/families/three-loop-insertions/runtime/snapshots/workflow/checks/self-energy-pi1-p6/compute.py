"""Pi1 P6: source factorization, complete normalized I and evaluated Hlog series."""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib
import sympy as s
from sympy.printing.mathematica import mathematica_code
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha
e,z=s.symbols('eps z',positive=True)
x,y,u,v=s.symbols('x y u v',positive=True)
T0,T1,T2,L=s.symbols('T0 T1 T2 L')
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def eq(name,value):
    assert s.simplify(s.expand(value))==0,(name,value)
    checks.append({'name':name,'passed':True})
def main():
    package=BASE/'prepared/4e-Pi1.json';manifest=BASE/'prepared/manifest.json'
    meta=json.loads(manifest.read_text(encoding='utf-8'))
    for rec in meta['inputs']+meta['outputs']+[meta['normalization'],meta['preparer'],meta['color_backend']]:assert sha(ROOT/rec['path'])==rec['sha256']
    o=json.loads(package.read_text(encoding='utf-8'))['orbits'][-1];assert o['orbit_id']=='4e-Pi1-11'
    assert o['I_target']['a1']=='-1/2' and o['I_target']['sym']==2
    j=o['J'];assert j['iterated_order_outer_to_inner']==['s11','s21','s31','s41']
    assert all(b==['0','1'] for b in j['bounds'].values())
    parse=lambda t:s.sympify(t.replace('^','**'),locals={'eps':e,'z1':z,'s11':x,'s31':y,'s21':u,'s41':v})
    F=z*(1-x)*(1-y)+x*y;D=(1-u)*(1-v)+z*u*v
    N=-2*(1+z)*F+(1-2*e)*z
    eq('SE distance',parse(j['geometry_polynomials'][0]['f'])-F)
    eq('bare distance',parse(j['geometry_polynomials'][1]['f'])-D)
    eq('complete source numerator',parse(j['se_numerator'])-N)
    eq('actual source product',parse(j['integrand'])-N*F**(-2+2*e)/(2*(1+e))*(-1-z)*D**(-1+e))
    eq('SE crossing distance',F-z*((1-x)*(1-y)+x*y/z))
    eq('SE crossing numerator',N-z**2*(-2*(1+1/z)*((1-x)*(1-y)+x*y/z)+(1-2*e)/z))
    A,B=s.symbols('H0_eps H0_2eps')
    J=z**(3*e/2)*A*B+z**(e/2)*A*(1+z**(2*e))/(4*e*(1+e))
    G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
    symbolic=G**3*A*(B+s.cosh(e*s.log(z))/(2*e*(1+e)))/(64*e*(1-2*e)*s.cosh(3*e*s.log(z)/2))
    assembled=G**3/(64*e*(1-2*e))*(J+z**(3*e)*J.subs(z,1/z))/(1+z**(3*e))
    eq('complete normalization and sym=2',s.powsimp((assembled/symbolic).rewrite(s.exp),force=True)-1)
    eq('J self-crossing',s.powsimp(J-z**(3*e)*J.subs(z,1/z),force=True))
    delta=s.Symbol('delta')
    H=-(1+z)*z**(-delta/2)*s.meijerg([delta,0,0],[],[0,delta,delta],[],z)/(s.gamma(1-delta)*s.gamma(1+delta)**2)
    full=symbolic.subs({A:H.subs(delta,e),B:H.subs(delta,2*e)})
    regular=1+s.Rational(3,2)*s.zeta(2)*e**2
    series=s.series(regular*(T0+e*T1+e**2*T2)*(T0+2*e*T1+s.cosh(e*L)/(2*e*(1+e)))/(64*e*(1-2*e)*s.cosh(3*e*L/2)),e,0,1).removeO().expand()
    expected={-2:T0/128,-1:(T0+T1+2*T0**2)/128,0:(T2+T1+(3+s.Rational(3,2)*s.zeta(2)-5*L**2/8)*T0+4*T0**2+6*T0*T1)/128}
    basis_path=ROOT/'benchmarks/2ne/families/single-opposite/basis-functions.json';basis=json.loads(basis_path.read_text(encoding='utf-8'))
    vector=lambda key:{(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in basis[key]['hlog_terms']}
    vectors={T0:vector('tau0q0'),T1:vector('tau0q1'),T2:vector('tau0q2'),L:{((0,),(0,0,0,0,0,0)):Fraction(1)}}
    # Match log vector's zeta-coordinate length to the established algebra.
    sample=next(iter(vectors[T0]));vectors[L]={((0,),tuple(0 for _ in sample[1])):Fraction(1)}
    coeff={}
    for p in range(-6,1):
        value=series.coeff(e,p);eq('coefficient '+str(p),value-expected.get(p,0))
        hlog=ha.from_sympy(value,vectors)
        coeff[str(p)]={'expression':str(s.simplify(value)),'language':'SymPy/str','latex':s.latex(expected.get(p,s.S.Zero)), 'hlog_terms':ha.serialize(hlog)}
    dependencies=[basis_path,ROOT/'benchmarks/2ne/families/single-opposite/basis.md',BASE/'results/2e-Pi1-02/r001.json']
    result={'format':'wilson-computational-result/1','orbit_id':'4e-Pi1-11','revision':'r001','family':'Pi1-P6',
      'input':{'package':package.relative_to(ROOT).as_posix(),'package_sha256':sha(package),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},
      'normalization':{'scheme':'pi1-color-euler-normalization-1.0','loop_order':3,'a1':'-1/2','sym':2,'Euler_count':1,'external_color':'CA*CF*((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)'},
      'I':{'expression':mathematica_code(full),'language':'WolframLanguage/InputForm','sympy':str(full),'latex':s.latex(full)},
      'laurent':{'from_power':-6,'through_power':0,'leading_power':-2,'remainder':'O(eps)','coefficients':coeff},
      'functions':{'T0':'tau0q0 = M01','T1':'tau0q1 = M11/2','T2':'tau0q2 = [delta^2]H0(delta,z)','L':'log(z)','H0':mathematica_code(H)},
      'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in dependencies],
      'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/2; meromorphic near zero','complex_kinematics':'not verified'},'formal':'not_run'}
    p=BASE/'results/4e-Pi1-11/r001.json';p.parent.mkdir(parents=True,exist_ok=True);data=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if p.exists():assert p.read_text(encoding='utf-8')==data
    else:p.write_text(data,encoding='utf-8')
    (HERE/'exact-checks.json').write_text(json.dumps({'count':len(checks),'checks':checks,'source_domain_and_sym_checked':True,'new_positions':7,'nonzero':3,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)} exact checks; P6 double pole; 7 slots / 3 nonzero')
if __name__=='__main__':main()
