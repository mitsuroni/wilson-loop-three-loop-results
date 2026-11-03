from pathlib import Path
from fractions import Fraction
import sys,json,hashlib
import sympy as s
from sympy.printing.mathematica import mathematica_code
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha
e,z,L=s.symbols('eps z L',positive=True);checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def eq(n,v):
    assert s.simplify(s.expand(v))==0,(n,v)
    checks.append({'name':n,'passed':True})
def main():
    package=BASE/'prepared/4e-Pi1.json';manifest=BASE/'prepared/manifest.json';meta=json.loads(manifest.read_text(encoding='utf-8'))
    for d in meta['inputs']+meta['outputs']+[meta['normalization'],meta['preparer'],meta['color_backend']]:assert sha(ROOT/d['path'])==d['sha256']
    selected=[o for o in json.loads(package.read_text(encoding='utf-8'))['orbits'] if o['orbit_id'] in ['4e-Pi1-04','4e-Pi1-09']]
    symbols={k:s.Symbol(k,positive=True) for k in ['s21','s22','s31','s41','s42']};symbols.update(eps=e,z1=z)
    parse=lambda t:s.sympify(t.replace('^','**'),locals=symbols)
    for o in selected:
        j=o['J'];se,bare=j['geometry_polynomials'];f=parse(se['f']);d=parse(se['d'])
        eq(o['orbit_id']+' complete numerator',parse(j['se_numerator'])-2*d*f)
        eq(o['orbit_id']+' complete kernel',parse(j['integrand'])-parse(j['se_numerator'])*f**(-2+2*e)/(2*(1+e))*parse(bare['d'])*parse(bare['f'])**(-1+e))
        assert o['I_target']['sym']==4 and o['I_target']['a1']=='-1/2'
        expected={'s21':['0','1'],'s31':['0','1'],'s41':['0','s42'],'s42':['0','1']} if o['orbit_id'].endswith('04') else {'s21':['0','s22'],'s22':['0','1'],'s31':['0','1'],'s41':['0','1']}
        assert j['bounds']==expected
    t=s.Symbol('t',positive=True)
    eq('04 upper SE primitive',s.diff((1-t**(2*e))/(2*e),t)+t**(2*e-1))
    eq('09 lower SE primitive',s.diff((1-(1-t)**(2*e))/(2*e),t)-(1-t)**(2*e-1))
    A,B,C=s.symbols('H0 H2 H2cross');G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
    j4=z**(5*e/2)*(A-B)/(4*e**2*(1+e));j9=z**(e/2)*(A-C)/(4*e**2*(1+e))
    crossing=lambda f:f.xreplace({z:1/z,B:C,C:B})
    eq('raw J reflection',s.powsimp(j9-z**(3*e)*crossing(j4),force=True))
    numerator=s.cosh(e*s.log(z))*(A-(B+C)/2)-s.sinh(e*s.log(z))*(B-C)/2
    compact=G**3*numerator/(128*e**3*(1-2*e)*(1+e)*s.cosh(3*e*s.log(z)/2))
    for key,j in [('04',j4),('09',j9)]:
        assembled=G**3/(32*e*(1-2*e))*(j+z**(3*e)*crossing(j))/(1+z**(3*e))
        eq(key+' normalized I',(assembled-compact).rewrite(s.exp))
    def H(c,d,zz):return -(1+zz)*zz**(-e/2)*s.meijerg([e,0,-c*e],[],[0,e,(d+1)*e],[],zz)/(s.gamma(1-e)*s.gamma(1+e)*s.gamma(1+(c+d+1)*e))
    full=compact.subs({A:H(0,0,z),B:H(2,0,z),C:H(2,0,1/z)})
    basis_path=ROOT/'benchmarks/2ne/families/single-opposite/basis-functions.json';basis=json.loads(basis_path.read_text(encoding='utf-8'))
    vector=lambda key:{(tuple(v['word']),tuple(v['zeta_powers'])):Fraction(v['coefficient']) for v in basis[key]['hlog_terms']}
    assert ha.add(vector('tau0q0'),ha.scale(vector('tau2q0'),-1))=={} and vector('eta2q0')=={}
    names=['tau0q1','tau0q2','tau0q3','tau2q1','tau2q2','tau2q3','eta2q1','eta2q2']
    syms={n:s.Symbol(n) for n in names};D=[syms['tau0q'+str(k)]-syms['tau2q'+str(k)] for k in [1,2,3]]
    E1,E2=syms['eta2q1'],syms['eta2q2']
    num=s.cosh(e*L)*sum(D[k-1]*e**k for k in [1,2,3])-s.sinh(e*L)/L*(e*E1+e**2*E2)
    series=s.series((1+3*s.zeta(2)*e**2/2)*num/(128*e**3*(1-2*e)*(1+e)*s.cosh(3*e*L/2)),e,0,1).removeO().expand()
    expected={-2:D[0]/128,-1:(D[0]+D[1]-E1)/128,0:(D[2]+D[1]-E1-E2+(3+3*s.zeta(2)/2-5*L**2/8)*D[0])/128}
    vectors={syms[n]:vector(n) for n in names};sample=next(iter(vector('tau0q0')));vectors[L]={((0,),tuple(0 for _ in sample[1])):Fraction(1)}
    coeff={}
    for power in range(-6,1):
        value=series.coeff(e,power);eq('coefficient '+str(power),value-expected.get(power,0))
        coeff[str(power)]={'expression':str(value),'language':'SymPy/str','latex':s.latex(expected.get(power,s.S.Zero)),'hlog_terms':ha.serialize(ha.from_sympy(value,vectors))}
    for o in selected:
        result={'format':'wilson-computational-result/1','orbit_id':o['orbit_id'],'revision':'r001','family':'Pi1-P4',
          'input':{'package':package.relative_to(ROOT).as_posix(),'package_sha256':sha(package),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},
          'normalization':{'scheme':'pi1-color-euler-normalization-1.0','loop_order':3,'sym':4,'a1':'-1/2','Euler_count':1,'external_color':'CA*CF*((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)'},
          'I':{'language':'WolframLanguage/InputForm','expression':mathematica_code(full),'sympy':str(full)},
          'laurent':{'from_power':-6,'through_power':0,'leading_power':-2,'remainder':'O(eps)','coefficients':coeff},
          'functions':{'L':'log(z)','tau0qk':'[eps^k]H0','tau2qk':'[eps^k](H2(z)+H2(1/z))/2','eta2qk':'log(z)*[eps^k](H2(z)-H2(1/z))/2'},
          'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in [basis_path,ROOT/'benchmarks/2ne/families/single-opposite/basis.md']],
          'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/2; meromorphic near zero'},'formal':'not_run'}
        p=BASE/f"results/{o['orbit_id']}/r001.json";p.parent.mkdir(parents=True,exist_ok=True);data=json.dumps(result,indent=2,ensure_ascii=False)+'\n'
        if p.exists():assert p.read_text(encoding='utf-8')==data
        else:p.write_text(data,encoding='utf-8')
    (HERE/'exact-checks.json').write_text(json.dumps({'count':len(checks),'checks':checks,'constant_Hlog_cancellation':True,'positions':14,'nonzero':6,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)} exact checks; two I equal; double poles; 14 slots / 6 nonzero')
if __name__=='__main__':main()
