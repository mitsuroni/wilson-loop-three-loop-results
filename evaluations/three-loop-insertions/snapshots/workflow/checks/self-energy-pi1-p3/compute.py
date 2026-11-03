from pathlib import Path
import json,sys,hashlib
import sympy as s
from sympy.printing.mathematica import mathematica_code
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi1'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha
e,z,L=s.symbols('eps z L',positive=True)
M={(k,n):s.Symbol(f'M{k}{n}') for k in range(4) for n in range(1,5-k)}
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def eq(name,value):
    assert s.simplify(s.expand(value))==0,(name,value)
    checks.append({'name':name,'passed':True})
def hseries(c):
    W,T=s.symbols('W T');avg=(c+1)*W/2+c*(T-L)/2
    out=0
    for q in range(4):
        poly=s.Poly(s.expand(((avg+T/2)**(q+1)-(avg-T/2)**(q+1))/s.factorial(q+1)),W,T)
        value=sum(co*M[k,n] for (k,n),co in poly.terms())
        out+=(2*e)**q*value
    return s.expand(out)
def cross(f):return f.xreplace({L:-L,**{m:(-1)**(n+1)*m for (k,n),m in M.items()}})
def main():
    package=BASE/'prepared/4e-Pi1.json';manifest=BASE/'prepared/manifest.json';meta=json.loads(manifest.read_text(encoding='utf-8'))
    for d in meta['inputs']+meta['outputs']+[meta['normalization'],meta['preparer'],meta['color_backend']]:assert sha(ROOT/d['path'])==d['sha256']
    selected=[o for o in json.loads(package.read_text(encoding='utf-8'))['orbits'] if o['orbit_id'] in ['4e-Pi1-03','4e-Pi1-10']]
    x,y=s.symbols('x y',positive=True);D=(1-x)*(1-y)+z*x*y;N=-2*(1+z)*D+(1-2*e)*z
    for o in selected:
        j=o['J'];local={'eps':e,'z1':z,'s21':x if o['orbit_id'].endswith('03') else s.Symbol('a',positive=True),'s22':x,'s41':y,'s31':s.Symbol('b',positive=True),'s42':s.Symbol('t',positive=True)}
        parse=lambda t:s.sympify(t.replace('^','**'),locals=local)
        eq(o['orbit_id']+' SE distance',parse(j['geometry_polynomials'][0]['f'])-D)
        eq(o['orbit_id']+' numerator',parse(j['se_numerator'])-N)
        bare=j['geometry_polynomials'][1]
        eq(o['orbit_id']+' complete source kernel',parse(j['integrand'])-N*D**(-2+2*e)/(2*(1+e))*parse(bare['d'])*parse(bare['f'])**(-1+e))
        assert o['I_target']['sym']==4 and o['I_target']['a1']=='-1/2'
        expected={'s21':['0','1'],'s31':['0','1'],'s41':['0','s42'],'s42':['0','1']} if o['orbit_id'].endswith('03') else {'s21':['0','s22'],'s22':['0','1'],'s31':['0','1'],'s41':['0','1']}
        assert j['bounds']==expected
    mixed=2*e*(s.diff(D,x,y)*D+(2*e-1)*s.diff(D,x)*s.diff(D,y))
    eq('complete-kernel derivative identity',N/(2*(1+e))+(1+z)*D-mixed/(4*e*(1+e)))
    Q=s.gamma(1+e)*s.gamma(1+2*e)/s.gamma(1+3*e)
    boundary=(z**(2*e)+3*(1-Q))/(12*e*(1+e))
    integrated=(z**(2*e)*(1/(2*e)-1/(3*e))+1/(2*e)-s.gamma(1+e)*s.gamma(2*e)/s.gamma(1+3*e))/(2*(1+e))
    eq('all derivative boundary terms',integrated-boundary)
    G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
    A,B,C=s.symbols('H0_2eps Hhalf_2eps Hhalf_cross')
    j3=z**(2*e)*(A-B)/e**2+z**e*boundary/e**2
    jcross=j3.xreplace({z:1/z,B:C,C:B})
    bulk=G**3*(s.cosh(e*s.log(z)/2)*(A-(B+C)/2)-s.sinh(e*s.log(z)/2)*(B-C)/2)/(32*e**3*(1-2*e)*s.cosh(3*e*s.log(z)/2))
    edge=G**3*(1+3*(1-Q)*s.cosh(e*s.log(z)/2)/s.cosh(3*e*s.log(z)/2))/(384*e**4*(1-2*e)*(1+e))
    eq('complete C4 and boundary assembly',(G**3*(j3+z**(3*e)*jcross)/(32*e*(1-2*e)*(1+z**(3*e)))-bulk-edge).rewrite(s.exp))
    def H(c,zz):return -(1+zz)*zz**(-e)*s.meijerg([2*e,0,-2*c*e],[],[0,2*e,2*e],[],zz)/(s.gamma(1-2*e)*s.gamma(1+2*e)*s.gamma(1+2*(c+1)*e))
    full=(bulk+edge).subs({A:H(0,z),B:H(s.Rational(1,2),z),C:H(s.Rational(1,2),1/z)})
    h0=hseries(0);hp=hseries(s.Rational(1,2));hm=cross(hp)
    gs=s.series(s.exp(sum(3*s.zeta(k)*e**k/k for k in range(2,5))),e,0,5).removeO()
    qs=s.series(s.exp(sum((-1)**k*s.zeta(k)*(1+2**k-3**k)*e**k/k for k in range(2,5))),e,0,5).removeO()
    bs=gs*(s.cosh(e*L/2)*(h0-(hp+hm)/2)-s.sinh(e*L/2)*(hp-hm)/2)/(32*e**3*(1-2*e)*s.cosh(3*e*L/2))
    es=gs*(1+3*(1-qs)*s.cosh(e*L/2)/s.cosh(3*e*L/2))/(384*e**4*(1-2*e)*(1+e))
    series=s.series(bs+es,e,0,1).removeO().expand()
    run=ROOT/json.loads((ROOT/'workflow/checks/single-opposite/current-moments.json').read_text())['run']
    files=[run/f'M_{k}_{n}.txt' for k,n in M];vectors={M[k,n]:ha.parse(run/f'M_{k}_{n}.txt') for k,n in M};vectors[L]=ha.LOG
    basis=json.loads((ROOT/'benchmarks/2ne/families/single-opposite/basis-functions.json').read_text())
    from fractions import Fraction
    vector=lambda key:{(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in basis[key]['hlog_terms']}
    # Check the general moment formula against all existing H0 and H2 even/odd coefficients through q=3.
    h2=hseries(s.Integer(2))
    for q in range(4):
        for name,value in [('tau0q'+str(q),h0.coeff(e,q)/2**q),('tau2q'+str(q),s.expand((h2+cross(h2))/2).coeff(e,q)/2**q),('eta2q'+str(q),L*s.expand((h2-cross(h2))/2).coeff(e,q)/2**q)]:
            assert ha.add(ha.from_sympy(value,vectors),ha.scale(vector(name),-1))=={},name
            checks.append({'name':'general M moment matches '+name,'passed':True})
    coeff={}
    for power in range(-6,1):
        value=s.factor(series.coeff(e,power));eq('crossing coefficient '+str(power),value-cross(value))
        coeff[str(power)]={'expression':str(value),'language':'SymPy/str','latex':s.latex(value),'hlog_terms':ha.serialize(ha.from_sympy(value,vectors))}
        print(power,value)
    dependencies=[ROOT/'benchmarks/2ne/families/single-opposite/basis.md',ROOT/'benchmarks/2ne/families/single-opposite/basis-functions.json',*files]
    for o in selected:
        result={'format':'wilson-computational-result/1','orbit_id':o['orbit_id'],'revision':'r001','family':'Pi1-P3','input':{'package':package.relative_to(ROOT).as_posix(),'package_sha256':sha(package),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},'normalization':{'scheme':'pi1-color-euler-normalization-1.0','loop_order':3,'sym':4,'a1':'-1/2','Euler_count':1,'external_color':'CA*CF*((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)'},'I':{'expression':mathematica_code(full),'language':'WolframLanguage/InputForm','sympy':str(full)},'laurent':{'from_power':-6,'through_power':0,'leading_power':-4,'remainder':'O(eps)','coefficients':coeff},'functions':{'L':'log(z)','Mkn':'existing M_{k,n}; definitions and evaluated files in dependencies'},'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in dependencies],'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/2; meromorphic near zero'},'formal':'not_run'}
        p=BASE/f"results/{o['orbit_id']}/r001.json";p.parent.mkdir(parents=True,exist_ok=True);data=json.dumps(result,indent=2,ensure_ascii=False)+'\n'
        if p.exists():assert p.read_text(encoding='utf-8')==data
        else:p.write_text(data,encoding='utf-8')
    (HERE/'exact-checks.json').write_text(json.dumps({'count':len(checks),'checks':checks,'positions':14,'nonzero':10,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)} exact checks; 03/10 complete; fourth-order poles')
if __name__=='__main__':main()
