"""Source-bound Pi1 P1/P2 reduction and Laurent coefficients; no formal proof."""
from pathlib import Path
import json, hashlib
import sympy as s
from sympy.printing.mathematica import mathematica_code

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'benchmarks/self-energy/pi1'
e,z,L=s.symbols('eps z L',positive=True)
a,b,c,d=s.symbols('s31 s32 s41 s42',positive=True)
x=s.Symbol('s21',positive=True)
G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
Q=s.gamma(1+e)*s.gamma(1+2*e)/s.gamma(1+3*e)
checks=[]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,ok):
    assert ok,name
    checks.append({'name':name,'passed':True})
def eq(name,v):check(name,s.simplify(s.expand(v))==0)
def freeze(p,obj):
    data=(json.dumps(obj,ensure_ascii=False,indent=2)+'\n').encode('utf-8')
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==data, p
    else:p.write_bytes(data)
def parse(t):return s.sympify(t.replace('^','**'),locals={'eps':e,'z1':z,'s21':x,'s31':a,'s32':b,'s41':c,'s42':d})

def main():
    package=BASE/'prepared/4e-Pi1.json';manifest=BASE/'prepared/manifest.json'
    meta=json.loads(manifest.read_text(encoding='utf-8'))
    for item in meta['outputs']+[meta['normalization'],meta['preparer'],meta['color_backend']]:
        check('bound source '+item['path'],sha(ROOT/item['path'])==item['sha256'])
    selected={o['orbit_id'][-2:]:o for o in json.loads(package.read_text(encoding='utf-8'))['orbits'] if o['orbit_id'][-2:] in ['01','02','05','06']}
    simple={
      '01':z**(3*e)/(1+e)*(1-a)**(2*e-1)*c**(2*e-1)*(1-b)**(e-1)*d**(e-1),
      '02':z**(3*e)/(1+e)*(1-b)**(2*e-1)*d**(2*e-1)*(1-a)**(e-1)*c**(e-1),
      '05':z**e/(1+e)*(1-x)**(2*e-1)*b**(2*e-1)*(1-a)**(e-1)*c**(e-1),
      '06':z**(2*e)/(1+e)*(1-a)**(2*e-1)*c**(2*e-1)*(1-x)**(e-1)*b**(e-1)}
    for key,o in selected.items():
        j=o['J'];se,bare=j['geometry_polynomials']
        f=parse(se['f']);f0=parse(bare['f']);ds=parse(se['d']);db=parse(bare['d'])
        eq(key+' full SE numerator',parse(j['se_numerator'])-2*ds*f)
        expected=parse(j['se_numerator'])*f**(-2+2*e)/(2*(1+e))*db*f0**(-1+e)
        eq(key+' actual integrand assembly',parse(j['integrand'])-expected)
        # Positive factor coordinates: 1-s > 0 is justified by the open source domain.
        A,B,X=s.symbols('A B X',positive=True)
        sub={a:1-A,b:1-B,x:1-X}
        ratio=s.powsimp(s.expand_power_base(s.factor(expected.subs(sub))/simple[key].subs(sub),force=True),force=True)
        eq(key+' positive-domain monomial kernel',ratio-1)
        target=o['I_target'];check(key+' normalization',target['sym']==4 and target['a1']=='-1/2' and target['total_loop_order']==3)
        bounds={'s31':['0','s32'],'s32':['0','1'],'s41':['0','s42'],'s42':['0','1']} if key in ['01','02'] else {'s21':['0','1'],'s31':['0','s32'],'s32':['0','1'],'s41':['0','1']}
        check(key+' bounds',j['bounds']==bounds)
        order=['s32','s31','s42','s41'] if key in ['01','02'] else ['s21','s32','s31','s41']
        check(key+' iterated order',j['iterated_order_outer_to_inner']==order)
    # Elementary primitives; endpoint powers vanish for positive real p,q.
    p,q,t=s.symbols('p q t',positive=True)
    eq('BB inner primitive',s.diff(t**p/p,t)-t**(p-1))
    eq('AA inner primitive',s.diff((1-(1-t)**p)/p,t)-(1-t)**(p-1))
    AA=1/(q*(p+q));BB=1/(p*(p+q))
    eq('AA integrated moments',(1/q-1/(p+q))/p-AA)
    eq('BB integrated moments',1/p/(p+q)-BB)
    eq('P1 ordered factors',AA.subs({p:2*e,q:e})*BB.subs({p:2*e,q:e})-1/(18*e**4))
    eq('P1 reversed ordered factors',AA.subs({p:e,q:2*e})*BB.subs({p:e,q:2*e})-1/(18*e**4))
    qgeneric=s.gamma(1+p)*s.gamma(1+q)/s.gamma(1+p+q)
    eq('BA Beta to Gamma',(1/q-s.gamma(q)*s.gamma(1+p)/s.gamma(1+p+q))/p-(1-qgeneric)/(p*q))
    eq('BA exponent interchange',qgeneric-qgeneric.xreplace({p:q,q:p}))
    J1=z**(3*e)/(18*e**4*(1+e))
    F2=(1-Q)/(4*e**4*(1+e))
    Js={'01':J1,'02':J1,'05':z**e*F2,'06':z**(2*e)*F2}
    I1=G**3/(576*e**5*(1-2*e)*(1+e))
    I2=G**3*(1-Q)*s.cosh(e*s.log(z)/2)/(128*e**5*(1-2*e)*(1+e)*s.cosh(3*e*s.log(z)/2))
    for key,j in Js.items():
        assembled=G**3/(32*e*(1-2*e))*(j+z**(3*e)*j.subs(z,1/z))/(1+z**(3*e))
        goal=I1 if key in ['01','02'] else I2
        eq(key+' complete C4 assembly',s.powsimp((assembled/goal).rewrite(s.exp),force=True)-1)
    # Analytic log-Gamma expansion with a controlled Taylor remainder.
    glog=sum(3*s.zeta(k)*e**k/k for k in range(2,6))
    qlog=sum((-1)**k*s.zeta(k)*(1+2**k-3**k)*e**k/k for k in range(2,6))
    trunc=lambda v:s.series(v,e,0,6).removeO().expand()
    gp=trunc(s.exp(glog));qp=trunc(s.exp(qlog))
    rp=trunc(1/((1-2*e)*(1+e)))
    cp=trunc(s.cosh(e*L/2)/s.cosh(3*e*L/2))
    ser1=trunc(gp*rp)/576/e**5
    ser2=trunc(gp*(1-qp)*rp*cp)/128/e**5
    Z2,Z3,Z4,Z5=s.zeta(2),s.zeta(3),s.zeta(4),s.zeta(5)
    expected1={-5:s.Rational(1,576),-4:s.Rational(1,576),-3:(3+3*Z2/2)/576,-2:(5+3*Z2/2+Z3)/576,
      -1:(11+9*Z2/2+Z3+57*Z4/16)/576,
      0:(21+15*Z2/2+3*Z3+57*Z4/16+3*Z2*Z3/2+3*Z5/5)/576}
    expected2={-3:Z2/64,-2:(2*Z2-6*Z3)/128,
      -1:(2*(3-L**2)*Z2-6*Z3+s.Rational(37,2)*Z4)/128,
      0:(2*(5-L**2)*Z2-6*(3-L**2)*Z3+s.Rational(37,2)*Z4+5*Z2*Z3-42*Z5)/128}
    for key,o in selected.items():
        full,ser,expected,pole=(I1,ser1,expected1,5) if key in ['01','02'] else (I2,ser2,expected2,3)
        coeff={}
        for power in range(-6,1):
            v=s.expand(ser).coeff(e,power);eq(key+' Laurent '+str(power),v-expected.get(power,0))
            coeff[str(power)]={'expression':str(s.simplify(v.subs(L,s.log(z)))),'language':'SymPy/str','latex':s.latex(expected.get(power,s.S.Zero))}
        item={'format':'wilson-computational-result/1','family':'Pi1-P1' if pole==5 else 'Pi1-P2','orbit_id':o['orbit_id'],'revision':'r001',
          'input':{'package':package.relative_to(ROOT).as_posix(),'package_sha256':sha(package),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},
          'normalization':{'scheme':'pi1-color-euler-normalization-1.0','loop_order':3,'sym':4,'a1':'-1/2','G':'exp(-EulerGamma*eps)*Gamma(1-eps)','external_color':'CA*CF*((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)','Euler_count':1},
          'J':{'expression':str(Js[key]),'language':'SymPy/str'},
          'I':{'kind':'Gamma_closed_form','expression':mathematica_code(full),'language':'WolframLanguage/InputForm','sympy':str(full),'latex':s.latex(full)},
          'laurent':{'variable':'eps','from_power':-6,'through_power':0,'leading_power':-pole,'remainder':'O(eps)','coefficients':coeff},
          'domain':{'z':'positive real','epsilon':'initial 0<Re(eps)<1/2; meromorphic near zero','complex_kinematics':'not verified'},
          'verification':{'exact':(HERE/'exact-checks.json').relative_to(ROOT).as_posix(),'numerical':(HERE/'numeric-checks.json').relative_to(ROOT).as_posix()},'formal':'not_run'}
        freeze(BASE/f"results/{o['orbit_id']}/r001.json",item)
    (HERE/'exact-checks.json').write_text(json.dumps({'checks':checks,'count':len(checks),'new_results':4,'positions':28,'nonzero':20,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)} exact checks; four results; 28 positions / 20 nonzero')

if __name__=='__main__':main()
