"""Independent source/operator and scoped finite-coefficient controls for P5.

Imports no P5 evaluator and reads no historical final-result table. The
polynomial tests concern J, not the singular Gamma prefactor at integer eps.
"""
from __future__ import annotations
import argparse
import ast
import copy
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sp

PACK=Path(__file__).resolve().parents[1]
SNAP=PACK/'snapshots'
EPS=sp.Symbol('eps');Z=sp.Symbol('z1',positive=True)
NAMES={'eps':EPS,'z':Z,'z1':Z,'Pi':sp.pi,'pi':sp.pi,'EulerGamma':sp.EulerGamma}
NAMES.update({n:sp.Symbol(n) for n in ['CA','CF','TF','nf','alpha_s','v14','I0']})
ONE={((),(0,0,0)):Q(1)}
PATHS={'pi1':'benchmarks/current/prepared/3loop/4e-Pi1.json',
       'chain':'benchmarks/current/prepared/3loop/2e-Pi1+Pi1.json',
       'pi2':'benchmarks/current/prepared/3loop/2e-Pi2.json'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def require(v,m):
    if not v:raise ValueError(m)

def parse(text,names=None):
    names={**NAMES,**(names or {})}
    def visit(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return sp.Integer(n.value)
        if isinstance(n,ast.Name) and n.id in names:return names[n.id]
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):
            return visit(n.operand) if isinstance(n.op,ast.UAdd) else -visit(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=visit(n.left),visit(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
            if isinstance(n.op,ast.Pow):return a**b
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and not n.keywords:
            functions={'gamma':sp.gamma,'Gamma':sp.gamma,'exp':sp.exp,'Exp':sp.exp,'zeta':sp.zeta,'log':sp.log}
            if n.func.id in functions:return functions[n.func.id](*(visit(a) for a in n.args))
        raise ValueError('Input outside declared scalar grammar: '+ast.dump(n))
    return visit(ast.parse(text.replace('^','**').replace('[','(').replace(']',')'),mode='eval').body)

def source_check(orbit):
    j=orbit['J'];order=j['iterated_order_outer_to_inner']
    require(set(order)==set(j['variables'])==set(j['bounds']),'Variable/order inventory')
    names={s:sp.Symbol(s,positive=True) for s in j['variables']}
    available=set()
    for s in order:
        lo,hi=(parse(v,names) for v in j['bounds'][s])
        require(lo==0 and hi.free_symbols<={names[t] for t in available},'Unavailable dependent bound')
        available.add(s)
    norm=orbit['normalization'];prep=orbit['preparation']
    require(norm['euler_compensation_count']==1 and norm['assembly_sign']==1,'Euler/sign normalization')
    require(prep['native_sign_conversion_count']==0 and prep['Euler_compensation_count']==1,'Unexpected preparation/sign switch')
    require(parse(orbit['C']['components']['euler'])==sp.exp(-3*sp.EulerGamma*EPS),'Wrong Euler weight')
    core=orbit['integral_structure']['core'];operator=core['operator'];se=j['geometry_polynomials'][0]
    ca,tf,nf=(NAMES[k] for k in ['CA','TF','nf'])
    B1=((10-6*EPS)*ca-8*(1-EPS)*tf*nf)/(3-2*EPS)
    if operator in ['Pi1','Pi1+Pi1']:
        expected_color=B1 if operator=='Pi1' else B1**2
        require(sp.cancel(parse(orbit['C']['components']['self_energy_color'])-expected_color)==0,'Actual B1/full chain colour operator')
    elif operator=='Pi2':
        expected_scalar=sum(parse(row['operator_monomial'])*parse(row['p']) for row in orbit['operator']['channels'].values())
        require(sp.expand(parse(orbit['C']['components']['self_energy_scalar'])-expected_scalar)==0,'Full Pi2 scalar/channels connection')
    f,d,u,w=(parse(se[k],names) for k in ['f','d','u','w'])
    numerator=parse(j['se_numerator'],names)
    if operator=='Pi1':
        expected=(1+2*EPS)*d*f+(1-2*EPS)*u*w;rate=2
        metric=parse(j['metric_numerator'],names);longitudinal=parse(j['derivative_numerator'],names)
        require(sp.expand(metric-2*(1+EPS)*d*f)==0,'Metric tensor contraction')
        require(sp.expand(longitudinal-(-d*f+(1-2*EPS)*u*w))==0,'Longitudinal tensor contraction')
        require(sp.expand(numerator-metric-longitudinal)==0,'Dropped longitudinal numerator')
    elif operator in ['Pi1+Pi1','Pi2']:
        expected=(1+4*EPS)*d*f+(1-3*EPS)*u*w;rate=3
        require(core['insertion_count']==(2 if operator=='Pi1+Pi1' else 1),'Insertion count/operator identity')
        if operator=='Pi2':require(orbit['operator']['kind']=='1pi' and orbit['operator']['not_Pi1_squared'],'Pi2 is not a chain')
    else:raise ValueError('Unknown insertion operator')
    require(sp.expand(numerator-expected)==0,'Complete source tensor numerator')
    assembled=numerator*f**(-2+rate*EPS)/(2*(1+(rate-1)*EPS))
    if rate==2:
        bare=j['geometry_polynomials'][1];assembled*=parse(bare['d'],names)*parse(bare['f'],names)**(-1+EPS)
    raw=parse(j['integrand'],names)
    require(sp.simplify(raw-assembled)==0,'Original integrand assembly')
    return raw,names,operator

def integrate_original(expr,j,names):
    expr=sp.expand(expr)
    for name in reversed(j['iterated_order_outer_to_inner']):
        v=names[name];lo,hi=(parse(t,names) for t in j['bounds'][name]);out=0
        for (degree,),coefficient in sp.Poly(expr,v).terms():out+=coefficient*(hi**(degree+1)-lo**(degree+1))/(degree+1)
        expr=sp.expand(out)
    return sp.cancel(expr)

def integrate_cube(expr,j,names):
    order=j['iterated_order_outer_to_inner'];xs=sp.symbols('x0:'+str(len(order)),positive=True);mapping={}
    for x,name in zip(xs,order):mapping[names[name]]=x*parse(j['bounds'][name][1],names).subs(mapping,simultaneous=True)
    determinant=sp.factor(sp.Matrix([mapping[names[n]] for n in order]).jacobian(xs).det())
    require(determinant.is_nonnegative,'Positive cube orientation')
    transformed=sp.Poly(sp.expand(expr.subs(mapping,simultaneous=True)*determinant),*xs)
    integral=sum(coefficient/sp.prod(k+1 for k in degrees) for degrees,coefficient in transformed.terms())
    return sp.cancel(integral),str(determinant)

def polynomial_checks(orbit):
    raw,names,_=source_check(orbit);rows=[]
    for e in [1,2]:
        expr=raw.subs({EPS:e,Z:sp.Rational(1,2)})
        a=integrate_original(expr,orbit['J'],names);b,jac=integrate_cube(expr,orbit['J'],names)
        require(sp.cancel(a-b)==0,'Original/cube integral mismatch '+orbit['orbit_id'])
        rows.append({'eps':e,'z':'1/2','original_order_integral':str(a),'cube_integral':str(b),'absolute_jacobian':jac,'residual':'0'})
    return {'orbit_id':orbit['orbit_id'],'paired_checks':rows,'scope':'Original J only. Gamma prefactor is not evaluated at these regulator values.'}

def vector(rows):
    out={}
    for t in rows:
        key=(tuple(t['word']),tuple(t['zeta_powers']))
        require(all(a in (0,-1) for a in key[0]) and len(key[1])==3,'HLOG coordinate conventions')
        out[key]=out.get(key,Q(0))+Q(t['coefficient'])
    return {k:v for k,v in out.items() if v}

def scalar_vector(expr):
    z2,z3,z5=sp.symbols('_z2 _z3 _z5')
    expr=sp.expand(expr).subs({sp.zeta(3):z3,sp.zeta(5):z5}).subs(sp.pi**2,6*z2)
    poly=sp.Poly(sp.expand(expr),z2,z3,z5)
    out={}
    for powers,c in poly.terms():
        require(c.is_Rational,'Nonperiod scalar coefficient')
        if c:out[((),powers)]=Q(int(c.p),int(c.q))
    return out

def multiply_periods(a,b):
    out={}
    for (w,pa),ca in a.items():
        for (empty,pb),cb in b.items():
            require(not empty,'Convolution multiplier must be a period')
            k=(w,tuple(x+y for x,y in zip(pa,pb)));out[k]=out.get(k,Q(0))+ca*cb
    return {k:v for k,v in out.items() if v}

def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,Q(0))+v
    return {k:v for k,v in out.items() if v}

def geometry_rows(result):
    return result.get('I_geom_laurent',result.get('laurent'))['coefficients']

def color_rows(channel):
    return channel.get('coefficients',channel.get('laurent',{}).get('coefficients'))

def rational_jet(expr,lower,upper):
    numerator,denominator=map(lambda v:sp.Poly(v,EPS),sp.fraction(sp.cancel(expr)))
    npowers={k[0]:v for k,v in numerator.terms() if v};dpowers={k[0]:v for k,v in denominator.terms() if v}
    if not npowers:return {k:sp.S.Zero for k in range(lower,upper+1)}
    na,da=min(npowers),min(dpowers);shift=na-da
    ps={k-na:v for k,v in npowers.items()};qs={k-da:v for k,v in dpowers.items()};coefficients={}
    for k in range(max(0,upper-shift)+1):
        coefficients[k]=sp.cancel((ps.get(k,0)-sum(v*coefficients.get(k-j,0) for j,v in qs.items() if j>0))/qs[0])
    return {k:coefficients.get(k-shift,sp.S.Zero) for k in range(lower,upper+1)}

def multipliers(orbit):
    operator=orbit['integral_structure']['core']['operator']
    if operator=='Pi1':
        return {'CA^2*CF':(10-6*EPS)/(3-2*EPS),'CA*CF*TF*nf':-8*(1-EPS)/(3-2*EPS)}
    if operator=='Pi1+Pi1':
        a=(10-6*EPS)/(3-2*EPS);b=-8*(1-EPS)/(3-2*EPS)
        return {'CA^2*CF':a*a,'CA*CF*TF*nf':2*a*b,'CF*(TF*nf)^2':b*b}
    # U is expanded by the defining Gamma logarithm, independently of a
    # previously saved operator expansion or the evaluator's coefficient code.
    slopes=[(1,2),(-1,1),(-3,1),(-2,-2),(2,-1)]
    loginv={k:-sp.zeta(k)/k*sum(mult*(-slope)**k for slope,mult in slopes) for k in range(2,6)}
    inv={0:sp.S.One}
    for k in range(1,6):inv[k]=sp.expand(sum(j*v*inv[k-j] for j,v in loginv.items() if j<=k)/k)
    out={}
    for name,entry in orbit['operator']['channels'].items():
        a,b=(parse(entry[k]) for k in ['A','B'])
        aa,bb=(rational_jet(v,-3,2) for v in [a,b])
        coefficients={p:sp.expand(4*(aa.get(p-2,0)+sum(bb.get(p-2-k,0)*v for k,v in inv.items()))) for p in range(-1,5)}
        require(coefficients[-1]==0,'Spurious Pi2 operator pole')
        out[entry['full_monomial']]=sum(v*EPS**p for p,v in coefficients.items() if p>=0)
    return out

def convolution_check(orbit,result):
    geometry={int(p):vector(row['hlog_terms']) for p,row in geometry_rows(result).items()}
    require(set(geometry)==set(range(-6,1)),'Seven geometry powers')
    bound=min(p for p,v in geometry.items() if v)
    expected=multipliers(orbit);seen=set();rows=[]
    for channel in result['color_channels']:
        monomial=channel.get('monomial',channel.get('color_monomial'))
        actual=color_rows(channel);require(actual is not None,'Missing full colored Laurent table')
        if monomial not in expected:
            require(all(not vector(row['hlog_terms']) for row in actual.values()),'Unproved nonzero extra color')
            continue
        require(monomial not in seen,'Duplicate color channel');seen.add(monomial)
        phi=sp.series(expected[monomial],EPS,0,1-bound).removeO().expand()
        phi_rows={k:scalar_vector(phi.coeff(EPS,k)) for k in range(0,1-bound)}
        for p in range(-6,1):
            assembled={}
            for k,f in phi_rows.items():assembled=add(assembled,multiply_periods(geometry.get(p-k,{}),f))
            require(assembled==vector(actual[str(p)]['hlog_terms']),'Wrong full-color convolution '+orbit['orbit_id']+monomial+str(p))
        rows.append({'monomial':monomial,'positions':7,'operator_Taylor_through':-bound,'residual':'0'})
    require(seen==set(expected),'Missing nonzero physical color channel')
    return {'orbit_id':orbit['orbit_id'],'geometry_leading_power':bound,'channels':rows,'scope':'Independent exact operator Taylor and HLOG/period convolution, not independent geometry integration.'}

class HLog:
    """Independent log-polynomial Taylor evaluator at 0<z<1."""
    def __init__(self,z,terms=160):self.z=mp.mpf(z);self.terms=terms;self.cache={():{0:[mp.mpf(1)]}}
    def expansion(self,word):
        if word in self.cache:return self.cache[word]
        old=self.expansion(word[1:]);letter=word[0];out={}
        if letter==-1:
            coefficients={};previous=[]
            for n in range(self.terms):
                a=old.get(n,[]);v=[(a[k] if k<len(a) else 0)-(previous[k] if k<len(previous) else 0) for k in range(max(len(a),len(previous)))]
                coefficients[n+1]=v;previous=v
        elif letter==0:coefficients=old
        else:raise ValueError('Unsupported HLOG letter')
        for n,poly in coefficients.items():
            if n>self.terms:continue
            if n==0:out[0]=[mp.mpf(0)]+[v/(k+1) for k,v in enumerate(poly)];continue
            target=[mp.mpf(0)]*len(poly)
            for k,c in enumerate(poly):
                for j in range(k+1):target[k-j]+=c*(-1)**j*mp.factorial(k)/mp.factorial(k-j)/mp.mpf(n)**(j+1)
            out[n]=target
        self.cache[word]=out;return out
    def value(self,word):
        l=mp.log(self.z)
        return mp.fsum(self.z**n*mp.fsum(c*l**k for k,c in enumerate(poly)) for n,poly in self.expansion(tuple(word)).items())
    def coefficient(self,row):
        return mp.fsum(mp.mpf(t['coefficient'])*self.value(t['word'])*mp.fprod(mp.zeta(k)**p for k,p in zip([2,3,5],t['zeta_powers'])) for t in row['hlog_terms'])

def gamma_cauchy(orbits,results):
    evaluator=HLog('.5');z=mp.mpf('.5');records=[]
    require(abs(evaluator.value((0,-1))+mp.polylog(2,-z))<mp.mpf('1e-42'),'HLOG independent evaluator sanity')
    for oid,pole in [('4e-Pi1-01',5),('4e-Pi1-02',5),('4e-Pi1-05',3),('4e-Pi1-06',3),('2e-Pi1+Pi1-01',4),('2e-Pi2-01',4)]:
        orbit=orbits[oid];op=orbit['integral_structure']['core']['operator'];comp=orbit['C']['components']
        rational=parse(comp['rational'])
        if op=='Pi1':
            wilson=parse(comp['wilson_color']);CA,CF=NAMES['CA'],NAMES['CF']
            a,b=sp.symbols('_a _b')
            numerator=sp.fraction(sp.cancel((wilson-a*CA*CF-b*CF**2).subs(CF,(CA**2-1)/(2*CA))))[0]
            solutions=sp.solve(sp.Poly(numerator,CA).all_coeffs(),[a,b],dict=True)
            require(len(solutions)==1 and solutions[0][a].is_Rational,'Source fixed-basis color reconstruction')
            a1=solutions[0][a]
            rational*=sp.Rational(orbit['c4_symmetry']['orbit_multiplicity'],2)*a1
        else:
            rational*=sp.Rational(orbit['c4_symmetry']['orbit_multiplicity'],2)
            if op=='Pi2':rational/=4*EPS**2
        fn=sp.lambdify(EPS,rational,'mpmath')
        def full(e):
            G=mp.exp(-mp.euler*e)*mp.gamma(1-e)
            if op=='Pi1':
                if oid[-2:] in ['01','02']:
                    j=z**(3*e)/(18*e**4*(1+e));jr=z**(-3*e)/(18*e**4*(1+e))
                else:
                    q=mp.gamma(1+e)*mp.gamma(1+2*e)/mp.gamma(1+3*e)
                    degree=1 if oid.endswith('05') else 2
                    j=z**(degree*e)*(1-q)/(4*e**4*(1+e));jr=z**(-degree*e)*(1-q)/(4*e**4*(1+e))
                return fn(e)*G**3*(j+z**(3*e)*jr)/(1+z**(3*e))
            u=mp.gamma(1+e)**2*mp.gamma(1-e)*mp.gamma(1-3*e)/(mp.gamma(1-2*e)**2*mp.gamma(1+2*e))
            return fn(e)*G**3*u*(2+e)/(18*e**2*(1+2*e))
        actual={p:evaluator.coefficient(row) for p,row in ((int(k),v) for k,v in geometry_rows(results[oid]).items())}
        for p in range(-6,-pole):require(not vector(geometry_rows(results[oid])[str(p)]['hlog_terms']),'Spurious pole below Gamma bound')
        max_error=mp.mpf(0);radii=[]
        for radius in [mp.mpf(1)/20,mp.mpf(1)/40]:
            n=96;values=[]
            for j in range(n):
                angle=2*mp.pi*j/n;e=radius*mp.exp(1j*angle);values.append(e**pole*full(e))
            errors=[]
            for p in range(-pole,1):
                k=p+pole;c=mp.fsum(values[j]*mp.exp(-2j*mp.pi*j*k/n) for j in range(n))/n/radius**k
                error=abs(c-actual[p]);errors.append(error);max_error=max(max_error,error)
            radii.append({'radius':str(radius),'nodes':n,'max_observed_residual':mp.nstr(max(errors),12)})
        require(max_error<mp.mpf('1e-34'),'Gamma finite-coefficient control failed '+oid+' '+str(max_error))
        records.append({'orbit_id':oid,'positions':7,'z':'1/2','source_C_prefactor':str(rational),'radii':radii,
            'scope':'Independent Cauchy extraction of source-derived adjacent Gamma formula versus separately evaluated HLOG coefficients. No certified numerical enclosure.'})
    return records

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--results',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    require(not args.out.exists(),'Choose a NEW output');started=time.perf_counter();mp.mp.dps=75
    source_files={k:SNAP/v for k,v in PATHS.items()};orbits={o['orbit_id']:o for p in source_files.values() for o in read(p)['orbits']}
    results=read(args.results);require(set(orbits)==set(results) and len(orbits)==15,'Fifteen-orbit inventory')
    polynomial=[polynomial_checks(o) for o in orbits.values()]
    convolution=[convolution_check(o,results[oid]) for oid,o in orbits.items()]
    cauchy=gamma_cauchy(orbits,results)
    # Whole-source tensors imply the canonical opposite derivative relation;
    # test both rates without evaluating endpoint distributions numerically.
    x,y=sp.symbols('x y');D=(1-x)*(1-y)+Z*x*y;identities=[]
    for rate in [2,3]:
        N=(-(2 if rate==2 else 2+EPS)*(1+Z)*D+(1-rate*EPS)*Z)
        differentiated=rate*EPS*(sp.diff(D,x,y)*D+(rate*EPS-1)*sp.diff(D,x)*sp.diff(D,y))
        residual=sp.expand(N/(2*(1+(rate-1)*EPS))+(1+Z)*D-differentiated/(2*rate*EPS*(1+(rate-1)*EPS)))
        require(sp.cancel(residual)==0,'Opposite full numerator derivative identity')
        identities.append({'exponent_rate':rate,'residual':'0','all_four_corners':{'00':'1','01':'0','10':'0','11':'z'}})
    for index in ['01','02']:
        a,b=(orbits[p+index]['J'] for p in ['2e-Pi1+Pi1-','2e-Pi2-'])
        require(a['integrand']==b['integrand'] and a['bounds']==b['bounds'],'Pi2/chain shared J geometry')
    negative=[]
    def reject(label,operation):
        try:operation()
        except (ValueError,AssertionError,KeyError):negative.append({'case':label,'rejected':True})
        else:raise ValueError('Negative control accepted '+label)
    bad=copy.deepcopy(orbits['4e-Pi1-03']);bad['J']['se_numerator']=bad['J']['metric_numerator']
    reject('omit longitudinal tensor',lambda:source_check(bad))
    bad=copy.deepcopy(orbits['4e-Pi1-01']);bad['C']['components']['euler']='exp(-6*EulerGamma*eps)'
    reject('double Euler compensation',lambda:source_check(bad))
    bad=copy.deepcopy(orbits['2e-Pi2-01']);bad['integral_structure']['core']['insertion_count']=2
    reject('Pi2 treated as a two-insertion chain',lambda:source_check(bad))
    bad=copy.deepcopy(orbits['4e-Pi1-01']);bad['J']['iterated_order_outer_to_inner'].reverse()
    reject('illegal dependent integration order',lambda:source_check(bad))
    bad=copy.deepcopy(results['2e-Pi2-01']);bad['color_channels'][0]['coefficients']['0']['hlog_terms'].append({'word':[],'zeta_powers':[0,0,0],'coefficient':'1'})
    reject('alter physical finite coefficient',lambda:convolution_check(orbits['2e-Pi2-01'],bad))
    bad=copy.deepcopy(results['2e-Pi1+Pi1-01']);bad['color_channels']=bad['color_channels'][:2]
    reject('drop full chain flavor-squared channel',lambda:convolution_check(orbits['2e-Pi1+Pi1-01'],bad))
    report={'format':'wilson-P5-independent-controls/1','status':'passed','results_sha256':sha(args.results),
        'script_sha256':sha(Path(__file__)),'source_files':{k:{'path':str(p.relative_to(PACK)),'sha256':sha(p)} for k,p in source_files.items()},
        'source_orbits':15,'paired_polynomial_integrals':30,'polynomial_controls':polynomial,
        'operator_and_color_convolution':convolution,'Gamma_finite_coefficient_controls':cauchy,
        'opposite_IBP_identities':identities,'negative_controls':negative,'elapsed_seconds':time.perf_counter()-started,
        'coverage':{'all15':'Source tensor, Euler/operator identity, original ordered J polynomial integration and full physical color convolution.',
            '4e-Pi1-01/02/05/06 and adjacent chain/Pi2':'Independent source-derived Gamma Laurent extraction including finite terms.',
            'other nine':'No new original multidimensional finite-integral numerical integration; analytic moment/dictionary reuse remains explicit.'},
        'limitations':['Integer-regulator checks do not determine Laurent finite terms of divergent J.',
            'Cauchy/HLOG truncation errors are observed numerical residuals, not rigorous error enclosures.',
            'Pi2 internal full-D literature operator is adopted, not independently integrated.',
            'No formal original-integral proof is claimed.'],'formal':'not_run'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
    print(json.dumps({'status':'passed','source_orbits':15,'paired_polynomial_integrals':30,'Gamma_finite_controls':len(cauchy),'negative_controls':len(negative),'elapsed_seconds':report['elapsed_seconds']}))

if __name__=='__main__':main()
