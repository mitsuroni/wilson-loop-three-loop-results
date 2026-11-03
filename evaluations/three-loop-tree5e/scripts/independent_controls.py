"""Scoped independent controls for P6; no historical final-I table is read.

The source algebra is reconstructed from the four-edge Gram matrix and SU(N)
Fierz identity.  Polynomial probes concern S/P and the ordered domain, not the
singular full J at integer epsilon.  Gamma/Cauchy controls below do reach finite
I for four explicitly identified orbits.  Replaying raw HyperInt primitives is
a different activity, performed by the generation adapters.
"""
from __future__ import annotations

import argparse
import ast
import copy
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import re
import time

import mpmath as mp
import sympy as sp

PACK = Path(__file__).resolve().parents[1]
E, Z, NC = sp.symbols('eps z1 Nc')
AL = sp.symbols('alpha1:4')
GRAM = sp.Matrix([[0,Z,-1-Z,1],[Z,0,1,-1-Z],[-1-Z,1,0,Z],[1,-1-Z,Z,0]])
NAMES = {'eps':E,'z1':Z,'z':Z,'I':sp.I,'Pi':sp.pi,'pi':sp.pi,'EulerGamma':sp.EulerGamma}
FAMILIES = {'t1-empty','t2-empty','t3-empty','t1-opposite','t2-opposite','t3-opposite','t4-empty','t4-opposite'}
EXPECTED_COUNTS = {'t1-empty':7,'t2-empty':7,'t3-empty':6,'t1-opposite':3,'t2-opposite':3,'t3-opposite':5,'t4-empty':8,'t4-opposite':5}
CORRECTED = {'13','14','15','33','34','35'}


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


@lru_cache(maxsize=None)
def parse(text):
    """Closed scalar grammar; source expressions are data, never Python code."""
    def visit(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:
            return sp.Integer(n.value)
        if isinstance(n,ast.Name):
            return NAMES.get(n.id,sp.Symbol(n.id))
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):
            return -visit(n.operand) if isinstance(n.op,ast.USub) else visit(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=visit(n.left),visit(n.right)
            for kind,operation in [(ast.Add,lambda:a+b),(ast.Sub,lambda:a-b),(ast.Mult,lambda:a*b),(ast.Div,lambda:a/b),(ast.Pow,lambda:a**b)]:
                if isinstance(n.op,kind):return operation()
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and not n.keywords:
            functions={'Gamma':sp.gamma,'gamma':sp.gamma,'Exp':sp.exp,'exp':sp.exp,'Log':sp.log,'log':sp.log,'zeta':sp.zeta}
            if n.func.id in functions:return functions[n.func.id](*(visit(a) for a in n.args))
        raise ValueError('Outside declared scalar grammar: '+ast.dump(n))
    return visit(ast.parse(text.replace('^','**').replace('[','(').replace(']',')'),mode='eval').body)


def equal(a,b,message):
    require(sp.expand(a-b)==0,message)


def pos(edge,name):
    return sp.Matrix([1 if j<edge-1 else sp.Symbol(name) if j==edge-1 else 0 for j in range(4)])


def dot(a,b):
    return sp.expand((a.T*GRAM*b)[0])


def kernel(topology):
    legs=topology['tgv_legs'];xs=[pos(t['edge'],t['s_var']) for t in legs]
    S=sum(AL[i]*AL[j]*dot(xs[i]-xs[j],xs[i]-xs[j])/2 for i in range(3) for j in range(i+1,3))
    P=0
    for a,b,c in [(0,1,2),(1,2,0),(2,0,1)]:
        rho=GRAM[legs[a]['edge']-1,legs[b]['edge']-1];v=sp.eye(4)[:,legs[c]['edge']-1]
        for k,sign in [(b,1),(a,-1)]:
            for j in range(3):
                if j!=k:P+=sign*rho*AL[k]*AL[j]*dot(v,xs[k]-xs[j])/4
    return sp.expand(S),sp.expand(P)


def traces_sum(traces):
    """Close fundamental indices using T^a_ij T^a_kl=(delta_il delta_jk-delta_ij delta_kl/N)/2."""
    occurrences={};n=0
    for trace in traces:
        indices=list(range(n,n+len(trace)));n+=len(trace)
        for k,color in enumerate(trace):occurrences.setdefault(color,[]).append((indices[k],indices[(k+1)%len(trace)]))
    require(all(len(v)==2 for v in occurrences.values()),'Each color must occur twice in the Fierz closure')
    total=0
    for bits in itertools.product((0,1),repeat=len(occurrences)):
        parents=list(range(n))
        def find(a):
            while parents[a]!=a:a=parents[a]
            return a
        def union(a,b):parents[find(a)]=find(b)
        for bit,pair in zip(bits,occurrences.values()):
            (i,j),(k,l)=pair
            for a,b in ([(i,l),(k,j)] if bit==0 else [(i,j),(k,l)]):union(a,b)
        loops=len({find(i) for i in range(n)})
        total+=sp.Rational(1,2)**len(occurrences)*(-1)**sum(bits)*NC**(loops-sum(bits))
    return sp.expand(total)


def color(topology):
    vertex_colors={v['color'] for v in topology['vertices']};rename={}
    loops=[t['loop'] for t in topology['tgv_legs']]
    for propagator in topology['propagators']:
        if propagator['loop'] in loops:
            exterior=next(c for c in propagator['color'] if c in vertex_colors)
            rename[exterior]='t'+str(loops.index(propagator['loop']))
        else:
            for c in propagator['color']:rename[c]=propagator['loop']
    # Slots are early-to-late, while U'=ig A U puts the later generator
    # on the left.  This is the actual current path-ordered color word.
    word=[rename[v['color']] for v in reversed(topology['vertices'])]
    return sp.cancel(-2*sp.I*(traces_sum([word,['t0','t1','t2']])-traces_sum([word,['t1','t0','t2']]))/NC)


def order_check(j):
    order=j['iterated_order_outer_to_inner'];require(set(order)==set(j['variables'])==set(j['bounds']),'Variable/order inventory')
    require(len(order)==len(set(order))==5,'Five distinct ordered endpoints')
    available=set()
    for name in order:
        lo,hi=map(parse,j['bounds'][name]);require(lo==0 and hi.free_symbols<=available,'Bound depends on an unavailable endpoint')
        available.add(sp.Symbol(name))
    measure=j['alpha_measure']
    require(measure['kind']=='simplex' and measure['substitution']=={'alpha3':'1-alpha1-alpha2'},'Explicit alpha simplex')
    require(measure['iterated_order_outer_to_inner']==['alpha1','alpha2'] and measure['bounds']=={'alpha1':['0','1'],'alpha2':['0','1-alpha1']},'Ordered alpha measure')
    return order


def normalization_check(native,raw,legacy):
    norm=native['normalization'];prep=native['preparation'];comp=native['C']['components']
    require(norm['euler_compensation_count']==1 and norm['assembly_sign']==1,'Single Euler factor and assembly sign')
    require(prep['native_sign_conversion_count']==0 and prep['Euler_compensation_count']==1,'No hidden sign conversion / duplicate Euler factor')
    require(native['original_C']==raw['C'],'Native original-C preservation')
    require({k:v for k,v in comp.items() if k!='euler'}==raw['C']['components'],'C changes only by declared Euler compensation')
    require(parse(comp['euler'])==sp.exp(-3*sp.EulerGamma*E),'Exactly one weight-three Euler factor')
    require(native['C']['factor_order'].count('euler')==1,'Euler component count')
    actual=sp.prod(parse(comp[k]) for k in native['C']['factor_order'])
    require(sp.simplify(parse(native['C']['combined'])/actual)==1,'C multiplicative assembly')
    CA,CF,alphaS,v14,I0=sp.symbols('CA CF alphaS v14 I0')
    common=CF/16*(alphaS/sp.pi)**3*sp.exp(-3*sp.EulerGamma*E)*sp.gamma(1-E)*sp.gamma(2-2*E)*sp.pi**(3*E)*(I0-2*v14)**(3*E)
    require(sp.simplify(actual/common)==1,'Native full C common normalization')
    require(sp.simplify(parse(legacy['C']['expression'])/(CA**2*common))==1,'Legacy a1 C bridge')
    require(legacy['C']['fixed_basis_coordinates']==['-1/4','1/2','0'],'Full historical fixed-basis color coordinates')
    require(sp.expand(CF-(CA**2*CF-2*CA*CF**2)).subs(CF,(CA**2-1)/(2*CA)).expand()==0,'Exact finite-N color relation')
    require(native['c4_symmetry']['orbit_multiplicity']==legacy['c4_symmetry']['orbit_multiplicity']==4,'Orbit weight')


def source_check(native,raw,legacy):
    oid=native['orbit_id'];normalization_check(native,raw,legacy)
    for field in ['J','J_raw_s','topology','c4_symmetry']:
        require(native[field]==raw[field],oid+' native/raw '+field)
    require(native['J_raw_s']==legacy['J_raw_s'],oid+' legacy/native ordered root')
    j=native['J_raw_s'];order=order_check(j);top=native['topology'];S,P=kernel(top)
    equal(S,parse(j['tgv']['S_polynomial']),oid+' Gram S')
    equal(P,parse(j['tgv']['pnum']),oid+' cyclic-derivative P')
    require(j['tgv']['alpha_monomial']=='(alpha1*alpha2*alpha3)^(-eps)','Alpha epsilon measure')
    require(parse(j['tgv']['denom_power'])==2-2*E,'S exponent')
    pure=[]
    for k,pair in enumerate(top['pure_pairings'],1):
        xa,xb=pos(pair['edge_a'],pair['s_a']),pos(pair['edge_b'],pair['s_b'])
        f=dot(xa-xb,xa-xb)/2;d=GRAM[pair['edge_a']-1,pair['edge_b']-1]
        equal(f,parse(j['f_polynomials']['f'+str(k)]),oid+' pure distance')
        pure.append((sp.expand(f),d))
    require(len(pure)==1,'Exactly one pure propagator in 5e-3g')
    equal(pure[0][1],parse(j['pure']['lorentz_z']),oid+' pure tangent contraction')
    stored_color=parse(native['C']['components']['color']).subs({sp.Symbol('CA'):NC,sp.Symbol('CF'):(NC**2-1)/(2*NC)})
    require(sp.cancel(color(top)-stored_color)==0,oid+' independent Fierz color')
    cube=native['J'];mapping={sp.Symbol(k):parse(v) for k,v in native['representative']['hypercube']['s_to_x'].items()}
    xs=list(map(sp.Symbol,cube['variables']));jac=sp.factor(sp.Matrix([mapping[sp.Symbol(n)] for n in order]).jacobian(xs).det())
    expected_jac=parse(native['representative']['hypercube']['jacobian'])
    equal(jac**2,expected_jac**2,oid+' cube Jacobian absolute determinant')
    equal(S.subs(mapping,simultaneous=True),parse(cube['tgv']['S_polynomial']),oid+' cube S')
    equal(P.subs(mapping,simultaneous=True),parse(cube['tgv']['pnum']),oid+' cube P')
    equal(pure[0][0].subs(mapping,simultaneous=True),parse(cube['pure']['f_polynomials']['f1']),oid+' cube pure distance')
    # Exact pointwise algebra at integer exponents is legal; integration at these
    # epsilon values would be divergent because alpha^(-epsilon) remains.
    sample={Z:sp.Rational(2,3),AL[0]:sp.Rational(1,5),AL[1]:sp.Rational(3,10),AL[2]:sp.Rational(1,2)}
    sample.update({x:sp.Rational(k+2,k+5) for k,x in enumerate(xs)})
    expected=expected_jac*P.subs(mapping,simultaneous=True)*sp.prod(AL)**(-E)*S.subs(mapping,simultaneous=True)**(-2+2*E)*pure[0][1]*pure[0][0].subs(mapping,simultaneous=True)**(-1+E)
    got=parse(cube['integrand']).subs(sp.Symbol('f1'),parse(cube['pure']['f_polynomials']['f1']))
    ordered=parse(j['integrand']).subs(sp.Symbol('f1'),pure[0][0]).subs(mapping,simultaneous=True)*expected_jac
    for ep in [0,1,2]:
        require(sp.cancel((got-expected).subs(E,ep).subs(sample))==0,oid+' reconstructed cube integrand')
        require(sp.cancel((got-ordered).subs(E,ep).subs(sample))==0,oid+' ordered/cube pointwise integrand')
    return {'orbit_id':oid,'Gram_S':True,'cyclic_P':True,'Fierz_color':str(sp.factor(stored_color)),'ordered_source_and_cube':True,'C_bridge':True,'orbit_multiplicity':4,'absolute_cube_jacobian':str(expected_jac),'pointwise_integer_exponents':[0,1,2]},(S,P,mapping,expected_jac)


def integrate_poly(expr,order,bounds):
    expr=sp.expand(expr)
    for name in reversed(order):
        v=sp.Symbol(name);lo,hi=map(parse,bounds[name]);out=0
        for (power,),coefficient in sp.Poly(expr,v).terms():out+=coefficient*(hi**(power+1)-lo**(power+1))/(power+1)
        expr=sp.expand(out)
    return sp.cancel(expr)


def polynomial_probe(native,data):
    S,P,mapping,jac=data;j=native['J_raw_s'];a,b=sp.symbols('_a _b')
    simplex={AL[0]:a,AL[1]:(1-a)*b,AL[2]:(1-a)*(1-b)}
    equal(sp.det(sp.Matrix([simplex[AL[0]],simplex[AL[1]]]).jacobian([a,b])),1-a,'Simplex projective Jacobian')
    rows=[]
    for label,probe in [('1+S',1+S),('P^2',P**2)]:
        expr=probe.subs(Z,sp.Rational(2,3));direct=integrate_poly(expr,j['iterated_order_outer_to_inner'],j['bounds'])
        direct=integrate_poly(direct.subs(AL[2],1-AL[0]-AL[1]),['alpha1','alpha2'],{'alpha1':['0','1'],'alpha2':['0','1-alpha1']})
        transformed=sp.Poly(sp.expand(expr.subs(mapping,simultaneous=True).subs(simplex,simultaneous=True)*jac*(1-a)),*map(sp.Symbol,native['J']['variables']),a,b)
        via_cube=sp.cancel(sum(c/sp.prod(k+1 for k in powers) for powers,c in transformed.terms()))
        equal(direct,via_cube,'Original-order/projective polynomial probe')
        rows.append({'probe':label,'z':'2/3','ordered_simplex_moment':str(direct),'cube_projective_moment':str(via_cube),'residual':'0'})
    return {'orbit_id':native['orbit_id'],'moments':rows,'scope':'Polynomial S/P/domain probes only; these are not full-J or finite-I values.'}


def chains(topology):
    output={i:[] for i in range(1,5)}
    for name in topology['bounds']:output[int(name[1])].append(name)
    for names in output.values():names.sort(key=lambda name:int(name[2:]))
    return output


def reflection(target,source,rho,alpha_permutation=(2,1,0),invert=False,reverse=True,p_sign=1):
    tc,sc=chains(target['topology']),chains(source['topology']);sub={AL[k]:AL[v] for k,v in enumerate(alpha_permutation)}
    if invert:sub[Z]=1/Z
    for edge,vs in tc.items():
        ws=list(reversed(sc[rho[edge]])) if reverse else sc[rho[edge]]
        require(len(vs)==len(ws),'Reflected edge occupation')
        sub.update({sp.Symbol(v):1-sp.Symbol(w) for v,w in zip(vs,ws)})
        before=[sp.S.Zero]+list(map(sp.Symbol,vs))+[sp.S.One]
        after=[sp.S.Zero]+list(map(sp.Symbol,sc[rho[edge]]))+[sp.S.One]
        require([sp.expand((before[k+1]-before[k]).subs(sub,simultaneous=True)) for k in range(len(before)-1)]==[after[k+1]-after[k] for k in reversed(range(len(after)-1))],'Reflected ordered gaps')
    tj,sj=target['J_raw_s'],source['J_raw_s'];factor=Z if invert else 1
    equal(parse(tj['tgv']['S_polynomial']).subs(sub,simultaneous=True)*factor,parse(sj['tgv']['S_polynomial']),'Reflection S')
    equal(p_sign*parse(tj['tgv']['pnum']).subs(sub,simultaneous=True)*factor**2,parse(sj['tgv']['pnum']),'Reflection P')
    equal(parse(tj['f_polynomials']['f1']).subs(sub,simultaneous=True)*factor,parse(sj['f_polynomials']['f1']),'Reflection pure f')
    equal(parse(tj['pure']['lorentz_z']).subs(sub,simultaneous=True)*factor,parse(sj['pure']['lorentz_z']),'Reflection tangent')
    endvars=list(map(sp.Symbol,target['J_raw_s']['variables']));other=list(map(sp.Symbol,source['J_raw_s']['variables']))
    determinant=sp.Matrix([sub[v] for v in endvars]).jacobian(other).det()
    alpha_jac=sp.Matrix([sub[v].subs(AL[2],1-AL[0]-AL[1]) for v in AL[:2]]).jacobian(AL[:2]).det()
    require(abs(determinant)==abs(alpha_jac)==1,'Reflection absolute endpoint/simplex Jacobians')
    require(target['C']==source['C'],'Reflection C equality')
    return {'target':target['orbit_id'],'source':source['orbit_id'],'z_inversion':invert,'absolute_jacobians':[1,1],'J_relation':'J_target(z)=z^(3eps)*J_source(1/z)' if invert else 'J_target(z)=J_source(z)','physical_I_equal':True}


def reflection_checks(native):
    pairs=[];checks=[]
    for a,b in {'05':'04','06':'03','07':'02','08':'01','20':'19','21':'18','31':'17','16':'32','36':'12','37':'11'}.items():
        checks.append(reflection(native[a],native[b],{1:2,2:1,3:4,4:3}));pairs.append((a,b))
    for a,b in {'39':'09','38':'10','30':'26','25':'22','24':'23'}.items():
        checks.append(reflection(native[a],native[b],{1:1,2:4,3:3,4:2},(0,2,1),True));pairs.append((a,b))
    for a,b in {'35':'13','34':'14','33':'15','43':'40','29':'27','42':'41'}.items():
        checks.append(reflection(native[a],native[b],{1:3,2:2,3:1,4:4},(2,1,0),True));pairs.append((a,b))
    return checks,pairs


def t4_structural_checks():
    x,y,w,z,a1,a2,a3=sp.symbols('x y w z a1 a2 a3')
    S=a1*a2*y*(1-x)+a1*a3*((1-x)*(1-w)+z*x*w)+z*a2*a3*w*(1-y)
    P=(a1*a2*(-(1+z)*(1-x)+z*y)+a1*a3*(1+3*z-(1+z)*x+z*(1+z)*w)+z*a2*a3*(1-y-(1+z)*w))/4
    recoil=(1+z)*a3*(a1*(1-x+z*w)-z*a2*w)/2
    derivative=(-z*sp.diff(S,x)-(1+z)*sp.diff(S,y)+sp.diff(S,w))/4
    equal(P,derivative+recoil,'T4 exact IBP generator with recoil')
    require(sp.expand(P-derivative)!=0,'Recoil cannot be dropped')
    moments=[]
    for n,m in itertools.product(range(4),repeat=2):
        lo=sp.Rational(1,(n+1)*(n+m+2));hi=sp.Rational(1,(m+1)*(n+m+2));full=sp.Rational(1,(n+1)*(m+1))
        equal(lo+hi,full,'41+42 complementary ordered-domain partition');moments.append({'n':n,'m':m,'sum':str(full)})
    u,lam,g3,g2,jcore,jpure=sp.symbols('u lambda G3 G2 Jcore Jpure',nonzero=True)
    from_source=lam*g3*jcore*jpure/(8*(1+u**6))
    tgv=lam*g2*jcore/(1+u**4);H0=jpure/u
    ratio=(u**2+u**-2)/(u**3+u**-3)
    factorized=(g3/g2)*ratio*tgv*H0/8
    require(sp.cancel(from_source-factorized)==0,'41/42 factorization weight 1/8')
    require(sp.cancel(from_source-factorized*8)!=0,'Missing 1/8 negative control')
    # Absolute Laurent indices, not number of terms relative to the pole.
    omitted={1:sp.Rational(1,6),2:sp.Rational(1,24)}
    require(sp.series(sp.exp(E)/E**2,E,0,3).removeO().coeff(E,2)==omitted[2],'Absolute Laurent budget witness')
    return {'exact_IBP_recoil':str(recoil),'recoil_required':True,'complementary_partition_moments':moments,'factorized_normalization':'1/8','absolute_Laurent_budget':{'Q_min':-2,'Q_max_required':2,'prefactor_max_required':4,'relative_three_terms_inadequate_for_pole_minus2':{'omitted_Q_epsilon1':'1/6','omitted_Q_epsilon2':'1/24'},'affected_current_r002_orbits':sorted(CORRECTED)},'scope':'Exact source and coefficient-bookkeeping identities; no independent evaluation of a T4 finite integral.'}


class HLog:
    """Independent convergent log-polynomial series for the 0,-1 alphabet."""
    def __init__(self,z,terms=220):self.z=mp.mpf(z);self.terms=terms;self.cache={():{0:[mp.mpf(1)]}}
    def expansion(self,word):
        if word in self.cache:return self.cache[word]
        old=self.expansion(word[1:]);letter=word[0];out={}
        if letter==-1:
            coefficients={};previous=[]
            for n in range(self.terms):
                a=old.get(n,[]);v=[(a[k] if k<len(a) else 0)-(previous[k] if k<len(previous) else 0) for k in range(max(len(a),len(previous)))]
                coefficients[n+1]=v;previous=v
        elif letter==0:coefficients=old
        else:raise ValueError('Gamma control expects only 0,-1 HLOG letters')
        for n,poly in coefficients.items():
            if n>self.terms:continue
            if n==0:out[0]=[mp.mpf(0)]+[v/(k+1) for k,v in enumerate(poly)];continue
            target=[mp.mpf(0)]*len(poly)
            for k,c in enumerate(poly):
                for j in range(k+1):target[k-j]+=c*(-1)**j*mp.factorial(k)/mp.factorial(k-j)/mp.mpf(n)**(j+1)
            out[n]=target
        self.cache[word]=out;return out
    def value(self,word):
        L=mp.log(self.z)
        return mp.fsum(self.z**n*mp.fsum(c*L**k for k,c in enumerate(poly)) for n,poly in self.expansion(tuple(word)).items())
    def coefficient(self,row):
        out=[]
        for term in row['hlog_terms']:
            powers=term.get('zeta_powers',term.get('period_powers'))
            require(powers is not None and not any(powers[3:]),'Non-zeta period in Gamma control')
            c=Fraction(term['coefficient']);out.append(mp.mpf(c.numerator)/c.denominator*self.value(term['word'])*mp.fprod(mp.zeta(k)**p for k,p in zip([2,3,5],powers)))
        return mp.fsum(out)


def cauchy(function,degree,radius,nodes=192):
    r=mp.mpf(radius);values=[]
    for j in range(nodes):
        angle=2*mp.pi*j/nodes;phase=mp.exp(mp.j*angle);values.append((phase,function(r*phase)))
    return [mp.fsum(value/phase**k for phase,value in values)/(nodes*r**k) for k in range(degree+1)]


def gamma_controls(pack,results=None):
    with mp.workdps(80):
        z=mp.mpf('.5');evaluator=HLog(z)
        require(abs(evaluator.value((0,-1))+mp.polylog(2,-z))<mp.mpf('1e-62'),'HLOG independent sanity Li2')
        def factors(e):
            G=mp.exp(-3*mp.euler*e)*mp.gamma(1-e)*mp.gamma(1-2*e)
            A=mp.gamma(1+e)*mp.gamma(1-e);B=mp.gamma(1-e)**2/mp.gamma(2-2*e)
            V=mp.gamma(1+e)*mp.gamma(1+2*e)/mp.gamma(1+3*e)
            R=mp.cosh(e*mp.log(z)/2)/mp.cosh(3*e*mp.log(z)/2)
            return G,A,B,V,R
        controls=[]
        for representative,ids in [('01',['01','08']),('17',['17','31'])]:
            def scaled(e):
                G,A,B,V,R=factors(e);U=2*A-B
                return G*U/576 if representative=='01' else G*R*U*(1-V)/128
            a,b=(cauchy(scaled,6,r) for r in ['.065','.085'])
            require(max(abs(x-y) for x,y in zip(a,b))<mp.mpf('1e-60'),'Dual-radius Cauchy consistency')
            if results is not None:
                for oid in ids:
                    coefficients=results['5e-3g-'+oid]['laurent']['coefficients'];residuals=[]
                    for p in range(-6,1):
                        residual=abs(evaluator.coefficient(coefficients[str(p)])-a[p+6]);require(residual<mp.mpf('1e-52'),'Independent full-I Gamma/Cauchy '+oid+' '+str(p));residuals.append(mp.nstr(residual,8))
                    controls.append({'orbit_id':'5e-3g-'+oid,'powers':list(range(-6,1)),'z':'1/2','max_residual':mp.nstr(max(map(mp.mpf,residuals)),8),'residuals':residuals,'source_formula':'G*(2*A-B)/(576*eps^6)' if representative=='01' else 'G*R*(2*A-B)*(1-V)/(128*eps^6)'})
            else:controls.append({'representative':representative,'dual_radius_consistency_only':True,'frozen_candidate_comparison':'not run'})
        base=pack/'snapshots/workflow/checks/tgv-t1-empty/runs/20260917T115009990802Z'
        def bare(e):
            _,A,B,_,_=factors(e)
            return ((1-e)*B-A)/(e*(1-2*e))
        values=cauchy(bare,5,'.075');res=[]
        for k in range(6):
            expr=parse((base/('bareF'+str(k)+'.txt')).read_text().strip())
            expected=mp.mpf(str(sp.N(expr,85)));delta=abs(values[k]-expected)
            require(delta<mp.mpf('1e-60'),'Raw bareF Gamma coefficient '+str(k));res.append(mp.nstr(delta,8))
        return {'full_I':controls,'bareF':{'raw_moments':list(range(6)),'max_residual':mp.nstr(max(map(mp.mpf,res)),8),'Gamma_formula':'((1-eps)*Beta(1-eps,1-eps)-Gamma(1+eps)*Gamma(1-eps))/(eps*(1-2eps))'},'precision_decimal_digits':80,'Cauchy_nodes':192,'radii':['0.065','0.085'],'HLOG_series_terms':220,'tolerance_full_I':'1e-52','scope':'Full finite-I control applies only to 01,08,17,31; bareF control validates a shared primitive, not every remaining finite coefficient.'}


def vector(row):
    out={}
    for term in row['hlog_terms']:
        powers=term.get('period_powers',term.get('zeta_powers'))
        require(powers is not None and len(powers) in [3,6],'Explicit period coordinates')
        powers=tuple(powers)+(0,)*(6-len(powers));word=tuple(Fraction(str(letter)) for letter in term['word'])
        require(all(isinstance(p,int) and p>=0 for p in powers),'Nonnegative integer period powers')
        require(set(word)<={Fraction(0),Fraction(-1),Fraction(-1,2),Fraction(-2)},'Unexpected HLOG letters require an explicit review')
        key=(word,powers);out[key]=out.get(key,Fraction(0))+Fraction(term['coefficient'])
    return {k:v for k,v in out.items() if v}


def raw_vector(text):
    """Read the *original* linear Maple GPL output without the replay adapter."""
    text=text.split('=',1)[-1].strip().rstrip(':;');words=[]
    def hlog(match):
        word=tuple(Fraction(a.strip()) for a in match[1].split(',') if a.strip());words.append(word)
        return 'H'+str(len(words)-1)
    text=re.sub(r'Hlog\(z,\[([^]]*)\]\)',hlog,text)
    text=re.sub(r'zeta\[1,\s*-3\]','A4',text);text=re.sub(r'zeta\[1,\s*-5\]','A6',text)
    text=text.replace('ln(2)','L2')
    text=re.sub(r'zeta\[(\d+)\]',lambda m:{'2':'Z2','3':'Z3','4':'(2*Z2^2/5)','5':'Z5','6':'(8*Z2^3/35)'}[m[1]],text)
    hs=[sp.Symbol('H'+str(k)) for k in range(len(words))];ps=list(map(sp.Symbol,['Z2','Z3','Z5','L2','A4','A6']))
    poly=sp.Poly(sp.expand(parse(text)),*hs,*ps);out={}
    for powers,c in poly.terms():
        h,p=powers[:len(hs)],powers[len(hs):];require(sum(h)<=1 and c.is_Rational,'Original output is a rational linear HLOG combination')
        word=words[h.index(1)] if any(h) else ();key=(word,tuple(p));out[key]=out.get(key,Fraction(0))+Fraction(int(c.p),int(c.q))
    return {k:v for k,v in out.items() if v}


class PositiveGPL:
    """Independent rational-series evaluator in r=z/(1+z), for positive z.

    d/d r G(a,v;z(r)) = [1/(1-r)+1/(r-a/(1+a))] G(v;z(r));
    for a=-1 only the first term remains.  The transformed -1/2 and -2
    letters are -1 and 2, so this converges throughout 0<r<1.
    """
    def __init__(self,z,terms=380):self.r=mp.mpf(z)/(1+mp.mpf(z));self.terms=terms;self.cache={():{0:[mp.mpf(1)]}}
    def expansion(self,word):
        if word in self.cache:return self.cache[word]
        old=self.expansion(word[1:]);a=word[0];coefficients={};previous=[];previous_extra=[]
        for n in range(self.terms):
            v=old.get(n,[]);length=max(len(v),len(previous),len(previous_extra))
            main=[(v[k] if k<len(v) else 0)+(previous[k] if k<len(previous) else 0) for k in range(length)];previous=main
            extra=[]
            if a not in [0,-1]:
                pole=mp.mpf(a.numerator)/a.denominator/(1+mp.mpf(a.numerator)/a.denominator)
                extra=[(-(v[k] if k<len(v) else 0)+(previous_extra[k] if k<len(previous_extra) else 0))/pole for k in range(length)];previous_extra=extra
            coefficients[n+1]=[main[k]+(extra[k] if extra else 0) for k in range(length)]
        if a==0:
            for n,v in old.items():
                before=coefficients.get(n,[]);coefficients[n]=[(before[k] if k<len(before) else 0)+(v[k] if k<len(v) else 0) for k in range(max(len(before),len(v)))]
        out={}
        for n,poly in coefficients.items():
            if n==0:out[0]=[mp.mpf(0)]+[v/(k+1) for k,v in enumerate(poly)];continue
            target=[mp.mpf(0)]*len(poly)
            for k,c in enumerate(poly):
                for j in range(k+1):target[k-j]+=c*(-1)**j*mp.factorial(k)/mp.factorial(k-j)/mp.mpf(n)**(j+1)
            out[n]=target
        self.cache[word]=out;return out
    def value(self,word):
        L=mp.log(self.r)
        return mp.fsum(self.r**n*mp.fsum(c*L**k for k,c in enumerate(poly)) for n,poly in self.expansion(tuple(word)).items())
    def value_vector(self,v,periods):
        return mp.fsum(mp.mpf(c.numerator)/c.denominator*self.value(w)*mp.fprod(x**p for x,p in zip(periods,ps)) for (w,ps),c in v.items())
    def derivative_vector(self,v,periods):
        z=self.r/(1-self.r)
        return mp.fsum(mp.mpf(c.numerator)/c.denominator*self.value(w[1:])/(z-mp.mpf(w[0].numerator)/w[0].denominator)*mp.fprod(x**p for x,p in zip(periods,ps)) for (w,ps),c in v.items() if w)


def omega_controls(pack,results=None):
    base=pack/'snapshots/workflow/checks/tgv-t4-opposite/runs/20260918T034133130688Z-omega_basis'
    O,P,T=(raw_vector((base/(n+'.txt')).read_text()) for n in ['Omega6','Phi5','Tprimitive'])
    derivative={}
    for (word,powers),c in T.items():
        require(word and word[0]==Fraction(-1,2),'Tprimitive exact integration letter')
        key=(word[1:],powers);derivative[key]=derivative.get(key,Fraction(0))+c
    require(derivative==P,'Exact original Tprimitive derivative = Phi5/(z+1/2)')
    candidate_rows=[]
    if results is not None:
        for suffix,coefficient in [('27',Fraction(-1,64)),('28',Fraction(1,32)),('29',Fraction(-1,64))]:
            result=results['5e-3g-'+suffix];definition=result['function_definitions']['Omega6']
            require(vector({'hlog_terms':definition['hlog_terms']})==O,'Full original Omega coordinates retained')
            require(vector({'hlog_terms':definition['Phi5_hlog_terms']})==P and vector({'hlog_terms':definition['Tprimitive_hlog_terms']})==T,'Original Omega defining primitives retained')
            require(Fraction(result['known_function_view']['Omega6_coefficient'])==coefficient,'Individual Omega coefficient')
            finite=vector(result['laurent']['coefficients']['0']);wide=lambda v:{k:c for k,c in v.items() if any(a not in [0,-1] for a in k[0])}
            expected={k:c*coefficient for k,c in wide(O).items()}
            require(wide(finite)==expected,'Individual finite wide-letter sector '+suffix)
            candidate_rows.append({'orbit_id':'5e-3g-'+suffix,'Omega6_coefficient':str(coefficient),'complete_named_function_and_wide_sector':True})
    with mp.workdps(65):
        A4=mp.quad(lambda t:mp.log(1/t)*mp.log(1+t)**2/(2*t),[0,mp.mpf('.5'),1])
        A6=mp.quad(lambda t:mp.log(1/t)**3*mp.log(1+t)**2/(12*t),[0,mp.mpf('.5'),1])
        periods=[mp.zeta(2),mp.zeta(3),mp.zeta(5),mp.log(2),A4,A6]
        evaluators={s:PositiveGPL(s) for s in ['.5','1','2']};values={s:g.value_vector(O,periods) for s,g in evaluators.items()}
        require(abs(values['.5']-values['2'])<mp.mpf('1e-48'),'Original Omega numerical crossing')
        boundary=abs(values['1']-2*evaluators['1'].value_vector(T,periods));require(boundary<mp.mpf('1e-48'),'Omega(1)=2Tprimitive(1) boundary')
        differential=[]
        for z,zinv in [('.5','2'),('1','1'),('2','.5')]:
            g,gi=evaluators[z],evaluators[zinv];zz=mp.mpf(z)
            target=g.value_vector(P,periods)/(zz+mp.mpf('.5'))-gi.value_vector(P,periods)/(zz**2*(1/zz+mp.mpf('.5')))
            residual=abs(g.derivative_vector(O,periods)-target);require(residual<mp.mpf('1e-47'),'Omega defining derivative at '+z)
            differential.append({'z':z,'residual':mp.nstr(residual,8)})
        return {'raw_derivative_T_Phi_exact':True,'raw_Omega_coordinates':len(O),'candidate_individual_controls':candidate_rows,'numerical_definition':{'r_variable':'z/(1+z)','series_terms':380,'precision':65,'derivative_residuals':differential,'crossing_residual':mp.nstr(abs(values['.5']-values['2']),8),'boundary_at_1_residual':mp.nstr(boundary,8),'periods':'A4/A6 computed from their convergent real one-dimensional defining integrals'},'scope':'Independent representation check of the explicit Omega function and its per-orbit coordinates; not an independent evaluation of the complete original 27/28 finite integral.'}


def representation_checks(results,pairs,bindings):
    require(set(results)=={'5e-3g-'+str(k).zfill(2) for k in range(1,45)},'Exactly 44 actual candidate orbits')
    rows=[]
    for oid,result in results.items():
        coefficients=result['laurent']['coefficients'];require(set(coefficients)=={str(k) for k in range(-6,1)},'All seven Laurent positions')
        require('T4[]' not in json.dumps(result) and 'unevaluated_period' not in json.dumps(result),'No opaque unevaluated T4 period')
        vectors={int(k):vector(v) for k,v in coefficients.items()};require(any(vectors.values()),'Nonempty evaluated I')
        letter_set=sorted({str(a) for v in vectors.values() for w,p in v for a in w})
        if oid[-2:] in {'27','28','29'}:
            require(any(a not in [0,-1] for w,p in vectors[0] for a in w),'Individual wide-letter finite coordinate must survive for '+oid)
        if oid[-2:] in {'13','14','15','33','34','35','40','43'}:
            budget_check(result['epsilon_budget'])
        selected=bindings[oid]['selected_revision'];require(selected==('r002' if oid[-2:] in CORRECTED else 'r001'),'Current T4 budget-correction revision')
        rows.append({'orbit_id':oid,'Laurent_positions':7,'hlog_coordinates':sum(map(len,vectors.values())),'letters':letter_set,'selected_revision':selected})
    for a,b in pairs:
        for p in range(-6,1):require(vector(results['5e-3g-'+a]['laurent']['coefficients'][str(p)])==vector(results['5e-3g-'+b]['laurent']['coefficients'][str(p)]),'Source-derived reflection candidate equality '+a+'/'+b+' '+str(p))
    # A summed channel has a nontrivial kernel: +delta in one graph and -delta
    # in another is invisible to that sum but fails the individual checks above.
    delta=Fraction(1,7);require(delta-delta==0 and delta!=0,'Cancellation-kernel negative witness')
    return {'orbits':rows,'source_derived_reflection_pairs':len(pairs),'exact_pair_coefficient_equalities':7*len(pairs),'Omega_individual_letter_preservation':['27','28','29'],'sum_only_validation_rejected':True,'scope':'Representation completeness and source-derived equality controls; the nonstandard-letter check alone is not a numerical validation of Omega6.'}


def budget_check(budget):
    require(budget['Q_absolute_min']<=-2 and budget['Q_absolute_max']>=2 and budget['even_Q_absolute_max']>=2 and budget['prefactor_absolute_max']>=4 and budget['relative_order_used'] is False,'Absolute T4 Laurent budget required')


def negative_checks(native,raw,legacy,pack,results=None):
    rejected=[]
    for label,mutate in [('duplicate_Euler',lambda x:x['normalization'].__setitem__('euler_compensation_count',2)),('wrong_orbit_weight',lambda x:x['c4_symmetry'].__setitem__('orbit_multiplicity',8)),('wrong_assembly_sign',lambda x:x['normalization'].__setitem__('assembly_sign',-1))]:
        candidate=copy.deepcopy(native['01']);mutate(candidate)
        try:normalization_check(candidate,raw['01'],legacy['01'])
        except ValueError:rejected.append(label)
        else:raise ValueError('Negative was not rejected: '+label)
    candidate=copy.deepcopy(native['01']['J_raw_s']);candidate['iterated_order_outer_to_inner']=list(reversed(candidate['iterated_order_outer_to_inner']))
    try:order_check(candidate)
    except ValueError:rejected.append('unavailable_dependent_bound')
    else:raise ValueError('Unordered source was not rejected')
    for label,kwargs in [('missing_alpha_swap',{'alpha_permutation':(0,1,2)}),('missing_endpoint_reversal',{'reverse':False}),('wrong_numerator_sign',{'p_sign':-1})]:
        try:reflection(native['05'],native['04'],{1:2,2:1,3:4,4:3},**kwargs)
        except ValueError:rejected.append(label)
        else:raise ValueError('Reflection negative was not rejected: '+label)
    x,y,w,z,a,b,c=sp.symbols('x y w z a b c')
    S=a*b*y*(1-x)+a*c*((1-x)*(1-w)+z*x*w)+z*b*c*w*(1-y)
    P=(a*b*(-(1+z)*(1-x)+z*y)+a*c*(1+3*z-(1+z)*x+z*(1+z)*w)+z*b*c*(1-y-(1+z)*w))/4
    missing_recoil=(-z*sp.diff(S,x)-(1+z)*sp.diff(S,y)+sp.diff(S,w))/4
    for label,check in [
        ('missing_T4_recoil',lambda:equal(P,missing_recoil,'Missing recoil mutant')),
        ('missing_factorization_1_over_8',lambda:require(sp.Rational(1,8)==1,'Missing 1/8 mutant')),
        ('relative_T4_Laurent_budget',lambda:budget_check({'Q_absolute_min':-2,'Q_absolute_max':0,'even_Q_absolute_max':0,'prefactor_absolute_max':4,'relative_order_used':True})),
    ]:
        try:check()
        except ValueError:rejected.append(label)
        else:raise ValueError('Negative was not rejected: '+label)
    if results is not None:
        # The common Omega definition cannot be changed in one cancelling
        # direction while passing per-orbit checks, even when a sum is fixed.
        mutant=copy.deepcopy(results)
        for suffix,sign in [('27',1),('28',-1)]:
            rows=mutant['5e-3g-'+suffix]['laurent']['coefficients']['0']['hlog_terms']
            rows.append({'word':['-1/2'], 'period_powers':[0,0,0,0,0,0], 'coefficient':str(Fraction(sign,7))})
        try:omega_controls(pack,mutant)
        except ValueError:rejected.append('individual_Omega_change_hidden_by_channel_sum')
        else:raise ValueError('Individual cancelling Omega mutant was not rejected')
    return rejected


def run(package_root=None,results=None):
    started=time.perf_counter();pack=Path(package_root or PACK);snap=pack/'snapshots';inputs=read(pack/'inputs.json')
    bindings={o['orbit_id']:o for o in inputs['orbits']};require(len(bindings)==44,'Actual 44-orbit bindings')
    counts={family:sum(o['evaluation_family']==family for o in inputs['orbits']) for family in FAMILIES};require(counts==EXPECTED_COUNTS,'Eight-family source partition')
    native={o['orbit_id'][-2:]:o for o in read(snap/'benchmarks/current/prepared/3loop/5e-3g.json')['orbits']}
    raw={o['orbit_id'][-2:]:o for o in read(snap/'workflow/checks/shared-sign-production/generated/montopar/3loop/5e-3g/montopar-initial.json')['orbits']}
    legacy={o['orbit_id'][-2:]:o for o in read(snap/'benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json')['orbits']}
    source_rows=[];probe_rows=[];seen=set()
    for key in sorted(native):
        row,data=source_check(native[key],raw[key],legacy[key]);source_rows.append(row)
        family=bindings['5e-3g-'+key]['evaluation_family']
        if family not in seen:probe_rows.append({'family':family,**polynomial_probe(native[key],data)});seen.add(family)
    reflections,pairs=reflection_checks(native)
    candidate=read(results) if isinstance(results,(str,Path)) else results
    if candidate is not None:
        for wrapper in ['results','orbits']:
            if wrapper in candidate:candidate=candidate[wrapper];break
        if isinstance(candidate,list):candidate={o['orbit_id']:o for o in candidate}
    report={'format':'wilson-P6-independent-controls/1','passed':True,'source_count':44,'family_counts':counts,'source_controls':source_rows,'family_polynomial_probes':probe_rows,'source_reflections':reflections,'T4_source_and_budget':t4_structural_checks(),'Gamma_Cauchy':gamma_controls(pack,candidate),'Omega_definition':omega_controls(pack,candidate),'negative_cases_rejected':negative_checks(native,raw,legacy,pack,candidate),'representation':representation_checks(candidate,pairs,bindings) if candidate is not None else {'status':'not run; frozen candidate not supplied'},'independence_boundary':'No historical final-I table read. Original primitive reuse is not independent integration. Source controls cover 44; independent full finite-I Gamma/Cauchy controls cover four. No new original-domain finite quadrature for the other 40. Formal proof and RG were not run.','elapsed_seconds':round(time.perf_counter()-started,3)}
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--results');parser.add_argument('--output');args=parser.parse_args()
    report=run(results=args.results);text=json.dumps(report,ensure_ascii=False,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text,encoding='utf-8')
    else:print(text,end='')
