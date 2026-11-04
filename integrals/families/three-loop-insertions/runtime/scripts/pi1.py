"""Source-normalized eleven-orbit Pi1 evaluation using admitted moment data.

Only frozen source inputs and evaluated mother-function dictionaries are read.
No final integral result or selector is an input. Original computation scripts
are preserved separately; they are not executed as filesystem writers.
"""
from __future__ import annotations
import ast
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as s

PACK=Path(__file__).resolve().parents[1]
SNAP=PACK/'snapshots'
e,z=s.symbols('eps z',positive=True)
rho=s.Symbol('rho',real=True)
M={(k,j):s.Symbol(f'M_{k}_{j}') for k in range(4) for j in range(1,5-k)}
R0,R1=s.symbols('rho_id_0_0_0 rho_id_0_0_1')

def read(rel):return json.loads((SNAP/rel).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def check(value,message):
    if not value:raise ValueError(message)
def eq(a,b,message):check(s.cancel(s.expand(a-b))==0,message)

def parse(text,names):
    """Closed arithmetic reader for the accepted original source strings."""
    def visit(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return s.Integer(n.value)
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
        raise ValueError('Unsupported source arithmetic')
    return visit(ast.parse(text.replace('^','**'),mode='eval').body)

def log_gamma_jet(slopes,through):
    """Taylor recurrence for a product of Gamma(1+a*eps) after Euler removal."""
    logarithm={k:s.zeta(k)*sum(p*(-a)**k for a,p in slopes)/k for k in range(2,through+1)}
    f=[s.Integer(1)]
    for n in range(1,through+1):f.append(s.expand(sum(k*logarithm.get(k,0)*f[n-k] for k in range(1,n+1))/n))
    return sum(v*e**k for k,v in enumerate(f))

def hseries(c,nu_multiple,through):
    """H(nu,lambda=c*nu) after its original linear endpoint integration."""
    W,T=s.symbols('W T');average=(c+1)*W/2+c*(T-rho)/2
    out=0
    for q in range(through+1):
        p=s.Poly(s.expand(((average+T/2)**(q+1)-(average-T/2)**(q+1))/s.factorial(q+1)),W,T)
        out+=(nu_multiple*e)**q*sum(value*M[k,j] for (k,j),value in p.terms())
    return s.expand(out)

def cross(value):return value.xreplace({rho:-rho,**{m:(-1)**(j+1)*m for (k,j),m in M.items()}})

def paired_kernels():
    """Actual J crossing quotient, before current C and Gamma factors."""
    q=log_gamma_jet([(1,1),(2,1),(3,-1)],5)
    c3=s.cosh(3*e*rho/2)
    p1=1/(18*e**4*(1+e))
    p2=(1-q)*s.cosh(e*rho/2)/(4*e**4*(1+e)*c3)
    h0=hseries(s.Integer(0),2,3);hw=hseries(s.Rational(1,2),2,3);hb=cross(hw)
    p3=(s.cosh(e*rho/2)*(h0-(hw+hb)/2)-s.sinh(e*rho/2)*(hw-hb)/2)/(e**2*c3)
    p3+=(1+3*(1-q)*s.cosh(e*rho/2)/c3)/(12*e**3*(1+e))
    h0=hseries(s.Integer(0),1,3);h2=hseries(s.Integer(2),1,3);h2b=cross(h2)
    p4=(s.cosh(e*rho)*(h0-(h2+h2b)/2)-s.sinh(e*rho)*(h2-h2b)/2)/(4*e**2*(1+e)*c3)
    h03=hseries(s.Integer(0),3,2);h01=hseries(s.Integer(0),1,2)
    h2=hseries(s.Integer(2),1,2);h2b=cross(h2)
    boundary=h03+s.cosh(e*rho)*(h01-h2-h2b)-s.sinh(e*rho)*(h2-h2b)
    p5=(2*R0+3*e*R1+boundary/(2*e*(1+e)))/(2*c3)
    p6=hseries(s.Integer(0),1,2)*(hseries(s.Integer(0),2,2)+s.cosh(e*rho)/(2*e*(1+e)))/c3
    return {f'P{k}':v for k,v in enumerate([p1,p2,p3,p4,p5,p6],1)}

def source_check(native,ordered):
    check(native['orbit_id']==ordered['orbit_id'],'Native/ordered orbit identity')
    check(native['J']==ordered['J'],'Original ordered kernel changed')
    check(native['normalization']['assembly_sign']==1,'Unexpected native assembly sign')
    check(native['preparation']['native_sign_conversion_count']==0,'Native source sign conversion')
    check(native['preparation']['Euler_compensation_count']==1,'Euler compensation duplicated')
    check(native['preparation']['Gamma_Euler_weight']=='3','Gamma/Euler weight')
    check(native['C']['components']['euler']=='exp(-3*EulerGamma*eps)','Current Euler factor')
    names={k:s.Symbol(k,positive=True) for k in native['J']['variables']}
    CA,CF,TF,nf=s.symbols('CA CF TF nf');names.update(eps=e,z1=z,CA=CA,CF=CF,TF=TF,nf=nf,G=s.Symbol('G'))
    P=lambda t:parse(t,names)
    color=P(native['C']['components']['wilson_color'])
    a1=s.Poly(color,CA,CF).coeff_monomial(CA*CF)
    eq(color,CF**2+a1*CA*CF,'Full Wilson colour trace')
    eq(a1,s.Rational(ordered['I_target']['a1']),'Connected projected colour')
    rational=P(native['C']['components']['rational'])
    eq(rational,1/(32*e*(1-2*e)),'Native Fourier rational factor')
    b1=P(native['C']['components']['self_energy_color'])
    eq(b1,((10-6*e)*CA-8*(1-e)*TF*nf)/(3-2*e),'Full Pi1 external operator')
    target=ordered['I_target'];sym=native['c4_symmetry']['orbit_multiplicity']
    check(sym==target['sym'] and target['total_loop_order']==3,'C4 factor/loop order')
    J=native['J'];se,bare=J['geometry_polynomials'];F=P(se['f']);d=P(se['d']);N=P(J['se_numerator'])
    v,w=(names[k] for k in se['s_vars'])
    derivative=2*e*(s.diff(F,v,w)*F+(2*e-1)*s.diff(F,v)*s.diff(F,w))
    eq(N/(2*(1+e)),d*F+derivative/(4*e*(1+e)),'Complete transverse/derivative identity')
    expected=N*F**(-2+2*e)/(2*(1+e))*P(bare['d'])*P(bare['f'])**(-1+e)
    eq(P(J['integrand']),expected,'Full original Pi1 plus bare kernel')
    eq(P(J['integrand']),P(J['tensor_decomposition']['metric'])+P(J['tensor_decomposition']['derivative']),'Source tensor decomposition')
    seen=set()
    for variable in J['iterated_order_outer_to_inner']:
        for bound in J['bounds'][variable]:check(P(str(bound)).free_symbols<={names[k] for k in seen},'Unavailable dependent bound')
        seen.add(variable)
    return {'a1':str(a1),'sym':sym,'rational':str(rational),'full_wilson_color':str(color),'B1':str(b1),
            'Euler_count':1,'current_source_convention':'wilson-qcd-action-path-fourier/1.0','source_check_groups':13}

def primitives():
    path=SNAP/'workflow/checks/single-opposite/hlog_algebra.py'
    spec=importlib.util.spec_from_file_location('p5_pi1_hlog_math',path);ha=importlib.util.module_from_spec(spec);spec.loader.exec_module(ha)
    pointer=read('workflow/checks/single-opposite/current-moments.json')['run'];atoms={rho:ha.LOG}
    deps=[{'path':path.relative_to(PACK).as_posix(),'sha256':sha(path)}]
    for (k,j),symbol in M.items():
        p=SNAP/pointer/f'M_{k}_{j}.txt';atoms[symbol]=ha.parse(p);deps.append({'path':p.relative_to(PACK).as_posix(),'sha256':sha(p)})
    p=SNAP/'benchmarks/2ne/families/same-channel-two/structural-basis-functions.json'
    data=read(p.relative_to(SNAP))['functions']
    for symbol,key in [(R0,'rho_id_0_0_0'),(R1,'rho_id_0_0_1')]:
        atoms[symbol]={(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in data[key]['hlog_terms']}
    check(ha.add(atoms[R0],ha.scale(ha.mul(atoms[M[0,1]],atoms[M[0,1]]),-Fraction(1,6)))=={},'Ordered T00 moment connection')
    deps.append({'path':p.relative_to(PACK).as_posix(),'sha256':sha(p)})
    return ha,atoms,deps

def generate():
    native={o['orbit_id']:o for o in read('benchmarks/current/prepared/3loop/4e-Pi1.json')['orbits']}
    ordered={o['orbit_id']:o for o in read('benchmarks/self-energy/pi1/prepared/4e-Pi1-r002.json')['orbits']}
    check(set(native)==set(ordered) and len(native)==11,'Eleven Pi1 source roots')
    family={1:'P1',2:'P1',3:'P3',4:'P4',5:'P2',6:'P2',7:'P5',8:'P5',9:'P4',10:'P3',11:'P6'}
    poles={'P1':5,'P2':3,'P3':4,'P4':2,'P5':1,'P6':2}
    jet_orders={'P1':{'Gamma':5},'P2':{'Gamma':5,'Beta_ratio':5},
                'P3':{'Gamma':4,'Beta_ratio':4,'weighted_H':3},'P4':{'Gamma':2,'weighted_H':3},
                'P5':{'Gamma':0,'weighted_H':2,'ordered_mixed_T':1},'P6':{'Gamma':2,'H0':2}}
    kernels=paired_kernels();ha,atoms,deps=primitives();gamma=log_gamma_jet([(-1,3)],5)
    payloads={};answers={}
    # The regular rational external operator must be retained through the pole
    # order before convolution. Its 1/eps factor is already in source C.
    external={'CA^2*CF':(10-6*e)/(3-2*e),'CA*CF*TF*nf':-8*(1-e)/(3-2*e)}
    for oid,o in native.items():
        source=source_check(o,ordered[oid]);f=family[int(oid[-2:])]
        if f not in payloads:
            prefactor=s.Rational(source['sym'])*s.Rational(source['a1'])/2*parse(source['rational'],{'eps':e})
            series=s.series(gamma*prefactor*kernels[f],e,0,1).removeO().expand()
            coefficients={}
            for p in range(-6,1):
                value=s.factor(series.coeff(e,p));vector=ha.from_sympy(value,atoms)
                coefficients[str(p)]={'expression':s.sstr(value),'language':'SymPy/str','latex':s.latex(value),'hlog_terms':ha.serialize(vector)}
            colored=[]
            for monomial,phi in external.items():
                regular=s.series(phi,e,0,poles[f]+1).removeO().expand();cc={}
                for p in range(-6,1):
                    value=s.factor(sum(regular.coeff(e,k)*series.coeff(e,p-k) for k in range(max(0,p+poles[f])+1)))
                    cc[str(p)]={'expression':s.sstr(value),'language':'SymPy/str','hlog_terms':ha.serialize(ha.from_sympy(value,atoms))}
                colored.append({'monomial':monomial,'exact_scalar_multiplier':s.sstr(phi),'coefficients':cc})
            payloads[f]=(coefficients,colored)
        coefficients,colored=payloads[f]
        answers[oid]={'format':'wilson-p5-pi1-evaluation/1','orbit_id':oid,'family':'Pi1-'+f,
            'normalization':{'scheme':'wilson-self-energy-1pi-insertion/2.0','loop_order':3,**source,
                'external_color':'CA*CF*B1(eps)','formula':'I=(sym*a1/2)*Gcal^3*C_rational*(J(z)+z^(3eps)*J(1/z))/(1+z^(3eps))',
                'Ff_ff':'zero in the connected a1-projected Pi1 contribution; full raw trace retained separately'},
            'laurent':{'from_power':-6,'through_power':0,'leading_power':-poles[f],'remainder':'O(eps)','coefficients':coefficients},
            'color_channels':colored,'explicit_zero_color_channels':['CF^2*TF*nf','CF*(TF*nf)^2'],
            'subintegrals':{'raw_ordered_J':o['J'],'source_paired_J_jet_expression':s.sstr(kernels[f]),
                'required_jet_orders':jet_orders[f],'actual_conservative_Gamma_jet':5,
                'external_operator_Taylor_through':poles[f],
                'primitives':'Existing original M moments and rho_id_0_0_0/1 HLOG evaluations; no final I table read.'},
            'dependencies':deps,'scope':'Original ordered/native source checks, fresh finite-order reduction and assembly using evaluated moment data. No new internal gluon/ghost bubble or full CAS integration.',
            'domain':{'z':'positive real','epsilon':'initial small positive Re(eps); local meromorphic continuation at zero','complex_kinematics':'not newly established'},
            'functions':{'rho':'log(z), not the physical scale L','M_k_j':'Original no-factorial moments; complete frozen HLOG data',
                'rho_id_0_0_0':'[eps^0] T(eps,eps,z) = M_0_1^2/6',
                'rho_id_0_0_1':'[eps^1] T(eps,eps,z); existing complete weight5 HLOG dictionary'},
            'formal':'not_run'}
    return answers

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    if a.output.exists():p.error('Choose a NEW output file')
    result=generate();a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'orbits':len(result),'geometric_positions':77,'colored_positions':154,'new_CAS':0}))
