"""Compute each Pi2 color coefficient through eps^0, retaining full epsilon forms."""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib,re
import sympy as s
from sympy.printing.mathematica import mathematica_code
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'benchmarks/self-energy/pi2'
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha
e,z,x,y,L=s.symbols('eps z x y L',positive=True)
M01,M11,M21,M03=s.symbols('M01 M11 M21 M03')
G=s.exp(-s.EulerGamma*e)*s.gamma(1-e)
U=s.gamma(1+e)**2*s.gamma(1-e)*s.gamma(1-3*e)/(s.gamma(1-2*e)**2*s.gamma(1+2*e))
R=-G**3*U/(128*e**2*(1-2*e)**2)
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def freeze(p,obj):
    raw=(json.dumps(obj,ensure_ascii=False,indent=2)+'\n').encode('utf-8');p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==raw,p
    else:p.write_bytes(raw)
def check(name,cond):
    assert cond,name
    checks.append({'name':name,'passed':True})
def eq(name,a,b):check(name,s.simplify(s.expand_func(a-b))==0)
def parse(v):
    v=re.sub(r'Gamma\[([^\]]+)\]',r'gamma(\1)',v).replace('^','**')
    return s.sympify(v,locals={'eps':e,'z':z,'z1':z,'s31':x,'s21':x,'s41':y,'L':L,'M01':M01,'M11':M11,'M21':M21,'M03':M03,'gamma':s.gamma,'Pi':s.pi})
def coeffs(expr,lo,hi):
    pol=s.series(expr,e,0,hi+1).removeO().expand()
    return {k:s.simplify(pol.coeff(e,k)) for k in range(lo,hi+1)}
def multiply(a,b,lo=-6,hi=0):
    return {k:s.expand(sum(v*b.get(k-j,0) for j,v in a.items())) for k in range(lo,hi+1)}
def latex(expr):
    # Display even zetas compactly without inventing new basis functions.
    expr=s.expand(expr).subs(s.pi**4,90*s.Symbol('Z4')).subs(s.pi**2,6*s.Symbol('Z2'))
    return s.latex(expr).replace('M_{01}',r'M_{0,1}').replace('M_{11}',r'M_{1,1}').replace('M_{21}',r'M_{2,1}').replace('M_{03}',r'M_{0,3}').replace('Z_{2}',r'\zeta_2').replace('Z_{4}',r'\zeta_4')
def main():
    manifest=BASE/'prepared/manifest.json';meta=load(manifest)
    for r in meta['inputs']+meta['outputs']+meta['code']:check('prepared dependency '+r['path'],sha(ROOT/r['path'])==r['sha256'])
    package_path=BASE/'prepared/2e-Pi2.json';package=load(package_path)
    # Exact Gamma logarithms, sufficient before cancellation of the cubic poles.
    logu=sum(s.zeta(k)*e**k/s.Integer(k)*sum(mult*(-slope)**k for slope,mult in [(1,2),(-1,1),(-3,1),(-2,-2),(2,-1)]) for k in range(2,6))
    useries=coeffs(s.exp(logu),0,5);rseries=coeffs(s.exp(-logu),0,5)
    eq('U cubic coefficient',useries[3],6*s.zeta(3));eq('U quartic coefficient',useries[4],9*s.zeta(4));eq('U fifth coefficient',useries[5],42*s.zeta(5))
    ur=multiply(useries,rseries,0,5)
    for k,v in ur.items():eq('U R expansion '+str(k),v,int(k==0))
    source=package['orbits'][0]['operator'];phi={};psi={};p_exact={}
    for name in ['AA','Af','Ff']:
        aa,bb=[parse(source['channels'][name][k]) for k in ['A','B']]
        p_exact[name]=parse(source['channels'][name]['p'])
        eq(name+' actual scalar source',p_exact[name],aa+bb/U)
        # Keeping A+B together makes cancellation of eps^-3 explicit.
        base=4*e**2*(aa+bb)
        phi_expr=base+sum(4*e**(2+k)*bb*rseries[k] for k in range(3,6))
        psi_expr=base+sum(4*e**(2+k)*aa*useries[k] for k in range(3,6))
        phi[name]=coeffs(phi_expr,-1,4);psi[name]=coeffs(psi_expr,-1,4)
        eq(name+' no spurious phi pole',phi[name].pop(-1),0);eq(name+' no spurious Uphi pole',psi[name].pop(-1),0)
        for k,value in multiply(phi[name],useries,0,4).items():eq(name+' source connection Uphi '+str(k),value,psi[name][k])
    eq('AA leading phi',phi['AA'][0],-s.Rational(25,3));eq('Af leading phi',phi['Af'][0],s.Rational(20,3))
    eq('Ff phi constant',phi['Ff'][0],0);eq('Ff phi linear',phi['Ff'][1],8)
    op_expansion={name:{str(k):{'expression':str(v),'latex':latex(v)} for k,v in row.items()} for name,row in phi.items()}
    freeze(BASE/'operator-expansion-r001.json',{'definition':'phi_c=4*eps^2*p_c','through_power':4,'remainder':'O(eps^5)',
        'p_c_recovery':'p_c=phi_c/(4*eps^2), valid through eps^2 with O(eps^3)',
        'Gamma_ratio_before_cancellation_through':5,'channels':op_expansion,'source_sha256':sha(package_path)})
    moments=load(ROOT/'workflow/checks/single-opposite/current-moments.json');run=ROOT/moments['run']
    files=[run/f'M_{k}_{n}.txt' for k,n in [(0,1),(1,1),(2,1),(0,3)]]
    vectors={atom:ha.parse(p) for atom,p in zip([M01,M11,M21,M03],files)};vectors[L]=ha.LOG
    gseries=s.series(s.exp(3*sum(s.zeta(k)*e**k/s.Integer(k) for k in range(2,5))),e,0,5).removeO()
    hseries=M01+3*e*M11/2+9*e**2*(M21/8+M03/24)
    direct_geom=[-gseries*(2+e)/(1152*e**4*(1-2*e)**2*(1+2*e)),
        -gseries/(384*e**3*(1-2*e)**2*(1+2*e))-gseries*hseries/(128*e**2*(1-2*e)**2*s.cosh(3*e*L/2))]
    result_counts=[]
    for index,o in enumerate(package['orbits']):
        oid=o['orbit_id'];chain_id='2e-Pi1+Pi1-'+oid[-2:]
        old_path=ROOT/f'benchmarks/self-energy/pi1-chain/results/{chain_id}/r001.json';old=load(old_path)
        for dep in old['dependencies']:check(oid+' reused coefficient dependency '+dep['path'],sha(ROOT/dep['path'])==dep['sha256'])
        old_package=ROOT/old['input']['package'];check(oid+' old geometry root hash',sha(old_package)==old['input']['package_sha256'])
        chain=next(v for v in load(old_package)['orbits'] if v['orbit_id']==chain_id)
        eq(oid+' original full integrand equals solved chain',parse(o['J']['integrand']),parse(chain['J']['integrand']))
        check(oid+' bounds/order/multiplicity unchanged',o['J']['bounds']==chain['J']['bounds'] and o['J']['iterated_order_outer_to_inner']==chain['J']['iterated_order_outer_to_inner'] and o['I_target']['sym']==chain['I_target']['sym'])
        check(oid+' one Euler compensation',o['normalization']['euler_compensation_count']==1)
        geoc={int(k):parse(v['expression']) for k,v in old['laurent']['coefficients'].items()}
        for k,v in geoc.items():
            expected={(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in old['laurent']['coefficients'][str(k)]['hlog_terms']}
            check(oid+' old evaluated Hlog coefficient '+str(k),ha.from_sympy(v,vectors)==expected)
        # Direct C route: use G^3(A U+B), without the R/U cancellation route.
        plain=coeffs(direct_geom[index],-6,0)
        channels=[];nonzero=0
        for name,monomial in [('AA','CA^2*CF'),('Af','CA*CF*TF*nf'),('Ff','CF^2*TF*nf')]:
            final=multiply(phi[name],geoc);direct=multiply(psi[name],plain)
            rows={}
            for k,v in final.items():
                v=s.simplify(v);eq(oid+' '+name+' direct-C versus reused-I coefficient '+str(k),v,direct[k])
                check(oid+' '+name+' fully expanded coefficient '+str(k),v.free_symbols <= {L,M01,M11,M21,M03} and not v.has(s.gamma,s.EulerGamma,e))
                if v!=0:nonzero+=1
                rows[str(k)]={'expression':str(v),'language':'SymPy/str','latex':latex(v),'hlog_terms':ha.serialize(ha.from_sympy(v,vectors))}
            leading=min(k for k,v in final.items() if s.simplify(v)!=0)
            shared=parse(old['I']['sympy']);full=4*e**2*p_exact[name]*shared
            channels.append({'channel':name,'color_monomial':monomial,'definition':'I_c=phi_c*I_geom',
                'exact_multiplier':str(4*e**2*p_exact[name]),'full_epsilon':{'sympy':str(full),'expression':mathematica_code(full),'language':'WolframLanguage/InputForm'},
                'laurent':{'from_power':-6,'through_power':0,'leading_power':leading,'remainder':'O(eps)','coefficients':rows}})
        eq(oid+' leading AA control',parse(channels[0]['laurent']['coefficients'][str(-4+index)]['expression']),s.Rational(25,1728) if index==0 else s.Rational(25,1152))
        deps=[old_path,old_package,ROOT/'workflow/checks/self-energy-pi1-chain-solve/compute.py',ROOT/'workflow/checks/self-energy-pi1-chain-solve/exact-checks.json',*files,BASE/'prepared/normalization.json',BASE/'theory.md',BASE/'operator-expansion-r001.json']
        result={'format':'wilson-computational-result/1','orbit_id':oid,'revision':'r001','family':'Pi2-'+('adjacent' if index==0 else 'opposite'),
            'input':{'package':package_path.relative_to(ROOT).as_posix(),'package_sha256':sha(package_path),'manifest':manifest.relative_to(ROOT).as_posix(),'manifest_sha256':sha(manifest)},
            'normalization':{'scheme':'pi2-full-color-euler-1.0','Euler_count':1,'sym':o['I_target']['sym'],
                'external':'(alpha_s/Pi)^3*(S14+S12)*sum(T_c*I_c)','shared_geometry':'I_geom equals theory I; I_c=4 eps^2 p_c I_geom',
                'a1_selection':'full result includes Ff; a1-selected view uses AA/Af only'},
            'J':old['J'],'I_geom':old['I'],'I_geom_laurent':old['laurent'],'color_channels':channels,
            'functions':old['functions'],'dependencies':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in deps],
            'domain':old['domain'],'verification':{'exact':(HERE/'exact-checks.json').relative_to(ROOT).as_posix(),'numeric':(HERE/'numeric-checks.json').relative_to(ROOT).as_posix()},'formal':'not_run'}
        freeze(BASE/f'results/{oid}/r001.json',result);result_counts.append({'orbit_id':oid,'positions':21,'nonzero':nonzero})
        print(oid+f': all 3 colors solved; {nonzero}/21 nonzero Laurent positions.',flush=True)
    freeze(HERE/'exact-checks.json',{'count':len(checks),'checks':checks,'orbits':result_counts,'positions':42,'nonzero':sum(o['nonzero'] for o in result_counts),'formal':'not_run'})
    print(f'{len(checks)} exact checks passed.',flush=True)
if __name__=='__main__':main()
