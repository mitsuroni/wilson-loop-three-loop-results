"""Exact per-orbit source, normalization, coefficient and function checks for P2."""
from __future__ import annotations
import argparse
import ast
from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import time
import sympy as sp

PACK=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACK/'engine'))
from wilson_calc.algebra import ast_to_sympy
from wilson_calc.assembly import prepared_normalization

def require(v,m):
    if not v: raise ValueError(m)
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fsha(v): return hashlib.sha256(json.dumps(v,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def pointer(v,p):
    for t in p.strip('/').split('/') if p else []:
        t=t.replace('~1','/').replace('~0','~');v=v[int(t)] if isinstance(v,list) else v[t]
    return v

def expression(text,names):
    def visit(n):
        if isinstance(n,ast.Constant) and type(n.value) is int: return sp.Integer(n.value)
        if isinstance(n,ast.Name) and n.id in names: return names[n.id]
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):
            return visit(n.operand) if isinstance(n.op,ast.UAdd) else -visit(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=visit(n.left),visit(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
            if isinstance(n.op,ast.Pow):return a**b
        raise ValueError('Expression outside arithmetic grammar')
    return visit(ast.parse(text.replace('^','**'),mode='eval').body)

def canonical(payload):
    answer={}
    for t in payload['terms']:
        require(not t['periods'],'Unexpected explicit period container')
        w=tuple(int(ast_to_sympy(a)) for a in t['gpl']['word'])
        require(all(a in (0,-1) for a in w),'Unexpected alphabet')
        require(t['gpl']['branch']=='positive-real' and ast_to_sympy(t['gpl']['basepoint'])==0,'GPL branch/base point')
        answer[w]=answer.get(w,0)+ast_to_sympy(t['coefficient'])
    return {w:sp.cancel(sp.expand(v)) for w,v in answer.items() if sp.cancel(sp.expand(v))!=0}

def from_old(terms):
    out={}
    for t in terms:
        a,b,c=t['zeta_powers'];w=tuple(t['word'])
        out[w]=out.get(w,0)+sp.Rational(t['coefficient'])*sp.zeta(2)**a*sp.zeta(3)**b*sp.zeta(5)**c
    return {w:sp.expand(v) for w,v in out.items() if sp.expand(v)!=0}

def equal(a,b): return all(sp.cancel(sp.expand(a.get(w,0)-b.get(w,0)))==0 for w in set(a)|set(b))

def vector(ts):
    out={}
    for t in ts:
        k=(tuple(t['word']),tuple(t['zeta_powers']))
        out[k]=out.get(k,Q(0))+Q(t['coefficient'])
    return {k:v for k,v in out.items() if v}

@lru_cache(None)
def shuffle(a,b):
    if not a:return {b:1}
    if not b:return {a:1}
    out=Counter()
    for w,n in shuffle(a[1:],b).items():out[(a[0],)+w]+=n
    for w,n in shuffle(a,b[1:]).items():out[(b[0],)+w]+=n
    return dict(out)

def multiply(a,b):
    out={}
    for (wa,pa),ca in a.items():
        for (wb,pb),cb in b.items():
            p=tuple(x+y for x,y in zip(pa,pb))
            for w,n in shuffle(wa,wb).items():out[(w,p)]=out.get((w,p),Q(0))+ca*cb*n
    return {k:v for k,v in out.items() if v}

def reconstruct(conversion,generators):
    out={}
    for t in conversion:
        product={((),(0,0,0)):Q(t['coefficient'])}
        for f in t['factors']: product=multiply(product,generators[f])
        for k,v in product.items():out[k]=out.get(k,Q(0))+v
    return {k:v for k,v in out.items() if v}

def source_bridge(orbit,binding,normalization):
    values={}
    for key,ref in binding['source_bindings'].items():
        p=PACK/ref['path'];require(sha(p)==ref['sha256'],'Changed source '+key)
        v=pointer(read(p),ref['json_pointer']);require(fsha(v)==ref['field_sha256'],'Changed source field')
        values[key]=v
    graph,raw,native,original=(values[k] for k in ['graph','native_raw','native_prepared','ordered_original'])
    for key in ['vertices','propagators','phases','bounds']:
        require(graph[key]==native['topology'][key]==orbit['topology'][key],'Topology mismatch '+key)
    require(original['J']==orbit['J'],'Original ordered root changed')
    require(native['original_C']==raw['C'] and native['J']==raw['J'],'Native C/J preservation')
    require(native['normalization']['euler_compensation_count']==1 and native['normalization']['assembly_sign']==1,'Native normalization')
    require(native['preparation']['native_sign_conversion_count']==0,'Unexpected sign conversion')
    require(orbit['c4_symmetry']==native['c4_symmetry'],'C4 mismatch')
    require(orbit['euler_preparation']['gamma_power']==3,'Euler Gamma weight')
    ep=sp.Symbol('eps');z=sp.Symbol('z1',positive=True)
    names={k:sp.Symbol(k,positive=True) for k in orbit['J']['variables']}
    names.update({f'x{i}':sp.Symbol(f'x{i}',positive=True) for i in range(1,7)})
    names.update(eps=ep,z1=z,CA=sp.Symbol('CA',positive=True),CF=sp.Symbol('CF'))
    parse=lambda t:expression(t,names)
    mapping={names[k]:parse(v) for k,v in orbit['input_export']['s_to_x'].items()}
    xs=[names[f'x{i}'] for i in range(1,7)]
    jac=sp.factor(sp.Matrix([mapping[names[k]] for k in orbit['J']['variables']]).jacobian(xs).det())
    jac=sp.simplify(abs(jac));require(sp.expand(jac-parse(orbit['input_export']['jacobian']))==0,'Jacobian mismatch')
    native_integrand=parse(native['J']['expanded'])
    powers=[t for t in sp.Mul.make_args(native_integrand) if isinstance(t,sp.Pow) and sp.expand(t.exp-(ep-1))==0]
    require(len(powers)==3,'Native power inventory')
    bases=[sp.expand(t.base) for t in powers]
    ordered=[sp.expand(parse(v).subs(mapping,simultaneous=True)) for v in orbit['J']['propagator_polynomials'].values()]
    require(sorted(map(sp.srepr,bases))==sorted(map(sp.srepr,ordered)),'Three base polynomials differ')
    require(sp.expand(native_integrand.xreplace({p:sp.Integer(1) for p in powers})-parse(orbit['J']['numerator'])*jac)==0,'Numerator/Jacobian')
    color=parse(raw['C']['components']['color']);CA,CF=names['CA'],names['CF']
    coordinates=list(map(sp.Rational,orbit['color_preparation']['coordinates']))
    reconstructed=sum(c*CA**(2-i)*CF**(i+1) for i,c in enumerate(coordinates))
    require(sp.cancel((color-reconstructed).subs(CF,(CA**2-1)/(2*CA)))==0,'Full color reconstruction')
    # Also determine each coordinate anew from the independent rational relation.
    u=sp.symbols('u0:3');basis=[CA**2*CF,CA*CF**2,CF**3]
    num=sp.fraction(sp.cancel((color-sum(a*b for a,b in zip(u,basis))).subs(CF,(CA**2-1)/(2*CA))))[0]
    solutions=sp.solve(sp.Poly(num,CA).all_coeffs(),u,dict=True)
    require(len(solutions)==1 and [solutions[0][x] for x in u]==coordinates,'Color coordinates nonunique/changed')
    expected={'n':3,'a1':str(coordinates[0]),'sym':orbit['c4_symmetry']['orbit_multiplicity'],
              'kappa':str(sp.Rational(orbit['c4_symmetry']['orbit_multiplicity'])*coordinates[0]/(2*(-2)**3))}
    for k in ['n','kappa']: require(normalization[k]==expected[k],'Wrong normalization '+k)
    return {'status':'passed','check_groups':14,'jacobian':str(jac),'color_coordinates':list(map(str,coordinates)),
            'normalization':expected,'domain':'Positive real z and open cube; inverse maps use parent ratios. Polynomial/Jacobian source bridge only; convergence and endpoint continuation are supplied by family derivations.'}

def verify(results):
    started=time.perf_counter();rows=read(results);manifest=read(PACK/'inputs/source-manifest.json')
    original={o['orbit_id']:o for o in read(PACK/'engine/benchmarks/2ne/prepared/6e.json')['orbits']}
    require(set(rows)==set(original) and len(rows)==39,'39 orbit inventory')
    target=read(PACK/'references/basis-and-conversions.json');targets={(t['orbit_id'],t['eps_power']):t for t in target['targets']}
    generators={g['name']:vector(g['hlog_terms']) for g in target['generators']}
    outputs={};positions=0;backsubs=0
    for bind in manifest['orbits']:
        oid=bind['orbit_id'];r=rows[oid];o=original[oid]
        require(r['status']=='computed' and not r['residuals'],'Incomplete '+oid)
        n=r['provenance']['normalization'];bridge=source_bridge(o,bind,n)
        reference=read(PACK/f'references/results/{oid}.json')
        require(reference['orbit_id']==oid,'Reference identity')
        for key in ['n','a1','sym','kappa']:
            if key in reference['normalization']:require(reference['normalization'][key]==bridge['normalization'][key],'Reference C/I mismatch')
        if 'loop_order' in reference['normalization']:require(reference['normalization']['loop_order']==3,'Reference loop order')
        powers=[]
        for p in range(-6,1):
            if p in r['complete_powers']:actual=canonical(r['coefficients'][str(p)])
            else:
                require(r['lower_bound'] is not None and p<r['lower_bound'],'Missing power '+oid+str(p));actual={}
            t=targets[oid,p];require(equal(actual,from_old(t['hlog_terms'])),'HLOG residual '+oid+str(p))
            require(reconstruct(t['conversion'],generators)==vector(t['hlog_terms']),'Function conversion residual')
            # The common target is additionally tied to the actual currently selected revision.
            old=reference['laurent']['coefficients'].get(str(p))
            if old is not None and 'hlog_terms' in old:require(equal(actual,from_old(old['hlog_terms'])),'Selected result mismatch')
            elif old is None:require(not actual,'Nonzero structural zero')
            elif 'expression' in old:
                from sympy.parsing.mathematica import parse_mathematica
                expr=parse_mathematica(old['expression'])
                expr=expr.replace(lambda x:x.is_Function and x.func.__name__=='Zeta',lambda x:sp.zeta(*x.args))
                ell=sp.Symbol('ell');poly=sp.Poly(sp.expand(expr).subs(sp.log(sp.Symbol('z')),ell),ell)
                scalar_words={(0,)*degree[0]:coefficient*sp.factorial(degree[0]) for degree,coefficient in poly.terms() if coefficient!=0}
                require(equal(actual,scalar_words),'Selected scalar result mismatch '+oid+str(p))
            else:raise ValueError('Selected reference coefficient lacks supported complete representation')
            powers.append({'power':p,'coefficient_match':True,'function_back_substitution':True,'HLOG_word_count':len(actual),
                           'conversion':t['conversion']})
            positions+=1;backsubs+=1
        outputs[oid]={'family':bind['family'],'source_bridge':bridge,'powers':powers,'formal':'not_run'}
    alias_checks=[]
    for a,b,scope in [('6e-10','6e-23','full I'),('6e-11','6e-25','full I'),('6e-12','6e-24','full I'),
                      ('6e-21','6e-22','full I'),('6e-32','6e-34','full J and I'),('6e-37','6e-38','finite coefficient only')]:
        ps=[0] if 'only' in scope else range(-6,1)
        for p in ps:
            aa=canonical(rows[a]['coefficients'][str(p)]) if p in rows[a]['complete_powers'] else {}
            bb=canonical(rows[b]['coefficients'][str(p)]) if p in rows[b]['complete_powers'] else {}
            require(equal(aa,bb),'Alias mismatch')
        alias_checks.append({'a':a,'b':b,'scope':scope,'coefficient_checks':'passed','full_epsilon_scope':'Reuses original family derivation; finite truncation alone does not prove full-epsilon identity.'})
    base=canonical(rows['6e-01']['coefficients']['0']);wrong=dict(base);wrong[()]=wrong.get((),0)+1
    require(not equal(base,wrong),'Altered finite coefficient accepted')
    broken=deepcopy(targets['6e-04',0]['conversion']);broken[0]['coefficient']=str(Q(broken[0]['coefficient'])+1)
    require(reconstruct(broken,generators)!=vector(targets['6e-04',0]['hlog_terms']),'Altered function conversion accepted')
    return {'format':'wilson-6e-batch-verification/1','status':'passed','result_sha256':sha(results),
        'orbits':outputs,'orbit_count':39,'coefficient_positions':positions,'function_back_substitutions':backsubs,
        'alias_checks':alias_checks,'negative_controls':{'changed_finite_coefficient_rejected':True,'changed_function_conversion_rejected':True},
        'elapsed_seconds':time.perf_counter()-started,'scope':'Exact source bridge, independent full-color reconstruction, normalization, HLOG comparison and independent rational shuffle conversion. No new CAS integration or formal proof.',
        'formal':'not_run'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    require(not a.out.exists(),'Output must be new');d=verify(a.results)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:d[k] for k in ['status','orbit_count','coefficient_positions','function_back_substitutions','elapsed_seconds']}))
