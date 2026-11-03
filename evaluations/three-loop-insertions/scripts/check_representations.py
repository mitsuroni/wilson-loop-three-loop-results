"""Back-substitute every named-function coefficient with rational shuffles."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from collections import Counter
import argparse,ast,hashlib,importlib.util,json,math
import sympy as sp
PACK=Path(__file__).resolve().parents[1]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def vector(ts):
    out={}
    for t in ts:
        k=tuple(t['word']),tuple(t['zeta_powers']);out[k]=out.get(k,Q(0))+Q(t['coefficient'])
    return {k:v for k,v in out.items() if v}
def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,Q(0))+v
    return {k:v for k,v in out.items() if v}
def scale(a,q):return {k:v*q for k,v in a.items() if v*q}
ONE={((),(0,0,0)):Q(1)}
@lru_cache(None)
def shuffle(a,b):
    if not a:return {b:1}
    if not b:return {a:1}
    out=Counter()
    for w,n in shuffle(a[1:],b).items():out[(a[0],)+w]+=n
    for w,n in shuffle(a,b[1:]).items():out[(b[0],)+w]+=n
    return dict(out)
def mul(a,b):
    out={}
    for (wa,pa),ca in a.items():
        for (wb,pb),cb in b.items():
            p=tuple(x+y for x,y in zip(pa,pb))
            for w,n in shuffle(wa,wb).items():out[w,p]=out.get((w,p),Q(0))+ca*cb*n
    return {k:v for k,v in out.items() if v}
def constant(expr):
    out={}
    for (pi,z3,z5),c in sp.Poly(sp.expand(expr),sp.pi,sp.zeta(3),sp.zeta(5)).terms():
        if pi%2:raise ValueError('Odd pi constant outside dictionary')
        out[(),(pi//2,z3,z5)]=Q(str(c))*6**(pi//2)
    return {k:v for k,v in out.items() if v}
def parse(text,atoms):
    def visit(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return scale(ONE,Q(n.value))
        if isinstance(n,ast.Name) and n.id in atoms:return atoms[n.id]
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):return scale(visit(n.operand),-1 if isinstance(n.op,ast.USub) else 1)
        if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Pow):
            power=n.right.value if isinstance(n.right,ast.Constant) and type(n.right.value) is int else None
            if power is None or not 0<=power<=8:raise ValueError('Unsupported function power')
            if isinstance(n.left,ast.Name) and n.left.id=='pi':return constant(sp.pi**power)
            base=visit(n.left);out=ONE
            for _ in range(power):out=mul(out,base)
            return out
        if isinstance(n,ast.BinOp):
            a,b=visit(n.left),visit(n.right)
            if isinstance(n.op,ast.Add):return add(a,b)
            if isinstance(n.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(n.op,ast.Mult):return mul(a,b)
            if isinstance(n.op,ast.Div):
                if set(b)!={((),(0,0,0))}:raise ValueError('Non-rational function denominator')
                return scale(a,1/b[(),(0,0,0)])
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='zeta' and len(n.args)==1 and isinstance(n.args[0],ast.Constant):return constant(sp.zeta(n.args[0].value))
        raise ValueError('Unsupported named coefficient expression '+ast.dump(n))
    return visit(ast.parse(text,mode='eval').body)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    if a.out.exists():raise ValueError('Use a new output file')
    src=PACK/'snapshots/workflow/checks/single-opposite/hlog_algebra.py'
    spec=importlib.util.spec_from_file_location('p5_dictionary_parser',src);ha=importlib.util.module_from_spec(spec);spec.loader.exec_module(ha)
    atoms={'rho':{((0,),(0,0,0)):Q(1)}};bindings={}
    for f in (PACK/'snapshots/workflow/checks/single-opposite/runs/20260916T090538777022Z').glob('M_*.txt'):
        atoms[f.stem]=ha.parse(f);bindings[f.relative_to(PACK).as_posix()]=hashlib.sha256(f.read_bytes()).hexdigest()
    for old,new in [('M01','M_0_1'),('M11','M_1_1'),('M21','M_2_1'),('M03','M_0_3')]:atoms[old]=atoms[new]
    f=PACK/'snapshots/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json';dictionary=read(f)['functions']
    for name in ['rho_id_0_0_0','rho_id_0_0_1']:atoms[name]=vector(dictionary[name]['hlog_terms'])
    assert atoms['rho_id_0_0_0']==scale(mul(atoms['M01'],atoms['M01']),Q(1,6))
    rows=read(a.results);checks=[]
    for oid,result in rows.items():
        sets=[('I',result['laurent']['coefficients'])]+[(c.get('monomial',c.get('color_monomial')),c['coefficients']) for c in result['color_channels']]
        for name,cs in sets:
            for p,row in cs.items():
                if parse(row['expression'],atoms)!=vector(row['hlog_terms']):raise ValueError('Representation residual '+oid+name+p)
                checks.append({'orbit_id':oid,'quantity':name,'power':int(p),'residual':'0'})
    bad='M_0_1+1';assert parse(bad,atoms)!=atoms['M_0_1']
    report={'status':'passed','positions':len(checks),'checks':checks,'dictionary_bindings':bindings,
        'changed_named_coefficient_rejected':True,'scope':'Existing evaluated function dictionaries are read with the admitted parser; arithmetic/products/shuffles are reconstructed independently. This is representation verification, not original integration.','formal':'not_run'}
    a.out.write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode());print(json.dumps({'status':'passed','positions':len(checks)}))
if __name__=='__main__':main()
