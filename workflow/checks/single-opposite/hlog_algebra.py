"""Exact arithmetic in HyperInt's Hlog/period output; adapted from the gf42 reference."""
import ast,re,math
from fractions import Fraction as Q
from collections import defaultdict,Counter
from functools import lru_cache
from pathlib import Path
import sympy as sp

ZW=(2,3,5); ZERO=(0,0,0)
Z2,Z3,Z5=sp.symbols("Z2 Z3 Z5")
ZS=(Z2,Z3,Z5)
ONE={((),ZERO):Q(1)}
LOG={((0,),ZERO):Q(1)}
def clean(v): return {k:Q(c) for k,c in v.items() if c}
def add(*vs):
    out=defaultdict(Q)
    for v in vs:
        for k,c in v.items():out[k]+=c
    return clean(out)
def scale(v,c):return clean({k:Q(c)*x for k,x in v.items()})
@lru_cache(None)
def shuffles(a,b):
    if not a:return {b:1}
    if not b:return {a:1}
    out=Counter()
    for w,c in shuffles(a[1:],b).items():out[(a[0],)+w]+=c
    for w,c in shuffles(a,b[1:]).items():out[(b[0],)+w]+=c
    return dict(out)
def mul(a,b):
    out=defaultdict(Q)
    for (wa,za),ca in a.items():
        for (wb,zb),cb in b.items():
            zp=tuple(x+y for x,y in zip(za,zb))
            for word,m in shuffles(wa,wb).items():out[(word,zp)]+=ca*cb*m
    return clean(out)
def power(v,n):
    out=ONE
    for _ in range(n):out=mul(out,v)
    return out
def period(k):
    if k in ZW:
        p=list(ZERO);p[ZW.index(k)]=1
        return {((),tuple(p)):Q(1)}
    if k==4:return scale(power(period(2),2),Q(2,5))
    if k==6:return scale(power(period(2),3),Q(8,35))
    raise ValueError(k)
def parse(path):
    text=Path(path).read_text(encoding="utf-8").split("=",1)[1].strip()
    words=[]
    def replace(m):
        word=tuple(int(a) for a in m[1].split(",") if a.strip())
        assert all(a in (0,-1) for a in word),word
        if word not in words:words.append(word)
        return "g"+str(words.index(word))
    text=re.sub(r"Hlog\(z,\[([^]]*)\]\)",replace,text)
    text=re.sub(r"zeta\[(\d+)\]",r"Z\1",text).replace("^","**")
    gs=sp.symbols("g:"+str(len(words)))
    names={str(s):s for s in (*gs,*ZS)}
    def read(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:return sp.Integer(node.value)
        if isinstance(node,ast.Name) and node.id in names:return names[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -read(node.operand)
        if isinstance(node,ast.BinOp):
            a,b=read(node.left),read(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if isinstance(node.op,ast.Pow) and b.is_Integer and b>=0:return a**b
        raise ValueError("unsupported HyperInt output")
    expr=sp.expand(read(ast.parse(text,mode="eval").body))
    if gs:assert sp.Poly(expr,*gs).total_degree()<=1
    out={}
    for g,w in [*zip(gs,words),(None,())]:
        c=expr.subs(dict.fromkeys(gs,0)) if g is None else expr.coeff(g)
        for powers,val in sp.Poly(c,*ZS).terms():
            assert val.is_Rational
            if val:out[(w,powers)]=Q(int(val.p),int(val.q))
    return out
def from_sympy(expr,atoms):
    if expr.is_Rational:return scale(ONE,Q(int(expr.p),int(expr.q)))
    if expr in atoms:return atoms[expr]
    if expr.func==sp.zeta:return period(int(expr.args[0]))
    if expr==sp.pi:raise ValueError("odd Pi is outside the real period basis")
    if expr.func==sp.Add:return add(*(from_sympy(a,atoms) for a in expr.args))
    if expr.func==sp.Mul:
        out=ONE
        for a in expr.args:out=mul(out,from_sympy(a,atoms))
        return out
    if expr.func==sp.Pow and expr.args[1].is_Integer and expr.args[1]>=0:
        base,exponent=expr.args
        if base==sp.pi:
            assert int(exponent)%2==0
            return scale(power(period(2),int(exponent)//2),6**(int(exponent)//2))
        return power(from_sympy(base,atoms),int(exponent))
    raise ValueError(str(expr))
def serialize(v):
    return [{"word":list(w),"zeta_powers":list(p),"coefficient":str(c)}
            for (w,p),c in sorted(v.items())]
def expression(v,variable="z"):
    z=sp.Symbol(variable,positive=True)
    H=sp.Function("Hlog")
    return sum(sp.Rational(c.numerator,c.denominator)*
               sp.prod(sp.zeta(k)**m for k,m in zip(ZW,p))*
               (H(z,sp.Tuple(*w)) if w else 1) for (w,p),c in v.items())
