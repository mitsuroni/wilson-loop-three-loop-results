"""Map exact generic word exports to complete physical source numerators.

The W construction uses the exported Horn parent for the coupled J term, plus
the three explicitly integrated corner terms. Lambda uses its complete sector
export. All endpoint-log powers are retained until their exact coefficients
cancel. The physical outputs are generated without old target coefficients.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import hashlib
import importlib
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=ROOT/'workflow/checks/pure-three-generators-20261015'
sys.dont_write_bytecode=True
sys.path.insert(0,str(OLD))
al=importlib.import_module('w_source_check')
la=importlib.import_module('w_lambda_gpl_check')
z=sp.Symbol('z',positive=True)
VALUES={'1':sp.Integer(1),'Z':-z,'IZ':-1/z}
ENDPOINT_CHECKS=[]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))

@lru_cache(None)
def token(expr):
    expr=sp.cancel(expr)
    for name,value in VALUES.items():
        if sp.cancel(expr-value)==0:return name
    if expr==0:return '0'
    raise ValueError('unexpected physical letter '+str(expr))

def canon(words):
    out=[]
    for word in words:
        word=tuple(word)
        if word and all(x=='0' for x in word):return None
        if word:out.append(word)
    return tuple(sorted(out))

def add_term(out,words,c):
    key=canon(words)
    if key is not None:out[key]+=Q(c)

def cleaned(v):return {k:c for k,c in v.items() if c}
def add(*vectors):
    out=defaultdict(Q)
    for v in vectors:
        for k,c in v.items():out[k]+=c
    return cleaned(out)
def scale(v,c):return {k:x*Q(c) for k,x in v.items() if x*c}
def mul(a,b):
    out=defaultdict(Q)
    for wa,ca in a.items():
        for wb,cb in b.items():add_term(out,wa+wb,ca*cb)
    return cleaned(out)

EXPORT=HERE/'runs/r001'
fd=read(EXPORT/'W-parent-words.json')
ld=read(EXPORT/'Lambda-sector-words.json')
F={tuple(r['powers']):r['terms'] for r in fd['coefficients']}
LS={r['source_degree']:r['terms'] for r in ld['coefficients']}

def specialize(rows,C,U,V,boundary=False,tag=''):
    cc,uu,vv=(VALUES[x] for x in (C,U,V))
    mapping={'0':'0','A':token(1/uu),'PX':token(1/cc),
             'PXY':token(uu/cc),'PZ':token(1/vv)}
    total=defaultdict(Q)
    for row in rows:
        cw=tuple(row['constant_word']);vw=tuple(mapping[x] for x in row['outer_word'])
        pieces=la.leading_one_regularization(cw) if U=='1' else {(0,cw):1}
        cd=row.get('boundary_J_degree',0)
        for (power,word),coef in pieces.items():
            words=[tuple(mapping[x] for x in word),vw]
            if boundary and cd:words.append(('0',)+('Z',)*cd)
            key=canon(words)
            if key is not None:total[(power,key)]+=Q(row['coefficient'])*Q(str(coef))
    total={key:c for key,c in total.items() if c}
    divergent={key:c for key,c in total.items() if key[0]}
    assert not divergent,(tag,divergent)
    if U=='1':ENDPOINT_CHECKS.append({'tag':tag,'exact_log_cancellation':True,'remaining_terms':len(total)})
    return {words:c for (power,words),c in total.items()}

@lru_cache(None)
def parent(i,j,k,C,U,V):
    return specialize(F[(i,j,k)],C,U,V,tag=f'F_{i}_{j}_{k}_{C}_{U}_{V}')

def primitive(degree,param):
    return {((('0',)+(token(1/VALUES[param]),)*degree),):Q(1)}

def primitive_difference(degree,X,A):
    return add(primitive(degree,X),scale(primitive(degree,A),-1))

def inserted_primitive(i,k,X,A):
    """Integral du/u [G((1/X)^i;u)-G((1/A)^i;u)] G(0,(1/A)^k;u)."""
    out=defaultdict(Q)
    for sign,param in ((1,X),(-1,A)):
        pw=(token(1/VALUES[param]),)*i
        jw=('0',)+(token(1/VALUES[A]),)*k
        for word,multiplicity in al.shuffle(pw,jw).items():
            add_term(out,[('0',)+word],sign*multiplicity)
    return cleaned(out)

def w_sector(total,X,Y,A,B,C):
    out={}
    for i in range(1,total-1):
        for j in range(1,total-i):
            k=total-i-j
            coupled={}
            for signu,U in ((1,X),(-1,A)):
                for signv,V in ((1,Y),(-1,B)):
                    coupled=add(coupled,scale(parent(i,j,k,C,U,V),signu*signv))
            edgeu=primitive_difference(i,X,A)
            edgev=primitive_difference(j,Y,B)
            minus_a=scale(mul(inserted_primitive(i,k,X,A),edgev),-1)
            minus_b=scale(mul(edgeu,inserted_primitive(j,k,Y,B)),-1)
            corner=mul(mul(edgeu,edgev),primitive(k,'IZ'))
            out=add(out,coupled,minus_a,minus_b,corner)
    return out

def physical_numerators(total):
    # Parameter tuples agree with the actual oriented p/q changes of variable.
    ws=[w_sector(total,'Z','Z','1','1','Z'),
        w_sector(total,'1','1','IZ','IZ','IZ'),
        scale(w_sector(total,'Z','1','1','IZ','1'),-2)]
    ls=[specialize(LS[total],'Z','1','1',True,f'L{total}_pp'),
        specialize(LS[total],'IZ','IZ','IZ',True,f'L{total}_mm'),
        scale(specialize(LS[total],'1','1','IZ',True,f'L{total}_pm'),-2)]
    return add(*ws),add(*ls)

MAPLE={'0':'0','1':'1','Z':'-z','IZ':'-1/z'}
def maple_expression(v):
    parts=[]
    for words,c in sorted(v.items()):
        factors=[f'Hlog(1,[{",".join(MAPLE[x] for x in word)}])' for word in words]
        parts.append('('+str(c)+')'+''.join('*'+f for f in factors))
    return '+'.join(parts) or '0'

def serialize(v):
    return [{'coefficient':str(c),'gpl_factors':[list(word) for word in words]}
            for words,c in sorted(v.items())]

def build(out,max_order=6):
    records={};lines=[]
    for n in range(3,max_order+1):
        wn,ln=physical_numerators(n)
        for kind,vec in (('W',wn),('L',ln)):
            assert all(sum(len(word) for word in words)==n+3 for words in vec)
            records[f'{kind}_num{n}']={'source_degree':n,'weight':n+3,'terms':serialize(vec)}
            lines.append(f'{kind}n{n}:='+maple_expression(vec)+':')
        print(f'Physical source words N={n}: W={len(wn)}, L={len(ln)}',flush=True)
    (out/'physical-numerators.json').write_text(json.dumps({'letter_map':MAPLE,'records':records,
        'endpoint_checks':ENDPOINT_CHECKS,'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    (out/'source-inputs.mpl').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return records
