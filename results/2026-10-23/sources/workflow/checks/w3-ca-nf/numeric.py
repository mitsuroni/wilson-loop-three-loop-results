"""Independent real GPL series comparison: 31-graph sum vs final expressions."""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import importlib.util
import hashlib
import json
import re
import mpmath as mp
import sympy as sp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
dep=HERE.parent/'gamma1-03-z/check_numeric_z.py'
spec=importlib.util.spec_from_file_location('z_gpl',dep)
gpl=importlib.util.module_from_spec(spec);spec.loader.exec_module(gpl)
mp.mp.dps=100
words={'z':set(),'r':set()};syms={};z,L,r=sp.symbols('z L r')
def parse(txt):
    def replace(m):
        arg=m[1];w=tuple(int(x) for x in m[2].split(',') if x.strip());words[arg].add(w)
        name='H'+arg+''.join('m' if a==-1 else str(a) for a in w);syms[arg,w]=sp.Symbol(name)
        return name
    txt=re.sub(r'Hlog\(([zr]),\s*\[([^]]*)\]\)',replace,txt)
    txt=re.sub(r'zeta\[(\d+)\]',r'Z\1',txt).replace('Pi','pi').replace('^','**').replace('ln(','log(')
    return sp.sympify(txt)
data=json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
bridge=json.loads((HERE/'physical-bridge.json').read_text(encoding='utf-8'));assert bridge['passed']
source=[(row['part'],{int(p):parse(e) for p,e in row['coefficients'].items()}) for row in data['parts']]
target={}
for line in (HERE/'physical-coefficients.txt').read_text().splitlines():
    kind,p,e=line.split('\t')
    if kind in ['raw','bare','post','shape']:target[kind,int(p)]=parse(e)
ct=json.loads((ROOT/'benchmarks/renormalization/three-loop/counterterms-r001.json').read_text())
moment_defs={}
for name in ['M_0_1','M_1_1','M_2_1','M_1_3','M_1_2']:
    txt=(ROOT/f'workflow/checks/single-opposite/runs/20260916T090538777022Z/{name}.txt').read_text().split('=',1)[1]
    moment_defs[sp.Symbol(name)]=parse(txt)
moment_defs[sp.Symbol('N_1_2')]=moment_defs[sp.Symbol('M_1_3')]-sp.log(z)*moment_defs[sp.Symbol('M_1_2')]
counter={int(p):sp.sympify(row['sympy']).subs(moment_defs).subs(sp.Symbol('rho'),sp.log(z)) for p,row in ct['color_channels']['CF*CA*TF*nf'].items()}
for n in [2,3,5]:
    counter={p:e.subs(sp.zeta(n),sp.Symbol('Z'+str(n))) for p,e in counter.items()}
gpl.WORDS=sorted(words['z'])

@lru_cache(None)
def rcoeff(w):
    assert w[-1]==1
    if len(w)==1:return tuple([mp.mpf(0)]+[-mp.mpf(1)/n for n in range(1,261)])
    tail=rcoeff(w[1:])
    if w[0]==0:return tuple([mp.mpf(0)]+[tail[n]/n for n in range(1,261)])
    total=mp.mpf(0);out=[mp.mpf(0)]
    for n in range(1,261):total+=tail[n-1];out.append(-total/n)
    return tuple(out)
def rvalues(point,depth):
    cache={():mp.mpf(1)}
    def get(w):
        if w in cache:return cache[w]
        if not any(w):v=mp.log(point)**len(w)/mp.factorial(len(w))
        elif w[-1]==1:v=mp.polyval(tuple(reversed(rcoeff(w)[:depth+1])),point)
        else:
            prefix=w[:-1];sh=Counter(prefix[:j]+(0,)+prefix[j:] for j in range(len(w)));m=sh.pop(w)
            v=(get(prefix)*mp.log(point)-mp.fsum(n*get(q) for q,n in sh.items()))/m
        cache[w]=v;return v
    return {w:get(w) for w in words['r']}
def number(ex,subs):
    v=ex.evalf(85,subs=subs);assert not v.free_symbols and abs(sp.im(v))<sp.Rational(1,10**70),v
    return mp.mpf(str(sp.re(v)))
records=[];samples={}
for zv in [mp.mpf('0.5'),mp.mpf(1),mp.mpf(2)]:
    for lv in [mp.mpf('-0.37'),mp.mpf(0),mp.mpf('0.43')]:
        allvalues={}
        for depth in [220,260]:
            zh=dict(zip(gpl.WORDS,gpl.values(zv,depth)));rh=rvalues(1/(1+zv),depth)
            subs={z:zv,L:lv,r:1/(1+zv),**{sp.Symbol('Z'+str(n)):mp.zeta(n) for n in [2,3,5]}}
            subs.update({syms['z',w]:v for w,v in zh.items()});subs.update({syms['r',w]:v for w,v in rh.items()})
            raw={p:mp.fsum(bridge['weights'][part]*number(cs[p],subs) for part,cs in source) for p in range(-6,1)}
            scaled={q:mp.fsum(raw[p]*3**(q-p)*((lv+mp.log(zv)/2)**(q-p)+(lv-mp.log(zv)/2)**(q-p))/mp.factorial(q-p) for p in range(-6,q+1)) for q in range(-6,1)}
            post={q:scaled[q]+number(counter[q],subs) for q in range(-6,1)}
            allvalues[depth]=post
            for q in range(-6,1):
                expected=number(target['post',q],subs);err=abs(post[q]-expected)
                assert err<mp.mpf('1e-28'),(zv,lv,q,err)
        for q in range(-6,1):
            change=abs(allvalues[220][q]-allvalues[260][q]);assert change<mp.mpf('1e-28')
            records.append({'z':str(zv),'L':str(lv),'power':q,'value':mp.nstr(allvalues[260][q],55),'series_220_260_change':mp.nstr(change,8),'pass':True})
        if lv==0:samples[str(zv)]={str(q):mp.nstr(allvalues[260][q],55) for q in range(-6,1)}
        print(json.dumps({'z':str(zv),'L':str(lv),'passed_coefficients':7}),flush=True)
out={'scope':'63 Laurent coefficient comparisons at three positive z and three physical scales; 31 source graphs evaluated separately, real r/z GPL series; no new multidimensional integral computation','rows':records,'samples_at_L0':samples,'all_pass':True,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'inputs.json',HERE/'physical-coefficients.txt',HERE/'physical-bridge.json',dep]},'formal':'not_run','limitations':['Series refinement is not a rigorous interval bound','Uses the current scalar integral results, not an independent re-integration of all original graphs']}
(HERE/'numeric-checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_pass':True,'checks':len(records),'z1_L0_finite':samples['1.0']['0']}))
