"""Bind the 31 current graph contributions to the CA*TF*nf exponent channel."""
from pathlib import Path
import hashlib
import json
import re
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
eps=s.Symbol('eps'); sources={}; rows=[]; hlogs={}

def read(rel):
    p=ROOT/rel; sources[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text(encoding='utf-8-sig'))

def parse(expr):
    if isinstance(expr,dict):
        if 'hlog_terms' in expr:
            out=0
            for t in expr['hlog_terms']:
                h=hl('z',t['word']) if t['word'] else s.Integer(1)
                out+=s.sympify(t['coefficient'])*h*s.prod(s.zeta(n)**p for n,p in zip([2,3,5],t['zeta_powers']))
            return s.expand(out)
        expr=expr.get('sympy',expr['expression'])
    expr=re.sub(r'Hlog\(([zr]),\s*\[([^]]*)\]\)',lambda m:str(hl(m[1],[int(x) for x in m[2].split(',') if x.strip()])),expr)
    expr=re.sub(r'zeta\[(\d+)\]',r'zeta(\1)',expr)
    expr=re.sub(r'\b(?:Z|zeta_)([235])\b',r'zeta(\1)',expr)
    expr=expr.replace('Pi','pi').replace('Zeta(','zeta(').replace('^','**')
    out=s.sympify(expr)
    assert not out.atoms(s.Float),expr
    return out

def hl(arg,word):
    assert set(word)<={0,1,-1}
    name='H_'+arg+'_'+('_'.join('m1' if x==-1 else str(x) for x in word) or 'empty')
    hlogs[name]='Hlog('+arg+',['+','.join(map(str,word))+'])'
    return s.Symbol(name)

def maple(ex):
    txt=str(s.expand(ex)).replace('**','^')
    txt=re.sub(r'\bpi\b','Pi',txt)
    txt=re.sub(r'zeta\((\d+)\)',r'zeta[\1]',txt)
    txt=re.sub(r'\bM_?(\d)_?(\d)\b',r'Mom[\1,\2]',txt)
    txt=re.sub(r'\bN_?(\d)_?(\d)\b',r'Nom[\1,\2]',txt)
    for name,val in hlogs.items():txt=re.sub(r'\b'+name+r'\b',val,txt)
    return txt

def add(oid,part,rel,cs,multiplier=1):
    values={int(q):parse(v) for q,v in cs.items()}
    lo=min(values); mult=s.series(multiplier,eps,0,-lo+1).removeO() if multiplier!=1 else s.Integer(1)
    # Extract only the analytic multiplier's coefficients. Expanding a full
    # Laurent expression with rational z denominators can put eps inside an
    # additive denominator, where .coeff(eps,q) would silently miss that term.
    mc={j:s.expand(mult).coeff(eps,j) for j in range(-lo+1)}
    out={q:s.expand(sum(v*mc.get(q-p,0) for p,v in values.items())) for q in range(-6,1)}
    rows.append({'orbit_id':oid,'part':part,'result':rel,'sha256':sources[rel],'scalar_multiplier':str(multiplier),'coefficients':{str(q):maple(v) for q,v in out.items()}})

bf=-8*(1-eps)/(3-2*eps)
cur=read('benchmarks/self-energy/current-r002.json')
for entry in cur['entries']:
    oid=entry['orbit_id']; rel=entry['result']['path']
    if oid.startswith('4e-Pi1-'):
        d=read(rel);assert d['normalization']['assembly_sign']==1
        add(oid,'Pi1',rel,d['laurent']['coefficients'],bf)
    elif oid.startswith('2e-Pi1+Pi1-'):
        d=read(rel);ch=next(c for c in d['color_channels'] if c['monomial']=='CA*CF*TF*nf')
        add(oid,'chain',rel,ch['coefficients'])
    elif oid.startswith('2e-Pi2-'):
        d=read(rel);ch=next(c for c in d['color_channels'] if c['channel']=='Af')
        add(oid,'Pi2',rel,ch['laurent']['coefficients'])
tc=read('benchmarks/self-energy/tgv-pi1/current.json')
assert tc['complete_count']==12
for oid,entry in tc['orbits'].items():
    rel='benchmarks/self-energy/tgv-pi1/'+entry['result'];d=read(rel)
    assert d['normalization']['euler_compensation_count']==1 and d['normalization']['c4_multiplicity']==4
    cs=d.get('laurent_Hlog_coefficients',d['laurent_coefficients'])
    add(oid,'legs',rel,cs,bf)
gc=read('benchmarks/self-energy/gamma1/current.json')
assert gc['complete_count']==4
for oid,entry in gc['orbits'].items():
    selected=entry.get('preferred_representation',entry)
    rel='benchmarks/self-energy/gamma1/'+selected['result'];d=read(rel)
    assert sources[rel]==selected['sha256']
    assert d['normalization']['euler_compensation_count']==1 and not d['normalization']['external_leg_self_energy_included']
    add(oid,'Gamma',rel,d['channels']['q']['I_coefficients'])
assert len(rows)==31
assert {p:sum(x['part']==p for x in rows) for p in ['Pi1','chain','Pi2','legs','Gamma']}==dict(Pi1=11,chain=2,Pi2=2,legs=12,Gamma=4)
lower=read('benchmarks/renormalization/three-loop/lower-order-inputs-r001.json')
ct=read('benchmarks/renormalization/three-loop/counterterms-r001.json')
two=read('benchmarks/renormalization/two-loop/result-r002.json')
lines=['diagrams:=table():parts:=table():CT:=table():']
for row in rows:
    for q,v in row['coefficients'].items():lines.append(f'diagrams["{row["orbit_id"]}",{q}]:={v}:')
for part in ['Pi1','chain','Pi2','legs','Gamma']:
    for q in range(-6,1):
        names=[f'diagrams["{x["orbit_id"]}",{q}]' for x in rows if x['part']==part]
        lines.append(f'parts["{part}",{q}]:='+ '+'.join(names)+':')
for q,rec in ct['color_channels']['CF*CA*TF*nf'].items():lines.append(f'CT[{q}]:={maple(parse(rec["sympy"]))}:')
CA,CF,T=s.symbols('CA CF T')
for label,key,monomial in [('F2A','F2_over_CF',CA),('F2f','F2_over_CF',T),('k2A','k2',CA),('k2f','k2',T),('h2A','h2',CA),('h2f','h2',T)]:
    ex=parse(two[key]['sympy'])
    if key in ['k2','h2']:ex=s.expand(ex/CF)
    ex=s.expand(ex).coeff(monomial)
    lines.append(label+':='+maple(ex)+':')
lines.append('F1:='+maple(parse(two['F1_over_CF']['sympy']))+':')
moment_dir='workflow/checks/single-opposite/runs/20260916T090538777022Z/'
for name in ['M_0_1','M_1_1','M_2_1','M_3_1','M_1_3','M_1_2']:
    rel=moment_dir+name+'.txt';sources[rel]=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
    k,j=name.split('_')[1:];text=(ROOT/rel).read_text(encoding='utf-8-sig').strip().split('=',1)[1].rstrip(';')
    lines.append(f'Mom[{k},{j}]:={text}:')
lines.append('Nom[1,2]:=Mom[1,3]-ln(z)*Mom[1,2]:')
(HERE/'source_inputs.mpl').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(HERE/'inputs.json').write_text(json.dumps({'date':'2026-10-06','channel':'CF*CA*TF*nf','parts':rows,'source_sha256':sources,'external_scale':'S3=exp(3 eps Ls)+exp(3 eps Lt)','coupling':'lambda=alpha_s/pi','counterterm_scale':'retain the existing S2 and S1, already expanded in L,rho','double_counting':'existing I includes C4/web and Euler; do not multiply those again','formal':'not_run'},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'diagrams':len(rows),'parts':{p:sum(x['part']==p for x in rows) for p in ['Pi1','chain','Pi2','legs','Gamma']},'source_bindings':len(sources),'unscaled_leading_minus6':str(sum(parse(x['coefficients']['-6']) for x in rows))},ensure_ascii=False))
