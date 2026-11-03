"""Generate the existing Euler P and Mellin-P families from the Gauss ODE.

F={}_2F_1(a,a;1+b;-z) obeys [theta(theta+b)+z(theta+a)^2]F=0.
V_P=(F-1)/(a*g), g=Gamma(1+b)/(Gamma(1+a)*Gamma(1+b-a)).
No target coefficient is read or fitted by this construction.
"""
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import math
import json
import hlog9 as h

def prefix(v,letter,times=1):
    if letter==0:
        assert all(word and set(word)!={0} for word,ps in v)
    return {((letter,)*times+word,ps):c for (word,ps),c in v.items()}

def poly_add(*polys):
    out={}
    for poly in polys:
        for key,v in poly.items():out[key]=h.add(out.get(key,{}),v)
    return {key:v for key,v in out.items() if v}

def poly_mul(a,b):
    out={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            key=(i+k,j+l);out[key]=h.add(out.get(key,{}),h.mul(v,w))
    return {key:v for key,v in out.items() if v}

def build(max_weight=9):
    assert max_weight<=9
    f={(0,0):h.ONE};q={(0,0):{}}
    for total in range(1,max_weight+1):
        for r in range(total+1):
            s=total-r
            prev=q.get((r,s-1),{})
            value=h.add(h.scale(prefix(prev,0),-1),prefix(prev,-1),
                        h.scale(prefix(q.get((r-1,s),{}),-1),-2),
                        h.scale(prefix(f.get((r-2,s),{}),-1),-1))
            q[(r,s)]=value;f[(r,s)]=prefix(value,0)
    # Gamma ratio inverse: exp sum (-1)^n*zeta(n)/n*(a^n+(b-a)^n-b^n).
    exponent={}
    for total in range(2,max_weight):
        row={}
        for i in range(1,total+1):
            coefficient=math.comb(total,i)*(-1)**i+(1 if i==total else 0)
            if coefficient:row[(i,total-i)]=h.scale(h.period(total),Q((-1)**total*coefficient,total))
        exponent[total]=row
    inv_gamma={0:{(0,0):h.ONE},1:{}}
    for total in range(2,max_weight):
        row={}
        for n in range(2,total+1):
            piece=poly_mul(exponent[n],inv_gamma[total-n])
            row=poly_add(row,{key:h.scale(v,Q(n,total)) for key,v in piece.items()})
        inv_gamma[total]=row
    gamma={key:v for row in inv_gamma.values() for key,v in row.items()}
    p={}
    for degree in range(1,max_weight):
        for k in range(1,degree+1):
            j=degree-k;value={}
            for (r,s),fv in f.items():
                if r+s<2 or r>k+1 or s>j:continue
                gv=gamma.get((k+1-r,j-s),{})
                value=h.add(value,h.mul(fv,gv))
            value=h.scale(value,math.factorial(k)*math.factorial(j))
            assert all(h.weight(key)==degree+1 for key in value)
            assert all(word and set(word)!={0} for word,ps in value)
            p[(k,j)]=value
    assert p[(1,0)]=={((0,-1),h.ZERO):Q(-1)}
    assert p[(2,0)]=={((0,-1,-1),h.ZERO):Q(4)}
    assert p[(1,1)]=={((0,0,-1),h.ZERO):Q(1),((0,-1,-1),h.ZERO):Q(-1)}
    phat={}
    for (k,j),value in p.items():
        for m in range(max_weight-k-j-1):
            v=h.scale(prefix(value,0,m+1),(-1)**m*math.factorial(m))
            assert all(h.weight(key)==k+j+m+2 for key in v)
            phat[(k,j,m)]=v
    return p,phat

def maple(v):
    periods=['zeta[2]','zeta[3]','zeta[5]','zeta[7]','zeta[3,5]','zeta[9]']
    terms=[]
    for (word,ps),c in sorted(v.items()):
        factors=[f'({atom})^{power}' for atom,power in zip(periods,ps) if power]
        if word:factors.append('Hlog(z,['+','.join(map(str,word))+'])')
        terms.append('('+str(c)+')'+''.join('*'+f for f in factors))
    return '+'.join(terms) or '0'

def export(out):
    p,ph=build();records={}
    for (k,j),v in p.items():records[f'P_{k}_{j}']={'weight':k+j+1,'terms':h.serialize(v)}
    for (k,j,m),v in ph.items():records[f'Phat_{k}_{j}_{m}']={'weight':k+j+m+2,'terms':h.serialize(v)}
    data={'definition':'Existing P and Mellin-P of single-opposite three-loop.md RG8/RG9/PM6.',
          'construction':'Exact triangular Gauss ODE and Gamma-ratio expansion; Mellin moments prepend zeros.',
          'max_weight':9,'records':records,'low_identity_checks':3,'formal':'not_run'}
    (out/'gauss-families.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    lines=[]
    for name,v in [(f'P_{k}_{j}',v) for (k,j),v in p.items()]+[(f'Phat_{k}_{j}_{m}',v) for (k,j,m),v in ph.items()]:
        lines.append(name+':='+maple(v)+':')
    (out/'gauss-inputs.mpl').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return records
