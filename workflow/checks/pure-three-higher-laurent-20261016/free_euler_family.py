"""Three independent Euler indices and one independent Mellin index.

E(a,b,c;z)=int t^(b-1)(1-t)^(c-b)*[(1+zt)^(-a)-1] dt.
F=2F1(a,b;1+c;-z)=1+b*g(b,c)*E, g=Gamma(1+c)/(Gamma(1+b)Gamma(1+c-b)).
"""
from pathlib import Path
from fractions import Fraction as Q
import math
import json
import gauss_family as diagonal
import hlog9 as h

prefix=diagonal.prefix

def inverse_gamma(max_degree):
    expo={};parts={0:{(0,0):h.ONE},1:{}}
    for n in range(2,max_degree+1):
        expo[n]={}
        for i in range(1,n+1):
            co=math.comb(n,i)*(-1)**i+(i==n)
            if co:expo[n][(i,n-i)]=h.scale(h.period(n),Q((-1)**n*co,n))
        row={}
        for m in range(2,n+1):
            piece=diagonal.poly_mul(expo[m],parts[n-m])
            row=diagonal.poly_add(row,{key:h.scale(v,Q(m,n)) for key,v in piece.items()})
        parts[n]=row
    return {key:value for row in parts.values() for key,value in row.items()}

def build(max_weight=9):
    f={(0,0,0):h.ONE};q={(0,0,0):{}}
    for n in range(1,max_weight+1):
        for r in range(n+1):
            for s in range(n-r+1):
                t=n-r-s;previous=q.get((r,s,t-1),{})
                value=h.add(h.scale(prefix(previous,0),-1),prefix(previous,-1),
                            h.scale(prefix(q.get((r-1,s,t),{}),-1),-1),
                            h.scale(prefix(q.get((r,s-1,t),{}),-1),-1),
                            h.scale(prefix(f.get((r-1,s-1,t),{}),-1),-1))
                q[(r,s,t)]=value;f[(r,s,t)]=prefix(value,0)
    for (r,s,t),v in f.items():assert v==f[(s,r,t)]
    gamma=inverse_gamma(max_weight-1)
    e={}
    for degree in range(1,max_weight):
        for r in range(1,degree+1):
            for s in range(degree-r+1):
                t=degree-r-s;value={}
                for u in range(1,s+2):
                    for v in range(t+1):
                        value=h.add(value,h.mul(f.get((r,u,v),{}),gamma.get((s+1-u,t-v),{})))
                value=h.scale(value,math.factorial(r)*math.factorial(s)*math.factorial(t))
                assert all(h.weight(key)==degree+1 for key in value)
                assert all(word and set(word)!={0} for word,ps in value)
                e[(r,s,t)]=value
    eh={}
    for (r,s,t),value in e.items():
        for m in range(max_weight-r-s-t-1):
            eh[(r,s,t,m)]=h.scale(prefix(value,0,m+1),(-1)**m*math.factorial(m))
    old,oldhat=diagonal.build(max_weight)
    checks=0
    for (k,j),v in old.items():
        rebuilt=h.add(*(h.scale(e[(r,k-r,j)],math.comb(k,r)) for r in range(1,k+1)))
        assert rebuilt==v,(k,j);checks+=1
    for (k,j,m),v in oldhat.items():
        rebuilt=h.add(*(h.scale(eh[(r,k-r,j,m)],math.comb(k,r)) for r in range(1,k+1)))
        assert rebuilt==v,(k,j,m);checks+=1
    return e,eh,checks

def export(out):
    e,eh,checks=build();records={}
    for (r,s,t),value in e.items():records[f'E_{r}_{s}_{t}']={'weight':r+s+t+1,'terms':h.serialize(value)}
    for (r,s,t,m),value in eh.items():records[f'Ehat_{r}_{s}_{t}_{m}']={'weight':r+s+t+m+2,'terms':h.serialize(value)}
    data={'definition':'Independent-index Euler source E(a,b,c;z) and a Mellin extension with independent exponent.',
          'moment':'E_rst=(-1)^r int C_z(t)^r*A(t)^s*B(t)^t dt/t; r>=1. Factorials are included in the moment.',
          'construction':'Triangular general Gauss ODE, inverse Gamma normalization, Mellin zero-prefix rule.',
          'max_weight':9,'records':records,'diagonal_regressions':checks,'formal':'not_run'}
    (out/'free-euler-families.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    lines=[]
    for name,row in records.items():lines.append(name+':='+diagonal.maple(h.lift(row['terms']))+':')
    (out/'euler-inputs.mpl').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return records
