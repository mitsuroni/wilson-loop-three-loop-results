"""Build full Laurent coefficients from freshly computed universal moment coefficients."""
from pathlib import Path
import json
import sympy as sp
from functools import lru_cache
from reduce_inputs import ROOT,HPARAMS,H,HB,eps,z,ell,R2,G,exp_series,load_inputs,reduce
import hlog_algebra as ha

HERE=Path(__file__).resolve().parent
def moment_run():
    return ROOT/json.loads((HERE/"current-moments.json").read_text())["run"]
@lru_cache(None)
def coefficient_atoms():
    run=moment_run()
    atoms={};meta={}
    for key in HPARAMS:
        for q in range(5):
            for ch in ("plus","minus"):
                symbol=sp.Symbol(f"{key}_{q}_{ch}")
                v=ha.parse(run/f"{key}_{q}_{ch}.txt")
                atoms[symbol]=v
                meta[symbol]={"moment":key,"q":q,"channel":ch,"weight":q+2,"parity":1 if ch=="plus" else -1}
    return atoms,meta

def coefficients(row,a):
    P=2*row["loop_order"]-2
    atoms,meta=coefficient_atoms()
    substitutions={}
    for key in HPARAMS:
        plus=sum(sp.Symbol(f"{key}_{q}_plus")*eps**q for q in range(P+1))
        minus=sum(sp.Symbol(f"{key}_{q}_minus")*eps**q for q in range(P+1))
        substitutions[H[key]]=plus+minus
        if HB[key]!=H[key]:substitutions[HB[key]]=plus-minus
    reg=sp.cancel(eps**P*a["I_compact"]).subs(substitutions)
    gs=exp_series({k:sp.zeta(k)/k for k in range(2,P+1)},P)
    rs=exp_series({k:(-1)**k*sp.zeta(k)*(2-2**k)/k for k in range(2,P+1)},P)
    # Expand universal prefactors before multiplying the coefficient polynomials.
    reg=reg.subs({G:sp.Symbol("g"),R2:sp.Symbol("r")})
    reg=sp.series(reg,eps,0,P+1).removeO().subs({sp.Symbol("g"):gs,sp.Symbol("r"):rs})
    reg=sp.Poly(sp.expand(reg),eps)
    out={};raw={}
    mapping={**atoms,ell:ha.LOG}
    for p in range(-2*row["loop_order"],1):
        expr=reg.coeff_monomial(eps**(p+P)) if p+P>=0 else sp.Integer(0)
        raw[p]=expr
        out[p]=ha.from_sympy(expr,mapping)
    return out,raw,atoms,meta

def primary_generators(atoms,meta):
    """Old tau first; other even coefficients and ell times odd coefficients follow."""
    gens=[]
    for k in ha.ZW:
        gens.append({"name":f"zeta_{k}","latex":r"\zeta("+str(k)+")","weight":k,"v":ha.period(k)})
    for key in ("H0","H1","H2","HW","HD","HX"):
        for q in range(5):
            symbol=sp.Symbol(f"{key}_{q}_plus")
            if not atoms[symbol]:continue
            label=(r"\tau_{"+str({"H0":0,"H1":1,"H2":2}[key])+","+str(q)+"}") if key in ("H0","H1","H2") else (r"\tau^{+}_{"+key[1:]+","+str(q)+"}")
            gens.append({"name":f"tau_{key}_{q}","latex":label,"weight":q+2,"v":atoms[symbol]})
    for key in ("H1","H2","HD"):
        for q in range(4):
            symbol=sp.Symbol(f"{key}_{q}_minus")
            if atoms[symbol]:
                gens.append({"name":f"eta_{key}_{q}","latex":r"\eta_{"+key[1:]+","+str(q)+"}",
                             "weight":q+3,"v":ha.mul(ha.LOG,atoms[symbol])})
    return gens

def m_generators(with_log=False):
    run=moment_run();gens=[]
    for k in ha.ZW:
        gens.append({"name":f"zeta_{k}","weight":k,"parity":1,"v":ha.period(k)})
    for wt in range(2,7):
        for nn in range(1,wt):
            kk=wt-nn-1
            gens.append({"name":f"M_{kk}_{nn}","weight":wt,"parity":(-1)**(nn+1),"v":ha.parse(run/f"M_{kk}_{nn}.txt")})
    if with_log:gens.append({"name":"ell","weight":1,"parity":-1,"v":ha.LOG})
    return gens

def candidates(gens,weight):
    result=[]
    def visit(start,left,parity,inds):
        if left==0:
            if parity==1:result.append(tuple(inds))
            return
        for i in range(start,len(gens)):
            g=gens[i]
            if g["weight"]>left:continue
            if g["name"]=="ell" and i in inds:continue
            visit(i,left-g["weight"],parity*g.get("parity",1),inds+[i])
    visit(0,weight,1,[])
    result.sort(key=lambda x:(len(x),tuple(x)))
    return result

def fit(target,gens,weight,witness=False):
    if not target:return {"rank":0,"augmented_rank":0,"fit":[],"candidate_count":0}
    inds=candidates(gens,weight)
    vectors=[]
    for ind in inds:
        v=ha.ONE
        for i in ind:v=ha.mul(v,gens[i]["v"])
        vectors.append(v)
    keys=sorted(set(target).union(*(v for v in vectors)))
    A=sp.Matrix([[v.get(k,0) for v in vectors] for k in keys])
    b=sp.Matrix([target.get(k,0) for k in keys])
    red,piv=A.row_join(b).rref()
    rank=sum(p<len(inds) for p in piv)
    row={"rank":rank,"augmented_rank":len(piv),"candidate_count":len(inds),
         "candidate_terms":[[gens[i]["name"] for i in ind] for ind in inds]}
    if rank==len(piv):
        co=[sp.Integer(0)]*len(inds)
        for i,p in enumerate(piv):co[p]=red[i,-1]
        assert A*sp.Matrix(co)==b
        row["fit"]=[{"factors":[gens[i]["name"] for i in ind],"coefficient":str(v)}
                    for ind,v in zip(inds,co) if v]
    elif witness:
        for w in A.T.nullspace():
            value=(w.T*b)[0]
            if value:
                assert w.T*A==sp.zeros(1,A.cols)
                row["witness"]=[{"coordinate":ha.serialize({k:1})[0],"coefficient":str(v)}
                                for k,v in zip(keys,w) if v]
                row["obstruction"]=str(value);break
        assert "obstruction" in row
    return row

if __name__=="__main__":
    for row in load_inputs():
        a=reduce(row);cs,raw,atoms,meta=coefficients(row,a)
        print(row["orbit_id"],{k:len(v) for k,v in cs.items()},flush=True)
