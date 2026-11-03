"""Exact Hlog coefficients and comparison with the existing weight-six basis."""
from pathlib import Path
from fractions import Fraction as Q
import json,sys,math
import sympy as sp
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"workflow/checks/single-opposite"))
import hlog_algebra as ha
def read(p):return json.loads(p.read_text(encoding="utf-8"))
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
def run():return ROOT/read(HERE/"current-run.json")["run"]
def vec(terms):return {(tuple(t["word"]),tuple(t["zeta_powers"])):Q(t["coefficient"]) for t in terms}
def weight(key):return len(key[0])+sum(a*b for a,b in zip(key[1],ha.ZW))
def conv(a,b,n):return [ha.add(*(ha.mul(a[j],b[k-j]) for j in range(k+1) if j<len(a) and k-j<len(b))) for k in range(n+1)]
def coefficients():
    ell=sp.Symbol("ell")
    atom=lambda x:ha.from_sympy(sp.sympify(x),{ell:ha.LOG})
    G=list(map(atom,[1,0,5*sp.zeta(2)/2,3*sp.zeta(3)]))
    U=list(map(atom,[1,-2,3*sp.zeta(2)-4,2*sp.zeta(3)+2*sp.zeta(2)-8]))
    sech=list(map(atom,[1,0,-9*ell**2/8,0]))
    GU=conv(conv(G,U,3),sech,3);Gs=conv(G,sech,3)
    ET=[{}]+[ha.parse(run()/f"ET{q}.txt") for q in range(1,5)]
    EK=[{}]+[ha.parse(run()/f"EK{q}.txt") for q in range(1,5)]
    assert EK[1]=={}
    c32=conv(GU,ET,4);c12=conv(Gs,EK,4)
    out={oid:{} for oid in ("11","12","32")}
    for k in range(-6,1):
        a=ha.scale(c32[k+4],Q(-1,128)) if k>=-4 else {}
        b=ha.scale(c12[k+4],Q(-1,128)) if k>=-4 else {}
        out["32"][k]=a;out["12"][k]=b;out["11"][k]=ha.add(a,ha.scale(b,-1))
    return out
def check_reused_T():
    prior=ROOT/read(ROOT/"workflow/checks/single-opposite/current-moments.json")["run"]
    for q in range(1,5):
        from_old={}
        for j in range(q+1):
            diff=ha.add(ha.parse(prior/f"H0_{q-j}_plus.txt"),ha.parse(prior/f"H0_{q-j}_minus.txt"),
                        ha.scale(ha.parse(prior/f"H2_{q-j}_plus.txt"),-1),ha.scale(ha.parse(prior/f"H2_{q-j}_minus.txt"),-1))
            from_old=ha.add(from_old,ha.scale(ha.mul(ha.power(ha.LOG,j),diff),-Q(1,2**j*math.factorial(j))))
        assert from_old==ha.parse(run()/f"T{q}.txt"),q
    return prior
def fit(v,basis):
    if not v:return []
    keys=sorted(set(v).union(*(vec(b["hlog_terms"]) for b in basis)))
    vs=[vec(b["hlog_terms"]) for b in basis]
    A=sp.Matrix(len(keys),len(vs),lambda i,j:vs[j].get(keys[i],0));b=sp.Matrix([v.get(k,0) for k in keys])
    try:x,params=A.gauss_jordan_solve(b)
    except ValueError:return None
    assert not params.rows
    assert A*x==b
    return [{"coefficient":str(c),"factors":basis[i]["factors"]} for i,c in enumerate(x) if c]
def main():
    assert "T1_OPPOSITE_COMPLETE" in (run()/"maple.log").read_text()
    prior=check_reused_T()
    basispath=ROOT/read(ROOT/"workflow/checks/unified-weight-six/current-run.json")["run"]/"basis-and-conversions.json"
    basis=read(basispath);cs=coefficients();report={}
    for oid,coeff in cs.items():
        report[oid]={}
        for k,v in coeff.items():
            by_weight={}
            for w in range(7):
                part={key:c for key,c in v.items() if weight(key)==w}
                if part:by_weight[str(w)]=fit(part,basis["weight_bases"][str(w)])
            report[oid][str(k)]={"hlog_terms":ha.serialize(v),"existing_basis_fit":by_weight}
            print(oid,k,"terms",len(v),"fit",by_weight,flush=True)
    dump(run()/"coefficients.json",report)
    dump(run()/"representation.json",{"prior_T_source":prior.relative_to(ROOT).as_posix(),
      "basis_path":basispath.relative_to(ROOT).as_posix(),"coefficients":report,
      "formal":"not_run","all_fit":all(None not in x["existing_basis_fit"].values() for c in report.values() for x in c.values())})
if __name__=="__main__":main()
