"""T3[] Hlog coefficients and existing-basis comparisons."""
from pathlib import Path
from fractions import Fraction as Q
import sys,json,importlib.util
import sympy as sp
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"workflow/checks/single-opposite"))
import hlog_algebra as ha
def read(p):return json.loads(p.read_text(encoding="utf-8"))
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
def run():return ROOT/read(HERE/"current-run.json")["run"]
def module(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def vec(terms):return {(tuple(t["word"]),tuple(t["zeta_powers"])):Q(t["coefficient"]) for t in terms}
def weight(key):return len(key[0])+sum(a*b for a,b in zip(key[1],ha.ZW))
def fit(v,basis):
 if not v:return []
 vs=[vec(b["hlog_terms"]) for b in basis];keys=sorted(set(v).union(*vs))
 A=sp.Matrix(len(keys),len(vs),lambda i,j:vs[j].get(keys[i],0));rhs=sp.Matrix([v.get(k,0) for k in keys])
 try:x,params=A.gauss_jordan_solve(rhs)
 except ValueError:return None
 assert params.rows==0 and A*x==rhs
 return [{"coefficient":str(c),"factors":basis[i]["factors"]} for i,c in enumerate(x) if c]
def coefficients():
 out={}
 ell=sp.Symbol("ell");g2=ha.from_sympy(5*sp.zeta(2)/2-9*ell**2/8,{ell:ha.LOG})
 for oid,base in [("09",9),("39",9),("10",10),("38",10),("26",26),("30",26)]:
  h={q:ha.parse(run()/f"E{base}q{q}.txt") for q in range(2,5)}
  out[oid]={k:{} for k in range(-6,1)}
  out[oid][-2]=ha.scale(h[2],Q(-1,32))
  out[oid][-1]=ha.scale(h[3],Q(-1,32))
  out[oid][0]=ha.scale(ha.add(h[4],ha.mul(g2,h[2])),Q(-1,32))
 return out
def main():
 assert "T3_EMPTY_COMPLETE" in (run()/"maple.log").read_text()
 bp=ROOT/read(ROOT/"workflow/checks/unified-weight-six/current-run.json")["run"]/"basis-and-conversions.json"
 basis=read(bp);report={}
 for oid,cs in coefficients().items():
  report[oid]={}
  for k,v in cs.items():
   fits={}
   for w in range(7):
    part={key:c for key,c in v.items() if weight(key)==w}
    if part:fits[str(w)]=fit(part,basis["weight_bases"][str(w)])
   report[oid][str(k)]={"hlog_terms":ha.serialize(v),"existing_basis_fit":fits}
   print(oid,k,len(v),fits,flush=True)
 dump(run()/"representation.json",{"basis_path":bp.relative_to(ROOT).as_posix(),"coefficients":report,
  "all_fit":all(None not in x["existing_basis_fit"].values() for c in report.values() for x in c.values()),"formal":"not_run"})
if __name__=="__main__":main()
