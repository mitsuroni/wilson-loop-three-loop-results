"""T3[1] full Laurent data and exact existing-basis comparison."""
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
def basispath():return ROOT/read(ROOT/"workflow/checks/unified-weight-six/current-run.json")["run"]/"basis-and-conversions.json"
def fit(v,basis):
 if not v:return []
 vs=[vec(b["hlog_terms"]) for b in basis];keys=sorted(set(v).union(*vs))
 A=sp.Matrix(len(keys),len(vs),lambda i,j:vs[j].get(keys[i],0));rhs=sp.Matrix([v.get(k,0) for k in keys])
 try:x,params=A.gauss_jordan_solve(rhs)
 except ValueError:return None
 assert params.rows==0 and A*x==rhs
 return [{"coefficient":str(c),"factors":basis[i]["factors"]} for i,c in enumerate(x) if c]
def coefficients():
 out={oid:{k:{} for k in range(-6,1)} for oid in ("22","23","24","25","44")}
 for oid,base,factor in (("22",22,32),("25",22,32),("23",23,64),("24",23,64)):
  out[oid][-1]=ha.scale(ha.parse(run()/f"E{base}q1.txt"),Q(1,factor))
  out[oid][0]=ha.scale(ha.parse(run()/f"E{base}q2.txt"),Q(1,factor))
 atoms={sp.Symbol(g["name"]):vec(g["hlog_terms"]) for g in read(basispath())["generators"]}
 M0,M1,M2,Z2=sp.symbols("M_0_1 M_1_1 M_2_1 zeta_2")
 out["44"][-1]=ha.from_sympy(-M0*(2*M0+M1)/32,atoms)
 out["44"][0]=ha.from_sympy(-((4-Z2)*M0**2+3*M0*M1+(M0*M2+M1**2)/2)/32,atoms)
 return out
def main():
 assert "T3_OPPOSITE_COMPLETE" in (run()/"maple.log").read_text()
 basis=read(basispath());report={}
 for oid,cs in coefficients().items():
  report[oid]={}
  for k,v in cs.items():
   fits={}
   for w in range(7):
    part={key:c for key,c in v.items() if weight(key)==w}
    if part:fits[str(w)]=fit(part,basis["weight_bases"][str(w)])
   report[oid][str(k)]={"hlog_terms":ha.serialize(v),"existing_basis_fit":fits}
   if v:print(oid,k,len(v),fits,flush=True)
 dump(run()/"representation.json",{"basis_path":basispath().relative_to(ROOT).as_posix(),"coefficients":report,
  "all_fit":all(None not in x["existing_basis_fit"].values() for c in report.values() for x in c.values()),"formal":"not_run"})
if __name__=="__main__":main()
