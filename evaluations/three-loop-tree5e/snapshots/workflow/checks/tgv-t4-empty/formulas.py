"""Exact Hlog coordinates and fits to the already defined weight-six basis."""
from pathlib import Path
from fractions import Fraction as Q
import sys,json
import sympy as sp
from check_math import ROOT,HERE,read,dump
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as ha

def vec(ts):return {(tuple(t['word']),tuple(t['zeta_powers'])):Q(t['coefficient']) for t in ts}
def weight(k):return len(k[0])+sum(a*b for a,b in zip(k[1],ha.ZW))
def fit(v,basis):
 if not v:return []
 vs=[vec(b['hlog_terms']) for b in basis];keys=sorted(set(v).union(*vs))
 A=sp.Matrix(len(keys),len(vs),lambda i,j:vs[j].get(keys[i],0));rhs=sp.Matrix([v.get(k,0) for k in keys])
 try:sol,pars=A.gauss_jordan_solve(rhs)
 except ValueError:return None
 assert pars.rows==0 and A*sol==rhs
 return [{'coefficient':str(c),'factors':basis[i]['factors']} for i,c in enumerate(sol) if c]

def main():
 names=sys.argv[1:] or ['assemble_13_14_15_40']
 bp=ROOT/read(ROOT/'workflow/checks/unified-weight-six/current-run.json')['run']/'basis-and-conversions.json';basis=read(bp)
 out={}
 for name in names:
  run=ROOT/read(HERE/('current-'+name+'.json'))['run']
  assert 'ASSEMBLY_COMPLETE' in (run/'maple.log').read_text()
  for p in run.glob('I*q*.txt'):
   oid,q=p.stem[1:].split('q');v=ha.parse(p);fits={}
   for w in range(7):
    part={k:c for k,c in v.items() if weight(k)==w}
    if part:fits[str(w)]=fit(part,basis['weight_bases'][str(w)])
   out.setdefault(oid,{})[q]={'hlog_terms':ha.serialize(v),'existing_basis_fit':fits}
   print(oid,q,len(v),fits,flush=True)
 report={'coefficients':out,'basis_path':bp.relative_to(ROOT).as_posix(),'all_fit':all(None not in v['existing_basis_fit'].values() for cs in out.values() for v in cs.values()),'formal':'not_run'}
 dump(HERE/('representation-'+('-'.join(names))+'.json'),report)
if __name__=='__main__':main()