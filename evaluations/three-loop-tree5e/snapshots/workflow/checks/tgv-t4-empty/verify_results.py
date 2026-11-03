"""Read-only, offline validation of the frozen eight-member T4[] delivery."""
import hashlib,json
import sympy as sp
from formulas import ROOT,HERE,read,vec,ha,weight
from publish import FULL,DEFS,FINAL
from check_math import IDS,PAIRS,main as check_roots

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 delivery=ROOT/read(HERE/'current-delivery.json')['delivery'];manifest=read(delivery/'manifest.json')
 for b in manifest['bindings']:
  assert sha(ROOT/b['path'])==b['sha256'],b['path']
  assert sha(delivery/b['snapshot'])==b['sha256'],b['snapshot']
 assert check_roots(False)==read(HERE/'math-checks.json')
 run=ROOT/manifest['evidence_run'];rep=read(HERE/('representation-'+FINAL+'.json'));basis=read(ROOT/rep['basis_path'])
 definitions={g['name']:g for g in basis['generators']};atoms={sp.Symbol(n):vec(g['hlog_terms']) for n,g in definitions.items()}
 rows=read(ROOT/'benchmarks/tgv/families/t4-empty/results-r001.json')['results'];assert {r['orbit_id'][-2:] for r in rows}==set(IDS)
 count=nonzero=0;allcs={}
 for row in rows:
  p=ROOT/row['path'];assert sha(p)==row['sha256'];d=read(p);oid=d['orbit_id'][-2:];base=PAIRS.get(oid,oid)
  assert d['full_epsilon']['expression']==FULL and d['full_epsilon']['definitions']==DEFS and d['full_epsilon']['g']==oid
  assert d['basis_definitions']==definitions and d['formal']=='not_run'
  assert sha(ROOT/d['basis_source']['path'])==d['basis_source']['sha256']
  for label in ('prepared','canonical'):assert sha(ROOT/d['source'][label+'_path'])==d['source'][label+'_sha256']
  assert set(d['laurent']['coefficients'])=={str(k) for k in range(-6,1)}
  allcs[oid]={}
  for k,s in d['laurent']['coefficients'].items():
   v=ha.from_sympy(sp.sympify(s,locals={str(a):a for a in atoms}),atoms)
   assert v==vec(d['hlog_coefficients'][k])==ha.parse(run/f'I{base}q{k}.txt')
   assert all(weight(key)==6+int(k) for key in v)
   allcs[oid][k]=v;count+=1;nonzero+=bool(v)
  assert d['pole_order']==-min(int(k) for k,v in allcs[oid].items() if v)
 for target,source in PAIRS.items():assert allcs[target]==allcs[source]
 num=read(HERE/'numeric.json');reg=read(HERE/'regularizer-numeric.json');blocks=read(HERE/'series-checks.json');cons=read(HERE/'consistency-checks.json')
 assert num['passed'] and len(num['raw'])==16 and len(num['weighted_IBP'])==16
 assert all(r['passed'] and abs(r['difference'])<=r['tolerance'] for r in num['raw']+num['weighted_IBP'])
 assert reg['passed'] and len(reg['checks'])==8 and all(r['passed'] and abs(r['difference'])<=r['tolerance'] for r in reg['checks'])
 assert len(reg['independent_operator_point_certificates'])==6
 assert blocks['passed'] and len(blocks['blocks'])==12
 assert cons['passed'] and len(cons['checks'])==11
 assert count==56 and nonzero==32
 report={'passed':True,'integrals':8,'Laurent_coefficients':count,'nonzero_coefficients':nonzero,'exact_basis_Hlog_roundtrips':count,'uniform_weight_nonzero_coefficients':nonzero,'full_epsilon_equalities':4,'coefficient_equalities':28,'new_basis_generators':0,'raw_QMC_cases':16,'full_weighted_IBP_QMC_cases':16,'regularizer_QMC_cases':8,'exact_operator_point_checks':6,'finite_series_blocks':12,'order_encoding_assembly_consistency_checks':11,'source_and_snapshot_hash_checks':2*len(manifest['bindings']),'I_remainder_check':'not_run','formal':'not_run'}
 print(json.dumps(report,indent=2));return report
if __name__=='__main__':main()