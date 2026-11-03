"""Read-only verification of frozen T4[1] inputs, exact coefficients and evidence."""
import hashlib,json
import sympy as sp
from check_math import ROOT,HERE,read,check,IDS
from formulas import coefficients,basis_data
from publish import FULL_ORDERED,FULL_FACTOR,DEFS,PERIOD_DEFS
import wide_hlog as wh

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 delivery=ROOT/read(HERE/'current-delivery.json')['delivery'];manifest=read(delivery/'manifest.json')
 for b in manifest['bindings']:
  assert sha(ROOT/b['path'])==b['sha256'],b['path']
  assert sha(delivery/b['snapshot'])==b['sha256'],b['snapshot']
 assert check()==read(HERE/'math-checks.json')
 assert (ROOT/read(HERE/'current-global_ibp.json')['run']).relative_to(ROOT).as_posix()==manifest['evidence_run']
 rows=read(ROOT/'benchmarks/tgv/families/t4-opposite/results-r001.json')['results'];assert {r['orbit_id'][-2:] for r in rows}==set(IDS)
 cs=coefficients();rep=read(HERE/'representation.json');bp,basis=basis_data();count=nz=0;allcs={}
 for row in rows:
  p=ROOT/row['path'];assert sha(p)==row['sha256'];d=read(p);oid=d['orbit_id'][-2:]
  assert d['format']=='wilson-tgv-integral-result/2' and d['formal']=='not_run'
  assert d['full_epsilon']['definitions']==DEFS and d['full_epsilon']['g']==oid
  assert d['full_epsilon']['expression']==(FULL_FACTOR if oid in ('41','42') else FULL_ORDERED)
  assert d['period_definitions']==PERIOD_DEFS and d['hlog_encoding']['period_labels']==list(wh.PERIODS)
  defs=d['basis_definitions'];atoms={sp.Symbol(n):wh.vec(g['hlog_terms']) for n,g in defs.items()}
  for label in ['prepared','canonical']:assert sha(ROOT/d['source'][label+'_path'])==d['source'][label+'_sha256']
  assert sha(ROOT/d['basis_source']['path'])==d['basis_source']['sha256']
  if 'two_loop_source' in d:assert sha(ROOT/d['two_loop_source']['path'])==d['two_loop_source']['sha256']
  assert set(d['laurent']['coefficients'])=={str(k) for k in range(-6,1)}
  allcs[oid]={}
  for k,expr in d['laurent']['coefficients'].items():
   v=wh.from_sympy(sp.sympify(expr,locals={str(a):a for a in atoms}),atoms)
   assert v==wh.vec(d['hlog_coefficients'][k])==cs[oid][int(k)]
   assert all(wh.weight(key)==6+int(k) for key in v)
   allcs[oid][k]=v;count+=1;nz+=bool(v)
  assert d['pole_order']==-min(int(k) for k,v in allcs[oid].items() if v)
 assert allcs['27']==allcs['29'] and allcs['41']==allcs['42']
 om=wh.vec(rep['Omega6']['hlog_terms']);olds=[wh.vec(b['hlog_terms']) for b in basis['weight_bases']['6']]
 witness=wh.vec([rep['new_letter_witness']]);key=next(iter(witness));assert om[key]==witness[key] and all(v.get(key,0)==0 for v in olds)
 vectors=olds+[om];keys=sorted(set().union(*vectors));rank=sp.Matrix(len(keys),len(vectors),lambda i,j:vectors[j].get(keys[i],0)).rank();assert rank==15
 assert rep['new_quotient_dimension_for_this_batch']==1
 combo=wh.add(wh.scale(allcs['27']['0'],2),allcs['28']['0']);wh.lower(combo)
 # Successful Maple sessions must have no error, a completion marker, and the required files.
 for tag,mark in [('partial','PARTIAL_COMPLETE'),('weights','WEIGHTS_COMPLETE'),('global_ibp','GLOBAL_COMPLETE'),('outer28','OUTER_COMPLETE'),('omega_basis','OMEGA_COMPLETE'),('basis_reduce','BASIS_REDUCE_COMPLETE'),('cusp_export','CUSP_EXPORT_COMPLETE')]:
  r=ROOT/read(HERE/('current-'+tag+'.json'))['run'];log=(r/'maple.log').read_text(errors='replace')
  assert mark in log and 'Error,' not in log,(tag,log[-1000:]);assert read(r/'exit.json')['process_exit_code']==0
 outer=ROOT/read(HERE/'current-outer28.json')['run'];glob=ROOT/read(HERE/'current-global_ibp.json')['run']
 for tag in ['m1','zero']:assert wh.parse(outer/f'I28{tag}.txt')==wh.parse(glob/f'I28{tag}.txt')
 num=read(HERE/'numeric.json');assert num['passed'] and num['raw_root_sample_cases']==12
 for key,n in [('ordered_IBP',8),('original_root_comparisons',6),('factorization',2),('reflection',2)]:
  assert len(num[key])==n and all(r['passed'] and abs(r['difference'])<=r['tolerance'] for r in num[key])
 series=read(HERE/'series-checks.json');assert series['passed'] and len(series['finite_blocks'])==12
 omega=read(HERE/'omega-numeric.json');assert omega['passed'] and len(omega['Omega_definition_checks'])==2 and len(omega['period_checks'])==2
 cusp=read(HERE/'cusp-numeric.json');assert cusp['passed'] and len(cusp['raw_boundary_checks'])==2 and len(cusp['cusp_remainder_samples'])==3
 assert all(r['passed'] and abs(r['difference'])<=r['tolerance'] for r in cusp['raw_boundary_checks'])
 assert all(r['returncode']==0 for r in read(HERE/'reuse.json').values())
 assert count==35 and nz==12
 report={'passed':True,'integrals':5,'Laurent_positions':count,'nonzero_coefficients':nz,'exact_basis_Hlog_roundtrips':count,'coefficient_equalities':14,'raw_root_sample_cases':12,'ordered_IBP_comparisons':8,'root_reduction_comparisons':6,'factorization_comparisons':2,'reflection_comparisons':2,'finite_series_blocks':12,'cusp_boundary_comparisons':2,'cusp_remainder_samples':3,'independent_Omega_definition_checks':2,'independent_period_checks':2,'new_weight_six_generators':1,'augmented_weight_six_rank':rank,'source_and_snapshot_hash_checks':2*len(manifest['bindings']),'full_I_remainder_check':'not_run','formal':'not_run'}
 print(json.dumps(report,indent=2));return report
if __name__=='__main__':main()