"""Hlog coefficients, factorized channel, and existing-basis comparisons."""
import sys
import sympy as sp
from check_math import ROOT,HERE,read,dump,module
import wide_hlog as wh
from fractions import Fraction as Q
f=module('t4_base_formulas',ROOT/'workflow/checks/tgv-t4-empty/formulas.py')
ha,vec,weight,fit=f.ha,f.vec,f.weight,f.fit

def basis_data():
 p=ROOT/read(ROOT/'workflow/checks/unified-weight-six/current-run.json')['run']/'basis-and-conversions.json'
 return p,read(p)

def factorized():
 bp,b=basis_data();atoms={sp.Symbol(g['name']):vec(g['hlog_terms']) for g in b['generators']}
 M0,M1,M2,Z2,Z3=sp.symbols('M_0_1 M_1_1 M_2_1 zeta_2 zeta_3');e=sp.Symbol('eps')
 ib={-2:Z2/8,-1:(Z3-M1)/8,0:-M2/8+M0**2/48+sp.Rational(9,10)*Z2**2}
 old=read(ROOT/'benchmarks/tgv/3e-3g/results/3e-3g-04/r001.json')['laurent']['coefficients']
 for k,v in ib.items():assert ha.from_sympy(v,atoms)==vec(old[str(k)]['hlog_terms'])
 H=M0+M1*e/2+(M2/8-M0**2/24)*e**2
 pref=1+(sp.Rational(5,4)*M0+sp.Rational(17,4)*Z2)*e**2
 prod=sp.expand(H*pref*sum(v*e**k for k,v in ib.items())/8)
 exprs={k:prod.coeff(e,k) if k>=-2 else sp.Integer(0) for k in range(-6,1)}
 return {k:ha.from_sympy(v,atoms) for k,v in exprs.items()},exprs

def coefficients(preview=False):
 out={}
 tag='outer28' if preview else 'global_ibp';run=ROOT/read(HERE/('current-'+tag+'.json'))['run']
 ids=('28',) if preview else ('27','28')
 for oid in ids:
  cs={k:{} for k in range(-6,1)};cs[-1]=wh.parse(run/f'I{oid}m1.txt');cs[0]=wh.parse(run/f'I{oid}zero.txt');out[oid]=cs
 if not preview:out['29']=out['27']
 c41,exprs=factorized();out['41']={k:wh.lift(v) for k,v in c41.items()};out['42']=out['41']
 return out

def main(preview=False):
 bp,b=basis_data();cs=coefficients(preview);out={}
 omega_run=ROOT/read(HERE/'current-omega_basis.json')['run'];omega=wh.parse(omega_run/'Omega6.txt')
 phi=ha.parse(omega_run/'Phi5.txt');assert all(weight(k)==5 for k in phi)
 assert all(wh.weight(k)==6 for k in omega)
 omega_coeff={'27':Q(-1,64),'29':Q(-1,64),'28':Q(1,32),'41':Q(0),'42':Q(0)}
 for oid,cols in cs.items():
  out[oid]={}
  for k,v in cols.items():
   oc=omega_coeff[oid] if k==0 else Q(0)
   residual=wh.add(v,wh.scale(omega,-oc));old=wh.lower(residual);fits={}
   for w in range(7):
    part={key:c for key,c in old.items() if weight(key)==w}
    if part:fits[str(w)]=fit(part,b['weight_bases'][str(w)])
   assert None not in fits.values(),(oid,k)
   if oc:fits.setdefault('6',[]).append({'coefficient':str(oc),'factors':['Omega6']})
   out[oid][str(k)]={'hlog_terms':wh.serialize(v),'basis_fit':fits,'old_basis_spans':not bool(oc)}
   print(oid,k,len(v),fits,flush=True)
 witness=next((key,c) for key,c in omega.items() if Q(-1,2) in key[0] and key[1]==wh.ZERO)
 oldvs=[wh.vec(x['hlog_terms']) for x in b['weight_bases']['6']];vs=oldvs+[omega];keys=sorted(set().union(*vs))
 A=sp.Matrix(len(keys),len(vs),lambda i,j:vs[j].get(keys[i],0));assert A.rank()==15
 assert all(v.get(witness[0],0)==0 for v in oldvs)
 candidates=list(b['weight_bases']['5'])
 for item in b['weight_bases']['4']:candidates.append({'factors':['ell']+item['factors'],'hlog_terms':ha.serialize(ha.mul(ha.LOG,vec(item['hlog_terms'])))})
 report={'coefficients':out,'basis_path':bp.relative_to(ROOT).as_posix(),'all_fit_with_Omega6':True,'old_weight_six_rank':14,'augmented_weight_six_rank':15,'new_quotient_dimension_for_this_batch':1,
 'new_letter_witness':wh.serialize({witness[0]:witness[1]})[0],'period_labels':list(wh.PERIODS),'period_weights':list(wh.WEIGHTS),
 'Omega6':{'definition':'T(z)+T(1/z), T(z)=integral_0^z Phi5(t)/(t+1/2) dt','hlog_terms':wh.serialize(omega),'Phi5_hlog_terms':ha.serialize(phi),'Phi5_old_even_and_log_times_weight4_fit':fit(phi,candidates)},
 'scope':'exact rational span of the stated Hlog words/formal period coordinates; not a full-epsilon IBP minimality theorem','formal':'not_run'}
 dump(HERE/('representation-preview.json' if preview else 'representation.json'),report)
if __name__=='__main__':main('--preview28' in sys.argv)
