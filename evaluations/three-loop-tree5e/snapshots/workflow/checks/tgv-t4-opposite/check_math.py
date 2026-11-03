"""Exact T4[1] roots, reflected domains, factorization and C/I bridges."""
from pathlib import Path
import importlib.util,json,hashlib,copy
import sympy as sp
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent

def module(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
old=module('t4_empty_math',ROOT/'workflow/checks/tgv-t4-empty/check_math.py')
parse=old.parse;eq=old.eq;S,P=old.S,old.P;x,y,w,z,e=old.x,old.y,old.w,old.z,old.e;a1,a2,a3=old.a1,old.a2,old.a3
IDS=('27','28','29','41','42');PAIR={'29':'27','42':'41'}
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def dump(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def check():
 rp=ROOT/'benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json';srcp=ROOT/'benchmarks/tgv/5e-3g/sources/canonical-r001.json'
 roots={g['orbit_id'][-2:]:g for g in read(rp)['orbits']};src={g['orbit_id'][-2:]:g for g in read(srcp)['orbits']}
 cat=read(ROOT/'benchmarks/tgv/families/catalog-r001.json')
 members=set().union(*(set(g['mixed_members']) for g in cat['groups'] if g['group_id'] in ['t4-opposite-13','t4-opposite-24']))
 assert members=={'5e-3g-'+i for i in IDS}
 expected={'27':{1:'t1 p1a',2:'t2',3:'t3 p1b'},'28':{1:'t1 p1a',2:'t2',3:'p1b t3'},'29':{1:'p1a t1',2:'t2',3:'p1b t3'},'41':{1:'t1',2:'t2 p1b',3:'t3',4:'p1a'},'42':{1:'t1',2:'p1b t2',3:'t3',4:'p1a'}}
 zz=sp.Symbol('z1');al=sp.symbols('alpha1:4');pp,qq=sp.symbols('p q')
 for oid in IDS:
  g=roots[oid];j=g['J_raw_s'];top=g['representative']['canonical_topology'];legs=top['tgv_legs'];pair=top['pure_pairings'][0]
  assert j==src[oid]['J_raw_s'] and j['bounds']==top['bounds']
  order=j['iterated_order_outer_to_inner']
  for k,v in enumerate(order):assert parse(j['bounds'][v][1]).free_symbols<={sp.Symbol(n) for n in order[:k]}
  tok={leg['s_var']:'t'+str(k+1) for k,leg in enumerate(legs)};tok.update({pair['s_a']:'p1a',pair['s_b']:'p1b'})
  assert {ed:' '.join(tok[n] for n in vs) for ed,vs in old.bridge.chains(top).items() if vs}==expected[oid]
  sub={sp.Symbol(leg['s_var']):v for leg,v in zip(legs,(x,y,w))};sub.update(dict(zip(al,(a1,a2,a3))));sub[zz]=z
  eq(parse(j['tgv']['S_polynomial']).subs(sub),z*S.subs(z,1/z));eq(parse(j['tgv']['pnum']).subs(sub),z*z*P.subs(z,1/z))
  sub.update({sp.Symbol(pair['s_a']):pp,sp.Symbol(pair['s_b']):qq})
  wanted=z*(1-pp)*(1-qq)+pp*qq if oid in ('27','28','29') else (1-pp)*(1-qq)+z*pp*qq
  eq(parse(j['f_polynomials']['f1']).subs(sub),wanted);eq(parse(j['pure']['lorentz_z']).subs(sub),-(1+z))
  assert g['C']['fixed_basis_coordinates']==['-1/4','1/2','0'] and g['C']['euler_compensation_count']==1 and g['c4_symmetry']['orbit_multiplicity']==4
  CA,CF,as_,v14,I0=sp.symbols('CA CF alphaS v14 I0')
  common=CA**2*CF/16*(as_/sp.pi)**3*sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(2-2*e)*sp.pi**(3*e)*(I0-2*v14)**(3*e)
  assert sp.simplify(parse(g['C']['expression']).subs(zz,z)/common)==1
 refs=[old.reflect(src[t],src[s],roots) for t,s in [('29','27'),('28','28'),('42','41')]]
 negative=[]
 for name,kwargs in [('alpha_swap_missing',{'swap':False}),('order_reversal_missing',{'reverse':False}),('wrong_P_sign',{'sgn':-1})]:
  try:old.reflect(src['29'],src['27'],roots,**kwargs)
  except AssertionError:negative.append(name)
  else:raise AssertionError(name)
 bad=copy.deepcopy(roots);bad['29']['C']['expression']='2*('+bad['29']['C']['expression']+')'
 try:old.reflect(src['29'],src['27'],bad)
 except AssertionError:negative.append('wrong_C_factor')
 else:raise AssertionError('wrong_C_factor')
 # Complementary y<q and q<y domains; independent monomial moments.
 partition=[]
 for n in range(4):
  for m in range(4):
   lo=sp.Rational(1,(n+1)*(n+m+2));hi=sp.Rational(1,(m+1)*(n+m+2));prod=sp.Rational(1,(n+1)*(m+1));assert lo+hi==prod
   partition.append({'y_power':n,'q_power':m,'difference':'0'})
 # Match the reused two-loop T4 input, rather than just its result label.
 bare_path=ROOT/'benchmarks/tgv/3e-3g/prepared/3e-3g.json'
 bare=next(g for g in read(bare_path)['orbits'] if g['orbit_id']=='3e-3g-04')
 bs={sp.Symbol('s21'):x,sp.Symbol('s31'):y,sp.Symbol('s41'):w,sp.Symbol('z1'):z,**dict(zip(al,(a1,a2,a3)))}
 eq(parse(bare['J']['S']).subs(bs),S);eq(parse(bare['J']['numerator']).subs(bs),P)
 assert all(v==['0','1'] for v in bare['J']['bounds'].values())
 old_result=read(ROOT/'benchmarks/tgv/3e-3g/results/3e-3g-04/r001.json')
 assert sha(bare_path)==old_result['input']['package_sha256']
 # Gamma and crossing factors in the factorized physical I.
 g2=sp.exp(-2*sp.EulerGamma*e)*sp.gamma(1-2*e);g3=sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(1-2*e)
 assert sp.simplify(g3/g2-sp.exp(-sp.EulerGamma*e)*sp.gamma(1-e))==0
 u=sp.Symbol('u',positive=True)
 assert sp.factor(u*(1+u**4)/(1+u**6)-(u**2+u**-2)/(u**3+u**-3))==0
 # Compare the physical-root normalization to the factored expression.
 jc,jp,la,gg2,gg3=sp.symbols('jc jp la gg2 gg3',nonzero=True)
 from_root=la*gg3*jc*jp/(8*(1+u**6))
 old_I=la*gg2*jc/(1+u**4);H0=jp/u
 from_product=gg3/gg2*(u**2+u**-2)/(u**3+u**-3)*old_I*H0/8
 assert sp.factor(from_root-from_product)==0
 assert sp.factor(from_root-8*from_product)!=0
 negative.append('factorized_I_missing_1_over_8')
 # Scalar H endpoint reflection and optimized recoil polynomial identity.
 eq(S.subs({x:1-w,y:1-y,w:1-x,a1:a3,a3:a1,z:1/z},simultaneous=True),S/z)
 eq(P,(-z*sp.diff(S,x)-(1+z)*sp.diff(S,y)+sp.diff(S,w))/4+old.R)
 return {'passed':True,'roots':list(IDS),'prepared_sha256':sha(rp),'canonical_sha256':sha(srcp),'reflections':refs,'negative_examples':negative,'partition_monomial_checks':partition,'factorization':'J41+J42=J_T4*J_pure; I41=I42','I_factorized':'exp(-EulerGamma*e)*Gamma(1-e)*cosh(e*log(z))/(8*cosh(3*e*log(z)/2))*I_T4(z,e)*H0(e,z)','full_D_to_F_bridge':'J_F(z)=z^(3*e)*J_D(1/z), for 27/28/29','bare_T4_input_sha256':sha(bare_path),'bare_T4_SP_match':True,'formal':'not_run'}
if __name__=='__main__':d=check();dump(HERE/'math-checks.json',d);print(json.dumps(d,ensure_ascii=False))