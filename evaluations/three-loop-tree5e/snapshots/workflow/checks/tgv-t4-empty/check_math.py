"""Exact input/IBP checks for T4[], independent of the computed coefficients."""
from pathlib import Path
import importlib.util,json,hashlib,copy
import sympy as sp
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent

def loadmod(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
parse=loadmod('audit',ROOT/'workflow/checks/tgv-audit/audit.py').parse
bridge=loadmod('t2bridge',ROOT/'workflow/checks/tgv-t2-empty/check_math.py')
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def eq(a,b):assert sp.expand(a-b)==0,sp.factor(a-b)
IDS=('13','14','15','33','34','35','40','43')
PAIRS={'35':'13','34':'14','33':'15','43':'40'}
x,y,w,z,e=sp.symbols('x y w z eps');a1,a2,a3=sp.symbols('a1 a2 a3')
S=a1*a2*y*(1-x)+a1*a3*((1-x)*(1-w)+z*x*w)+z*a2*a3*w*(1-y)
P=(a1*a2*(-(1+z)*(1-x)+z*y)+a1*a3*(1+3*z-(1+z)*x+z*(1+z)*w)+z*a2*a3*(1-y-(1+z)*w))/4
R=(1+z)*a3*(a1*(1-x+z*w)-z*a2*w)/2

def reflect(target,source,roots,swap=True,reverse=True,sgn=1):
 tt=target['representative']['canonical_topology'];st=source['representative']['canonical_topology']
 tc=bridge.chains(tt);sc=bridge.chains(st);rho={1:3,2:2,3:1,4:4}
 zz=sp.Symbol('z1');als=sp.symbols('alpha1:4');sub={zz:1/zz}
 sub.update({als[0]:als[2] if swap else als[0],als[1]:als[1],als[2]:als[0] if swap else als[2]})
 for edge,vs in tc.items():
  ws=list(reversed(sc[rho[edge]])) if reverse else sc[rho[edge]]
  assert len(vs)==len(ws)
  sub.update({sp.Symbol(v):1-sp.Symbol(wv) for v,wv in zip(vs,ws)})
  gaps=[sp.Integer(0)]+[sp.Symbol(v) for v in vs]+[sp.Integer(1)]
  orig=[sp.Integer(0)]+[sp.Symbol(v) for v in sc[rho[edge]]]+[sp.Integer(1)]
  assert [sp.expand((gaps[i+1]-gaps[i]).subs(sub,simultaneous=True)) for i in range(len(gaps)-1)]==[orig[i+1]-orig[i] for i in reversed(range(len(orig)-1))]
 tj,sj=target['J_raw_s'],source['J_raw_s']
 for key,power in [('S_polynomial',1),('pnum',2)]:
  eq((sgn if key=='pnum' else 1)*parse(tj['tgv'][key]).subs(sub,simultaneous=True)*zz**power,parse(sj['tgv'][key]))
 eq(parse(tj['f_polynomials']['f1']).subs(sub,simultaneous=True)*zz,parse(sj['f_polynomials']['f1']))
 eq(parse(tj['pure']['lorentz_z']).subs(sub,simultaneous=True)*zz,parse(sj['pure']['lorentz_z']))
 assert sp.simplify(parse(roots[target['orbit_id'][-2:]]['C']['expression'])/parse(roots[source['orbit_id'][-2:]]['C']['expression']))==1
 jac=sp.det(sp.Matrix([sub[sp.Symbol(v)] for v in tt['bounds']]).jacobian([sp.Symbol(v) for v in st['bounds']]))
 ajac=sp.det(sp.Matrix([sub[a].subs(als[2],1-als[0]-als[1]) for a in als[:2]]).jacobian(als[:2]))
 assert abs(jac)==abs(ajac)==1
 return {'target':target['orbit_id'],'source':source['orbit_id'],'substitution':{str(k):str(v) for k,v in sub.items()},'S_ratio':'1/z','P_ratio':'1/z^2','J_relation':'J_target(z)=z^(3*eps)*J_source(1/z)','I_equal':True}

def main(write=True):
 path=ROOT/'benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json'
 roots={g['orbit_id'][-2:]:g for g in read(path)['orbits']}
 src={g['orbit_id'][-2:]:g for g in read(ROOT/'benchmarks/tgv/5e-3g/sources/canonical-r001.json')['orbits']}
 al=sp.symbols('alpha1:4');zs=sp.Symbol('z1')
 expected={
 '13':{1:'t1',2:'t2 p1a',3:'t3 p1b'},'14':{1:'t1',2:'p1a t2',3:'t3 p1b'},'15':{1:'t1',2:'p1a t2',3:'p1b t3'},
 '33':{1:'t1 p1a',2:'t2 p1b',3:'t3'},'34':{1:'p1a t1',2:'t2 p1b',3:'t3'},'35':{1:'p1a t1',2:'p1b t2',3:'t3'},
 '40':{1:'t1 p1b',2:'t2',3:'t3',4:'p1a'},'43':{1:'t1',2:'t2',3:'p1b t3',4:'p1a'}}
 cat=read(ROOT/'benchmarks/tgv/families/catalog-r001.json')
 assert set(next(q['mixed_members'] for q in cat['groups'] if q['group_id']=='t4-empty'))=={'5e-3g-'+i for i in IDS}
 for oid in IDS:
  g=roots[oid];j=g['J_raw_s'];top=g['representative']['canonical_topology'];legs=top['tgv_legs'];pair=top['pure_pairings'][0]
  assert j==src[oid]['J_raw_s'];assert j['bounds']==top['bounds']
  order=j['iterated_order_outer_to_inner']
  for k,vv in enumerate(order):assert parse(j['bounds'][vv][1]).free_symbols<={sp.Symbol(v) for v in order[:k]}
  tok={leg['s_var']:'t'+str(k+1) for k,leg in enumerate(legs)};tok.update({pair['s_a']:'p1a',pair['s_b']:'p1b'})
  words={ed:' '.join(tok[v] for v in vs) for ed,vs in bridge.chains(top).items() if vs};assert words==expected[oid]
  sub={sp.Symbol(leg['s_var']):v for leg,v in zip(legs,(x,y,w))};sub.update(dict(zip(al,(a1,a2,a3))));sub[zs]=z
  eq(parse(j['tgv']['S_polynomial']).subs(sub),z*S.subs(z,1/z))
  eq(parse(j['tgv']['pnum']).subs(sub),z*z*P.subs(z,1/z))
  assert g['C']['fixed_basis_coordinates']==['-1/4','1/2','0'];assert g['C']['euler_compensation_count']==1
  assert g['c4_symmetry']['orbit_multiplicity']==4
  CA,CF,as_,v14,I0=sp.symbols('CA CF alphaS v14 I0')
  common=CA**2*CF/16*(as_/sp.pi)**3*sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(2-2*e)*sp.pi**(3*e)*(I0-2*v14)**(3*e)
  assert sp.simplify(parse(g['C']['expression']).subs(zs,z)/common)==1
 reflections=[reflect(src[a],src[b],roots) for a,b in PAIRS.items()]
 negative=[]
 for label,kw in [('missing_alpha_swap',{'swap':False}),('missing_order_reversal',{'reverse':False}),('wrong_P_sign',{'sgn':-1})]:
  try:reflect(src['35'],src['13'],roots,**kw)
  except AssertionError:negative.append(label)
  else:raise AssertionError(label)
 bad=copy.deepcopy(roots);bad['35']['C']['expression']='2*('+bad['35']['C']['expression']+')'
 try:reflect(src['35'],src['13'],bad)
 except AssertionError:negative.append('doubled_C')
 else:raise AssertionError('doubled_C')
 eq(P,(-z*sp.diff(S,x)-(1+z)*sp.diff(S,y)+sp.diff(S,w))/4+R)
 assert sp.expand(R)!=0
 lam=1-2*e
 for v in (x,y,w):assert sp.simplify(sp.diff(S**(-lam),v)+lam*sp.diff(S,v)*S**(-lam-1))==0
 u,b=sp.symbols('u b',positive=True)
 sub={a1:u,a2:(1-u)*b,a3:(1-u)*(1-b)}
 eq(S.subs(w,0).subs(sub,simultaneous=True),u*(1-u)*(1-x)*(1-b+b*y))
 sub={a3:u,a1:(1-u)*(1-b),a2:(1-u)*b}
 eq(S.subs(x,1).subs(sub,simultaneous=True),z*u*(1-u)*w*(1-b*y))
 # Primitive certificates retain endpoint differences; the positive-e limits give the stated weights.
 assert sp.simplify(sp.diff(w**e/e,w)-w**(e-1))==0
 assert sp.simplify(sp.diff(-(1-y)**e/e,y)-(1-y)**(e-1))==0
 eq(S.subs({x:1-w,y:1-y,w:1-x,a1:a3,a3:a1,z:1/z},simultaneous=True),S/z)
 eq(P.subs({x:1-w,y:1-y,w:1-x,a1:a3,a3:a1,z:1/z},simultaneous=True),P/z**2)
 report={'passed':True,'roots':list(IDS),'root_hash':hashlib.sha256(path.read_bytes()).hexdigest(),'reflections':reflections,'negative_examples':negative,'IBP_P_residual':'0','recoil_nonzero':True,'D_to_F_S':'S_F(z)=z*S_D(1/z)','D_to_F_P':'P_F(z)=z^2*P_D(1/z)','formal':'not_run'}
 if write:dump(HERE/'math-checks.json',report);print(json.dumps(report,ensure_ascii=False))
 return report
if __name__=='__main__':main()