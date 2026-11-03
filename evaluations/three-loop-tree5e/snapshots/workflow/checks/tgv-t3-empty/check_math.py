"""Check original roots, reflections, density reduction and integrated polynomial blocks."""
import copy,hashlib,math
import sympy as sp
from formulas import ROOT,HERE,read,dump,run,module,coefficients,ha,weight
audit=module("t3_audit",ROOT/"workflow/checks/tgv-audit/audit.py");parse=audit.parse
bridge=module("t3_bridge",ROOT/"workflow/checks/tgv-t2-empty/check_math.py")
PAIR={"39":"09","38":"10","30":"26"}
IDS=("09","10","26","30","38","39")
def eq(a,b):assert sp.expand(a-b)==0,(a,b)
def records(path):return {d["orbit_id"][-2:]:d for d in read(path)["orbits"]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def reflection(target,source,rt,rs,swap=True,reverse=True,sign=1):
 tt=target["representative"]["canonical_topology"];st=source["representative"]["canonical_topology"]
 tc=bridge.chains(tt);sc=bridge.chains(st);rho={1:1,2:4,3:3,4:2}
 z=sp.Symbol("z1");als=sp.symbols("alpha1:4");sub={z:1/z}
 sub.update({als[0]:als[0],als[1]:als[2] if swap else als[1],als[2]:als[1] if swap else als[2]})
 for edge,vs in tc.items():
  ws=list(reversed(sc[rho[edge]])) if reverse else sc[rho[edge]]
  assert len(vs)==len(ws)
  sub.update({sp.Symbol(v):1-sp.Symbol(w) for v,w in zip(vs,ws)})
  gaps=[sp.Symbol(v) for v in vs]
  gaps=[sp.Integer(0)]+gaps+[sp.Integer(1)]
  orig=[sp.Integer(0)]+[sp.Symbol(w) for w in sc[rho[edge]]]+[sp.Integer(1)]
  assert [sp.expand((gaps[i+1]-gaps[i]).subs(sub,simultaneous=True)) for i in range(len(gaps)-1)]==[orig[i+1]-orig[i] for i in reversed(range(len(orig)-1))]
 tj,sj=target["J_raw_s"],source["J_raw_s"]
 eq(parse(tj["tgv"]["S_polynomial"]).subs(sub,simultaneous=True)*z,parse(sj["tgv"]["S_polynomial"]))
 eq(sign*parse(tj["tgv"]["pnum"]).subs(sub,simultaneous=True)*z*z,parse(sj["tgv"]["pnum"]))
 for expr in ("f1",):
  eq(parse(tj["f_polynomials"][expr]).subs(sub,simultaneous=True)*z,parse(sj["f_polynomials"][expr]))
 eq(parse(tj["pure"]["lorentz_z"]).subs(sub,simultaneous=True)*z,parse(sj["pure"]["lorentz_z"]))
 assert sp.simplify(parse(rt["C"]["expression"])/parse(rs["C"]["expression"]))==1
 jac=sp.det(sp.Matrix([sub[sp.Symbol(v)] for v in tt["bounds"]]).jacobian([sp.Symbol(v) for v in st["bounds"]]))
 ajac=sp.det(sp.Matrix([sub[a].subs(als[2],1-als[0]-als[1]) for a in als[:2]]).jacobian(als[:2]))
 assert abs(jac)==abs(ajac)==1
 return {"target":target["orbit_id"],"source":source["orbit_id"],"substitution":{str(k):str(v) for k,v in sub.items()},"S_ratio":"1/z","P_ratio":"1/z^2","pure_f_ratio":"1/z","pure_lorentz_ratio":"1/z","J_relation":"J_target(z)=z^(3*eps)*J_source(1/z)","I_equal":True,"absolute_Jacobians":[1,1]}
def check():
 src=records(ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json");roots=records(ROOT/"benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json")
 cat=read(ROOT/"benchmarks/tgv/families/catalog-r001.json")
 assert set(next(g["mixed_members"] for g in cat["groups"] if g["group_id"]=="t3-empty"))=={"5e-3g-"+v for v in IDS}
 e,z=sp.symbols("eps z1");t,x,y,b,a,v,w=sp.symbols("t x y b a v w");als=sp.symbols("alpha1:4")
 F=lambda u:z*(1-t)*(1-u)+t*u
 expected={"09":{1:"t1",2:"p1a",3:"t2 t3 p1b"},"10":{1:"t1",2:"p1a",3:"t2 p1b t3"},
 "26":{1:"t1 p1b",3:"t2 t3",4:"p1a"},"30":{1:"p1a t1",2:"p1b",3:"t2 t3"},
 "38":{1:"t1",3:"t2 p1a t3",4:"p1b"},"39":{1:"t1",3:"p1a t2 t3",4:"p1b"}}
 for oid in IDS:
  raw=src[oid];d=roots[oid];top=raw["representative"]["canonical_topology"];legs=top["tgv_legs"];pair=top["pure_pairings"][0]
  assert d["J_raw_s"]==raw["J_raw_s"];j=d["J_raw_s"];assert j["bounds"]==top["bounds"]
  order=j["iterated_order_outer_to_inner"]
  for i,name in enumerate(order):assert parse(j["bounds"][name][1]).free_symbols<={sp.Symbol(q) for q in order[:i]}
  token={leg["s_var"]:"t"+str(i+1) for i,leg in enumerate(legs)}
  token.update({pair["s_a"]:"p1a",pair["s_b"]:"p1b"})
  words={edge:" ".join(token[n] for n in names) for edge,names in bridge.chains(top).items() if names}
  assert words==expected[oid]
  sub={sp.Symbol(leg["s_var"]):var for leg,var in zip(legs,(t,x,y))}
  eq(parse(j["tgv"]["S_polynomial"]).subs(sub),als[0]*(als[1]*F(x)+als[2]*F(y)))
  eq(parse(j["tgv"]["pnum"]).subs(sub),-(1+z)*((1+z)*t-z)*als[0]*(als[1]-als[2])/4)
  ep={pair["edge_a"]:sp.Symbol(pair["s_a"]),pair["edge_b"]:sp.Symbol(pair["s_b"])}
  if oid in ("09","10"):ps={ep[2]:v,ep[3]:w};pf=(1-v)*w;lor=1
  elif oid=="26":ps={ep[4]:v,ep[1]:w};pf=(1-v)*w;lor=1
  elif oid=="30":ps={ep[1]:w,ep[2]:v};pf=z*(1-w)*v;lor=z
  else:ps={ep[3]:w,ep[4]:v};pf=z*(1-w)*v;lor=z
  eq(parse(j["f_polynomials"]["f1"]).subs(ps),pf);eq(parse(j["pure"]["lorentz_z"]),lor)
  assert d["C"]["fixed_basis_coordinates"]==["-1/4","1/2","0"] and d["C"]["euler_compensation_count"]==1
  assert d["c4_symmetry"]["orbit_multiplicity"]==4
  CA,CF,as_,v14,I0=sp.symbols("CA CF alphaS v14 I0")
  common=CA**2*CF/16*(as_/sp.pi)**3*sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(2-2*e)*sp.pi**(3*e)*(I0-2*v14)**(3*e)
  assert sp.simplify(parse(d["C"]["expression"])/common)==1
 reflected=[reflection(src[k],src[j],roots[k],roots[j]) for k,j in PAIR.items()]
 negatives=[]
 for label,kwargs in [("missing_alpha_swap",{"swap":False}),("missing_order_reversal",{"reverse":False}),("flipped_P",{"sign":-1})]:
  try:reflection(src["39"],src["09"],roots["39"],roots["09"],**kwargs)
  except AssertionError:negatives.append(label)
  else:raise AssertionError(label)
 altered=copy.deepcopy(roots["39"]);altered["C"]["expression"]="2*("+altered["C"]["expression"]+")"
 try:reflection(src["39"],src["09"],altered,roots["09"])
 except AssertionError:negatives.append("doubled_C")
 else:raise AssertionError("doubled_C")
 # Alpha scale and ordered IBP generator.
 alpha={als[0]:a,als[1]:(1-a)*b,als[2]:(1-a)*(1-b)}
 eq(sp.det(sp.Matrix([alpha[als[0]],alpha[als[1]]]).jacobian([a,b])),1-a)
 assert sp.simplify(sp.diff(a**e/e,a)/a**(e-1))==1
 Q=F(b*x+(1-b)*y);lam=1-2*e
 assert sp.simplify((sp.diff(Q**(-lam),x)-sp.diff(Q**(-lam),y))/Q**(-lam-1)+lam*(2*b-1)*((1+z)*t-z))==0
 # Density change u=b*x+(1-b)*y, with exact Jacobians.
 u=sp.Symbol("u")
 assert sp.diff((u-(1-b)*y)/b,u)==1/b
 assert sp.diff((u-b*x)/(1-b),u)==1/(1-b)
 # Differentiate moving density limits; resulting terms give the stated incomplete beta / 2F1 kernels.
 assert sp.simplify(sp.diff(((u-1+b)/b)**(1+e)/(1+e),u)/(((u-1+b)/b)**e/b))==1
 # Incomplete beta primitive for D_y'.
 bb=b**(-e)*(1-b)**(-2*e)
 assert sp.factor(sp.diff(bb,b)/bb-e*(3*b-1)/(b*(1-b)))==0
 # Endpoint monomial primitive and integrated-series moments.
 assert sp.simplify(sp.diff(w**e/e,w)/w**(e-1))==1
 density_moments=[]
 B0=sp.Symbol("B0")
 def beta_ratio(i,j):return sp.rf(1-e,i)*sp.rf(1-e,j)/sp.rf(2-2*e,i+j)
 for n in range(1,7):
  bs=sum(sp.binomial(n-1,j)*(2*beta_ratio(j+1,n-1-j)-beta_ratio(j,n-1-j))/(j+1) for j in range(n))
  for me in (0,e):
   direct=-n*bs/(n+1+me)
   target=(beta_ratio(n,0)-1/sp.Integer(n+1))/e if me==0 else ((n+1)*beta_ratio(n,0)-1)/(e*(n+1+e))
   assert sp.factor(direct-target)==0
   density_moments.append({"n":n,"weight":"1" if me==0 else "y^eps","passed":True})
 # No relation inferred from the taxonomy: three independent leading Hlog vectors.
 cs=coefficients();vecs=[cs[x][-2] for x in ("09","10","26")];keys=sorted(set().union(*vecs))
 rank=sp.Matrix(len(keys),3,lambda i,j:vecs[j].get(keys[i],0)).rank();assert rank==3
 for k,j in PAIR.items():assert cs[k]==cs[j]
 assert all(all(weight(key)==6+k for key in v) for k,v in cs["10"].items())
 return {"passed":True,"roots":list(IDS),"source_and_C_checks":6,"reflections":reflected,
 "negative_cases_rejected":negatives,"alpha_Jacobian":"1-a","alpha_scale_integral":"1/eps",
 "ordered_IBP_generator":True,"density_moment_checks":density_moments,
 "full_epsilon_pair_equalities":3,"coefficient_equalities":21,
 "leading_coefficient_Q_rank":rank,"10_uniform_weight":True,
 "pole_orders":{oid:-min(k for k,v in c.items() if v) for oid,c in cs.items()},
 "formal":"not_run"}
if __name__=="__main__":
 r=check();dump(run()/"math-checks.json",r);print(__import__("json").dumps(r,indent=2))
