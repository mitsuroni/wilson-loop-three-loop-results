"""Connect five ordered source roots to T3 density, reflection and factorization formulas."""
import copy,hashlib
import sympy as sp
from formulas import ROOT,HERE,read,dump,run,module,coefficients,ha,weight
audit=module("t3opp_audit",ROOT/"workflow/checks/tgv-audit/audit.py");parse=audit.parse
prior=module("t3opp_prior",ROOT/"workflow/checks/tgv-t3-empty/check_math.py")
IDS=("22","23","24","25","44")
def eq(a,b):assert sp.expand(a-b)==0,(a,b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def records(path):return {d["orbit_id"][-2:]:d for d in read(path)["orbits"]}
def check():
 src=records(ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json");roots=records(ROOT/"benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json")
 cat=read(ROOT/"benchmarks/tgv/families/catalog-r001.json")
 members=sum([g["mixed_members"] for g in cat["groups"] if g["group_id"] in ("t3-opposite-13","t3-opposite-24")],[])
 assert set(members)=={"5e-3g-"+i for i in IDS}
 e,z=sp.symbols("eps z1");t,x,y,p,q,b,u,a=sp.symbols("t x y p q b u a");als=sp.symbols("alpha1:4")
 F=lambda t,u:z*(1-t)*(1-u)+t*u
 expected={"22":{1:"t1 p1a",3:"t2 t3 p1b"},"23":{1:"t1 p1a",3:"t2 p1b t3"},
 "24":{1:"p1a t1",3:"t2 p1b t3"},"25":{1:"p1a t1",3:"p1b t2 t3"},"44":{1:"t1",2:"p1b",3:"t2 t3",4:"p1a"}}
 for oid in IDS:
  raw=src[oid];d=roots[oid];top=raw["representative"]["canonical_topology"];legs=top["tgv_legs"];pair=top["pure_pairings"][0];j=d["J_raw_s"]
  assert j==raw["J_raw_s"] and j["bounds"]==top["bounds"]
  order=j["iterated_order_outer_to_inner"];assert set(order)==set(j["bounds"])
  for i,name in enumerate(order):assert parse(j["bounds"][name][1]).free_symbols<={sp.Symbol(q) for q in order[:i]}
  token={leg["s_var"]:"t"+str(i+1) for i,leg in enumerate(legs)};token.update({pair["s_a"]:"p1a",pair["s_b"]:"p1b"})
  words={edge:" ".join(token[n] for n in names) for edge,names in prior.bridge.chains(top).items() if names};assert words==expected[oid]
  sub={sp.Symbol(leg["s_var"]):var for leg,var in zip(legs,(t,x,y))}
  eq(parse(j["tgv"]["S_polynomial"]).subs(sub),als[0]*(als[1]*F(t,x)+als[2]*F(t,y)))
  eq(parse(j["tgv"]["pnum"]).subs(sub),-(1+z)*((1+z)*t-z)*als[0]*(als[1]-als[2])/4)
  ps={sp.Symbol(pair["s_a"]):p,sp.Symbol(pair["s_b"]):q}
  pure=F(p,q) if oid!="44" else (1-p)*(1-q)+z*p*q
  eq(parse(j["f_polynomials"]["f1"]).subs(ps),pure);eq(parse(j["pure"]["lorentz_z"]),-1-z)
  assert d["C"]["fixed_basis_coordinates"]==["-1/4","1/2","0"] and d["C"]["euler_compensation_count"]==1 and d["c4_symmetry"]["orbit_multiplicity"]==4
  CA,CF,AS,v14,I0=sp.symbols("CA CF alphaS v14 I0")
  common=CA**2*CF/16*(AS/sp.pi)**3*sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(2-2*e)*sp.pi**(3*e)*(I0-2*v14)**(3*e)
  assert sp.simplify(parse(d["C"]["expression"])/common)==1
 pairs=[prior.reflection(src[a],src[b],roots[a],roots[b]) for a,b in (("25","22"),("24","23"))]
 negative=[]
 for label,kwargs in [("missing_alpha_swap",{"swap":False}),("missing_order_reversal",{"reverse":False}),("flipped_P",{"sign":-1})]:
  try:prior.reflection(src["25"],src["22"],roots["25"],roots["22"],**kwargs)
  except AssertionError:negative.append(label)
  else:raise AssertionError(label)
 bad=copy.deepcopy(roots["25"]);bad["C"]["expression"]="2*("+bad["C"]["expression"]+")"
 try:prior.reflection(src["25"],src["22"],bad,roots["22"])
 except AssertionError:negative.append("doubled_C")
 else:raise AssertionError("doubled_C")
 # Two-loop T3 input in its D frame, mapped into the current F frame.
 jold=next(d["J"] for d in read(ROOT/"benchmarks/tgv/3e-3g/prepared/3e-3g.json")["orbits"] if d["orbit_id"]=="3e-3g-03")
 oldsub={sp.Symbol("s21"):t,sp.Symbol("s41"):x,sp.Symbol("s42"):y,z:1/z}
 eq(z*parse(jold["S"]).subs(oldsub,simultaneous=True),als[0]*(als[1]*F(t,x)+als[2]*F(t,y)))
 eq(z*z*parse(jold["numerator"]).subs(oldsub,simultaneous=True),-(1+z)*((1+z)*t-z)*als[0]*(als[1]-als[2])/4)
 top44=src["44"]["representative"]["canonical_topology"];purevars={top44["pure_pairings"][0][k] for k in ("s_a","s_b")}
 tgvvars={v["s_var"] for v in top44["tgv_legs"]}
 assert purevars.isdisjoint(tgvvars)
 for v,bounds in top44["bounds"].items():
  symbols=set().union(*(parse(q).free_symbols for q in bounds))
  assert symbols<={sp.Symbol(w) for w in (purevars if v in purevars else tgvvars)}
 # Alpha scale, density differentiation and endpoint-order transformations.
 alpha={als[0]:a,als[1]:(1-a)*b,als[2]:(1-a)*(1-b)}
 eq(sp.det(sp.Matrix([alpha[als[0]],alpha[als[1]]]).jacobian([a,b])),1-a)
 lam=1-2*e;kernel=F(t,b*x+(1-b)*y)
 assert sp.simplify((sp.diff(kernel**(-lam),x)-sp.diff(kernel**(-lam),y))/kernel**(-lam-1)+lam*(2*b-1)*((1+z)*t-z))==0
 def br(i,j):return sp.rf(1-e,i)*sp.rf(1-e,j)/sp.rf(2-2*e,i+j)
 def moment(n):return (br(n,0)-1/sp.Integer(n+1))/e
 density=[]
 for n in range(1,5):
  low=0;mid=0
  for j in range(n):
   beta=sp.binomial(n-1,j)*(2*br(j+1,n-1-j)-br(j,n-1-j))
   low+=-n*beta*q**(n+1)/((j+1)*(n+1))
   mid+=-n*beta*q**(j+1)*(1-q**(n-j))/((j+1)*(n-j))
  upper=sum(sp.binomial(n,k)*q**(n-k)*(1-q)**(k+1)*moment(k) for k in range(n+1))
  assert sp.factor(low-q**(n+1)*moment(n))==0
  assert sp.factor(mid-(moment(n)-q**(n+1)*moment(n)-upper))==0
  density += [{"n":n,"low_triangle":True,"middle_rectangle":True}]
 # Pure partial endpoint integral, including its removable internal pole.
 primitive=F(p,q)**e/(e*((1+z)*q-z))
 assert sp.simplify(sp.diff(primitive,p)/F(p,q)**(e-1))==1
 # The t<p and p<t triangles exactly cover the middle-domain t,p square.
 probes=0
 for n in range(3):
  test=(t+2*p)**n
  assert sp.integrate(test,(t,0,p),(p,0,1))+sp.integrate(test,(p,0,t),(t,0,1))==sp.integrate(test,(t,0,1),(p,0,1));probes+=1
 cs=coefficients();assert cs["22"]==cs["25"] and cs["23"]==cs["24"]
 assert all(all(weight(key)==6+k for key in v) for k,v in cs["23"].items())
 # Factorized I44 expansion checked against generic multiplication of the known H0 and T3 series.
 M0,M1,M2,Z2=sp.symbols("M_0_1 M_1_1 M_2_1 zeta_2")
 H0=M0+e*M1/2
 CminusB=-e*(2*M0+M1)-e**2*((4-Z2)*M0+2*M1+M2/2)
 product=sp.expand(H0*CminusB/32)
 atoms={sp.Symbol(g["name"]):__import__("formulas").vec(g["hlog_terms"]) for g in read(__import__("formulas").basispath())["generators"]}
 assert ha.from_sympy(product.coeff(e,1),atoms)==cs["44"][-1]
 assert ha.from_sympy(product.coeff(e,2),atoms)==cs["44"][0]
 return {"passed":True,"original_roots":list(IDS),"normalizations":5,"reflections":pairs,"negative_cases_rejected":negative,
 "factorized_44":{"disjoint_domains":True,"two_loop_S_P_scale_bridge":True,"unique_I_expansion":True},
 "alpha_Jacobian":"1-a","density_moments":density,"ordered_t_p_partition_probes":probes,"pure_endpoint_primitive":True,
 "coefficient_equalities":14,"23_uniform_weight":True,"pole_orders":{oid:-min(k for k,v in c.items() if v) for oid,c in cs.items()},"formal":"not_run"}
if __name__=="__main__":
 r=check();dump(run()/"math-checks.json",r);print(__import__("json").dumps(r,indent=2))
