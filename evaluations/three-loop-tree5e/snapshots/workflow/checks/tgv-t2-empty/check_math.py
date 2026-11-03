"""Exact reflection bridge from the seven original T2[] roots to T1[]."""
from pathlib import Path
import json,sys,hashlib,datetime,importlib.util,copy
import sympy as sp
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
REFLECT={1:2,2:1,3:4,4:3}
PAIR={"05":"04","06":"03","07":"02","08":"01","20":"19","21":"18","31":"17"}
SOURCE=ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json"
PREPARED=ROOT/"benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json"
def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
audit=load_module("t2audit",ROOT/"workflow/checks/tgv-audit/audit.py")
parse=audit.parse
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding="utf-8"))
def dump(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
def run():return ROOT/read(HERE/"current-run.json")["run"]
def records(p):return {r["orbit_id"][-2:]:r for r in read(p)["orbits"]}
def eq(a,b):assert sp.expand(a-b)==0,(a,b)
def chains(top):
    out={}
    for edge in range(1,5):
        vs=sorted((v for v in top["bounds"] if v[1]==str(edge)),key=lambda v:int(v[2:]))
        for i,v in enumerate(vs):
            assert top["bounds"][v]==["0",vs[i+1] if i+1<len(vs) else "1"]
        out[edge]=vs
    return out
def check_pair(t,s,rt,rs,swap_alpha=True,reverse_order=True):
    tt=t["representative"]["canonical_topology"];st=s["representative"]["canonical_topology"]
    tc,sc=chains(tt),chains(st);sub={}
    for edge,vs in tc.items():
        ws=list(reversed(sc[REFLECT[edge]])) if reverse_order else sc[REFLECT[edge]]
        assert len(vs)==len(ws)
        sub.update({sp.Symbol(v):1-sp.Symbol(w) for v,w in zip(vs,ws)})
    al=sp.symbols("alpha1:4")
    sub.update(dict(zip(al,reversed(al) if swap_alpha else al)))
    # Endpoint identity includes the TGV labels and the unordered pure pair.
    for i in range(3):
        v=tt["tgv_legs"][i];w=st["tgv_legs"][2-i]
        assert REFLECT[v["edge"]]==w["edge"]
        eq(sub[sp.Symbol(v["s_var"])],1-sp.Symbol(w["s_var"]))
    tp=tt["pure_pairings"][0];spair=st["pure_pairings"][0]
    assert {str(1-sub[sp.Symbol(tp[k])]) for k in ("s_a","s_b")}=={spair["s_a"],spair["s_b"]}
    # All inequalities, not only topology names, map bijectively.
    for edge,vs in tc.items():
        original=[sp.Integer(0)]+[sp.Symbol(v) for v in vs]+[sp.Integer(1)]
        mapped=[x.subs(sub,simultaneous=True) for x in original]
        ws=sc[REFLECT[edge]]
        expected=[sp.Integer(0)]+[sp.Symbol(v) for v in ws]+[sp.Integer(1)]
        eqs=[sp.expand(mapped[i+1]-mapped[i]) for i in range(len(mapped)-1)]
        reverse_gaps=[sp.expand(expected[i+1]-expected[i]) for i in reversed(range(len(expected)-1))]
        assert eqs==reverse_gaps
    for raw,root,top in ((t,rt,tt),(s,rs,st)):
        j=root["J_raw_s"]
        assert j["bounds"]==top["bounds"]==raw["J_raw_s"]["bounds"]
        assert j["tgv"]==raw["J_raw_s"]["tgv"] and j["pure"]==raw["J_raw_s"]["pure"]
        assert j["f_polynomials"]==raw["J_raw_s"]["f_polynomials"]
        order=j["iterated_order_outer_to_inner"]
        assert set(order)==set(top["bounds"])
        for i,v in enumerate(order):
            assert parse(j["bounds"][v][1]).free_symbols<={sp.Symbol(q) for q in order[:i]}
        assert j["alpha_measure"]==raw["J_raw_s"]["alpha_measure"]
        assert j["alpha_measure"]["constraint"]=="alpha1+alpha2+alpha3=1"
        eq(parse(j["tgv"]["denom_power"]),2-2*sp.Symbol("eps"))
        assert root["C"]["fixed_basis_coordinates"]==["-1/4","1/2","0"]
        assert root["C"]["euler_compensation_count"]==1
        assert root["c4_symmetry"]["orbit_multiplicity"]==4
    tj,sj=t["J_raw_s"],s["J_raw_s"]
    for key in ("S_polynomial","pnum"):
        eq(parse(tj["tgv"][key]).subs(sub,simultaneous=True),parse(sj["tgv"][key]))
    eq(parse(tj["pure"]["lorentz_z"]).subs(sub,simultaneous=True),parse(sj["pure"]["lorentz_z"]))
    eq(parse(tj["f_polynomials"]["f1"]).subs(sub,simultaneous=True),parse(sj["f_polynomials"]["f1"]))
    assert tj["pure"]["f_product"]==sj["pure"]["f_product"]=="f1^(-1+eps)"
    assert tj["tgv"]["alpha_monomial"]==sj["tgv"]["alpha_monomial"]=="(alpha1*alpha2*alpha3)^(-eps)"
    eq(sp.prod(al).subs(sub,simultaneous=True),sp.prod(al))
    variables=[sp.Symbol(v) for v in tt["bounds"]]
    destination=[sp.Symbol(v) for v in st["bounds"]]
    sjac=sp.det(sp.Matrix([sub[v] for v in variables]).jacobian(destination))
    ajac=sp.det(sp.Matrix([sub[a].subs(al[2],1-al[0]-al[1]) for a in al[:2]]).jacobian(al[:2]))
    assert abs(sjac)==abs(ajac)==1
    assert sp.simplify(parse(rt["C"]["expression"])/parse(rs["C"]["expression"]))==1
    assert rt["original_C"]==rs["original_C"]
    return {"T2":t["orbit_id"],"T1":s["orbit_id"],"substitution_T2_in_T1_coordinates":{str(k):str(v) for k,v in sub.items()},
      "ordered_domain_bijection":True,"s_absolute_jacobian":str(abs(sjac)),"alpha_absolute_jacobian":str(abs(ajac)),
      "S_P_pure_equal":True,"colour_and_C_equal":True,"raw_J_ratio":1,"physical_I_ratio":1}
def compute():
    src=records(SOURCE);roots=records(PREPARED)
    cat=read(ROOT/"benchmarks/tgv/families/catalog-r001.json")
    members=next(g["mixed_members"] for g in cat["groups"] if g["group_id"]=="t2-empty")
    assert set(members)=={"5e-3g-"+x for x in PAIR}
    prep=read(PREPARED);assert prep["loop_order"]==3
    assert prep["I_definition"]=="I=sym/32*Exp[-3*EulerGamma*eps]*Gamma[1-eps]*Gamma[2-2*eps]*(J(z,eps)+z^(3*eps)*J(1/z,eps))/(1+z^(3*eps))"
    manifest=read(PREPARED.parent/"manifest.json");assert manifest["sha256"]==sha(PREPARED)
    rows=[check_pair(src[t],src[s],roots[t],roots[s]) for t,s in PAIR.items()]
    failures=[]
    for label,alter,kwargs in [
      ("missing_alpha_permutation",None,{"swap_alpha":False}),
      ("missing_order_reversal",None,{"reverse_order":False}),
      ("flipped_P",lambda d:d["J_raw_s"]["tgv"].update(pnum="-("+d["J_raw_s"]["tgv"]["pnum"]+")"),{}),
      ("changed_bound",lambda d:d["representative"]["canonical_topology"]["bounds"].update(s31=["0","1"]),{})]:
        damaged=copy.deepcopy(src["05"])
        if alter:alter(damaged)
        try:check_pair(damaged,src["04"],roots["05"],roots["04"],**kwargs)
        except AssertionError:failures.append(label)
        else:raise AssertionError("negative accepted: "+label)
    # Reflection is an isometry of the dimensionless Gram matrix, keeping z and v14.
    z=sp.Symbol("z1");gram=sp.Matrix([[0,z,-1-z,1],[z,0,1,-1-z],[-1-z,1,0,z],[1,-1-z,z,0]])
    assert gram==gram.extract([1,0,3,2],[1,0,3,2])
    return {"passed":True,"pairs":rows,"negative_cases_rejected":failures,"Gram_reflection":True,
      "z_transformation":"z -> z","initial_domain":"z>0; 1/2<Re(eps)<1; absolute convergence permits the recorded order change",
      "formal":"not_run","method":"exact change of variables, then meromorphic continuation; no new integration"}
def main():
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    r=HERE/"runs"/stamp;r.mkdir(parents=True);dump(HERE/"current-run.json",{"run":r.relative_to(ROOT).as_posix()})
    result=compute();dump(r/"reflection.json",result)
    dump(r/"environment.json",{"python":sys.version,"sympy":sp.__version__,"source_sha256":sha(SOURCE),"prepared_sha256":sha(PREPARED)})
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
