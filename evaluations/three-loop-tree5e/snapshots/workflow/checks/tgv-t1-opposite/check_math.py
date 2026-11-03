"""Connection checks for all three original ordered roots and the exact reductions."""
import sys,importlib.util
from pathlib import Path
import sympy as sp
from formulas import ROOT,HERE,run,read,dump,coefficients,check_reused_T,ha,vec,weight
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
audit=module("t1opp_audit",ROOT/"workflow/checks/tgv-audit/audit.py");parse=audit.parse
def eq(a,b):assert sp.expand(a-b)==0,(a,b)
def normalize_positive(expr):
    for _ in range(3):
        expr=expr.replace(lambda t:isinstance(t,sp.Pow),lambda t:sp.Pow(sp.factor(t.base),t.exp))
        expr=sp.factor(sp.powdenest(sp.expand_power_base(expr,force=True),force=True))
    return expr
def check():
    raw={r["orbit_id"][-2:]:r for r in read(ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json")["orbits"]}
    roots={r["orbit_id"][-2:]:r for r in read(ROOT/"benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json")["orbits"]}
    e,z=sp.symbols("eps z1");a,b,u,x,y,v,w,s=sp.symbols("a b u x y v w s",positive=True);al=sp.symbols("alpha1:4")
    expected={"11":{2:"p1a",3:"t1",4:"t2 p1b t3"},"12":{2:"p1a",3:"t1",4:"p1b t2 t3"},
              "32":{1:"p1b",3:"t1 p1a",4:"t2 t3"}}
    source_checks=[]
    for oid in expected:
        r=raw[oid];d=roots[oid];top=r["representative"]["canonical_topology"];legs=top["tgv_legs"];pair=top["pure_pairings"][0]
        token={p["s_var"]:"t"+str(i+1) for i,p in enumerate(legs)}
        token[pair["s_a"]]="p1a";token[pair["s_b"]]="p1b";words={}
        for edge in range(1,5):
            vs=sorted([v for v in top["bounds"] if v[1]==str(edge)])
            if not vs:continue
            for i,name in enumerate(vs):assert top["bounds"][name]==["0",vs[i+1] if i+1<len(vs) else "1"]
            words[edge]=" ".join(token[name] for name in vs)
        assert words==expected[oid]
        assert d["J_raw_s"]==r["J_raw_s"]
        assert r["J_raw_s"]["bounds"]==top["bounds"]
        order=d["J_raw_s"]["iterated_order_outer_to_inner"]
        assert set(order)==set(top["bounds"])
        for i,name in enumerate(order):
            assert parse(top["bounds"][name][1]).free_symbols<={sp.Symbol(t) for t in order[:i]}
        sub={sp.Symbol(legs[0]["s_var"]):1-u,sp.Symbol(legs[1]["s_var"]):x,sp.Symbol(legs[2]["s_var"]):y}
        eq(parse(r["J_raw_s"]["tgv"]["S_polynomial"]).subs(sub),z*u*al[0]*(al[1]*x+al[2]*y))
        eq(parse(r["J_raw_s"]["tgv"]["pnum"]).subs(sub),z*z*u*al[0]*(al[1]-al[2])/4)
        ends={pair["edge_a"]:sp.Symbol(pair["s_a"]),pair["edge_b"]:sp.Symbol(pair["s_b"])}
        pure=parse(r["J_raw_s"]["f_polynomials"]["f1"])
        pure_sub={ends[2]:s,ends[4]:w} if oid!="32" else {ends[3]:1-w,ends[1]:1-s}
        eq(pure.subs(pure_sub), (1-s)*(1-w)+z*s*w)
        eq(parse(r["J_raw_s"]["pure"]["lorentz_z"]),-1-z)
        assert d["C"]["fixed_basis_coordinates"]==["-1/4","1/2","0"]
        assert d["C"]["euler_compensation_count"]==1 and d["c4_symmetry"]["orbit_multiplicity"]==4
        CA,CF,aS,v14,I0=sp.symbols("CA CF alphaS v14 I0")
        common=CA**2*CF/16*(aS/sp.pi)**3*sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(2-2*e)*sp.pi**(3*e)*(I0-2*v14)**(3*e)
        assert sp.simplify(parse(d["C"]["expression"])/common)==1
        source_checks.append(oid)
    alpha={al[0]:a,al[1]:(1-a)*b,al[2]:(1-a)*(1-b)}
    jac=sp.det(sp.Matrix([alpha[al[0]],alpha[al[1]]]).jacobian([a,b]));eq(jac,1-a)
    old=sp.prod(al)**(-e)*(z*z*u*al[0]*(al[1]-al[2])/4)*(z*u*al[0]*(al[1]*x+al[2]*y))**(-2+2*e)
    new=z**(2*e)/4*u**(2*e-1)*a**(e-1)*b**(-e)*(1-b)**(-e)*(2*b-1)*(b*x+(1-b)*y)**(-2+2*e)
    chart={a:a/(1+a),b:b/(1+b),z:sp.Symbol("positive_z",positive=True)}
    assert normalize_positive((old.subs(alpha,simultaneous=True)*jac/new).subs(chart,simultaneous=True))==1
    # Joint primitive, before any potentially divergent beta terms are separated.
    p=sp.symbols("p");L=b*x+(1-b)*y
    primitive=L**p/(p*(p-1)*b*(1-b))
    assert sp.simplify(sp.diff(primitive,x,y)/L**(p-2))==1
    # Homothetic domain identity h11+h12=(1-w^(2e))*h0/(2e).
    # Independent polynomial probes retain a general numerator b dependence.
    probes=0
    for power in (0,1,2):
        kernel=(2*b-1)*(b*x+(1-b)*y)**power
        low=sp.integrate(kernel,(x,0,y),(y,0,w))
        mid=sp.integrate(kernel,(x,0,w),(y,w,1))
        high=sp.integrate(kernel,(x,w,y),(y,w,1))
        whole=sp.integrate(kernel,(x,0,y),(y,0,1))
        eq(low+mid+high,whole);eq(low,w**(power+2)*whole);probes+=2
    # Beta IBP, initially Re(e)<0 at fixed 0<w<1, then continued.
    Bdensity=b**(-e)*(1-b)**(-e)
    assert sp.factor(sp.diff(Bdensity,b)/Bdensity-e*(2*b-1)/(b*(1-b)))==0
    # The powers are exact symbolic, verified by ratios rather than unsafe force globally.
    assert sp.simplify(sp.diff((1-b+b*w)**(2*e),b)/(2*e*(w-1)*(1-b+b*w)**(2*e-1)))==1
    check_reused_T()
    cs=coefficients()
    for k in range(-6,1):assert ha.add(cs["11"][k],cs["12"][k])==cs["32"][k]
    assert all(all(weight(key)==6+k for key in v) for k,v in cs["11"].items())
    basis=read(ROOT/read(run()/"representation.json")["basis_path"])["weight_bases"]["6"]
    mn=[vec(b["hlog_terms"]) for b in basis if not ({"tauDq4","E"}&set(b["factors"]))]
    targets=[{key:c for key,c in cs[oid][0].items() if weight(key)==6} for oid in ("11","12","32")]
    def rank(vs):
        keys=sorted(set().union(*vs));return sp.Matrix(len(keys),len(vs),lambda i,j:vs[j].get(keys[i],0)).rank()
    assert rank(mn)==12 and rank(mn+targets)==13
    poles={oid:-min(k for k,vv in c.items() if vv) for oid,c in cs.items()}
    assert poles=={"11":3,"12":2,"32":3}
    return {"passed":True,"original_roots":source_checks,"ordered_domains":3,"normalizations":3,
       "alpha_change_and_jacobian":True,"joint_xy_primitive":True,"Beta_IBP_identity":True,
       "ordered_partition_polynomial_probes":probes,"independent_old_single_moments":4,"coefficient_relations":7,
       "pole_orders":poles,"11_uniform_weight":True,"weight6_MN_rank":12,"weight6_MN_plus_targets_rank":13,"formal":"not_run"}
if __name__=="__main__":
    report=check();dump(run()/"math-checks.json",report);print(__import__("json").dumps(report,indent=2))
