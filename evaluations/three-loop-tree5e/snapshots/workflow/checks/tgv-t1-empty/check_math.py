"""Check exact roots, changes of variables, C/I normalization and analytic moment controls."""
import json,sys,importlib.util
import sympy as sp
from formulas import ROOT,HERE,IDS,run,coefficients,gamma_series,moment,e,ell
spec=importlib.util.spec_from_file_location("t1audit",ROOT/"workflow/checks/tgv-audit/audit.py")
am=importlib.util.module_from_spec(spec);spec.loader.exec_module(am)
parse=am.parse
def eq(a,b):assert sp.expand(a-b)==0,(a,b)
def positive_power_normalize(expr):
    # Used only after mapping all power bases into the stated positive real domain.
    for _ in range(3):
        expr=expr.replace(lambda t:isinstance(t,sp.Pow),
                          lambda t:sp.Pow(sp.factor(t.base),t.exp))
        expr=sp.powdenest(sp.expand_power_base(expr,force=True),force=True)
        expr=sp.factor(expr)
    return expr
def main():
    source=json.loads((ROOT/"benchmarks/tgv/5e-3g/sources/canonical-r001.json").read_text(encoding="utf-8"))["orbits"]
    roots=json.loads((ROOT/"benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json").read_text(encoding="utf-8"))["orbits"]
    roots={r["orbit_id"]:r for r in roots}
    cat=json.loads((ROOT/"benchmarks/tgv/families/catalog-r001.json").read_text(encoding="utf-8"))
    members=next(g["mixed_members"] for g in cat["groups"] if g["group_id"]=="t1-empty")
    assert members==["5e-3g-"+s for s in IDS]
    a,b,u,x,y,v,w,t,rho=sp.symbols("a b u x y v w t rho",positive=True)
    z=sp.Symbol("z1");al=am.al
    alpha={al[0]:a,al[1]:(1-a)*b,al[2]:(1-a)*(1-b)}
    jac=sp.det(sp.Matrix([alpha[al[0]],alpha[al[1]]]).jacobian([a,b]))
    eq(jac,1-a)
    density=sp.prod(al)**(-e)*(z*z*u*al[0]*(al[1]-al[2])/4)*(z*u*al[0]*(al[1]*x+al[2]*y))**(-2+2*e)
    after=density.subs(alpha,simultaneous=True)*jac
    expected=z**(2*e)/4*u**(-1+2*e)*a**(e-1)*b**(-e)*(1-b)**(-e)*(2*b-1)*(b*x+(1-b)*y)**(-2+2*e)
    # a,b in (0,1), z>0: factor on an explicitly positive rational chart.
    chart={a:a/(1+a),b:b/(1+b),z:sp.Symbol("positive_z",positive=True)}
    assert positive_power_normalize((after/expected).subs(chart,simultaneous=True))==1
    expected_words={
      "01":{3:"t1 p1a",4:"t2 t3 p1b"},"02":{3:"t1 p1a",4:"t2 p1b t3"},
      "03":{3:"p1a t1",4:"t2 p1b t3"},"04":{3:"p1a t1",4:"p1b t2 t3"},
      "17":{2:"p1a",3:"t1 p1b",4:"t2 t3"},
      "18":{1:"p1b",3:"t1",4:"t2 p1a t3"},
      "19":{1:"p1b",3:"t1",4:"p1a t2 t3"}}
    checks=[]
    for raw in source:
        oid=raw["orbit_id"]
        if oid not in members:continue
        top=raw["representative"]["canonical_topology"];legs=top["tgv_legs"]
        token={t["s_var"]:"t"+str(i) for i,t in enumerate(legs,1)}
        pair=top["pure_pairings"][0];token[pair["s_a"]]="p1a";token[pair["s_b"]]="p1b"
        actual_words={}
        for edge in range(1,5):
            names=sorted([v for v in top["bounds"] if v[1]==str(edge)],key=lambda v:int(v[2:]))
            if not names:continue
            for i,vname in enumerate(names):
                assert top["bounds"][vname]==["0",names[i+1] if i+1<len(names) else "1"]
            actual_words[edge]=" ".join(token[vname] for vname in names)
        assert actual_words==expected_words[oid[-2:]]
        assert raw["J_raw_s"]["bounds"]==top["bounds"]
        assert raw["J_raw_s"]["alpha_measure"]["constraint"]=="alpha1+alpha2+alpha3=1"
        sub={sp.Symbol(legs[0]["s_var"]):1-u,sp.Symbol(legs[1]["s_var"]):x,sp.Symbol(legs[2]["s_var"]):y}
        eq(parse(raw["J_raw_s"]["tgv"]["S_polynomial"]).subs(sub),z*u*al[0]*(al[1]*x+al[2]*y))
        eq(parse(raw["J_raw_s"]["tgv"]["pnum"]).subs(sub),z*z*u*al[0]*(al[1]-al[2])/4)
        prepared=roots[oid]
        assert prepared["c4_symmetry"]["orbit_multiplicity"]==4
        assert prepared["C"]["euler_compensation_count"]==1
        assert prepared["C"]["fixed_basis_coordinates"]==["-1/4","1/2","0"]
        CA,CF,aS,v14,I0=sp.symbols("CA CF alphaS v14 I0")
        common=CA**2*CF/16*(aS/sp.pi)**3*sp.exp(-3*sp.EulerGamma*e)*sp.gamma(1-e)*sp.gamma(2-2*e)*sp.pi**(3*e)*(I0-2*v14)**(3*e)
        assert sp.simplify(parse(prepared["C"]["expression"])/common)==1
        p=top["pure_pairings"][0]
        ends={p["edge_a"]:sp.Symbol(p["s_a"]),p["edge_b"]:sp.Symbol(p["s_b"])}
        actual=parse(raw["J_raw_s"]["f_polynomials"]["f1"])
        expected_f=(1-ends[3])*ends[4]*z if oid[-2:] in IDS[:4] else (1-ends[2])*ends[3] if oid.endswith("17") else (1-ends[4])*ends[1]
        eq(actual,expected_f)
        checks.append({"orbit_id":oid,"source_kernel":True,"pure_kernel":True,"C_I_connection":True})
    # Projective beta substitution separates the only H_m pole explicitly.
    betaK=b**(-e)*(1-b)**(-e)*(2*b-1)*(1-b+b*rho)**(-2+2*e)
    t_density=betaK.subs(b,1/(1+t))/(1+t)**2
    expected_t=t**(-e)*(1-t)/(1+t)*(t+rho)**(-2+2*e)
    assert positive_power_normalize(t_density/expected_t)==1
    # Independent polynomial probes of the hand-derived endpoint-scale integrals.
    probes=0
    for n in (1,2,3):
        before=sp.integrate(u**(2*n-1)*sp.integrate(v**(n-1),(v,0,u)),(u,0,1))
        after=sp.integrate(v**(n-1)*sp.integrate(u**(2*n-1),(u,0,v)),(v,0,1))
        eq(before,sp.Rational(1,3*n*n));eq(after,sp.Rational(1,6*n*n));probes+=2
        eq(sp.integrate(y**(2*n-1)*sp.integrate(w**(n-1),(w,y,1)),(y,0,1)),sp.Rational(1,6*n*n))
        eq(sp.integrate(y**(2*n-1)*sp.integrate(w**(n-1),(w,rho*y,y)),(y,0,1)),(1-rho**n)/(3*n*n))
        eq(sp.integrate(y**(2*n-1)*sp.integrate(w**(n-1),(w,0,rho*y)),(y,0,1)),rho**n/(3*n*n));probes+=3
    eq(sp.det(sp.Matrix([rho*y,y]).jacobian([rho,y])),y)
    # Compare target integrations in reversed order and a bare-core Gamma control.
    G,A,B,V,W,R=gamma_series()
    Fbare=sp.expand(((1-e)*B-A)*sum((2*e)**k for k in range(7))/e)
    for q in range(6):
        eq(moment("alternateF",q),moment("F",q))
        eq(moment("bareF",q),sp.simplify(Fbare.coeff(e,q)))
    cs=coefficients()
    for k in range(-6,1):
        eq(cs["02"][k],2*cs["03"][k])
        eq(cs["01"][k],cs["03"][k]+cs["04"][k])
        eq(cs["17"][k],cs["18"][k]+cs["19"][k])
        for oid in IDS:
            eq(cs[oid][k].subs(ell,-ell),cs[oid][k])
    poles={oid:min(k for k,v in cs[oid].items() if v!=0) for oid in IDS}
    assert poles=={"01":-6,"02":-6,"03":-6,"04":-6,"17":-4,"18":-4,"19":-3}
    report={"passed":True,"source_checks":checks,"alpha_change_and_jacobian":True,
       "projective_change":True,"ordered_domains":7,"endpoint_scale_polynomial_probes":probes,"alternate_F_coefficients":6,"bare_Gamma_controls":6,
       "exact_linear_coefficient_relations":21,"crossing_coefficients":49,"pole_orders":poles,"formal":"not_run"}
    (run()/"math-checks.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
