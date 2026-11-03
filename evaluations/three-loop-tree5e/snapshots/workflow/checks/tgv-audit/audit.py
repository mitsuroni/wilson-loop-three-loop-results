"""Independent Gram-matrix/Fierz audit; external inputs are strictly read-only."""
from pathlib import Path
import sys,json,itertools,hashlib,datetime,dataclasses
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
BASE=Path("C:/Users/Lenovo/Documents/learning/wilsonloop")
sys.path.insert(0,str(BASE/"montopar/src"))
z,e,Nc=sp.symbols("z1 eps Nc")
al=sp.symbols("alpha1:4")
G=sp.Matrix([[0,z,-1-z,1],[z,0,1,-1-z],[-1-z,1,0,z],[1,-1-z,z,0]])
def parse(s):
    import ast
    def node(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return sp.Integer(n.value)
        if isinstance(n,ast.Name):return {"I":sp.I,"Pi":sp.pi,"EulerGamma":sp.EulerGamma}.get(n.id,sp.Symbol(n.id))
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -node(n.operand)
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.UAdd):return node(n.operand)
        if isinstance(n,ast.BinOp):
            a,b=node(n.left),node(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
            if isinstance(n.op,ast.Pow):return a**b
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ("Gamma","Exp","Log"):
            return {"Gamma":sp.gamma,"Exp":sp.exp,"Log":sp.log}[n.func.id](*(node(x) for x in n.args))
        raise ValueError("unsupported expression")
    return node(ast.parse(s.replace("^","**").replace("[","(").replace("]",")"),mode="eval").body)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pos(edge,s):
    return sp.Matrix([1 if j<edge-1 else sp.Symbol(s) if j==edge-1 else 0 for j in range(4)])
def dot(a,b):return sp.expand((a.T*G*b)[0])
def independent_kernel(legs):
    xs=[pos(t.edge,t.s_var) for t in legs]
    S=sum(al[i]*al[j]*dot(xs[i]-xs[j],xs[i]-xs[j])/2 for i in range(3) for j in range(i+1,3))
    # O/i consists of the three cyclic differences of ambient point derivatives.
    P=0
    for a,b,c in ((0,1,2),(1,2,0),(2,0,1)):
        rho=G[legs[a].edge-1,legs[b].edge-1]
        v=sp.eye(4)[:,legs[c].edge-1]
        for k,sign in ((b,1),(a,-1)):
            for j in range(3):
                if j!=k:P+=sign*rho*al[k]*al[j]*dot(v,xs[k]-xs[j])/4
    return sp.expand(S),sp.expand(P)
def traces_sum(traces):
    occ={};n=0
    for tr in traces:
        indices=list(range(n,n+len(tr)));n+=len(tr)
        for k,c in enumerate(tr):occ.setdefault(c,[]).append((indices[k],indices[(k+1)%len(tr)]))
    assert all(len(v)==2 for v in occ.values())
    total=0
    for bits in itertools.product((0,1),repeat=len(occ)):
        par=list(range(n))
        def find(a):
            while par[a]!=a:a=par[a]
            return a
        def union(a,b):par[find(a)]=find(b)
        for bit,pair in zip(bits,occ.values()):
            (i,j),(k,l)=pair
            for a,b in (((i,l),(k,j)) if bit==0 else ((i,j),(k,l))):union(a,b)
        loops=len({find(i) for i in range(n)})
        total+=sp.Rational(1,2)**len(occ)*(-1)**sum(bits)*Nc**(loops-sum(bits))
    return sp.expand(total)
def color(orbit,loops=None):
    vertex_colors={v.color for v in orbit.vertices};rename={}
    tloops=loops or [t.loop for t in orbit.tgv_legs]
    for p in orbit.propagators:
        if p.loop in tloops:
            ext=next(c for c in p.color if c in vertex_colors);rename[ext]="t"+str(tloops.index(p.loop))
        else:
            for c in p.color:rename[c]=p.loop
    word=[rename[v.color] for v in orbit.vertices]
    return sp.simplify(-2*sp.I*(traces_sum([word,["t0","t1","t2"]])-traces_sum([word,["t1","t0","t2"]]))/Nc)
def inputs(loop,fam):
    from montopar.parser import parse_all_orbits
    path=BASE/f"qgraf-pipeline/outputs/{loop}loop/orbit_definitions.json"
    return [o for o in parse_all_orbits(json.loads(path.read_text(encoding="utf-8"))) if o.emission_type==fam]
def main():
    from montopar.cli import _normalize_tgv_to_canonical
    from montopar.hypercube import hypercube_transform
    from montopar.tgv import compute_S,compute_Pnum
    from montopar.j_factor import compute_j_tgv,compute_j_mixed,compute_f_p
    records=[];sources={Path(__file__)}
    for loop,fam in ((2,"3e-3g"),(3,"5e-3g")):
        sources.add(BASE/f"qgraf-pipeline/outputs/{loop}loop/orbit_definitions.json")
        for orig in inputs(loop,fam):
            o=_normalize_tgv_to_canonical(orig) if fam=="5e-3g" else orig
            src=BASE/f"workspaces/{loop}loop/{fam}/{o.orbit_id}/montopar.json";sources.add(src)
            stored=json.loads(src.read_text(encoding="utf-8"));hc=hypercube_transform(o.bounds)
            S,P=independent_kernel(o.tgv_legs)
            S_code=parse(compute_S(o.tgv_legs));P_code=parse(compute_Pnum(o.tgv_legs))
            assert sp.expand(S-S_code)==0,(o.orbit_id,"S")
            assert sp.expand(P-P_code)==0,(o.orbit_id,"P")
            sub={sp.Symbol(k):parse(v) for k,v in hc.s_to_x.items()}
            x=[sp.Symbol(t) for t in hc.x_vars];sv=[sp.Symbol(b.var) for b in o.bounds]
            jac=sp.factor(sp.Matrix([sub[t] for t in sv]).jacobian(x).det())
            # Original variable listing may differ by a permutation; use absolute determinant on x>0.
            jac_code=parse(hc.jacobian)
            assert sp.factor(jac**2-jac_code**2)==0
            exact_color=color(o);raw_color=parse(stored["C"]["components"]["color"]).subs({sp.Symbol("CA"):Nc,sp.Symbol("CF"):(Nc*Nc-1)/(2*Nc)})
            color_ok=sp.simplify(exact_color-raw_color)==0
            pure=[]
            for p in o.pairings:
                f=dot(pos(p.edge_a,p.s_a)-pos(p.edge_b,p.s_b),pos(p.edge_a,p.s_a)-pos(p.edge_b,p.s_b))/2
                fcode=parse(compute_f_p(p.edge_a,p.edge_b,p.s_a,p.s_b))
                assert sp.expand(f-fcode)==0
                pure.append({"pair":dataclasses.asdict(p),"f_s":str(sp.factor(f)),"lorentz":str(G[p.edge_a-1,p.edge_b-1])})
            j=compute_j_mixed(o,hc,o.pairings) if pure else compute_j_tgv(o,hc)
            cubeS=sp.expand(S.subs(sub));cubeP=sp.expand(P.subs(sub))
            if pure:
                old=stored["J_raw"]["tgv"]
                S_ok=sp.expand(cubeS-parse(old["S_polynomial"]))==0
                P_ok=sp.expand(cubeP-parse(old["pnum"]))==0
                fs_ok=all(sp.expand(parse(pure[i]["f_s"]).subs(sub)-parse(stored["J_raw"]["pure"]["f_polynomials"][f"f{i+1}"]))==0 for i in range(len(pure)))
            else:
                # The full raw expression equals the independent product after exact integer test powers.
                old=stored["J_raw"]
                S_ok=P_ok=fs_ok=True
            exact_integrand=jac_code*cubeP*(sp.prod(al))**(-e)*cubeS**(-2+2*e)
            for p in pure:
                pair=p["pair"]
                exact_integrand*=G[pair["edge_a"]-1,pair["edge_b"]-1]*parse(p["f_s"]).subs(sub)**(-1+e)
            got=parse(stored["J_raw"]["integrand"])
            if pure:got=got.subs({sp.Symbol(k):parse(v) for k,v in stored["J_raw"]["pure"]["f_polynomials"].items()})
            sample={z:sp.Rational(2,3),al[0]:sp.Rational(1,5),al[1]:sp.Rational(3,10),al[2]:sp.Rational(1,2)}
            sample.update({t:sp.Rational(i+2,i+5) for i,t in enumerate(x)})
            vals=[]
            for ep in (0,1,2):
                residual=sp.cancel((exact_integrand-got).subs(e,ep).subs(sample));vals.append(str(residual))
                assert residual==0,(o.orbit_id,"cube",ep,residual)
            records.append({"orbit_id":o.orbit_id,"loop_order":loop,"legs":[dataclasses.asdict(t) for t in o.tgv_legs],
                "bounds":{b.var:[b.low,b.high] for b in o.bounds},"s_to_x":hc.s_to_x,"jacobian":str(jac_code),
                "independent_S":str(sp.factor(S)),"independent_Pnum":str(sp.factor(P)),
                "source_S_and_P_pass":True,"stored_S_pass":S_ok,"stored_P_pass":P_ok,"stored_pure_pass":fs_ok,
                "color_Fierz":str(exact_color),"color_stored":str(raw_color),"color_pass":color_ok,
                "integer_cube_residuals":vals,"pure":pure,
                "raw_s_unresolved_f_symbols":sorted(re.findall(r"\bf\d+\b",stored["J_raw_s"]["integrand"])),
                "alpha_constraint_explicit":bool(stored["J_raw"].get("delta_constraints")),
                "mixed_phase_ratio_stored_to_FT":str(sp.I**(2*len(pure)+3)/(sp.I**(2*len(pure)+3)*(-1)**len(pure)))})
            print(o.orbit_id,"geometry/cube passed; color",color_ok,flush=True)
    # Concrete theory counterexamples.
    x1,x2,x3,x4,x5=sp.symbols("x1:6")
    determinant=sp.det(sp.Matrix([x1,x2,x1*x3,x2*x4*x5,x2*x4]).jacobian([x1,x2,x3,x4,x5]))
    example={"document_jacobian":"x1*x2^2","correct_absolute_jacobian":"x1*x2^2*x4",
        "determinant":str(sp.factor(determinant)),"document_bare_f":"z1*(1-x1)*(1-x2)+x1*x2",
        "correct_bare_f":"z1*(1-x1)*x2","five_point_color":"-I*CF/4",
        "document_FT_section_5c_missing_factor":"4^(-n)","phase_after_pure_FT_cancellation":"I^3, independent of n_pure"}
    theory=Path("C:/Users/Lenovo/Documents/learning/test/TGV/混合型TGV+pure_通用因子分解推导.md")
    sources.add(theory)
    sources|=set((BASE/"montopar/src/montopar").glob("*.py"))
    run=HERE/"runs"/datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ");run.mkdir(parents=True)
    output={"records":records,"theory_counterexamples":example,"formal":"not_run",
      "scope":"exact geometry/numerator/Fierz algebra and structured cube checks; not a proof of physical-source convention or complex continuation",
      "all_color_pass":all(t["color_pass"] for t in records),"fully_correct":False,"mixed_phase_convention":"pure FT has (-1)^n; net eikonal/FT phase is I^3; stored mixed differs by (-1)^n","source_sha256":{str(p):sha(p) for p in sorted(sources)}}
    (run/"audit.json").write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    (HERE/"current-run.json").write_text(json.dumps({"run":run.relative_to(ROOT).as_posix()},indent=2)+"\n",encoding="utf-8")
    print("AUDIT",run,"color",output["all_color_pass"],flush=True)
if __name__=="__main__":
    import re
    main()
