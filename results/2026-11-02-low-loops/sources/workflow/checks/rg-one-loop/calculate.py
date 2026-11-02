"""One-loop physical assembly, exact checks and immutable computational result."""
from pathlib import Path
import hashlib
import itertools
import json
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEST = ROOT / "benchmarks/renormalization/one-loop"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "workflow/checks/single-opposite"))
from preprocessing.prepare import verify_prepared
import hlog_algebra as ha

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def data(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))

def encode(obj):
    return (json.dumps(obj, ensure_ascii=False, indent=2) + "\n").encode("utf-8")

def immutable(path, obj):
    blob = encode(obj)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_bytes() == blob, f"Immutable mismatch: {path}"
    else:
        path.write_bytes(blob)

def main():
    rows = []
    dependencies = set()
    def bind(path):
        path = Path(path)
        if not path.is_absolute():
            path = ROOT / path
        assert path.exists(), path
        dependencies.add(path)
        return data(path)
    def check(name, ok, detail=""):
        rows.append({"name": name, "passed": bool(ok), "detail": detail})
        assert ok, name
    def zero(name, expr):
        check(name, s.simplify(s.expand(expr)) == 0)

    source = bind(HERE / "sources/manifest.json")
    for record in source["files"] + source["modules"]:
        p = ROOT / record["snapshot"]
        dependencies.add(p)
        check("source_hash:" + record["snapshot"], sha(p) == record["sha256"])
    bind(HERE / "preprocessing.json")
    verify_prepared()
    package_path = ROOT / "benchmarks/2ne/prepared/2e.json"
    package = bind(package_path)
    manifest_path = ROOT / "benchmarks/2ne/prepared/manifest.json"
    bind(manifest_path)
    ordered_path = ROOT / "benchmarks/2ne/sources/montopar-ordered/2e/montopar-initial.json"
    ordered = bind(ordered_path)
    original = bind("benchmarks/2ne/sources/montopar/2e/montopar-initial.json")
    definitions = bind(HERE / "generated/1loop/orbit_definitions.json")
    metadata = bind(HERE / "generated/1loop/integrals.json")
    replay = bind(HERE / "generated/montopar-replay.json")
    checksums = []
    for root, src, raw, gen, meta, played in zip(package["orbits"], ordered["orbits"],
            original["orbits"], definitions["orbits"], metadata["orbits"], replay):
        oid = root["orbit_id"]
        check(oid + ":actual_producer_connection",
            root["topology"] == src["topology"] == raw["topology"] ==
            {k: gen[k] for k in root["topology"]})
        check(oid + ":ordered_J", root["J"] == src["J"])
        check(oid + ":source_C", root["color_preparation"]["original_C"] == src["C"] == raw["C"]
              and played["C_combined"] == src["C"]["combined"])
        check(oid + ":coupling_and_color", meta["g_power"] == 2 and root["color_preparation"]["coordinates"] == ["1"])
        check(oid + ":euler_once", root["C"]["combined"].count("Exp[-EulerGamma*eps]") == 1)
        result = bind(f"benchmarks/2ne/results/{oid}/r002.json")
        orbit_hash = hashlib.sha256(json.dumps(root, sort_keys=True, ensure_ascii=False,
                                               separators=(",", ":")).encode()).hexdigest()
        check(oid + ":result_binding", result["input"]["manifest_sha256"] == sha(manifest_path)
              and result["input"]["package_sha256"] == sha(package_path)
              and result["input"]["orbit_sha256"] == orbit_hash)
        ev = result["evidence"]
        for key in ("checks", "basis_snapshot", "derivation_snapshot"):
            if key in ev:
                p = ROOT / ev[key]
                dependencies.add(p)
                check(oid + ":" + key + "_hash", sha(p) == ev[key + "_sha256"])
        for file, digest in ev["code_sha256"].items():
            p = ROOT / file
            dependencies.add(p)
            check(oid + ":reused_code:" + file, sha(p) == digest)
        checksums.append(result)
    check("stored_adjacent_exact", checksums[0]["I"]["expression"] ==
          "-Exp[-EulerGamma*eps]*Gamma[1 - eps]/eps^2")
    check("stored_opposite_exact", checksums[1]["I_compact"]["expression"] ==
          "-1/2*G*H0/Cosh[(1/2)*ell*eps]")

    e = s.Symbol("eps")
    z, tau, mu2, gx2 = s.symbols("z tau mu2 mu_x2", positive=True)
    u, v = s.symbols("u v", real=True)
    G = s.exp(-s.EulerGamma*e) * s.gamma(1-e)
    a, b = -tau*z/2, -tau/2
    c = -a-b
    gram = s.Matrix([[0,a,c,b],[a,0,b,c],[c,b,0,a],[b,c,a,0]])
    check("lightlike_and_closure_Gram", all(gram[i,i] == 0 for i in range(4)) and
          gram*s.ones(4,1) == s.zeros(4,1))
    positions = [s.zeros(4,1)]
    for i in range(3):
        positions.append(positions[-1] + s.eye(4)[:, i])
    pairs = set(itertools.combinations(range(1,5),2))
    def orbit(i,j):
        return {tuple(sorted((((i-1+k)%4+1), ((j-1+k)%4+1)))) for k in range(4)}
    oa, ob = orbit(3,4), orbit(2,4)
    check("complete_six_edge_pairs", len(oa)==4 and len(ob)==2 and not oa & ob and oa | ob == pairs)
    for ix, (i,j) in enumerate([(3,4),(2,4)]):
        dx = positions[i-1]+u*s.eye(4)[:,i-1]-positions[j-1]-v*s.eye(4)[:,j-1]
        distance = s.expand((dx.T*gram*dx)[0]/(2*b))
        expected = z*(1-u)*v if ix == 0 else (1-u)*(1-v)+z*u*v
        zero(f"geometry_{i}{j}", distance-expected)
        zero(f"Lorentz_{i}{j}", gram[i-1,j-1]/b-(z if ix == 0 else -(1+z)))
        vs = s.symbols("v1 v2 v3 v4")
        phase = s.sympify(definitions["orbits"][ix]["phases"][0]["expression"],
                         locals=dict(zip([str(x) for x in vs],vs)) |
                         {f"s{i}1":u, f"s{j}1":v})
        signed = phase + sum(dx[k]*vs[k] for k in range(4))
        zero(f"source_phase_is_minus_X_{i}{j}", signed.subs(vs[3],-vs[0]-vs[1]-vs[2]))
        check(f"original_order_{i}{j}", package["orbits"][ix]["J"]["integration_order_outer_to_inner"] ==
              [f"s{i}1", f"s{j}1"])
    for i in range(4):
        zero(f"same_edge_numerator_{i+1}", gram[i,i])
    check("orbit_counts_not_n_diagrams", metadata["orbits"][0]["n_diagrams"] == 1
          and package["orbits"][0]["c4_symmetry"]["orbit_multiplicity"] == 4)

    # The inner opposite integration is a divided difference with a removable zero.
    A, B, q = s.symbols("A B q", positive=True)
    primitive = (A+B*v)**e/(e*B)
    zero("opposite_inner_primitive", s.diff(primitive,v)-(A+B*v)**(e-1))
    for ep in (1,2,3):
        ja = s.integrate(s.integrate(z*(z*(1-u)*v)**(ep-1),(v,0,1)),(u,0,1))
        jb = s.integrate(s.integrate(-(1+z)*((1-u)*(1-v)+z*u*v)**(ep-1),(v,0,1)),(u,0,1))
        zero(f"adjacent_original_order_eps_{ep}", ja-z**ep/ep**2)
        zero(f"opposite_crossing_eps_{ep}", jb-z**ep*jb.subs(z,1/z))
    H = s.Symbol("H0")
    k = s.Symbol("z_to_eps_half", positive=True)
    # Independent orbit assembly; k=z^(eps/2).
    jA, jArot = k**2/e**2, 1/(k**2*e**2)
    jB, jBrot = k*H, H/k
    pref = -G/2
    IA = -G/e**2
    IB = -G*H/(k+1/k)
    zero("adjacent_CJ_to_I", 2*pref*(jA+k**2*jArot)-(1+k**2)*IA)
    zero("opposite_CJ_to_I", pref*(jB+k**2*jBrot)-(1+k**2)*IB)
    check("reject_missing_rotation_scale", s.simplify(2*pref*(jA+jArot)-(1+k**2)*IA) != 0)
    check("reject_duplicate_orbit_multiplicity", s.simplify(4*(1+k**2)*IA-(1+k**2)*IA) != 0)

    # Raw alpha0 -> prepared coefficient: undo its single Euler normalization,
    # then make the physical MS-bar bare-coupling substitution.
    source_factor = s.pi**e*s.gamma(1-e) * mu2**e*s.exp(s.EulerGamma*e)/(4*s.pi)**e
    expected_factor = G*(mu2*s.exp(2*s.EulerGamma)/4)**e
    zero("MSbar_coordinate_scale_bridge",
         s.expand_power_base(source_factor-expected_factor,force=False))
    zero("prepared_scale_bridge",
         s.expand_power_base((s.pi*(mu2*s.exp(2*s.EulerGamma)/(4*s.pi)))**e -
                             (mu2*s.exp(2*s.EulerGamma)/4)**e,force=False))

    Ls, Lt, L, rho, CF, lam = s.symbols("L_s L_t L rho CF lambda", real=True)
    M01,M11,M21,M03 = s.symbols("M01 M11 M21 M03")
    moments = ROOT/"workflow/checks/single-opposite/runs/20260916T090538777022Z"
    mv = {}
    for tag in ("M_0_1","M_1_1","M_2_1","M_0_3","H0_0_plus","H0_1_plus","H0_2_plus"):
        p=moments/(tag+".txt")
        dependencies.add(p)
        mv[tag]=ha.parse(p)
    check("M01_elementary_hlog", mv["M_0_1"] == ha.add(ha.scale(ha.power(ha.LOG,2),s.Rational(-1,2)),
                                                        ha.scale(ha.period(2),-3)))
    check("H0_constant_M_binding", mv["H0_0_plus"] == mv["M_0_1"])
    check("H0_linear_M_binding", mv["H0_1_plus"] == ha.scale(mv["M_1_1"],s.Rational(1,2)))
    check("H0_quadratic_M_binding", mv["H0_2_plus"] == ha.add(
        ha.scale(mv["M_2_1"],s.Rational(1,8)),ha.scale(mv["M_0_3"],s.Rational(1,24))))
    # Also reproduce these weights directly from the moment generating function.
    h, r = s.symbols("h r")
    generator = s.series(s.exp(e*r/2)*2*s.sinh(e*h/2)/e,e,0,3).removeO()
    zero("moment_generator_coefficients", generator-h-e*r*h/2-e**2*(r**2*h/8+h**3/24))
    greg=s.series(G,e,0,5).removeO().expand()
    zero("Gamma_series",greg-(1+s.pi**2*e**2/12+s.zeta(3)*e**3/3+s.pi**4*e**4/160))
    hreg=M01+e*M11/2+e**2*(M21/8+M03/24)
    adj = -greg*(s.series(s.exp(e*Ls),e,0,5).removeO()+
                 s.series(s.exp(e*Lt),e,0,5).removeO())/e**2
    opp = -greg*s.series(s.exp(e*(Ls+Lt)/2),e,0,3).removeO()*hreg
    all_series=s.expand(adj+opp)
    coeff={j:s.expand(all_series.coeff(e,j)) for j in range(-2,3)}
    finite=s.simplify(coeff[0].subs(M01,-((Ls-Lt)**2+s.pi**2)/2))
    zero("finite_total",finite+Ls*Lt-s.pi**2/3)
    pole=CF*lam*(coeff[-2]/e**2+coeff[-1]/e)
    ren=CF*lam*finite
    def dmu(expr, dim=False):
        return s.expand(2*s.diff(expr,Ls)+2*s.diff(expr,Lt)-
                        (2*e*lam*s.diff(expr,lam) if dim else 0))
    rg=dmu(ren)
    from_finite=s.simplify(-s.diff(rg.subs({Ls:L+rho/2,Lt:L-rho/2}),L)/4)
    z_gamma=dmu(pole,True)
    from_Z=s.simplify(s.diff(z_gamma,Ls)/2)
    zero("finite_RG_equation",rg+2*CF*lam*(Ls+Lt))
    zero("Z_RG_equation",z_gamma-2*CF*lam*(Ls+Lt))
    zero("cusp_finite_vs_Z",from_finite-from_Z)
    zero("cusp_reference",from_Z-CF*lam)
    zero("noncusp_one_loop",rg+2*from_Z*(Ls+Lt))
    exact_unit=-G*s.exp(e*L)*(2*s.cosh(e*rho/2)/e**2+H)
    zero("exact_bare_D_dimensional_scale_invariance",2*s.diff(exact_unit,L)-2*e*exact_unit)
    dhks_div=-CF*lam*(s.exp(e*Ls)+s.exp(e*Lt))/e**2
    dhks_shape=CF*lam*((Ls-Lt)**2/2+s.pi**2/3)
    dhks_finite=s.series(dhks_div,e,0,1).removeO()-pole+dhks_shape
    zero("DHKS_finite_subtraction_conversion",dhks_finite-ren)
    zero("DHKS_regulated_evolution",s.diff(dhks_div,Ls)+s.diff(dhks_div,Lt)-
         (-CF*lam*(s.exp(e*Ls)+s.exp(e*Lt))/e))
    check("reject_DHKS_shape_as_MS_W",s.simplify(dmu(dhks_shape)+2*CF*lam*(Ls+Lt)) != 0)
    check("reject_opposite_sign",s.simplify((s.expand(adj-opp).coeff(e,0)-
          finite).subs(M01,-((Ls-Lt)**2+s.pi**2)/2)) != 0)
    check("reject_double_Euler", s.expand(s.series(-G*s.exp(-s.EulerGamma*e)/e**2,e,0,0).removeO()).coeff(e,-1) != 0)
    check("reject_Z_sign",s.simplify(-z_gamma-2*CF*lam*(Ls+Lt)) != 0)
    zero("crossing_finite",finite-finite.xreplace({Ls:Lt,Lt:Ls}))
    for j,expr in coeff.items():
        zero(f"crossing_coefficient_{j}",expr-expr.xreplace({Ls:Lt,Lt:Ls}))

    def obj(x):
        return {"sympy": str(x), "latex": s.latex(x)}
    output = {
        "format": "wilson-one-loop-rg-result/1", "revision": "r001",
        "scope": "fundamental QCD closed non-self-intersecting null quadrilateral; one loop",
        "domain": "s,t<0 as analytic invariants, z=s/t>0; fixed finite positive z",
        "formal": "not_run", "MCP": "not_available_in_session",
        "QFT_basis": "dimensional regularization and Feynman boundary-value propagator; not a Lean proof",
        "definitions": {"lambda": "alpha_s_MSbar/pi", "s": "2*v12", "t": "2*v14",
            "mu_x_squared": "exp(2*EulerGamma)*mu_MSbar_squared/4",
            "mu_prep_squared": "exp(2*EulerGamma)*mu_MSbar_squared/(4*pi)",
            "L_s": "log(mu_x_squared*(-s))", "L_t": "log(mu_x_squared*(-t))",
            "L": "(L_s+L_t)/2", "rho": "L_s-L_t=log(z)",
            "G": "exp(-EulerGamma*eps)*Gamma(1-eps)",
            "H0": "existing 2e-02/r002.json definitions.H0"},
        "exact_log_W_one_loop": {
            "prefactor": "CF*lambda",
            "expression": "-G*exp(eps*L)*(2*cosh(eps*rho/2)/eps^2+H0(eps,z))",
            "error_in_coupling": "O(lambda^2)"},
        "laurent_unit_coefficients": {str(j):obj(x) for j,x in coeff.items()},
        "epsilon_remainder": "O(eps^3) after poles; fixed geometry",
        "log_Z4_one_loop":obj(pole),
        "log_WR_one_loop":obj(ren),
        "RG_one_loop":obj(rg),
        "cusp_one_loop":obj(from_Z),
        "Gamma0_alpha_over_4pi":obj(4*CF),
        "gamma4_one_loop":obj(s.S.Zero),
        "gamma4_convention": "logs use mu_x_squared; a different constant log scale redistributes noncusp",
        "DHKS_shape":obj(dhks_shape),
        "ordinary_QCD_counterterms": "first enter at order lambda^2; none at this order",
        "dependencies": {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(dependencies)},
        "limitations": ["general complex z / physical i0 continuation not certified",
                        "no two-loop assembly or new formal receipts"]
    }
    immutable(DEST/"result-r001.json",output)
    report={"checks":rows,"passed":sum(x["passed"] for x in rows),"total":len(rows),
            "result_sha256":sha(DEST/"result-r001.json"),
            "code_sha256":sha(__file__),"formal":"not_run"}
    (HERE/"exact-checks.json").write_bytes(encode(report))
    print(json.dumps({"exact_checks":len(rows),"passed":True,"finite":str(ren),"cusp":str(from_Z)}))

if __name__=="__main__":
    main()
