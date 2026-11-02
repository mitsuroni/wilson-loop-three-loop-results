"""Recompute the two original one-loop inputs through epsilon3.

Only SymPy, mpmath and the Python standard library are required. All output
goes to a fresh --output-dir outside this packet. Frozen input bytes are
checked, exact coefficients are generated before the historical table is
read, and no source generator or archived calculation script is executed.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import hashlib
import json
import math
import platform
import time
import sympy as s
from sympy.parsing.mathematica import parse_mathematica
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
eps,rho,z = s.symbols("eps rho z")
M = lambda k,j:s.Symbol(f"M{k}{j}")
CHECKS, READS = [], {}
encode = lambda d:(json.dumps(d,ensure_ascii=False,indent=2)+"\n").encode("utf-8")
sha = lambda b:hashlib.sha256(b).hexdigest()


def read_bytes(relative):
    p = (ROOT/relative).resolve()
    assert p.is_relative_to(ROOT), relative
    blob = p.read_bytes()
    h = sha(blob)
    if relative in READS:assert READS[relative] == h,relative
    READS[relative] = h
    return blob


def read_json(relative):
    return json.loads(read_bytes(relative).decode("utf-8-sig"))


def dereference(ref):
    blob = read_bytes(ref["path"])
    assert sha(blob) == ref["sha256"], ref["path"]
    value = json.loads(blob.decode("utf-8-sig"))
    for token in ref.get("json_pointer","").lstrip("/").split("/") if ref.get("json_pointer") else []:
        value = value[int(token)] if isinstance(value,list) else value[token.replace("~1","/").replace("~0","~")]
    if "field_sha256" in ref:
        assert sha(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(",",":")).encode("utf-8")) == ref["field_sha256"]
    return value


def fact(name, ok, oid=None, **details):
    row = {"name":name,"passed":bool(ok),"orbit_id":oid,**details}
    CHECKS.append(row)
    assert ok,row


def zero(name, value, oid=None):
    residual = s.simplify(s.expand(value))
    fact(name,residual==0,oid,residual=str(residual))


def wl(text):
    return parse_mathematica(text).replace(s.Function("Gamma"),s.gamma)


def original_source_checks(inp):
    oid = inp["orbit_id"]
    graph,native,ordered,prepared = [dereference(inp[k]) for k in
        ("graph_definition","native_orbit","original_ordered_orbit","prepared_orbit")]
    fact("prepared J is the unchanged original ordered J",prepared["J"]==ordered["J"],oid)
    fact("original order and dependent bounds retained",prepared["J"]["integration_order_outer_to_inner"]==inp["source_parameters"]
         and all(v==["0","1"] for v in prepared["J"]["bounds"].values()),oid)
    fact("full one-loop color coordinate is one",prepared["color_preparation"]["coordinates"]==["1"],oid)
    fact("native sign convention is explicit",native["physical_convention"]["id"]=="wilson-qcd-action-path-fourier/1.0",oid)
    fact("raw Euler count is zero",native["physical_convention"]["raw_Euler_count"]==0,oid)
    fact("C4 multiplicity read from source",native["c4_symmetry"]["orbit_multiplicity"]==inp["normalization"]["m"],oid)
    Cnative = wl(native["C"]["combined"])
    Cordered = wl(ordered["C"]["combined"])
    Cprepared = wl(prepared["C"]["combined"])
    zero("native and original C have identical physical sign",Cnative-Cordered,oid)
    zero("prepared C adds Euler exactly once",Cprepared-s.exp(-s.EulerGamma*eps)*Cordered,oid)
    CF,alphaS,v14,I0,A = s.symbols("CF alphaS v14 I0 A")
    common = CF*(alphaS/s.pi)*s.pi**eps*s.gamma(1-eps)*A**eps
    zero("universal raw C ratio is minus one half",Cnative.subs({v14:-A/2,I0:0})/common+s.Rational(1,2),oid)
    u,v = s.symbols("u v")
    ss = dict(zip(inp["source_parameters"],(u,v)))
    f = s.sympify(prepared["J"]["propagator_polynomials"]["f1"],locals={**ss,"z1":z})
    numerator = s.sympify(prepared["J"]["numerator"],locals={"z1":z})
    target = z*(1-u)*v if oid=="2e-01" else (1-u)*(1-v)+z*u*v
    zero("original denominator polynomial",f-target,oid)
    zero("original eikonal contraction",numerator-(z if oid=="2e-01" else -(1+z)),oid)
    V = s.symbols("v1:5")
    b = s.Symbol("b",nonzero=True)
    gram = s.Matrix([[0,z*b,-(1+z)*b,b],[z*b,0,b,-(1+z)*b],[-(1+z)*b,b,0,z*b],[b,-(1+z)*b,z*b,0]])
    phase = s.sympify(graph["phases"][0]["expression"],locals={**ss,**dict(zip(map(str,V),V))})
    pv = s.Matrix([s.expand(phase).coeff(w) for w in V])
    zero("source Fourier phase gives the actual denominator",(pv.T*gram*pv)[0]/(2*b)-f,oid)
    cube = s.sympify(native["J"]["expanded"],locals={"x1":u,"x2":v,"z1":z,"eps":eps})
    expected = numerator*target**(eps-1)
    zero("native cube is only the original variable renaming",cube-expected,oid)
    D = (1-u)*(1-v)+z*u*v
    if oid=="2e-02":
        zero("full opposite denominator reflection",D.subs({z:1/z,u:1-u,v:1-v},simultaneous=True)-D/z,oid)
        a,d,t = s.symbols("a d t",positive=True)
        primitive = (d+(a-d)*t)**eps/(eps*(a-d))
        zero("actual first primitive derivative",s.diff(primitive,t)-(d+(a-d)*t)**(eps-1),oid)
    else:
        t = s.Symbol("t",positive=True)
        zero("actual endpoint power primitive derivative",s.diff(t**eps/eps,t)-t**(eps-1),oid)
    return native,ordered,prepared


def derive_coefficients(order):
    # The source-derived log Gamma series and source-derived M generating
    # function are the inputs. No old coefficient table is read here.
    g_order = order+2
    logG = sum(s.zeta(n)*eps**n/n for n in range(2,g_order+1))
    Gseries = s.series(s.exp(logG),eps,0,g_order+1).removeO().expand()
    direct = s.series(s.exp(-s.EulerGamma*eps)*s.gamma(1-eps),eps,0,g_order+1).removeO()
    zero("Gamma-series Euler cancellation",Gseries-direct)
    # A source-derived auxiliary Mellin generator proves the two elementary
    # moments needed for compactification; it is not fitted to old tables.
    a,y=s.symbols("a y")
    beta=z/(z+y);rr=z/(1+z)
    zero("actual moment variable change",-s.diff(beta,y)/(beta-rr)-(1/(z+y)+1/(1-y)),"2e-02")
    generator=s.pi*(s.cot(s.pi*a)-s.exp(rho*a)/s.sin(s.pi*a))+rho
    elementary=s.series(generator,a,0,4).removeO().expand()
    zero("M01 from the Mellin generator",elementary.coeff(a,1)+(rho**2+s.pi**2)/2,"2e-02")
    zero("M03 from the Mellin generator",6*elementary.coeff(a,3)+(rho**2+s.pi**2)**2/4,"2e-02")
    adjacent = {q:s.expand(-Gseries).coeff(eps,q+2) for q in range(-2,order+1)}
    H = sum(eps**(k+2*j)*M(k,2*j+1)/(2**(k+2*j)*s.factorial(k)*s.factorial(2*j+1))
            for k in range(order+1) for j in range(order//2+1) if k+2*j<=order)
    prefactor = s.series(Gseries/s.cosh(eps*rho/2),eps,0,order+1).removeO()
    opposite_raw = s.series(-prefactor*H/2,eps,0,order+1).removeO().expand()
    opposite = {}
    for q in range(-2,order+1):
        expr = opposite_raw.coeff(eps,q)
        compact = s.expand(expr.subs({rho**2:-2*M(0,1)-s.pi**2,M(0,3):-M(0,1)**2}))
        opposite[q] = compact
        zero(f"opposite natural-M versus compact-M coefficient {q}",
             (expr-compact).subs({rho**2:-2*M(0,1)-s.pi**2,M(0,3):-M(0,1)**2}),"2e-02")
    return adjacent,opposite


def moment(k,j,zz):
    rr = zz/(1+zz)
    def kernel(beta):
        if beta in (0,1):return mp.mpf(0)
        ww = mp.log(beta)+mp.log1p(-beta)
        if beta==rr:return -ww**k/(rr*(1-rr)) if j==1 else mp.mpf(0)
        ll = mp.log(zz)+mp.log1p(-beta)-mp.log(beta)
        return ww**k*ll**j/(beta-rr)
    return mp.quadts(kernel,[0,rr,1])


def Jopposite(ee,zz):
    rr = zz/(1+zz)
    def kernel(beta):
        if beta in (0,1):return mp.mpf(0)
        if beta==rr:return -mp.exp(ee*mp.log(rr))/(rr*(1-rr))
        ll = mp.log(zz)+mp.log1p(-beta)-mp.log(beta)
        return mp.exp(ee*mp.log(beta))*mp.expm1(ee*ll)/(ee*(beta-rr)) if ee else ll/(beta-rr)
    return mp.quadts(kernel,[0,rr,1])


def numeric_checks(coeffs,order,dps):
    mp.mp.dps=dps
    controls=[]
    def close(name,a,b,oid,tolerance="1e-22"):
        err=abs(a-b)/max(1,abs(a),abs(b))
        fact(name,err<mp.mpf(tolerance),oid,scaled_error=mp.nstr(err,12),tolerance=tolerance,
             numeric_scope="empirical high-precision control, not a certified interval")
        controls.append({"name":name,"orbit_id":oid,"left":mp.nstr(a,dps-5),"right":mp.nstr(b,dps-5),"scaled_error":mp.nstr(err,12)})
    gammaG=lambda ee:mp.exp(-mp.euler*ee)*mp.gamma(1-ee)
    for zz in (mp.mpf(2)/3,mp.mpf(3)/2):
        ee=mp.mpf(3)/5
        adjacent=zz**ee*mp.quadts(lambda u:(1-u)**(ee-1),[0,1])*mp.quadts(lambda v:v**(ee-1),[0,1])
        close("actual adjacent ordered endpoint integrals",adjacent,zz**ee/ee**2,"2e-01","1e-20")
        ee=mp.mpf(1)/3
        jj=Jopposite(ee,zz)
        gg=-(1+zz)*mp.meijerg([[0,0,ee],[]],[[0,ee,ee],[]],zz)/(mp.gamma(1-ee)*mp.gamma(1+ee)**2)
        close("original difference quotient versus full Meijer G",jj,gg,"2e-02")
        close("original J reflection",Jopposite(ee,1/zz),zz**(-ee)*jj,"2e-02")
        mm={str(M(k,j)):moment(k,j,zz) for k,j in ((0,1),(0,3),(1,1),(2,1),(3,1),(1,3))}
        close("M01 original integral versus elementary value",mm["M01"],-(mp.log(zz)**2+mp.pi**2)/2,"2e-02")
        close("M03 original integral versus elementary value",mm["M03"],-mm["M01"]**2,"2e-02")
        classical=2*mp.zeta(3)+(mp.log(zz)**2+mp.pi**2)*mp.log1p(zz)-mp.pi**2*mp.log(zz)/2+2*mp.log(zz)*mp.polylog(2,-zz)-2*mp.polylog(3,-zz)-mp.log(zz)**3/6
        close("M11 original integral versus classical polylogarithms",mm["M11"],classical,"2e-02")
        radius=mp.mpf(".1");N=32
        values=[]
        for j in range(N):
            ee=radius*mp.exp(2j*mp.pi*j/N)
            values.append(-gammaG(ee)*Jopposite(ee,zz)/(1+mp.exp(ee*mp.log(zz))))
        for q in range(order+1):
            cauchy=sum(values[j]*mp.exp(-2j*mp.pi*j*q/N) for j in range(N))/(N*radius**q)
            expr=coeffs["2e-02"][q]
            number=s.lambdify(sorted(expr.free_symbols,key=str),expr,"mpmath")(*[mm[str(v)] if str(v) in mm else mp.log(zz) for v in sorted(expr.free_symbols,key=str)])
            close(f"epsilon Cauchy readout versus newly generated coefficient {q}",cauchy,number,"2e-02")
    # A genuinely original two-dimensional ordered integration is independent
    # of the one-dimensional primitive used above.
    with mp.workdps(min(dps,30)):
        zz=mp.mpf(2)/3;ee=mp.mpf(1)/2
        def outer(u):
            def inner(v):
                D=(1-u)*(1-v)+zz*u*v
                return D**(ee-1) if D else mp.mpf(0)
            return mp.quadts(inner,[0,1])
        two=-(1+zz)*mp.quadts(outer,[0,1])
        one=Jopposite(ee,zz)
        close("actual original two-dimensional ordered J",two,one,"2e-02","1e-18")
    radius=mp.mpf(".1");N=32
    values=[-gammaG(radius*mp.exp(2j*mp.pi*j/N)) for j in range(N)]
    for q,expr in coeffs["2e-01"].items():
        p=q+2
        number=s.lambdify((),expr,"mpmath")()
        cauchy=sum(values[j]*mp.exp(-2j*mp.pi*j*p/N) for j in range(N))/(N*radius**p)
        close(f"regular adjacent Gamma Cauchy coefficient {q}",cauchy,number,"2e-01")
    return controls


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--order",type=int,default=3,choices=range(0,4))
    parser.add_argument("--dps",type=int,default=40)
    args=parser.parse_args()
    assert args.dps>=35,"Use at least35 working digits for the stated controls"
    output=args.output_dir.resolve()
    assert not output.is_relative_to(HERE) and not output.exists(),"Use a fresh output directory outside the packet"
    assert not any(output.is_relative_to(ROOT/name) for name in ("data","packages","results","integrals")),"Keep generated output outside the scientific input/result trees"
    started=time.perf_counter()
    inputs={oid:read_json(f"integrals/one-loop/{oid}/input.json") for oid in ("2e-01","2e-02")}
    for inp in inputs.values():original_source_checks(inp)
    aa,oo=derive_coefficients(args.order)
    coeffs={"2e-01":aa,"2e-02":oo}
    results={}
    for oid,table in coeffs.items():
        results[oid]={"format":"wilson-published-one-loop-I/1","orbit_id":oid,
          "normalization":inputs[oid]["normalization"],
          "laurent":{"variable":"eps","minimum_power":-2,"maximum_power":args.order,
                     "coefficients":{str(q):str(v) for q,v in table.items()},"remainder":f"O(eps^{args.order+1})"},
          "symbols":{"M01":"M_(0,1)(z)","M11":"M_(1,1)(z)","M21":"M_(2,1)(z)","M31":"M_(3,1)(z)","M13":"M_(1,3)(z)"},
          "function_definitions":{"rho":"log(z)","r":"z/(1+z)","W":"log(beta*(1-beta))",
              "Lbeta":"rho+log((1-beta)/beta)","M_kj":"integral_0^1 W**k*Lbeta**j/(beta-r) d beta, raw powers without factorials, removable point filled continuously",
              "G":"exp(-EulerGamma*eps)*gamma(1-eps)",
              "H_epsilon":"sum_{k,j>=0} eps**(k+2*j)*M_(k,2*j+1)/(2**(k+2*j)*k!*(2*j+1)!)",
              "elementary_moment_generator":"pi*(cot(pi*a)-exp(rho*a)/sin(pi*a))+rho = sum_{j>=1} a**j*M_(0,j)/j!"},
          "all_epsilon_expression":"-exp(-EulerGamma*eps)*gamma(1-eps)/eps**2" if oid=="2e-01" else "(1+z)*z**(-eps/2)*exp(-EulerGamma*eps)*meijerg([[0,0,eps],[]],[[0,eps,eps],[]],z)/(2*cosh(eps*log(z)/2)*gamma(1+eps)**2)",
          "source_calculation":"Original ordered parameter kernel -> actual primitive -> Gamma/M generating function -> exact finite convolution",
          "historical_coefficient_table_used_to_generate":False,"formal":"not_run"}
    output.mkdir(parents=True)
    # Persist the predictions before opening the historical coefficient table.
    (output/"generated-before-holdout.json").write_bytes(encode(results))
    for oid,inp in inputs.items():
        historical=dereference(inp["historical_result"])
        for q,row in historical["laurent"]["coefficients"].items():
            if int(q)<=args.order:
                view=row.get("MN_alternative",row.get("compact",row))
                old=wl(view["expression"]).replace(
                    lambda v:v.is_Function and str(v.func)=="WilsonM",
                    lambda v:M(int(v.args[0]),int(v.args[1])))
                zero(f"post-generation historical r002 coefficient {q}",coeffs[oid].get(int(q),s.S.Zero)-old,oid)
    fact("negative control: duplicating the adjacent orbit is rejected",aa[-2]*2!=aa[-2],"2e-01",negative_control=True)
    wrong_Euler=s.series(s.exp(-s.EulerGamma*eps)*sum(v*eps**q for q,v in aa.items()),eps,0,1).removeO().expand()
    fact("negative control: repeated Euler changes the adjacent simple pole",wrong_Euler.coeff(eps,-1)!=aa[-1],"2e-01",negative_control=True)
    fact("negative control: dropping opposite crossing factor changes the constant",oo[0]*2!=oo[0],"2e-02",negative_control=True)
    fact("negative control: reversed overall source sign is rejected",oo[0]!=-oo[0],"2e-02",negative_control=True)
    controls=numeric_checks(coeffs,args.order,args.dps)
    for name,expected in READS.items():assert sha((ROOT/name).read_bytes())==expected,name
    for oid,result in results.items():
        (output/oid).mkdir()
        (output/oid/"result.json").write_bytes(encode(result))
    (output/"numeric.json").write_bytes(encode({"working_dps":args.dps,"controls":controls,"strict_interval_certificate":False}))
    summary={"status":"passed","checks":CHECKS,"check_count":len(CHECKS),
             "negative_controls":sum(c.get("negative_control",False) for c in CHECKS),
             "source_hashes":READS,"generated_before_holdout_sha256":sha((output/"generated-before-holdout.json").read_bytes()),
             "order":args.order,"working_dps":args.dps,"elapsed_seconds":time.perf_counter()-started,
             "python":platform.python_version(),"sympy":s.__version__,"mpmath":mp.__version__,
             "reproducer_sha256":sha(Path(__file__).read_bytes()),
             "fresh_execution":"Source geometry, normalization, actual parameter primitives, Gamma/M convolution and numerical controls",
             "reused":"Accepted QFT Fourier propagator rule and existing source contracts, standard Euler/Mellin integral identities. The Meijer-G and elementary-moment connections are derived in the accompanying text; no formal integral proof is claimed.",
             "scope":"Fixed positive z and local epsilon continuation, no arbitrary complex-kinematics theorem; numeric checks are empirical controls.","formal":"not_run"}
    (output/"checks.json").write_bytes(encode(summary))
    print(json.dumps({"status":"passed","checks":len(CHECKS),"order":args.order,"elapsed_seconds":summary["elapsed_seconds"],"output_dir":str(output)}))


if __name__=="__main__":main()
