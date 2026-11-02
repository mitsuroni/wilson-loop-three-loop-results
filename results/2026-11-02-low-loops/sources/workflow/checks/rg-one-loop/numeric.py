"""Independent contour-integral/moment checks of the one-loop RG assembly."""
from pathlib import Path
import hashlib
import json
import sys
import mpmath as mp
import sympy as sp
from scipy.integrate import dblquad

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"workflow/checks/single-opposite"))
import hlog_algebra as ha
from numeric_hlog import HlogEvaluator

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def opposite_J(ep,z):
    """Actual s21 outer integral after analytically integrating s41."""
    def inner(u):
        a=1-u
        b=z*u
        if a==b: return a**(ep-1)
        if a==0: return b**(ep-1)/ep if ep else mp.inf
        if b==0: return a**(ep-1)/ep if ep else mp.inf
        logratio=mp.log(b/a)
        if not ep: return logratio/(b-a)
        return a**ep*mp.expm1(ep*logratio)/(ep*(b-a))
    u0=1/(1+z)
    return -(1+z)*mp.quadts(inner,[0,u0,1])

def moment(k,r,z):
    a=z/(1+z)
    def integrand(b):
        if b==0 or b==1: return mp.mpf(0)  # endpoints of measure zero only
        R=mp.log(b)+mp.log1p(-b)
        h=mp.log(z)+mp.log1p(-b)-mp.log(b)
        if b==a: return -R**k/(a*(1-a)) if r==1 else mp.mpf(0)
        return R**k*h**r/(b-a)
    return mp.quadts(integrand,[0,a,1])

def H0_g(ep,z):
    return -(1+z)*z**(-ep/2)*mp.meijerg([[ep,0,0],[]],[[0,ep,ep],[]],z)/(
        mp.gamma(1-ep)*mp.gamma(1+ep)**2)

def main():
    mp.mp.dps=70
    rows=[]
    def check(name,passed,**extra):
        row={"name":name,"passed":bool(passed),**extra}
        rows.append(row)
        if not passed: raise AssertionError(row)
    moments=ROOT/"workflow/checks/single-opposite/runs/20260916T090538777022Z"
    vectors={(k,r):ha.parse(moments/f"M_{k}_{r}.txt")
             for k,r in [(0,1),(1,1),(2,1),(0,3)]}
    cached={}
    for zstr in (".3",".7","1.4"):
        z=mp.mpf(zstr)
        ev=HlogEvaluator(z if z<1 else 1/z,terms=420)
        cached[zstr]={}
        for (k,r),vector in vectors.items():
            integ=moment(k,r,z)
            hlog=ev.vector(vector)  # all four r odd: exact z<->1/z symmetry
            err=abs(integ-hlog)/max(1,abs(integ))
            check(f"moment_vs_existing_Hlog_{k}_{r}_z{zstr}",err<mp.mpf("1e-42"),
                  relative_error=mp.nstr(err,8),working_dps=70)
            cached[zstr][(k,r)]=hlog
        elementary=-(mp.log(z)**2+mp.pi**2)/2
        err=abs(cached[zstr][(0,1)]-elementary)
        check("M01_elementary_z"+zstr,err<mp.mpf("1e-42"),absolute_error=mp.nstr(err,8))
    print("15 moment / elementary checks passed",flush=True)
    for es,zs in [(".21",".3"),(".37",".7"),(".17","1.4")]:
        ep,z=mp.mpf(es),mp.mpf(zs)
        j=opposite_J(ep,z)
        exact=z**(ep/2)*H0_g(ep,z)
        err=abs(j-exact)/max(1,abs(j))
        check("original_order_reduction_vs_MeijerG_"+es+"_"+zs,err<mp.mpf("1e-42"),
              relative_error=mp.nstr(err,8),working_dps=70)
        reflected=z**ep*opposite_J(ep,1/z)
        check("opposite_crossing_"+es+"_"+zs,abs(j-reflected)<mp.mpf("1e-42"),
              absolute_error=mp.nstr(abs(j-reflected),8))
    print("6 full-epsilon and crossing checks passed",flush=True)
    for ep,z in [(.4,.7),(.3,1.4)]:
        val,estimate=dblquad(lambda v,u:-(1+z)*((1-u)*(1-v)+z*u*v)**(ep-1),
                             0,1,lambda _:0,lambda _:1,epsabs=2e-9,epsrel=2e-9)
        ref=opposite_J(mp.mpf(str(ep)),mp.mpf(str(z)))
        err=abs(mp.mpf(val)-ref)
        check(f"original_two_dimensional_eps{ep}_z{z}",err<mp.mpf("2e-8"),
              absolute_error=mp.nstr(err,8),integrator_error_estimate=estimate,
              precision="SciPy double, not 70-digit or certified interval")
    result_path=ROOT/"benchmarks/renormalization/one-loop/result-r001.json"
    result=json.loads(result_path.read_text(encoding="utf-8"))
    names=sp.symbols("L_s L_t M01 M11 M21 M03")
    loc={str(x):x for x in names}
    coeffs={int(k):sp.lambdify(names,sp.sympify(v["sympy"],locals=loc),"mpmath")
            for k,v in result["laurent_unit_coefficients"].items()}
    for index,zs in enumerate((".3",".7","1.4")):
        z=mp.mpf(zs);Lt=mp.mpf(".2")+index/3;Ls=Lt+mp.log(z)
        m=cached[zs]
        vals=[Ls,Lt,m[(0,1)],m[(1,1)],m[(2,1)],m[(0,3)]]
        cs={k:f(*vals) for k,f in coeffs.items()}
        errors=[]
        for denominator in (32,64,128,256):
            ep=mp.mpf(1)/denominator
            G=mp.exp(-mp.euler*ep)*mp.gamma(1-ep)
            exact=-G*(mp.exp(ep*Ls)+mp.exp(ep*Lt))/ep**2-G*mp.exp(ep*Lt)*opposite_J(ep,z)
            expansion=sum(cs[k]*ep**k for k in cs)
            errors.append(abs(exact-expansion))
        ratios=[errors[i+1]/errors[i] for i in range(3)]
        check("full_original_I_assembly_Oeps3_z"+zs,
              all(abs(x-mp.mpf(1)/8)<mp.mpf(".025") for x in ratios),
              errors=[mp.nstr(x,16) for x in errors],ratios=[mp.nstr(x,14) for x in ratios])
        def finite(logmu):
            return -(Ls+2*logmu)*(Lt+2*logmu)+mp.pi**2/3
        delta=mp.mpf("1e-10")
        derivative=(finite(delta)-finite(-delta))/(2*delta)
        err=abs(derivative+2*(Ls+Lt))
        check("finite_mu_derivative_z"+zs,err<mp.mpf("1e-45"),
              absolute_error=mp.nstr(err,8))
    report={"checks":rows,"passed":sum(x["passed"] for x in rows),"total":len(rows),
            "code_sha256":sha(__file__),"result_sha256":sha(result_path),
            "limits":"Empirical numerical controls; no interval certification or physical complex continuation."}
    (HERE/"numeric-checks.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"numeric_checks":len(rows),"passed":True}),flush=True)

if __name__=="__main__": main()
