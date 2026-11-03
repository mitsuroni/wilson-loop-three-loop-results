"""Full-epsilon T1[] formulas and exact Laurent assembly from fresh HyperInt moments."""
from pathlib import Path
import sys,json,re
import sympy as sp
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
e,ell=sp.symbols("eps ell")
IDS=["01","02","03","04","17","18","19"]
def run():return ROOT/json.loads((HERE/"current-run.json").read_text())["run"]
def moment(name,q):
    text=(run()/f"{name}{q}.txt").read_text().strip()
    text=re.sub(r"zeta\[(\d+)\]",r"zeta(\1)",text)
    return sp.sympify(text.replace("^","**"),locals={"zeta":sp.zeta})
def chop(expr,n=7):return sp.series(expr,e,0,n).removeO().expand()
def exp_series(coeff,n=6):
    # exp(sum c_k e^k) via exact recurrence, avoiding unnecessary Gamma expansion.
    b=[sp.Integer(1)]
    for m in range(1,n+1):
        b.append(sp.expand(sum(k*coeff.get(k,0)*b[m-k] for k in range(1,m+1))/m))
    return sum(b[k]*e**k for k in range(n+1))
def gamma_series():
    G=exp_series({k:(1+2**k)*sp.zeta(k)/k for k in range(2,7)})
    A=exp_series({k:(1+(-1)**k)*sp.zeta(k)/k for k in range(2,7)})
    B=chop(exp_series({k:(2-2**k)*sp.zeta(k)/k for k in range(2,7)})*sum((2*e)**k for k in range(7)))
    V=exp_series({k:(-1)**k*(1+2**k-3**k)*sp.zeta(k)/k for k in range(2,7)})
    W=exp_series({k:(-1)**k*(2-2**k)*sp.zeta(k)/k for k in range(2,7)})
    R=chop(sp.cosh(e*ell/2)/sp.cosh(3*e*ell/2))
    return G,A,B,V,W,R
def coefficients():
    G,A,B,V,W,R=gamma_series();U=2*A-B
    F=sum(moment("F",q)*e**q for q in range(6))
    Z=sum(moment("Z",q)*e**q for q in range(1,5))
    H1=B/(2*e)-2*F
    # H0 is substituted directly through U to avoid series truncation in 1/lambda.
    I01=G*U/(576*e**6)
    I04=(1-2*e)*G*H1/(576*e**5)
    I03=I01-I04
    I17=G*R*U*(1-V)/(128*e**6)
    Q=B*(1-2*W+V)/(2*e**2)-2*Z
    I19=(1-2*e)*G*R*Q/(64*e**4)
    expressions={"01":I01,"02":2*I03,"03":I03,"04":I04,"17":I17,"18":I17-I19,"19":I19}
    return {oid:{k:sp.factor(sp.expand(f).coeff(e,k)) for k in range(-6,1)} for oid,f in expressions.items()}
def numerical_full(ep,z,dps=55):
    import mpmath as mp
    mp.mp.dps=dps;ep=mp.mpf(ep);z=mp.mpf(z)
    A=mp.gamma(1+ep)*mp.gamma(1-ep);B=mp.beta(1-ep,1-ep)
    V=mp.gamma(1+ep)*mp.gamma(1+2*ep)/mp.gamma(1+3*ep)
    W=mp.gamma(1+ep)**2/mp.gamma(1+2*ep)
    U=2*A-B;lam=1-2*ep
    G=mp.exp(-3*mp.euler*ep)*mp.gamma(1-ep)*mp.gamma(1-2*ep)
    R=mp.cosh(ep*mp.log(z)/2)/mp.cosh(3*ep*mp.log(z)/2)
    F=mp.beta(2-ep,1-ep)/(1+2*ep)*mp.hyper([1,1,2-ep],[2+2*ep,3-2*ep],1)
    def z_kernel(r):
        if not r:return mp.mpf(0)
        L=mp.beta(2-ep,1-ep)*mp.hyp2f1(1,2-ep,3-2*ep,1-r)
        Ydiff=(1-mp.hyp2f1(-ep,2*ep,1+2*ep,r))/(2*ep)
        return r**ep*L*Ydiff
    Z=mp.quad(z_kernel,[0,mp.mpf(".25"),mp.mpf(".75"),1])
    Q=B*(1-2*W+V)/(2*ep**2)-2*Z
    H0=U/(lam*ep);H1=B/(2*ep)-2*F
    js={"01":H0/(72*ep**5),"02":(H0-H1)/(36*ep**5),"03":(H0-H1)/(72*ep**5),
       "04":H1/(72*ep**5),"17":H0*(1-V)/(16*ep**5),"19":Q/(8*ep**4)}
    js["18"]=js["17"]-js["19"]
    I={k:lam*G/8*j*(1 if k in IDS[:4] else R) for k,j in js.items()}
    J={k:j*z**((3 if k in IDS[:4] else 2)*ep) for k,j in js.items()}
    return {"I":I,"J":J,"F":F,"Z":Z}
def main():
    c=coefficients();r=run()
    data={oid:{str(k):str(v) for k,v in cs.items()} for oid,cs in c.items()}
    (r/"coefficients.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    for oid,cs in c.items():
        print(oid, "pole",min(k for k,v in cs.items() if v!=0))
        for k,v in cs.items():print(k,sp.sstr(v))
if __name__=="__main__":main()
