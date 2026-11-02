"""Portable coefficient engines adapted from the preserved 2026-10-15 builders.

These expand evaluated Gamma/moment mother formulas, not old finite tables.
The integrated TGV04 HLOG values are an explicitly reused input.
"""
import math,re
import sympy as s

eps,ell=s.symbols('eps ell')
def M(k,j):return s.Symbol(f'M_{k}_{j}')
def N(k,j):return s.Symbol(f'N_{k}_{j}')

def parse(text):
    text=re.sub(r'\b([MN])(\d)_(\d)\b',r'\1_\2_\3',str(text))
    text=re.sub(r'\b([MN])(\d)(\d)\b',r'\1_\2_\3',text)
    text=re.sub(r'\bzeta_?(\d+)\b',r'zeta(\1)',text)
    return s.sympify(text,locals={'tauD4':s.Symbol('tau_D4')})

def canonical(x):
    x=s.expand(x.subs({M(0,3):-M(0,1)**2,N(0,2):M(0,1)**2/3+4*s.zeta(2)*M(0,1)}))
    return s.expand(sum(c*(-2*M(0,1)-s.pi**2)**(p[0]//2)*ell**(p[0]%2) for p,c in s.Poly(x,ell).terms()))

def exp_coeff(logs,n):
    out=[s.S.One]
    for q in range(1,n+1):out.append(s.expand(sum(k*logs.get(k,0)*out[q-k] for k in range(1,q+1))/q))
    return out

def at(a,q):return a[q] if 0<=q<len(a) else s.S.Zero
def conv(a,b,q):return sum(at(a,i)*at(b,q-i) for i in range(q+1)) if q>=0 else s.S.Zero
def mul(a,b,n):return [s.expand(conv(a,b,q)) for q in range(n+1)]
def sech_times(a,n):
    out=[]
    for q in range(n+1):out.append(s.expand(at(a,q)-sum(ell**(2*k)*out[q-2*k]/s.factorial(2*k) for k in range(1,q//2+1))))
    return out

def reconstruct(tgv04_compact):
    zsq=exp_coeff({k:2*s.zeta(k)/k for k in range(2,7)},6)
    r2=exp_coeff({k:(-1)**k*(2-2**k)*s.zeta(k)/k for k in range(2,7)},6)
    pref=sech_times(zsq,6)
    h=[sum(M(q-2*j,2*j+1)/s.Integer(2**q*math.factorial(q-2*j)*math.factorial(2*j+1)) for j in range(q//2+1)) for q in range(7)]
    nv=[sum(N(q-2*j,2*j+2)/s.Integer(math.factorial(q-2*j)*math.factorial(2*j+2)) for j in range(q//2+1)) for q in range(3)]
    d=[s.expand((2**q-1)*h[q]-sum((ell/2)**(2*k)*h[q-2*k]/s.factorial(2*k) for k in range(1,q//2+1))) for q in range(5)]
    qc=[s.expand(sum(h[i]*h[q-i] for i in range(q+1))/2-2*nv[q]+2**(q+2)*sum((-1)**k*s.zeta(k+2)*h[q-k] for k in range(q+1))) for q in range(3)]
    result={}
    for oid in [f'4e-{i:02}' for i in range(1,7)]:
        cs={}
        for q in range(-4,3):
            if oid=='4e-01':v=-at(zsq,q+4)/16
            elif oid in {'4e-02','4e-05'}:v=conv(pref,d,q+2)/4
            elif oid=='4e-03':v=conv(pref,[a-(1 if i==0 else 0) for i,a in enumerate(r2)],q+4)/4
            elif oid=='4e-04':v=-conv(pref,qc,q)/8
            else:v=-conv(pref,mul(h,h,2),q)/16
            cs[q]=canonical(v)
        result[oid]=cs
    G2=exp_coeff({k:s.Rational(2**k,k)*s.zeta(k) for k in range(2,7)},6)
    F=exp_coeff({k:s.Rational(2,k)*s.zeta(k) for k in range(2,7,2)},6)
    be=mul(exp_coeff({k:s.Rational(2-2**k,k)*s.zeta(k) for k in range(2,7)},6),[s.Integer(2)**q for q in range(7)],6)
    a01=mul(G2,[2*F[q]-be[q] for q in range(7)],6)
    result['3e-3g-01']={q:canonical(at(a01,q+4)/32) for q in range(-4,3)}
    result['3e-3g-02']=dict(result['3e-3g-01'])
    VM=[sum(M(q-2*j,2*j+1)/s.Integer(math.factorial(q-2*j)*math.factorial(2*j+1)) for j in range(q//2+1)) for q in range(5)]
    VM0=[M(0,q+1)/s.factorial(q+1) if q%2==0 else s.S.Zero for q in range(5)]
    diff=[s.expand(VM0[q]-conv(be,VM,q)) for q in range(5)]
    a03=mul(sech_times(G2,4),diff,4)
    result['3e-3g-03']={q:canonical(at(a03,q+2)/8) for q in range(-4,3)}
    result['3e-3g-04']={-4:s.S.Zero,-3:s.S.Zero,-2:s.zeta(2)/8,-1:(s.zeta(3)-M(1,1))/8,0:M(0,1)**2/48-M(2,1)/8+9*s.zeta(2)**2/10,1:parse(tgv04_compact['I1']),2:parse(tgv04_compact['I2'])}
    rat=[s.Integer(2)**q for q in range(6)]
    ratA=mul(rat,[s.Integer(-1)**q for q in range(6)],5)
    adj=mul(zsq,ratA,5)
    opposite=mul(sech_times(mul(zsq,rat,5),4),VM[:4],3)
    result['2e-Pi1-01']={q:canonical(-at(adj,q+3)/32) for q in range(-4,3)}
    result['2e-Pi1-02']={q:canonical(-at(opposite,q+1)/16-at(adj,q+2)/32) for q in range(-4,3)}
    return result
