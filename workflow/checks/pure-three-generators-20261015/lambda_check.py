"""Validate the all-epsilon Lambda single-index construction from its full kernel.

This file does not edit old results. Exact controls use positive integer epsilon,
where the two-variable defining integrand is polynomial after removable factors.
Noninteger controls compare endpoint moments and the complete sector source.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sympy as s
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
checks = []

def record(name, ok, **kw):
    checks.append(dict(name=name, passed=bool(ok), **kw))
    if not ok:
        raise AssertionError(name)

def binomial_coeffs(a, n):
    out = [s.Integer(1)]
    for k in range(1, n+1):
        out.append(s.factor(out[-1]*(k-1-a)/k))
    return out

def conv(a, b, n):
    return sum(a[j]*b[n-j] for j in range(n+1))

def polynomial_control(e, A, B, C, d):
    u, v = s.symbols('u v')
    nmax = 3*e
    be = binomial_coeffs(s.Integer(e), nmax)
    b2 = binomial_coeffs(s.Integer(2*e), nmax)
    je = [s.Integer(0)] + [be[n]/n for n in range(1, nmax+1)]
    he = [conv(be, je, n) for n in range(nmax+1)]
    h2 = [conv(b2, je, n) for n in range(nmax+1)]
    J = lambda x: sum(je[n]*x**n for n in range(1, e+1))
    h = lambda x: s.expand((1-x)**e*J(x))
    def moment(expr, X, k):
        return sum(co*X**n/(n+k) for (n,), co in s.Poly(expr(u),u).terms())
    JA, JB, JD = J(A), J(B), J(d)
    result = JA*JB*(JD-(JA+JB)/2)
    for k in range(1, nmax+1):
        BA = sum(be[n]*A**n/(n+k) for n in range(e+1))
        BB = sum(be[n]*B**n/(n+k) for n in range(e+1))
        JA_k = sum(je[n]*A**n/(n+k) for n in range(1,e+1))
        JB_k = sum(je[n]*B**n/(n+k) for n in range(1,e+1))
        HA = moment(h, A, k)
        HB = moment(h, B, k)
        sk = (h2[k]/k**2 - he[k]*(BA+BB)/k + je[k]*BA*BB
              - b2[k]*(JA_k+JB_k)/k + be[k]*(HA+HB)/k
              + be[k]*(JA_k*BB+BA*JB_k)
              + JD*(b2[k]/k**2-be[k]*(BA+BB)/k))
        result += C**k*sk
    K = (1-C*u*v)**e
    kernel = s.cancel((J(C*u*v)-J(A*u)-J(B*v)+JD)
                      *(K-(1-A*u)**e)*(K-(1-B*v)**e)/(u*v))
    poly = s.Poly(s.expand(kernel),u,v)
    direct = sum(co/s.Integer((i+1)*(j+1)) for (i,j),co in poly.terms())
    record(f'complete_sector_integer_e{e}_{A}_{B}_{C}', s.factor(result-direct)==0)
    # A single sign or missing zero-mode must actually be detectable.
    if e == 2 and A == s.Rational(1,3):
        record('negative_control_missing_zero_mode', s.factor(JA*JB*(JD-(JA+JB)/2)) != 0)

def exact_controls():
    v,a,e = s.symbols('v a e', positive=True)
    # Derivatives of J_r(1-v) are (1-v**r)/(1-v).
    primitive_derivative = (s.diff(v**(a+1),v)*s.Symbol('J')
        +v**(a+1)*(1-v**e)/(1-v)
        +(1-v**(a+1))/(1-v)-(1-v**(a+e+1))/(1-v))/(a+1)
    record('mixed_endpoint_primitive', s.simplify(primitive_derivative-v**a*s.Symbol('J'))==0)
    u,A,k=s.symbols('u A k', positive=True)
    r=s.symbols('r', positive=True)
    # Divide the differentiated weighted product by u^(k-1) r^e,
    # where r=1-A*u; J is retained as an independent symbol.
    JJ=s.Symbol('J')
    lhs=k*r*JJ-A*(e+1)*u*JJ+r*(r**e-1)
    rhs=(k-A*(k+e+1)*u)*JJ+r**(e+1)-r
    record('H_transport_total_derivative',s.expand(lhs-rhs).subs(r,1-A*u).simplify()==0)
    record('B_transport_total_derivative',s.expand(k*(1-A*u)-A*(a+1)*u-(k-A*(k+a+1)*u))==0)
    c = s.symbols('c')
    b_c = binomial_coeffs(c,12)
    b_e = binomial_coeffs(e,12)
    j_e = [s.Integer(0)]+[b_e[n]/n for n in range(1,13)]
    h = [conv(b_c,j_e,n) for n in range(13)]
    driver1=binomial_coeffs(c+e+1,12)
    driver0=binomial_coeffs(c+1,12)
    for n in range(1,13):
        record(f'h_coefficient_recurrence_{n}', s.expand(n*h[n]+(c-n+1)*h[n-1]-driver1[n]+driver0[n])==0)
    for ep in (1,2,3):
        polynomial_control(ep,s.Rational(1,3),-s.Rational(2,5),s.Rational(1,4),-s.Rational(3,7))
        for z in (s.Rational(1,2),s.Integer(1),s.Integer(2)):
            for A,B,C in ((1,1,-z),(-1/z,-1/z,-1/z),(1,-1/z,1)):
                polynomial_control(ep,s.sympify(A),s.sympify(B),s.sympify(C),-1/z)

def jfun(e, x):
    if x == 0:
        return mp.mpf(0)
    if x == 1:
        return -mp.euler-mp.digamma(1+e)
    if e == mp.mpf('0.5'):
        y=mp.sqrt(1-x)
        return 2*(y-1)-2*mp.log((1+y)/2)
    if x < -mp.mpf('.5') and 0 < e < 1:
        y=-x; r=y/(1+y)
        return jfun(-e,r)+((1+y)**e-1)/e-mp.log(1+y)
    if 0.5 < x < 1:
        q=1-x
        return -mp.euler-mp.digamma(1+e)-mp.log(x)-q**(e+1)*mp.hyper([1,e+1],[e+2],q)/(e+1)
    return -e*x*mp.hyper([1,1,1-e],[2,2],x)

def h1(e,A):
    return -((1-A)**(e+1)*jfun(e,A)+jfun(e+1,A)-jfun(2*e+1,A))/(A*(e+1))

def moment_arrays(e,A,N):
    exponents=(e,e+1,2*e+1)
    bs={q:[None,(1-(1-A)**(q+1))/(A*(q+1))] for q in exponents}
    hs=[None,h1(e,A)]
    J=jfun(e,A)
    for k in range(1,N):
        hs.append((k*hs[k]+bs[2*e+1][k]-bs[e+1][k]-(1-A)**(e+1)*J)/(A*(k+e+1)))
        for q in exponents:
            bs[q].append((k*bs[q][k]-(1-A)**(q+1))/(A*(k+q+1)))
    return bs[e],hs

def sector_series(e,A,B,C,d,N,return_coefficients=False):
    BA,HA=moment_arrays(e,A,N)
    BB,HB=moment_arrays(e,B,N)
    JA,JB,JD=jfun(e,A),jfun(e,B),jfun(e,d)
    result=JA*JB*(JD-(JA+JB)/2)
    exponents=[e,2*e,2*e+1,e+1,3*e+1]
    bs={q:mp.mpf(1) for q in exponents}
    he=h2=mp.mpf(0)
    tail=[]
    for k in range(1,N+1):
        for q in exponents:
            bs[q] *= (k-1-q)/k
        he=(bs[2*e+1]-bs[e+1]-(e-k+1)*he)/k
        h2=(bs[3*e+1]-bs[2*e+1]-(2*e-k+1)*h2)/k
        JA_k=JA/k-BA[k]/k+mp.mpf(1)/k**2
        JB_k=JB/k-BB[k]/k+mp.mpf(1)/k**2
        sk=(h2/k**2-he*(BA[k]+BB[k])/k+bs[e]*BA[k]*BB[k]/k
            -bs[2*e]*(JA_k+JB_k)/k+bs[e]*(HA[k]+HB[k])/k
            +bs[e]*(JA_k*BB[k]+BA[k]*JB_k)
            +JD*(bs[2*e]/k**2-bs[e]*(BA[k]+BB[k])/k))
        term=C**k*sk
        result+=term
        tail.append(term)
    return result,tail if return_coefficients else tail[-1]

def sector_quadrature(e,A,B,C,d):
    # Gauss-Legendre quadrature evaluates the complete removable integrand.
    # This independent route does not use moment or coefficient recurrences.
    JD=jfun(e,d)
    def fixed(n):
        nodes, weights=mp.gauss_quadrature(n,'legendre')
        nodes=[(x+1)/2 for x in nodes]; weights=[w/2 for w in weights]
        ju=[jfun(e,A*u) for u in nodes]; jv=[jfun(e,B*v) for v in nodes]
        pu=[(1-A*u)**e for u in nodes]; pv=[(1-B*v)**e for v in nodes]
        out=mp.mpf(0)
        for i,u in enumerate(nodes):
            for j,v in enumerate(nodes):
                K=(1-C*u*v)**e
                out+=weights[i]*weights[j]*(jfun(e,C*u*v)-ju[i]-jv[j]+JD)*(K-pu[i])*(K-pv[j])/(u*v)
        return out
    low,high=fixed(20),fixed(28)
    record('complete_sector_quadrature_refinement',abs(low-high)<mp.mpf('1e-40'),absolute_difference=str(abs(low-high)))
    return high

def numeric_controls():
    mp.mp.dps=70
    for ep,A in ((mp.mpf('.3'),mp.mpf('.2')),(mp.mpf('.3'),mp.mpf('-.7')),(mp.mpf('.75'),mp.mpf(1))):
        bs,hs=moment_arrays(ep,A,6)
        for k in (1,2,5):
            direct=mp.quad(lambda u:u**(k-1)*(1-A*u)**ep*jfun(ep,A*u),[0,1])
            err=abs(direct-hs[k])
            record(f'noninteger_H_{ep}_{A}_{k}',err<mp.mpf('1e-55'),absolute_error=str(err))
    print('endpoint controls passed; complete-sector series/quadrature',flush=True)
    # Strictly inside all three Taylor circles: independent full-sector check.
    ep=mp.mpf('.37'); A=mp.mpf('.2'); B=mp.mpf('-.3'); C=mp.mpf('.15'); d=mp.mpf('-.4')
    # Upward moment recurrences need working precision exceeding cancellation.
    series,last=sector_series(ep,A,B,C,d,75)
    direct=sector_quadrature(ep,A,B,C,d)
    error=abs(series-direct)
    record('noninteger_complete_sector_interior',error<mp.mpf('1e-40'),absolute_error=str(error),series=str(series),quadrature=str(direct),last_term=str(last))

if __name__=='__main__':
    print('exact controls',flush=True)
    exact_controls()
    print(f'{len(checks)} exact controls passed; noninteger controls',flush=True)
    numeric_controls()
    source_paths=['benchmarks/2ne/theory/generation/three-kernel-eps2.md','benchmarks/2ne/theory/generation/natural-beta-source.md','benchmarks/2ne/theory/generation/three-kernel-generators.md']
    report=dict(passed=all(c['passed'] for c in checks),checks=checks,formal='not_run',
                scope='complete Lambda sector coefficient sum and endpoint recurrence; not an original six-dimensional all-order numerical proof',
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths})
    (HERE/'lambda-verification-r001.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':report['passed'],'checks':len(checks)}),flush=True)
