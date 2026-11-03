"""Exact checks of the intrinsic endpoint/rectangle representation.

These are local CAS identities, not formal certificates.  Existing original
integrals are read-only.  Run once per immutable output revision.
"""
from pathlib import Path
import hashlib
import json
import math
import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    checks = []
    negatives = []

    def zero(name, expr):
        value = S.factor(S.cancel(S.expand(expr)))
        assert value == 0, (name, value)
        checks.append({'name': name, 'exact_zero': True})

    t, s, z, x, y = S.symbols('t s z x y', positive=True)
    h = (1-t)*(1-s)+z*t*s
    q = h/S.sqrt(z)
    scale = S.sqrt((1-t)*s)
    aa = S.sqrt((1-t)/(z*s))
    bb = 1/aa
    kap = (1+z)*scale/S.sqrt(z)
    # Work on the original real domain, 0<t,s<1. Squared identities avoid
    # asking SymPy to assume that arbitrary t,s are below one.
    inner = ((1-t)*x + z*s*y-(1+z)*(1-t)*s*x*y)/S.sqrt(z)
    zero('rectangle_original_change',
         (((1-(1-(1-t)*x))*(1-s*y)+z*(1-(1-t)*x)*s*y)/S.sqrt(z))-inner)
    zero('rectangle_normalized_A_squared', (scale*aa)**2-(1-t)**2/z)
    zero('rectangle_normalized_B_squared', (scale*bb)**2-z*s**2)
    zero('rectangle_normalized_cross_squared', (scale*kap)**2-(1+z)**2*(1-t)**2*s**2/z)
    zero('rectangle_positive_corner', h-((1-t)+z*s-(1+z)*(1-t)*s))
    zero('rectangle_J_first_argument', 1-(1+z)*s+(1+z)*s-1)
    ap = 1-(1+z)*s
    bp = 1-(1+z)*(1-t)/z
    zero('rectangle_J_product', ap*bp-(1-(1+z)*h/z))
    zero('two_line_denominator',
         (1+z)**2/(((1+z)*t-1)*((1+z)*s-1))-
         1/((t-1/(1+z))*(s-1/(1+z))))
    a0, b0 = S.sqrt(z)*t, (1-s)/S.sqrt(z)
    zero('opposite_endpoint_product_is_z_free', a0*b0-t*(1-s))
    zero('total_exponential_scale_squared',
         (t*(1-s))**2*((1-t)*s)*(q**2/(t*(1-s)))-
         t*(1-t)*s*(1-s)*q**2)
    for n in range(1,5):
        # a^e q^e + b^e q^e - q^(2e) -(ab)^e is the two-cosh form.
        zero(f'two_interval_product_integer_e_{n}',
             (a0**n-q**n)*(q**n-b0**n)-
             (a0**n*q**n+b0**n*q**n-q**(2*n)-(a0*b0)**n))
    ca = {z:1/z,t:1-s,s:1-t}
    def reflect(v):
        return S.powdenest(v.subs(ca, simultaneous=True), force=True)
    zero('crossing_q', reflect(q)-q)
    zero('crossing_weight_product', reflect(t*(1-t)*s*(1-s))-t*(1-t)*s*(1-s))
    zero('crossing_L_argument_inverse', reflect(z*t/(1-s))*z*t/(1-s)-1)
    zero('crossing_V_argument_squared', reflect(q**2/(t*(1-s)))-q**2/(t*(1-s)))
    zero('crossing_eta_argument_inverse', reflect((1-t)/(z*s))*(1-t)/(z*s)-1)
    zero('crossing_kappa_squared', reflect(kap**2)-kap**2)
    den = (t-1/(1+z))*(s-1/(1+z))
    zero('crossing_measure_denominator', reflect(den)-den)
    zero('removable_t_line_V_equals_L',
         (q**2/(t*(1-s))/(z*t/(1-s))-1).subs(t,1/(1+z)))
    zero('removable_s_line_V_equals_minus_L',
         (q**2/(t*(1-s))*z*t/(1-s)-1).subs(s,1/(1+z)))

    # Every ordinary epsilon coefficient of the universal rectangle is
    # a finite J coefficient convolution. Verify integer powers directly
    # against the original elementary square integral at rational shapes.
    def j_integer(n, w):
        return sum(S.binomial(n,k)*(-w)**k/S.Integer(k) for k in range(1,n+1))
    shape_rows = []
    for A,K in [(S.Rational(2), S.Rational(1,2)),
                (S.Rational(3,2), S.Rational(7,6)),
                (S.Rational(3,4), S.Rational(4,3))]:
        B=1/A
        assert K<A+B
        a,b=1-K/A,1-K/B
        for n in range(1,5):
            rhs=-K**(-n)*(j_integer(n,a)+j_integer(n,b)-j_integer(n,a*b)-j_integer(n,1))/n
            raw=-K*S.integrate(S.integrate(S.expand((A*x+B*y-K*x*y)**(n-1)),(x,0,1)),(y,0,1))
            zero(f'rectangle_J_integer_{A}_{K}_{n}', rhs-raw)
            shape_rows.append({'A':str(A),'kappa':str(K),'epsilon':n,'value':str(raw)})
        wrong=-K**(-1)*(j_integer(1,a)+j_integer(1,b)-j_integer(1,a*b))
        assert S.simplify(wrong+K) != 0
        negatives.append({'name':f'omit_rectangle_corner_{A}_{K}', 'rejected':True})

    # Independent Taylor multiplication before integrating. Moment indices
    # are raw derivatives, so factorials must remain in the coefficient rule.
    e,U,L,V=S.symbols('e U L V')
    bv=S.symbols('b0:5')
    exp_ser=sum((e*U)**m/S.factorial(m) for m in range(5))
    cos_ser=sum(2*e**(2*j-2)*(L**(2*j)-V**(2*j))/S.factorial(2*j) for j in range(1,4))
    rect_ser=sum(bv[n]*e**n/S.factorial(n) for n in range(5))
    lhs=S.expand(exp_ser*cos_ser*rect_ser)
    coefficient_terms={}
    for order in range(5):
        rhs=0
        rows=[]
        for j in range(1,(order+2)//2+1):
            rem=order+2-2*j
            for m in range(rem+1):
                n=rem-m
                coeff=S.Rational(2,math.factorial(m)*math.factorial(n)*math.factorial(2*j))
                rhs+=coeff*U**m*(L**(2*j)-V**(2*j))*bv[n]
                rows.append({'m':m,'n':n,'j':j,'coefficient':str(coeff),'weight':m+n+2*j+4})
        zero(f'coefficient_rule_D{order}',lhs.coeff(e,order)-rhs)
        coefficient_terms[str(order)]=rows
    assert S.expand(lhs.coeff(e,2)-(U**2*bv[0]+U*bv[1]+bv[2])*(L**2-V**2)) != 0
    negatives.append({'name':'omit_Taylor_factorials','rejected':True})
    wrong_cross=(a0*q+b0*q-q**2+a0*b0)
    assert S.factor(wrong_cross-(a0-q)*(q-b0)) != 0
    negatives.append({'name':'reverse_endpoint_product_constant','rejected':True})
    inputs=[Path(__file__),
            ROOT/'benchmarks/2ne/theory/generation/pure-three-cd-20261016/D/derivation.md',
            ROOT/'benchmarks/2ne/theory/generation/4e04-direct-mn.md']
    report={'passed':True,'formal':'not_run','exact_checks':checks,
            'negative_controls':negatives,'rectangle_integer_examples':shape_rows,
            'D_coefficient_terms':coefficient_terms,
            'assumptions':'z>0, 0<t,s<1; A>0, B=1/A, 0<kappa<A+B; real branches',
            'inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    out=HERE/'algebra-verification.json'
    with out.open('x',encoding='utf-8') as f:
        json.dump(report,f,ensure_ascii=False,indent=2)
    print(json.dumps({'passed':True,'exact_checks':len(checks),'negatives':len(negatives)}))


if __name__=='__main__':
    main()
