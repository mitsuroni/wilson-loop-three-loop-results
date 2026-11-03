"""Independent endpoint quadrature and integer polynomial checks for G closure."""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import time

import mpmath as mp

from moment_system import coefficients, initial, matvec, project, transfer


def polynomial_moment(n, eps, z, second=False):
    # Direct polynomial integration of p**(n-1)*(1+p)**eps or (1-p/z)**eps.
    return sum(F(math.comb(eps, j)) * ((-1/z)**j if second else 1)
               * (z**(n+j)-F(-1)**(n+j))/F(n+j)
               for j in range(eps+1))


def quad_moment(n, eps, z, second=False):
    return mp.quad(lambda p: p**(n-1)*(1-p/z if second else 1+p)**eps,
                   [-1, 0, z])


def closure_theta2(eps, z, x, count=140):
    return sum(n*n*project(c, eps, z)*x**n
               for n, c in enumerate(coefficients(eps, z, count), 1))


def direct_theta2(eps, z, x):
    def a(p):
        return (1+p)**eps-(1-p/z)**eps
    return eps*x/z*mp.quad(
        lambda p: a(p)*mp.quad(
            lambda q: a(q)*(1+x*p*q/z)**(eps-1), [-1, 0, z]),
        [-1, 0, z])


def run():
    started = time.time()
    records = []
    for z in (F(1,2),F(1),F(2)):
        rr=(0,0,0,1,1,1,1,2,2,2)
        dd=(1,-z,z*z,-z,1,z*z,-z,z*z,-z,1)
        for n in range(1,8):
            ss=transfer(n,F(0),z)
            for i in range(10):
                assert n*n*ss[i][i] == dd[i]*n**(4-rr[i])*(n+1)**rr[i]
        records.append({"type":"epsilon_zero_diagonal_factorization","z":str(z),
                        "monomial_comparisons":70,"status":"pass"})
    for e in (1, 2, 3, 4):
        eps = F(e)
        for z in (F(1, 2), F(1), F(2)):
            u = v = F(1)/(1+eps)
            aa = (1+z)**(eps+1)
            bb = aa*z**(-eps)
            y = initial(eps, z)
            fc = coefficients(eps, z, e+2)
            for n in range(1, 11):
                ud = polynomial_moment(n, e, z)
                vd = polynomial_moment(n, e, z, True)
                assert aa*u == ud and bb*v == vd
                expected = [u*u,u*v,v*v,u*z**n,u*(-1)**n,v*z**n,v*(-1)**n,
                            z**(2*n),(-z)**n,F(1)]
                assert y == expected
                if n <= e+2:
                    poch = math.prod(-eps+j for j in range(n))
                    direct = poch/F(math.factorial(n)*n)*(-1/z)**n*(ud-vd)**2
                    assert project(fc[n-1],eps,z) == direct
                d = n+eps+1
                y = [t/d**2 for t in matvec(transfer(n,eps,z), y)]
                u,v=(-n*u+z**n)/d,(z*n*v+(-1)**n)/d
            records.append({"type":"integer_polynomial", "eps":e,"z":str(z),
                            "moments_each":10,"status":"pass"})
    mp.mp.dps = 45
    for es,zs in (("0.23","2"),("-0.19","0.7"),("0.37","1")):
        eps,z=mp.mpf(es),mp.mpf(zs)
        aa=(1+z)**(eps+1)
        bb=aa*z**(-eps)
        u=v=1/(1+eps)
        errors=[]
        for n in range(1,8):
            if n in (1,2,4,7):
                ud=quad_moment(n,eps,z)
                vd=quad_moment(n,eps,z,True)
                errors.extend([abs(aa*u-ud)/max(1,abs(ud)),
                               abs(bb*v-vd)/max(1,abs(vd))])
            d=n+eps+1
            u,v=(-n*u+z**n)/d,(z*n*v+(-1)**n)/d
        err=max(errors)
        assert err < mp.mpf("1e-32")
        records.append({"type":"noninteger_moments","eps":es,"z":zs,
                        "comparisons":8,"max_scaled_error":mp.nstr(err,8),"status":"pass"})
    for es,zs,xs in (("0.23","2","0.1"),("0.37","0.7","0.2")):
        eps,z,x=map(mp.mpf,(es,zs,xs))
        direct=direct_theta2(eps,z,x)
        series=closure_theta2(eps,z,x)
        err=abs(direct-series)/max(1,abs(direct))
        assert err < mp.mpf("1e-32")
        records.append({"type":"full_theta2_G_quadrature","eps":es,"z":zs,"x":xs,
                        "direct":mp.nstr(direct,35),"series":mp.nstr(series,35),
                        "scaled_error":mp.nstr(err,8),"status":"pass"})
    # A wrong sign in the V boundary term changes a nonzero moment.
    eps,z=F(2),F(2)
    v1=1/(1+eps)
    wrong=(z*v1+1)/(eps+2)
    correct=polynomial_moment(2,2,z,True)/((1+z)**(eps+1)*z**(-eps))
    assert wrong != correct
    records.append({"type":"negative_control","mutation":"flip V recurrence endpoint sign at n=1",
                    "status":"rejected"})
    result={"status":"pass","formal":"not_run","digits":mp.mp.dps,
            "elapsed_seconds":round(time.time()-started,3),"records":records,
            "scope":"Full G moment closure. No Lambda closure, physical endpoint continuation, or full I claimed here."}
    out=Path(__file__).with_name("moment_report.json")
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!='records'}))


if __name__ == '__main__':
    run()
