"""Independent-endpoint Gauss solution and direct original-kernel checks."""
from functools import lru_cache
import math
import mpmath as mp
import numpy as np
from scipy.special import roots_legendre


def components(lam, a, b, z):
    lam, a, b, z = map(mp.mpf, (lam, a, b, z))
    phi = mp.hyp2f1(-a-lam, -b-lam, 1-lam, -z)
    xi = z**lam * mp.hyp2f1(-a, -b, 1+lam, -z)
    q = -(a+lam)*(b+lam)*z/(1-lam)*mp.hyper(
        [1, 1, 1-a-lam, 1-b-lam], [2, 2, 2-lam], -z)
    v = z**lam/lam * mp.hyper(
        [-a, -b, lam], [1+lam, 1+lam], -z)
    cap = mp.beta(1+a, lam)*mp.beta(1+b, lam)
    d = (mp.euler-mp.digamma(1-lam)
         +mp.digamma(1+a+lam)+mp.digamma(1+b+lam))
    return phi, xi, q, v, cap, d


def moment(lam, a, b, z):
    """F = integral t^a s^b D_z^(lam-1), generic nonresonant lam."""
    lam, a, b, z = map(mp.mpf, (lam, a, b, z))
    phi, xi, q, v, cap, d = components(lam, a, b, z)
    return (cap*phi+(xi*(mp.log(z)+d+q)-phi*v)/lam)/(1+z)**(1+a+b+lam)


def endpoint(lam, a, b, z):
    """The analytic finite difference, with its original normalization."""
    lam, a, b, z = map(mp.mpf, (lam, a, b, z))
    return -(1+z)*z**(-lam/2)/(a*b)*(
        moment(lam,a,b,z)-moment(lam,a,0,z)
        -moment(lam,0,b,z)+moment(lam,0,0,z))


def centered(lam, a, b, z):
    lam,a,b,z=map(mp.mpf,(lam,a,b,z))
    c=(a+b)/2
    return (z**c*endpoint(lam,a,b,z)+z**(-c)*endpoint(lam,a,b,1/z))/2


def physical(e, z):
    e,z=map(mp.mpf,(e,z))
    Z=mp.exp(-mp.euler*e)*mp.gamma(1-e)
    return -Z**3/(32*e**2*mp.cosh(3*e*mp.log(z)/2))*centered(e,e,e,z)


@lru_cache(None)
def grid(n):
    # The endpoint map used in numeric_physical.py, copied without imports.
    x,w=roots_legendre(n)
    x=(x.astype(np.longdouble)+1)/2
    d=x**3+(1-x)**3
    return x**3/d,w.astype(np.longdouble)/2*3*x*x*(1-x)**2/d**2


def original(lam,a,b,z,n):
    """Direct two-dimensional integration; no special-function solution."""
    lam,a,b,z=map(lambda x:np.longdouble(str(x)),(lam,a,b,z))
    x,w=grid(n);t=x[:,None];s=x[None,:]
    h=(1-t)*(1-s)+z*t*s
    weight=w[:,None]*w[None,:]*h**(lam-1)
    ua=-np.expm1(a*np.log(t))/a if a else -np.log(t)
    ub=-np.expm1(b*np.log(s))/b if b else -np.log(s)
    return {
        "F":float(np.sum(weight*t**a*s**b)),
        "E":float(-(1+z)*z**(-lam/2)*np.sum(weight*ua*ub)),
    }


def indexed(p,m,n,z,nodes):
    """Derivative-normalized centered moment, before Gamma and cosh."""
    z=np.longdouble(str(z));ell=np.log(z)
    x,w=grid(nodes);t=x[:,None];s=x[None,:]
    h=(1-t)*(1-s)+z*t*s;q=h/np.sqrt(z)
    def pm(v,ell,j):
        # Expanded difference avoids cancellation at v near one.
        L=np.log(v);c=ell/2
        return sum(math.comb(j+1,k)*c**(j+1-k)*L**k
                   for k in range(1,j+2))/(j+1)
    kernel=-(1+z)/h*np.log(q)**p/2
    pair=pm(t,ell,m)*pm(s,ell,n)+pm(1-t,-ell,m)*pm(1-s,-ell,n)
    return float(np.sum(w[:,None]*w[None,:]*kernel*pair))


def diagonal_coefficients(z,nodes,order=4):
    return [
        sum(indexed(p,m,q-p-m,z,nodes)/
            (math.factorial(p)*math.factorial(m)*math.factorial(q-p-m))
            for p in range(q+1) for m in range(q-p+1))
        for q in range(order+1)]
