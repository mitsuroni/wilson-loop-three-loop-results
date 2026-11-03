"""Intrinsic D candidate with no external rho power or Laurent convolution.

The physical scale is inside U and the universal dimensionless rectangle B.
Independent original-D comparison functions are not used in this module.
"""
from pathlib import Path
import sys
import math
import numpy as np
from scipy.special import spence

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT/'workflow/checks/pure-three-higher-laurent-20261016'))
from numeric_physical import grid
from numeric_shared_physical_generator import j


def sinhc(x):
    x = np.asarray(x)
    return np.divide(np.sinh(x), x, out=np.ones_like(x), where=x != 0)


def geometry(t, s, z):
    h = (1-t)*(1-s)+z*t*s
    H = h/np.sqrt(z)
    alpha = 1/(1+z)
    U = np.log(t*(1-t)*s*(1-s))/2+np.log(H)
    # The two differences vanish separately on s=alpha and t=alpha.
    ds = s-alpha
    dt = t-alpha
    rs = (1+z)*t*ds/(1-s)
    rt = (1+z)*(1-s)*dt/h
    ws = np.where(abs(rs) < .3, np.log1p(np.maximum(rs, -.5)), np.log(h)-np.log(1-s))
    wt = np.where(abs(rt) < .3, np.log1p(np.maximum(rt, -.5)), np.log(z*t)-np.log(h))
    shape = np.broadcast_shapes(np.shape(ws), np.shape(wt))
    qs = np.divide(ws, ds, out=np.broadcast_to((1+z)*t/(1-s),shape).copy(),
                   where=np.broadcast_to(abs(ds) > 1e-12, shape))
    qt = np.divide(wt, dt, out=np.broadcast_to((1+z)*(1-s)/h,shape).copy(),
                   where=np.broadcast_to(abs(dt) > 1e-12, shape))
    eta = np.log((1-t)/(z*s))/2
    kappa = (1+z)*np.sqrt((1-t)*s/z)
    return U, ws, wt, qs*qt, eta, kappa


def rectangle_j(e, eta, kappa):
    A = np.exp(eta)
    B = np.exp(-eta)
    p = 1-kappa/A
    q = 1-kappa/B
    if e == 0:
        return spence(np.asarray(1-p,dtype=float))+spence(np.asarray(1-q,dtype=float))-spence(np.asarray(1-p*q,dtype=float))-np.pi**2/6
    return -kappa**(-e)/e*(j(e,p)+j(e,q)-j(e,p*q)-float(j(e,np.array([1.]))[0]))


def rectangle_interval(eta, kappa, y):
    A = np.exp(eta)
    B = np.exp(-eta)
    slope = A-kappa*y
    h0 = B*y
    h1 = A+(B-kappa)*y
    r = slope/h0
    d = np.where(abs(r) < .3, np.log1p(np.maximum(r,-.5)), np.log(h1)-np.log(h0))
    value = np.divide(d, slope, out=np.broadcast_to(1/h0, np.shape(d)).copy(), where=abs(slope)>1e-12)
    center = (np.log(h1)+np.log(h0))/2
    return value, center, d


def full(e,z,n,wrong_scale=False):
    z=np.longdouble(str(z));x,w=grid(n)
    t=x[:,None];s=x[None,:]
    U,ws,wt,factor,eta,kappa=geometry(t,s,z)
    if wrong_scale:
        U = U-np.log((1-t)*s)/2
    B=rectangle_j(e,eta,kappa)
    result=w[:,None]*w[None,:]*factor*np.exp(e*U)*sinhc(e*ws/2)*sinhc(e*wt/2)*B
    return float(np.sum(result))


def coefficients(z,n,block=8,wrong_scale=False):
    """Extract D0..D4 analytically from the intrinsic source before quadrature."""
    z=np.longdouble(str(z));x,w=grid(n)
    s=x[None,:,None];y=x[None,None,:]
    total=np.zeros(5,dtype=np.longdouble)
    for lo in range(0,n,block):
        t=x[lo:lo+block,None,None]
        U,ws,wt,factor,eta,kappa=geometry(t,s,z)
        if wrong_scale:
            U = U-np.log((1-t)*s)/2
        value,center,d=rectangle_interval(eta,kappa,y)
        c=U+center
        v=(ws*ws+wt*wt+d*d)/24
        fourth=(ws**4+wt**4+d**4)/2880
        polynomials=(np.ones_like(c),c,c*c/2+v,c**3/6+c*v,
                     c**4/24+c*c*v/2+v*v/2-fourth)
        measure=w[lo:lo+block,None,None]*w[None,:,None]*w[None,None,:]
        base=-measure*factor*kappa*value
        for q,p in enumerate(polynomials):
            total[q]+=np.sum(base*p)
    return [float(v) for v in total]


def rectangle_controls(e,eta,kappa,n):
    y,w=grid(n)
    v,c,d=rectangle_interval(np.longdouble(str(eta)),np.longdouble(str(kappa)),y)
    integral=float(-kappa*np.sum(w*v*np.exp(e*c)*sinhc(e*d/2)))
    closed=float(rectangle_j(e,np.asarray(eta),np.asarray(kappa)))
    return {'quadrature':integral,'J_closed':closed,'absolute_difference':abs(integral-closed)}
