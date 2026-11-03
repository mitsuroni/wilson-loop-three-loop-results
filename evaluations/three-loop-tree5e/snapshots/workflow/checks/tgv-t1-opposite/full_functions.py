"""One-dimensional full-epsilon numerical reference and explicitly integrated series."""
import mpmath as mp
def density(e,w,z):
    A=1-w;B=z*w
    if A==B:return (1+z)*A**(e-1)
    return (1+z)*A**e*mp.expm1(e*mp.log(B/A))/(e*(B-A))
def N(e,w):
    Be=mp.beta(1-e,1-e)
    return Be*(mp.expm1(2*e*mp.log(w))+2*e*(1-w)*mp.hyp2f1(1-2*e,1-e,2-2*e,1-w))
def moments(e,z):
    e=mp.mpf(e);z=mp.mpf(z);a=1/(1+z)
    points=sorted(set([mp.mpf(0),mp.mpf(".25"),a,mp.mpf(".75"),mp.mpf(1)]))
    T=mp.quad(lambda w:density(e,w,z)*(-mp.expm1(2*e*mp.log(w))),points)
    K=mp.quad(lambda w:density(e,w,z)*N(e,w),points)
    return T,K
def all_values(ep,zval,dps=55):
    mp.mp.dps=dps;e=mp.mpf(ep);z=mp.mpf(zval)
    T,K=moments(e,z);Ti,Ki=moments(e,1/z)
    U=2*mp.gamma(1+e)*mp.gamma(1-e)-mp.beta(1-e,1-e)
    G=mp.exp(-3*mp.euler*e)*mp.gamma(1-e)*mp.gamma(1-2*e);lam=1-2*e
    J32=-z**(2*e)*U*T/(16*lam*e**4)
    J12=-z**(2*e)*K/(16*lam*e**4)
    I32=-G*U/(128*e**4)*(z**(2*e)*T+z**e*Ti)/(1+z**(3*e))
    I12=-G/(128*e**4)*(z**(2*e)*K+z**e*Ki)/(1+z**(3*e))
    return {"T":T,"K":K,"J":{"11":J32-J12,"12":J12,"32":J32},
      "I":{"11":I32-I12,"12":I12,"32":I32}}
def monomial_moments(e,j):
    j=mp.mpf(j);Be=mp.beta(1-e,1-e)
    R=1/(j+1)-1/(j+1+2*e)
    Nj=Be*(-R+2*e/((j+1)*(j+2))*mp.hyper([1-2*e,1-e,2],[2-2*e,j+3],1))
    return R,Nj
def integrated_series(e,z,Nmax=32):
    e=mp.mpf(e);z=mp.mpf(z);L=1+z
    # Polynomial coefficient arrays for a(w)=(z+w)/L, b(w)=1-z*w/L.
    def mul(a,b):
        c=[mp.mpf(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        return c
    ap=[mp.mpf(1)];bp=[mp.mpf(1)];symmetric=[mp.mpf(1)]
    momentsj=[monomial_moments(e,j) for j in range(Nmax+1)]
    T=K=mp.mpf(0);poch=mp.mpf(1);snapshots={}
    # h_n(a,b)=a*h_(n-1)+b^n, S_n=h_n/(n+1).
    for n in range(Nmax+1):
        if n:
            bp=mul(bp,[1,-z/L])
            symmetric=mul(symmetric,[z/L,1/L])
            symmetric=[x+y for x,y in zip(symmetric,bp)]
            poch*= (n-e)/n
        weights=[c/(n+1) for c in symmetric]
        T+=L**e*poch*mp.fsum(c*momentsj[j][0] for j,c in enumerate(weights))
        K+=L**e*poch*mp.fsum(c*momentsj[j][1] for j,c in enumerate(weights))
        if n in (8,16,32,64):snapshots[n]={"T":str(T),"K":str(K)}
    cnext=poch*(Nmax+1-e)/(Nmax+1)
    bound=L**e*cnext*(mp.log(z)**2+mp.pi**2)/2
    U=2*mp.gamma(1+e)*mp.gamma(1-e)-mp.beta(1-e,1-e)
    return {"T":T,"K":K,"T_tail_bound":bound,"K_tail_bound":abs(U)*bound,"partials":snapshots}
