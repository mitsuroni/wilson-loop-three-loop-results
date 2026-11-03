"""Physical-z transport of the existing complete W moment source.

The 40 x-jets plus one constant form a nonminimal rational z connection.
Four zero x-leading rows cause no division: their z theta_x coefficient is
also exactly zero. A simple integer gauge makes the connection Fuchsian.
"""
from pathlib import Path
import importlib.util
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OLD=ROOT/'workflow/checks/pure-three-generators-20261015/moment_system.py'
spec=importlib.util.spec_from_file_location('_old_complete_moment_source_z',OLD)
mom=importlib.util.module_from_spec(spec);spec.loader.exec_module(mom)
e,z,n=s.symbols('e z n')

def transport():
    alpha=-(1+e)*z/(1+z);beta=(e-z)/(1+z);c=1/(1+z)
    D0=s.zeros(10)
    D0[0,0]=2*alpha;D0[0,3]=2*c
    D0[1,1]=alpha+beta;D0[1,4]=c;D0[1,5]=c
    D0[2,2]=2*beta;D0[2,6]=2*c
    D0[3,3]=alpha;D0[3,7]=c
    D0[4,4]=alpha;D0[4,8]=c
    D0[5,5]=beta;D0[5,8]=c
    D0[6,6]=beta;D0[6,9]=c
    D1=s.diag(-1,0,1,0,-1,1,0,1,0,-1)
    return D0,D1

def auxiliary_operators():
    S=s.Matrix(mom.transfer(n,e,z));S0=S.subs(n,0);S1=S.diff(n).subs(n,0);S2=S.diff(n,2)/2
    I=s.eye(10)
    P4=z*I+S2;P3=2*z*e*I+S1-e*S2;P2=z*e**2*I+S0-e*S1;P1=-e*S0
    forcing=s.Matrix(mom.initial(e,z))*e*(1+e)**2
    return [M.applyfunc(s.cancel) for M in (P1,P2,P3,P4)],forcing.applyfunc(s.cancel)

def theta_connection():
    D0,D1=transport();(P1,P2,P3,P4),r=auxiliary_operators();M=s.zeros(41)
    for k in range(4):
        for i in range(10):
            for j in range(10):M[10*k+i,10*k+j]=D0[i,j]
            if k<3:M[10*k+i,10*(k+1)+i]=D1[i,i]
    for i in range(10):
        if D1[i,i]==0:
            assert s.cancel(P4[i,i])==0
            continue
        assert s.cancel(P4[i,i])!=0
        factor=D1[i,i]/P4[i,i]
        for j in range(10):
            M[30+i,10+j]-=factor*P1[i,j]
            M[30+i,20+j]-=factor*P2[i,j]
            M[30+i,30+j]-=factor*P3[i,j]
        M[30+i,40]=factor*r[i]
    return M.applyfunc(s.cancel)

def fuchsian_connection():
    M=theta_connection()
    # Scale uv, v^2 and vZ by z, at each of the four x-jet levels.
    powers=[1 if i in (1,2,5) else 0 for k in range(4) for i in range(10)]+[0]
    H=s.zeros(41)
    for i in range(41):
        for j in range(41):
            if M[i,j]!=0:H[i,j]=s.cancel(z**(powers[i]-powers[j])*M[i,j]/z)
        H[i,i]+=s.Rational(powers[i],1)/z
    A0=s.zeros(41);A1=s.zeros(41)
    for i in range(41):
        for j in range(41):
            if H[i,j]==0:continue
            A0[i,j]=s.limit(z*H[i,j],z,0)
            A1[i,j]=s.limit((z+1)*H[i,j],z,-1)
            assert s.cancel(H[i,j]-A0[i,j]/z-A1[i,j]/(z+1))==0,(i,j,H[i,j])
    return H.applyfunc(s.cancel),A0.applyfunc(s.cancel),A1.applyfunc(s.cancel),powers
