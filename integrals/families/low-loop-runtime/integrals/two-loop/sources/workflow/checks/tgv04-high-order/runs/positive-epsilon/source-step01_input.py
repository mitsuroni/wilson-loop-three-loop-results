"""Exact input certificates for 04A. Does not import any claimed Laurent result."""
import json
from pathlib import Path
import sympy as s

x1,x2,x3,a1,a2,a3,z,u,v = s.symbols('x1 x2 x3 a1 a2 a3 z u v')
xs=(x1,x2,x3)
D=lambda x,y:(1-x)*(1-y)+z*x*y
S=a1*a2*x2*(1-x1)+a1*a3*D(x1,x3)+z*a2*a3*x3*(1-x2)
P=(a1*a2*(-(1+z)*(1-x1)+z*x2)+a1*a3*(1+3*z-(1+z)*x1+z*(1+z)*x3)+z*a2*a3*(1-x2-(1+z)*x3))/4
PR=((1+z)*a1*a2*(x1-1)+z*a1*a3*(1+(1+z)*x3)+z*a2*a3*(1-x2))/2
faces={
 '21_0':a1*a3*(1-x3)+x2*a1*a2+z*x3*a2*a3*(1-x2),
 '21_1':z*x3*a3*(a1+a2*(1-x2)),
 '31_0':a1*a3*D(x1,x3)+z*x3*a2*a3,
 '31_1':a1*(a2*(1-x1)+a3*D(x1,x3)),
 '41_0':a1*(1-x1)*(a3+x2*a2),
 '41_1':x2*a1*a2*(1-x1)+z*x1*a1*a3+z*a2*a3*(1-x2),
}
checks={}
def check(name,expr):
    r=s.factor(expr)
    checks[name]=str(r)
    print(name, '=', r, flush=True)
    assert r==0, (name,r)
check('P_boundary_plus_recoil',P+(z*s.diff(S,x1)-(1+z)*s.diff(S,x2)+s.diff(S,x3))/4-PR)
for i,key in enumerate(('21','31','41')):
    for b in (0,1): check('face_'+key+'_'+str(b),S.subs(xs[i],b)-faces[key+'_'+str(b)])
mirror={x1:1-x3,x2:1-x2,x3:1-x1,a1:a3,a3:a1,z:1/z}
check('mirror_S',S.subs(mirror,simultaneous=True)-S/z)
check('mirror_P',P.subs(mirror,simultaneous=True)-P/z**2)
par1={a1:1-u,a2:u*(1-v),a3:u*v}
par3={a3:1-u,a1:u*(1-v),a2:u*v}
check('31_1_reduction',faces['31_1'].subs(par1,simultaneous=True)-u*(1-u)*D(x1,v*x3))
check('31_0_reduction',faces['31_0'].subs(par3,simultaneous=True)-u*(1-u)*D(v+(1-v)*x1,x3))
check('21_1_cusp',faces['21_1'].subs(par3,simultaneous=True)-u*(1-u)*z*x3*(1-v*x2))
check('41_0_cusp',faces['41_0'].subs(par1,simultaneous=True)-u*(1-u)*(1-x1)*(v+(1-v)*x2))
# A symmetry of the two hard boundary integrals, not of each integrand at fixed z.
face_mirror={x1:1-x3,x2:1-x2,a1:a3,a3:a1,z:1/z}
check('hard_face_pair',faces['41_1'].subs(face_mirror,simultaneous=True)-faces['21_0']/z)
out={'checks':checks,'S':str(s.expand(S)),'P':str(s.expand(P)),'PR':str(s.expand(PR)),
     'recoil_x2_A':str(s.factor(S.subs(x2,0))),
     'recoil_x2_B':str(s.factor(s.diff(S,x2))),
     'recoil_x2_p':str(s.factor(PR.subs(x2,0))),
     'recoil_x2_q':str(s.factor(s.diff(PR,x2)))}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('ALL CERTIFICATES PASSED')
