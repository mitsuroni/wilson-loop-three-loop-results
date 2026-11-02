"""Read optional author ancillary files as exact data; never execute Wolfram code.

The original downloads remain outside this repository. Their hashes and
attribution are in reference/provenance.json. Default replay uses the derived
mathematical fixture; --reference-dir repeats this translation from originals.
"""
from pathlib import Path
import hashlib, json, re
import sympy as s
from sympy.parsing.mathematica import parse_mathematica

HERE=Path(__file__).resolve().parent

def wl(text):
    return parse_mathematica(text.replace('\\[Mu]','mu')).replace(s.Function('Zeta'),s.zeta).replace(s.Function('PolyLog'),s.polylog)

def extract(directory):
    directory=Path(directory)
    metadata=json.loads((HERE/'reference/provenance.json').read_text(encoding='utf-8-sig'))
    for name in ['WRen.m','ADs.m']:
        assert hashlib.sha256((directory/name).read_bytes()).hexdigest()==metadata['files'][name]['sha256'],name
    checks=[]
    def zero(name,value):
        residual=s.cancel(s.expand(value))
        assert residual==0,(name,residual)
        checks.append({'name':name,'residual':'0'})
    r,L,B,CA,CF,nf,TR,a,T=s.symbols('rho L B CA CF nf TR as4pi T')
    ss,tt,mu=s.symbols('s t mu');z,U=s.symbols('z U',positive=True)
    Li2,Li3,Q=s.symbols('Li2 Li3 Q')
    x=wl((directory/'WRen.m').read_text(encoding='utf-8-sig')).subs({ss:s.sqrt(z),tt:1,mu:s.sqrt(U),TR:s.Rational(1,2)})
    logs={z:r,1/z:-r,1+z:B,1+1/z:B-r,-z:r+s.I*s.pi,-1/z:-r+s.I*s.pi,
        z*U:L+r/2,U:L-r/2,(z+1)/(z**2*U):B-L-3*r/2,(z+1)/U:B-L+r/2,
        z*U/(z+1):L+r/2-B,U/(z+1):L-r/2-B,(z+1)*U:B+L-r/2}
    for f in sorted(x.atoms(s.log),key=str):
        values=[v for k,v in logs.items() if s.cancel(f.args[0]-k)==0]
        assert len(values)==1,('Unmapped author log',f)
        x=x.xreplace({f:values[0]})
    sum3=-Li3-r**2*B/2+r**3/6-s.pi**2*r/3+s.pi**2*B/2+s.I*s.pi*(r*B-B**2-r**2/2)+s.zeta(3)
    polys={(2,-z):Li2,(2,-1/z):-Li2-s.pi**2/6-r**2/2,
        (2,z/(1+z)):-Li2-B**2/2,(2,1/(1+z)):s.pi**2/6+Li2+r*B-B**2/2,
        (3,-z):Li3,(3,-1/z):Li3+(r**3+s.pi**2*r)/6,(3,1+z):Q,(3,1+1/z):sum3-Q}
    for f in sorted(x.atoms(s.polylog),key=str):
        n,arg=f.args;values=[v for (nn,k),v in polys.items() if n==nn and s.cancel(arg-k)==0]
        assert len(values)==1,('Unmapped author polylog',f)
        x=x.xreplace({f:values[0]})
    def dr(f):return s.diff(f,r)+z*s.diff(f,z)+z/(1+z)*s.diff(f,B)-B*s.diff(f,Li2)+Li2*s.diff(f,Li3)
    zero('dilog inversion derivative',dr(polys[(2,-1/z)])-(B-r))
    zero('fraction dilog derivative',dr(polys[(2,z/(1+z))])-B/(1+z))
    zero('complement dilog derivative',dr(polys[(2,1/(1+z))])-z/(1+z)*(r-B))
    zero('trilog inversion derivative',dr(polys[(3,-1/z)])+polys[(2,-1/z)])
    dsum=(z*(s.pi**2/6-Li2-B*(r+s.I*s.pi))-(s.pi**2/6-polys[(2,-1/z)]-(B-r)*(-r+s.I*s.pi)))/(1+z)
    zero('above-cut trilog pair derivative',dr(sum3)-dsum)
    x=s.expand(x);W1=s.expand(x.coeff(a,1)/4);W2=s.expand(x.coeff(a,2)/16)
    F2=s.expand(W2-W1**2/2)
    zero('ordinary CF squared cancels in log W',F2.coeff(CF,2))
    assert not F2.has(s.I,Q)
    F1=s.expand(W1/CF).subs(nf,2*T)
    F2=s.expand(F2/CF).subs(nf,2*T)
    WR2=s.expand(W2/CF).subs(nf,2*T)
    ad=(directory/'ADs.m').read_text(encoding='utf-8-sig')
    assert ad.count('\\[Gamma]Thr=')==1
    body=ad.split('\\[Gamma]Thr=',1)[1].split('/.as->4*as4pi;',1)[0]
    h2=s.expand(wl(re.sub(r'\bas\b','coupling',body))).coeff(s.Symbol('coupling'),2)/CF
    h2=s.expand(h2.subs(nf,2*T))
    b0=(11*CA-4*T)/12
    d=s.expand(s.diff(F2,L)-b0*F1)
    k2=s.expand(-s.diff(d,L)/2)
    h2_finite=s.expand(-d.subs(L,0))
    zero('external finite RG is local',d+2*k2*L+h2_finite)
    zero('noncusp finite versus independent threshold AD',h2_finite-h2)
    zero('cusp finite versus adopted two-loop cusp',k2-(CA*(s.Rational(67,36)-s.pi**2/12)-5*T/9))
    assert k2.free_symbols<= {CA,T} and h2.free_symbols<={CA,T}
    return {'format':'wilson-derived-author-reference/1','quantity':'lambda=alpha_s/pi; all values divide out external CF; F2 is log W_R, WR2 is ordinary W_R','domain':'positive z and the inherited principal-polylog boundary used by the preserved comparison','symbol_dictionary':{'rho':'log(z)','B':'log(1+z)','Li2':'polylog(2,-z)','Li3':'polylog(3,-z)','T':'TF*nf, author TR=1/2'},'F1_over_CF':str(F1),'F2_over_CF':str(F2),'WR2_over_CF':str(WR2),'k2_over_CF':str(k2),'h2_over_CF':str(h2),'checks':checks,'source_provenance':'reference/provenance.json','external_Wolfram_code_executed':False,'formal':'not_run'}
