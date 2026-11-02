"""Generate report-local pure-exchange coefficients; preserve all old evidence."""
from pathlib import Path
import hashlib, importlib.util, json, math, sys
import sympy as s
from sympy.parsing.mathematica import parse_mathematica
sys.dont_write_bytecode = True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=ROOT/'benchmarks/renormalization/two-loop/report'
checks=[]; sources={}
def bind(path):
    p=ROOT/path if not isinstance(path,Path) else path
    sources[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p
def read(path): return json.loads(bind(path).read_text(encoding='utf8'))
def check(name,ok):
    checks.append({'name':name,'passed':bool(ok)})
    assert ok,name
def zero(name,value): check(name,s.expand(value)==0 or s.simplify(value)==0)
e,l=s.symbols('eps ell')
def M(k,j): return s.Symbol(f'M{k}_{j}')
def N(k,j): return s.Symbol(f'N{k}_{j}')
def ze(k): return s.Symbol(f'zeta{k}')
h=[sum(M(q-2*j,2*j+1)/s.Integer(2**q*math.factorial(q-2*j)*math.factorial(2*j+1)) for j in range(q//2+1)) for q in range(7)]
nv=[sum(N(q-2*j,2*j+2)/s.Integer(math.factorial(q-2*j)*math.factorial(2*j+2)) for j in range(q//2+1)) for q in range(3)]
zp=[s.Integer(1)]; pp=[s.Integer(1)]; rr=[s.Integer(1)]
for r in range(1,7):
    zp.append(s.expand(s.Rational(2,r)*sum(ze(k)*zp[r-k] for k in range(2,r+1))))
    pp.append(s.expand(zp[r]-sum(l**(2*k)/s.factorial(2*k)*pp[r-2*k] for k in range(1,r//2+1))))
    rr.append(s.expand(sum((-1)**k*(2-2**k)*ze(k)*rr[r-k] for k in range(2,r+1))/r))
bb=[s.Integer(1 if r==0 else 0)-rr[r] for r in range(7)]
d=[s.expand((2**q-1)*h[q]-sum((l/2)**(2*k)/s.factorial(2*k)*h[q-2*k] for k in range(1,q//2+1))) for q in range(5)]
qc=[s.expand(sum(h[i]*h[q-i] for i in range(q+1))/2-2*nv[q]+2**(q+2)*sum((-1)**k*ze(k+2)*h[q-k] for k in range(q+1))) for q in range(3)]
def get(a,p): return a[p] if 0<=p<len(a) else s.Integer(0)
def conv(a,b,p): return sum(get(a,i)*get(b,p-i) for i in range(p+1)) if p>=0 else s.Integer(0)
def coeff(g,p):
    if g=='4e-01': return -get(zp,p+4)/16
    if g in ('4e-02','4e-05'): return conv(pp,d,p+2)/4
    if g=='4e-03': return -conv(pp,bb,p+4)/4
    if g=='4e-04': return -conv(pp,qc,p)/8
    if g=='4e-06': return -sum(get(pp,r)*get(h,i)*get(h,p-r-i) for r in range(p+1) for i in range(p-r+1))/16 if p>=0 else s.Integer(0)
    raise ValueError(g)
def normalize(expr):
    expr=s.sympify(expr)
    for atom in list(expr.atoms(s.Function)):
        name=str(atom.func)
        if name in ('M','WilsonM'): expr=expr.subs(atom,M(*map(int,atom.args[:2])))
        elif name in ('N','WilsonN'): expr=expr.subs(atom,N(*map(int,atom.args[:2])))
        elif name in ('Zeta','zeta'): expr=expr.subs(atom,ze(int(atom.args[0])))
    expr=expr.subs({s.Symbol('Pi'):s.pi,s.log(s.Symbol('z')):l})
    expr=s.expand(expr).xreplace({s.pi**6:216*ze(2)**3,s.pi**4:36*ze(2)**2,s.pi**2:6*ze(2)})
    expr=expr.subs({ze(4):s.Rational(2,5)*ze(2)**2,ze(6):s.Rational(8,35)*ze(2)**3,M(0,3):-M(0,1)**2,N(0,2):M(0,1)**2/3+4*ze(2)*M(0,1)})
    return s.expand(sum(c*(-2*M(0,1)-6*ze(2))**(p[0]//2)*l**(p[0]%2) for p,c in s.Poly(s.expand(expr),l).terms()))
def parse_entry(entry):
    f=next((entry[k] for k in ('MN_alternative','M','compact') if k in entry),entry)
    return normalize(parse_mathematica(f['expression']) if f['language']=='WolframLanguage/InputForm' else s.sympify(f['expression']))
ids=[f'4e-{i:02}' for i in range(1,7)]
raw={g:{p:s.expand(coeff(g,p)) for p in range(-4,3)} for g in ids}
printed={g:{p:normalize(v) for p,v in row.items()} for g,row in raw.items()}
# Separate truncated-exponential construction, independent of recurrence code.
def trunc(expr,degree): return s.series(expr,e,0,degree+1).removeO().expand()
ZS=trunc(s.exp(sum(2*ze(k)*e**k/k for k in range(2,7))),6)
RS=trunc(s.exp(sum((-1)**k*(2-2**k)*ze(k)*e**k/k for k in range(2,7))),6)
PS=trunc(ZS/s.cosh(e*l),6)
HS=sum(h[q]*e**q for q in range(5))
DS=trunc(HS.subs(e,2*e)-s.cosh(e*l/2)*HS,4)
QS=trunc(HS**2/2-2*sum(nv[q]*e**q for q in range(3))+2*sum((-1)**k*2**(k+1)*ze(k+2)*e**k for k in range(3))*HS.subs(e,2*e),2)
series={'4e-01':-ZS/(16*e**4),'4e-02':trunc(PS*DS,4)/(4*e**2),'4e-03':trunc(PS*(RS-1),6)/(4*e**4),'4e-04':-trunc(PS*QS,2)/8,'4e-06':-trunc(PS*HS**2,2)/16}
series['4e-05']=series['4e-02']
for g in ids:
    for p in range(-4,3): zero(f'independent series {g} eps^{p}',raw[g][p]-s.expand(series[g]).coeff(e,p))
package=read('benchmarks/2ne/prepared/4e.json')
results={}
z=s.Symbol('z1')
s11,s21,s22,s31,s32,s41,s42=s.symbols('s11 s21 s22 s31 s32 s41 s42')
A=lambda u,v:(1-u)*v
D=lambda u,v:(1-u)*(1-v)+z*u*v
DB=lambda u,v:z*(1-u)*(1-v)+u*v
kernel_rows=[(z**2,z*A(s31,s41),z*A(s32,s42)),(-z*(1+z),D(s21,s41),z*A(s31,s42)),(z,A(s21,s32),z*A(s31,s41)),((1+z)**2,D(s21,s41),D(s22,s42)),(-(1+z),A(s21,s31),D(s22,s41)),((1+z)**2,DB(s11,s31),D(s21,s41))]
orders=[['s32','s42','s31','s41'],['s21','s31','s42','s41'],['s21','s32','s41','s31'],['s22','s42','s21','s41'],['s31','s22','s41','s21'],['s11','s31','s21','s41']]
upper=[{'s31':'s32','s41':'s42'},{'s41':'s42'},{'s31':'s32'},{'s21':'s22','s41':'s42'},{'s21':'s22'},{}]
for index,(g,orbit) in enumerate(zip(ids,package['orbits'])):
    jj=orbit['J']
    for key,expected in zip(['numerator','f1','f2'],kernel_rows[index]):
        text=jj[key] if key=='numerator' else jj['propagator_polynomials'][key]
        zero(f'printed original polynomial {g} {key}',parse_mathematica(text)-expected)
    check(f'printed original order {g}',jj['integration_order_outer_to_inner']==orders[index])
    check(f'printed original bounds {g}',jj['bounds']=={v:['0',upper[index].get(v,'1')] for v in orders[index]})
    check(f'common prepared C {g}',orbit['C']['combined']==package['orbits'][0]['C']['combined'])
for g,orbit in zip(ids,package['orbits']):
    check(f'ordered source ID {g}',orbit['orbit_id']==g)
    revision='r001' if g=='4e-05' else 'r002'
    record=read(f'benchmarks/2ne/results/{g}/{revision}.json'); results[g]=record
    norm=record['normalization']
    check(f'C4 multiplicity {g}',orbit['c4_symmetry']['orbit_multiplicity']==norm['sym'])
    zero(f'kappa {g}',s.Rational(norm['kappa'])-s.Rational(norm['sym'])*s.Rational(norm['a1'])/8)
    for p in range(-4,1): zero(f'old result {g} eps^{p}',printed[g][p]-parse_entry(record['laurent']['coefficients'][str(p)]))
mapping=read('benchmarks/2ne/theory/generation/mother-map-r002.json')
G,R2,V,V2,VN,chi,Up,Um=s.symbols('G R2 VM VM2 VN chi Up Um')
H0,H1,H1bar,Qid=s.symbols('H0 H1 H1bar Qid')
Q=V**2/2-2*VN+2*chi/e*V2
full={'4e-01':-G**2/(16*e**4),'4e-02':G**2*(V2-s.cosh(e*l/2)*V)/(4*e**2*s.cosh(e*l)),'4e-03':G**2*(R2-1)/(4*e**4*s.cosh(e*l)),'4e-04':-G**2*Q/(8*s.cosh(e*l)),'4e-06':-G**2*V**2/(16*s.cosh(e*l))}
full['4e-05']=full['4e-02']
for item in mapping['items']:
    g=item['orbit_id']
    if g not in ids: continue
    original=s.sympify(item['I_sympy'],locals={'G':G,'eps':e,'ell':l})
    original=original.subs({H0:V,H1:s.exp(-e*l/2)*Up/e,H1bar:-s.exp(e*l/2)*Um/e,Qid:Q})
    zero(f'full epsilon bridge {g}',s.simplify((original-full[g]).subs(V2,(Up-Um)/(2*e)).rewrite(s.exp)))
# Historical Hlog data are reused; no old Maple run is repeated.
spec=importlib.util.spec_from_file_location('pure_hlog8',bind('workflow/checks/4e04-high-order/hlog8.py'))
ha=importlib.util.module_from_spec(spec); spec.loader.exec_module(ha)
bind('workflow/checks/weight-seven/hlog7.py')
run=ROOT/read('workflow/checks/single-opposite/current-moments.json')['run']
moments={}
for k in range(5):
    for j in range(1,6-k): moments[M(k,j)]=ha.parse(bind(run/f'M_{k}_{j}.txt'))
def hlog(expr):
    expr=s.expand(expr)
    for atom in list(expr.free_symbols):
        name=str(atom)
        if name.startswith('N'):
            k,j=map(int,name[1:].split('_')); expr=expr.subs(atom,M(k,j+1)-l*M(k,j))
    atoms=moments|{l:ha.LOG}|{ze(k):ha.period(k) for k in range(2,7)}
    return ha.from_sympy(s.expand(expr),atoms)
high=read('benchmarks/2ne/theory/generation/4e04-high-order.json')
for q in range(3):
    rec=next(x for x in high['coefficients'] if x['epsilon_power']==q)
    expected=ha.lift(rec['I_hlog_terms'])
    check(f'historical 4e04 integrated Hlog eps^{q}',hlog(raw['4e-04'][q])==expected)
    alternative=normalize(s.sympify(rec['I_MN_expression'],locals={'M':s.Function('M'),'N':s.Function('N')}))
    check(f'historical 4e04 printed M/N eps^{q}',hlog(alternative)==expected)
    printed['4e-04'][q]=alternative
for g in ids:
    for p in range(-4,3): check(f'printed expression Hlog {g} eps^{p}',hlog(printed[g][p])==hlog(raw[g][p]))
check('negative: duplicate Euler changes 4e01 eps^-3',s.expand(-ZS*s.exp(-2*s.Symbol('gammaE')*e)).series(e,0,2).removeO().coeff(e,1)!=0)
check('negative: missing C4 changes 4e04 finite',hlog(printed['4e-04'][0]/2)!=hlog(printed['4e-04'][0]))
check('negative: missing cosh changes 4e03 finite',normalize(raw['4e-03'][0]-s.expand(trunc(ZS*(RS-1),6)/(4*e**4)).coeff(e,0))!=0)
for name in ['benchmarks/2ne/theory/generation/unified-language.md','benchmarks/2ne/theory/generation/4e04-high-order.md','benchmarks/2ne/theory/generation/4e04-direct-mn.md']: bind(name)
for g in ids: bind(f'benchmarks/2ne/workspaces/{g}/derivation.md')
latex_names={l:r'\ell'}|{M(k,j):rf'M_{{{k},{j}}}' for k in range(5) for j in range(1,6)}|{N(k,j):rf'N_{{{k},{j}}}' for k in range(5) for j in range(1,6)}|{ze(k):rf'\zeta_{{{k}}}' for k in range(2,7)}
def latex(x): return s.latex(x,symbol_names=latex_names)
def coeff_display(g):
    lines=[]
    for p,v in printed[g].items():
        if v==0: continue
        terms=s.Add.make_args(s.expand(v))
        chunks=[s.Add(*terms[i:i+3]) for i in range(0,len(terms),3)]
        lines.append(rf'a_{{\mathrm{{{g}}},{p}}} &= '+latex(chunks[0])+(r'\notag' if len(chunks)>1 else ''))
        for index,c in enumerate(chunks[1:],1):
            t=latex(c)
            lines.append(r'&\quad '+('+' if not t.startswith('-') else '')+t+(r'\notag' if index<len(chunks)-1 else ''))
    return '\\begin{align}\n'+' \\\\\n'.join(lines)+'\n\\end{align}\n'
intro=r'''% Generated by workflow/checks/two-loop-report-20261015/build_pure.py.
The six entries below include the colour projection and the complete cyclic
orbit. The Laurent convention is
\begin{equation}
I_g(\varepsilon,z)=\sum_{p=-4}^{2}a_{g,p}(z)\varepsilon^p+O(\varepsilon^3).
\end{equation}
Every coefficient in this range not displayed below is zero. All functions
are evaluated at the same fixed $z>0$; the remainder refers to the local
meromorphic expansion at $\varepsilon=0$.

\subsection{Montopar input and the original ordered domains}
For this family the prepared prefactor is common to all six representatives:
\begin{equation}
C_{4e}=-\frac{C_AC_F}{2}\,e^{-2\gamma_E\varepsilon}
(4\pi\alpha_s)^2
\frac{\Gamma(1-\varepsilon)^2}{(4\pi^{2-\varepsilon})^2}
\frac{(-2v_{14}+i0)^{2\varepsilon}}{4}.
\label{eq:pure-C}
\end{equation}
The Euler factor occurs once in this prepared $C$. With
\begin{equation}
A(u,v)=(1-u)v,\quad
D_z(u,v)=(1-u)(1-v)+zuv,\quad
\overline D_z(u,v)=z(1-u)(1-v)+uv,
\end{equation}
each original $J_g$ is the iterated integral of
$P_g f_{g,1}^{\varepsilon-1}f_{g,2}^{\varepsilon-1}$.
The following table specifies the entire integrand.
\begin{center}\small
\begin{tabular}{llll}
\toprule
ID & $P_g$ & $f_{g,1}$ & $f_{g,2}$\\
\midrule
4e-01 & $z^2$ & $zA(s_{31},s_{41})$ & $zA(s_{32},s_{42})$\\
4e-02 & $-z(1+z)$ & $D_z(s_{21},s_{41})$ & $zA(s_{31},s_{42})$\\
4e-03 & $z$ & $A(s_{21},s_{32})$ & $zA(s_{31},s_{41})$\\
4e-04 & $(1+z)^2$ & $D_z(s_{21},s_{41})$ & $D_z(s_{22},s_{42})$\\
4e-05 & $-(1+z)$ & $A(s_{21},s_{31})$ & $D_z(s_{22},s_{41})$\\
4e-06 & $(1+z)^2$ & $\overline D_z(s_{11},s_{31})$ & $D_z(s_{21},s_{41})$\\
\bottomrule
\end{tabular}
\end{center}
All lower limits are zero. Unless specified in the third column below,
all upper limits are one. The second column is the original outer-to-inner
integration order; dependent limits are part of the definition.
\begin{center}\small
\begin{tabular}{llll}
\toprule
ID & Outer $\to$ inner & Nontrivial upper limits & Cyclic multiplicity\\
\midrule
4e-01 & $s_{32},s_{42},s_{31},s_{41}$ & $s_{31}<s_{32},\ s_{41}<s_{42}$ & 4\\
4e-02 & $s_{21},s_{31},s_{42},s_{41}$ & $s_{41}<s_{42}$ & 4\\
4e-03 & $s_{21},s_{32},s_{41},s_{31}$ & $s_{31}<s_{32}$ & 4\\
4e-04 & $s_{22},s_{42},s_{21},s_{41}$ & $s_{21}<s_{22},\ s_{41}<s_{42}$ & 2\\
4e-05 & $s_{31},s_{22},s_{41},s_{21}$ & $s_{21}<s_{22}$ & 4\\
4e-06 & $s_{11},s_{31},s_{21},s_{41}$ & none & 1\\
\bottomrule
\end{tabular}
\end{center}
These are the original ordered-parameter roots, with $z=z_1$.
The fixed-basis colour coefficient is $a_1=-1/2$ for each representative.
After the cyclic sum and scale extraction the corresponding constants are
$\kappa_g=\mathrm{mult}_g a_1/[2(-2)^2]
=(-1/4,-1/4,-1/4,-1/8,-1/4,-1/16)$.
The $I_g$ below already contains this operation; no extra multiplicity or
Euler factor is to be applied when the six $I_g$ are added.

'''
parts=[intro]
sections={
'4e-01':(r'4e-01: two adjacent exchanges',r'''The ordered endpoint integrals give
\begin{equation}
I_{\mathrm{4e-01}}=-\frac{Z(\varepsilon)^2}{16\varepsilon^4}.
\end{equation}
There is no dependence on $z$ after the orbit normalization.
'''),
'4e-02':(r'4e-02 and 4e-05: one opposite-edge kernel',r'''The two labelled representatives are equal after the full cyclic
normalization, at unexpanded $\varepsilon$:
\begin{equation}
I_{\mathrm{4e-02}}=I_{\mathrm{4e-05}}
=\frac{Z(\varepsilon)^2}{4\varepsilon^2c_2}
\bigl[\mathcal V_M(2\varepsilon,z)-c_1\mathcal V_M(\varepsilon,z)\bigr].
\end{equation}
Only the kernel parameter is doubled in the first term: the external
factor remains $Z(\varepsilon)^2$. Since the constant term in the brackets
cancels, the leading pole is simple. The following coefficients apply to
both representatives, $a_{\mathrm{4e-05},p}=a_{\mathrm{4e-02},p}$.
'''),
'4e-03':(r'4e-03: the ordered adjacent-edge difference',r'''The two nested endpoints produce the beta ratio
$R_2(\varepsilon)=\Gamma(1+\varepsilon)^2/\Gamma(1+2\varepsilon)$:
\begin{equation}
I_{\mathrm{4e-03}}=\frac{Z(\varepsilon)^2[R_2(\varepsilon)-1]}
{4\varepsilon^4c_2}.
\end{equation}
The numerator begins at $\varepsilon^2$, leaving a double pole.
'''),
'4e-04':(r'4e-04: two ordered kernels in the same channel',r'''The evaluated same-channel kernel is
\begin{align}
I_{\mathrm{4e-04}}&=-\frac{Z(\varepsilon)^2}{8c_2}Q(\varepsilon,z),\\
Q&=\frac12\mathcal V_M(\varepsilon,z)^2-2\mathcal V_N(\varepsilon,z)
 +\frac{2\chi(\varepsilon)}{\varepsilon}\mathcal V_M(2\varepsilon,z),\\
\chi(\varepsilon)&=\gamma_E+\psi(1+2\varepsilon).
\end{align}
Here $\psi=\Gamma'/\Gamma$, and $\chi/\varepsilon$ is regular at zero.
This representative has no poles. The explicit $M/N$ coefficients below
also agree with the earlier direct integration through $\varepsilon^2$.
'''),
'4e-06':(r'4e-06: the product of two opposite-edge kernels',r'''The two channels factorize after the original domain is accounted for:
\begin{equation}
I_{\mathrm{4e-06}}=-\frac{Z(\varepsilon)^2}{16c_2}
\mathcal V_M(\varepsilon,z)^2.
\end{equation}
This representative is regular at $\varepsilon=0$.
''')}
for g,(title,note) in sections.items(): parts += ['\\subsection{'+title+'}\n',note,coeff_display(g)]
tex=''.join(parts)
check('six IDs retained in printed source',all(g in tex for g in ids))
check('no control characters in LaTeX',not any(ord(c)<32 and c not in '\r\n\t' for c in tex))
OUT.mkdir(parents=True,exist_ok=True)
(OUT/'pure-results.tex').write_text(tex,encoding='utf8',newline='\n')
report={'format':'two-loop-report-pure/1','domain':'fixed real z>0','powers':list(range(-4,3)),'remainder':'O(eps^3)','coefficients':{g:{str(p):{'sympy':str(v),'latex':latex(v)} for p,v in row.items()} for g,row in printed.items()},'normalization':{'Z':'exp(-gamma_E eps) Gamma(1-eps)','c_n':'cosh(n eps log(z)/2)','Euler':'one e^(-2 gamma_E eps) already included','C4':'already included; multiplicities 4,4,4,2,4,1','color':'CA CF; a1=-1/2'},'source_sha256':sources,'formal':'not_run'}
(HERE/'pure-data.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8',newline='\n')
bind(Path(__file__).resolve())
summary={'passed':True,'checks':checks,'counts':{'all_coefficient_positions':42,'old_coefficient_positions':30,'historical_high_order_positions':3,'negative_controls':3},'source_sha256':sources,'output_sha256':{'pure-results.tex':hashlib.sha256((OUT/'pure-results.tex').read_bytes()).hexdigest(),'pure-data.json':hashlib.sha256((HERE/'pure-data.json').read_bytes()).hexdigest()},'new_multidimensional_integrations':0,'formal':'not_run'}
(HERE/'pure-checks.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps({'passed':True,'checks':len(checks),'coefficient_positions':42,'old_positions':30,'historical_high_order_positions':3,'formal':'not_run'}))