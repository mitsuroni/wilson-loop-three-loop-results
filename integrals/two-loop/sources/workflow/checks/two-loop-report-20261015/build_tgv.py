"""Generate the TGV appendices from exact existing sources; no integral is rerun."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / 'benchmarks/renormalization/two-loop/report'
sys.path.insert(0, str(ROOT/'workflow/checks/single-opposite'))
import hlog_algebra as h

eps, ell = sp.symbols('eps ell')
Z = {n: sp.Symbol(f'zeta_{n}') for n in range(2, 7)}
M = lambda k,j: sp.Symbol(f'M_{k}_{j}')
N = lambda k,j: sp.Symbol(f'N_{k}_{j}')
sources = set()
checks = []

def read_json(rel):
    p = ROOT/rel
    sources.add(p)
    return json.loads(p.read_text(encoding='utf-8'))

def check(name, condition):
    checks.append({'name': name, 'passed': bool(condition)})
    assert condition, name

def eq(name, a, b=0):
    check(name, sp.expand(a-b)==0)

def exp_coeff(logs, n):
    coeff = [sp.Integer(1)]
    for q in range(1,n+1):
        coeff.append(sp.expand(sum(k*logs.get(k,0)*coeff[q-k] for k in range(1,q+1))/q))
    return coeff

def poly(cs):
    return sum(c*eps**q for q,c in enumerate(cs))

def norm(expr):
    return sp.expand(expr.subs({Z[4]:sp.Rational(2,5)*Z[2]**2,
        Z[6]:sp.Rational(8,35)*Z[2]**3}))

glog = {n:sp.Rational(2**n,n)*Z[n] for n in range(2,7)}
flog = {n:sp.Rational(2,n)*Z[n] for n in range(2,7,2)}
blog = {n:sp.Rational(2-2**n,n)*Z[n] for n in range(2,7)}
G = poly(exp_coeff(glog,6))
F = poly(exp_coeff(flog,6))
B = sp.Poly(sp.expand(poly(exp_coeff(blog,6))*sum((2*eps)**q for q in range(7))),eps)
B = sum(B.coeff_monomial(eps**q)*eps**q for q in range(7))
A01 = sp.expand(G*(2*F-B)/32)
I01 = {q:norm(A01.coeff(eps,q+4)) for q in range(-4,3)}

# Separate SymPy Taylor implementation checks the recurrence/product bookkeeping.
direct01 = sp.series(sp.exp(sum(v*eps**k for k,v in glog.items()))*
    (2*sp.exp(sum(v*eps**k for k,v in flog.items()))-
     sp.exp(sum(v*eps**k for k,v in blog.items()))/(1-2*eps)),eps,0,7).removeO()/32
for q in range(-4,3):
    eq(f'01 independent Taylor coefficient {q}', I01[q],norm(sp.expand(direct01).coeff(eps,q+4)))

VM = sum(M(k,2*j+1)*eps**(k+2*j)/sp.factorial(k)/sp.factorial(2*j+1)
         for k in range(5) for j in range(3) if k+2*j<=4)
VM0 = sum(M(0,2*j+1)*eps**(2*j)/sp.factorial(2*j+1) for j in range(3))
A03 = sp.expand(G*(1-ell**2*eps**2/2)*(VM0-B*VM)/8)
I03 = {q:norm(A03.coeff(eps,q+2).subs({ell**2:-2*M(0,1)-6*Z[2],M(0,3):-M(0,1)**2}))
       if q>=-2 else sp.Integer(0) for q in range(-4,3)}
m01,m11,m21,m31,m13,m41,m23 = [M(k,j) for k,j in [(0,1),(1,1),(2,1),(3,1),(1,3),(4,1),(2,3)]]
z2,z3,z5=Z[2],Z[3],Z[5]
pretty03 = {
 1:-m31/48-m13/48-m21/8-(m01+4*z2+4)*m11/8-sp.Rational(5,24)*m01**2+(z3/4-z2-1)*m01,
 2:-m41/192-m23/96-m31/24-m13/24-(m01/16+z2/4+sp.Rational(1,4))*m21-
   (m01/4+z2+z3/12+1)*m11+(sp.Rational(5,48)*z2-sp.Rational(5,12))*m01**2+
   (sp.Rational(59,80)*z2**2-2*z2-z3/6-2)*m01}
for q in (1,2):eq(f'03 displayed grouped coefficient {q}',I03[q],pretty03[q])

compact = read_json('workflow/checks/tgv04-high-order/compact.json')
I04 = {-4:sp.Integer(0),-3:sp.Integer(0),-2:z2/8,-1:(z3-m11)/8,
       0:m01**2/48-m21/8+sp.Rational(9,10)*z2**2,
       1:sp.sympify(compact['I1']),2:sp.sympify(compact['I2'])}
data = {'3e-3g-01':I01,'3e-3g-02':I01,'3e-3g-03':I03,'3e-3g-04':I04}
oldlocals={f'M{k}{j}':M(k,j) for k in range(5) for j in range(1,6)}
oldlocals.update({f'zeta_{n}':sp.zeta(n) for n in (2,3,5)})
oldlocals.update({f'zeta{n}':sp.zeta(n) for n in (2,3,5)})
for name, coeffs in data.items():
    old=read_json(f'benchmarks/tgv/3e-3g/results/{name}/r001.json')
    for q in range(-4,1):
        target=sp.sympify(old['laurent']['coefficients'][str(q)]['expression'],locals=oldlocals)
        actual=coeffs[q].subs({Z[n]:sp.zeta(n) for n in Z})
        eq(f'{name} preserved published coefficient {q}',actual,target)

# Reconstruct the compact high-order formulas in the complete stored Hlog algebra.
moment_index=read_json('workflow/checks/single-opposite/current-moments.json')
run=ROOT/moment_index['run']
tau_data=read_json('benchmarks/2ne/families/single-opposite/basis-functions.json')
tau_data=tau_data.get('functions',tau_data)['tauDq4']
tau={(tuple(t['word']),tuple(t['zeta_powers'])):Fraction(t['coefficient']) for t in tau_data['hlog_terms']}
atoms={str(Z[n]):h.period(n) for n in (2,3,5)}
atoms['tauD4']=tau
def atom(name):
    if name in atoms:return atoms[name]
    kind,k,j=name.split('_');k,j=int(k),int(j)
    if kind=='M':
        path=run/f'M_{k}_{j}.txt';sources.add(path);atoms[name]=h.parse(path)
    elif kind=='N':
        atoms[name]=h.add(atom(f'M_{k}_{j+1}'),h.scale(h.mul(h.LOG,atom(f'M_{k}_{j}')),-1))
    return atoms[name]
def evaluate(expr):
    parts=[]
    for term in sp.Add.make_args(sp.expand(expr)):
        c,f=term.as_coeff_Mul();v=h.scale(h.ONE,Fraction(c))
        for a,n in f.as_powers_dict().items():v=h.mul(v,h.power(atom(str(a)),int(n)))
        parts.append(v)
    return h.add(*parts)
for q in (1,2):
    path=ROOT/f'workflow/checks/tgv04-high-order/runs/positive-epsilon/I{q}.txt'
    sources.add(path)
    check(f'04 coefficient {q} exact full Hlog equality',not h.add(evaluate(I04[q]),h.scale(h.parse(path),-1)))
check('M03=-M01^2 exact Hlog identity',not h.add(atom('M_0_3'),h.power(atom('M_0_1'),2)))
check('ell squared reduction exact Hlog identity',not h.add(h.power(h.LOG,2),h.scale(atom('M_0_1'),2),h.scale(h.period(2),6)))
# An intentionally reversed tau sign must fail the complete coefficient comparison.
check('negative control rejects missing tau contribution',bool(h.add(evaluate(I04[2]-sp.Symbol('tauD4')/8),h.scale(h.parse(path),-1))))

# Match every displayed numerator, denominator, and ordered domain to current inputs.
prepared=read_json('benchmarks/tgv/3e-3g/prepared/3e-3g.json')
a,b,c,t,x,y,w,z=sp.symbols('a b c t x y w z')
D=lambda u,v:(1-u)*(1-v)+z*u*v
kernels=[
 (a*z*(1-t)*(b*x+c*y),a*z**2*(1-t)*(b-c)/4,{'s31':t,'s41':x,'s42':y},['s31','s42','s41']),
 (c*z*t*(a*(1-x)+b*(1-y)),c*z**2*t*(b-a)/4,{'s31':x,'s32':y,'s41':t},['s32','s41','s31']),
 (a*(b*D(t,x)+c*D(t,y)),(1+z)*a*(b-c)*(1-(1+z)*t)/4,{'s21':t,'s41':x,'s42':y},['s21','s42','s41']),
 (a*b*y*(1-x)+a*c*D(x,w)+z*b*c*w*(1-y),
 (a*b*(-(1+z)*(1-x)+z*y)+a*c*(1+3*z-(1+z)*x+z*(1+z)*w)+z*b*c*(1-y-(1+z)*w))/4,
 {'s21':x,'s31':y,'s41':w},['s21','s31','s41'])]
for row,(S,P,variables,order) in zip(prepared['orbits'],kernels):
    mapping={'alpha1':a,'alpha2':b,'alpha3':c,'z1':z,**variables}
    eq(row['orbit_id']+' S kernel',sp.sympify(row['J']['S'],locals=mapping),S)
    eq(row['orbit_id']+' P kernel',sp.sympify(row['J']['numerator'],locals=mapping),P)
    check(row['orbit_id']+' iterated order',row['J']['integration_order_outer_to_inner']==order)
    check(row['orbit_id']+' Euler factor once',row['C']['Gamma_Euler_order']==2 and row['C']['expression'].count('Exp[-2*EulerGamma*eps]')==1)

names={str(Z[n]):r'\zeta_{'+str(n)+'}' for n in Z}
names.update({f'M_{k}_{j}':rf'M_{{{k},{j}}}' for k in range(5) for j in range(1,6)})
names.update({f'N_{k}_{j}':rf'N_{{{k},{j}}}' for k in range(5) for j in range(1,6)})
names['tauD4']=r'\tau^+_{D,4}'
latex=lambda expr:sp.latex(expr,symbol_names={sp.Symbol(k):v for k,v in names.items()})
def aligned_coeffs(coeffs, max_terms=5):
    lines=[]
    for q,v in coeffs.items():
        if v==0:continue
        terms=sp.Add.make_args(sp.expand(v))
        # Order follows SymPy's reproducible polynomial order, then line-break terms.
        terms=sp.Add(*terms).as_ordered_terms()
        chunks=[terms[i:i+max_terms] for i in range(0,len(terms),max_terms)]
        for i,ch in enumerate(chunks):
            rhs=latex(sp.Add(*ch))
            if i and not rhs.startswith('-'):rhs='+ '+rhs
            lhs=rf't_{{{q}}}' if not i else ''
            lines.append(lhs+(' &={}' if not i else r' &\quad{}')+rhs)
    return ('\n'+r'\\[4pt]'+'\n').join(lines)

tex=r'''% Generated by build_tgv.py; source and arithmetic checks accompany the report.
The four representatives contain one tree-level three-gluon vertex. Their
attachment edges are $(3,4,4)$, $(3,3,4)$, $(2,4,4)$ and $(2,3,4)$, respectively.
Each orbit has multiplicity four. In this appendix write
$I_{3e\text{-}3g\text{-}r}=\sum_{k=-4}^{2}t^{(r)}_k\eps^k+O(\eps^3)$;
all coefficients below refer to the complete normalized $I$.

\subsection{The parameter integrals and normalization}
For $a,b,c\geq0$, $a+b+c=1$, let
\begin{equation}
 \dd\mu_\eps=\dd a\,\dd b\,\dd c\,\delta(1-a-b-c)(abc)^{-\eps},
 \qquad D_z(u,v)=(1-u)(1-v)+zuv.
\end{equation}
The original representative is $C_{\rm raw}J_r$, with
\begin{equation}
 C_{\rm raw}=\frac{C_AC_Fg_0^4\Gamma(2-2\eps)A_0^{2\eps}}
 {64\pi^{4-2\eps}},\qquad
 J_r=\int_{\mathcal D_r}\dd^3s\int\dd\mu_\eps\,
 \frac{P_r}{S_r^{2-2\eps}},\qquad A_0=-2v_{14}>0.
\end{equation}
Here $C=e^{-2\gamma_E\eps}C_{\rm raw}$ includes the Euler factor once.
After the four rotations and extraction of the two external scale powers,
\begin{equation}
 I_r=\frac{e^{-2\gamma_E\eps}\Gamma(2-2\eps)}{2}
 \frac{J_r(\eps,z)+z^{2\eps}J_r(\eps,1/z)}{1+z^{2\eps}}.
 \label{eq:tgv-normalization}
\end{equation}
The full parameter kernels are
\begin{align}
 S_{01}&=az(1-t)(bx+cy),&
 P_{01}&=\tfrac14az^2(1-t)(b-c),\\
 S_{02}&=czt[a(1-x)+b(1-y)],&
 P_{02}&=\tfrac14cz^2t(b-a),\\
 S_{03}&=a[bD_z(t,x)+cD_z(t,y)],&
 P_{03}&=\tfrac14(1+z)a(b-c)[1-(1+z)t],\\
 S_{04}&=ab\,y(1-x)+acD_z(x,w)+zbc\,w(1-y).
\end{align}
\begin{equation}
 \begin{aligned}
 4P_{04}={}&ab[-(1+z)(1-x)+zy]\,+ac[1+3z-(1+z)x+z(1+z)w]\\
 &+zbc[1-y-(1+z)w].
 \end{aligned}
\end{equation}
The original iterated path integrals, from outside to inside, are
\begin{equation}
\begin{array}{c|c|c}
 r&\text{variables}&\displaystyle\int_{\mathcal D_r}\dd^3s\\ \hline
 01&(t,x,y)=(s_{31},s_{41},s_{42})&\int_0^1\dd t\int_0^1\dd y\int_0^y\dd x\\[2pt]
 02&(x,y,t)=(s_{31},s_{32},s_{41})&\int_0^1\dd y\int_0^1\dd t\int_0^y\dd x\\[2pt]
 03&(t,x,y)=(s_{21},s_{41},s_{42})&\int_0^1\dd t\int_0^1\dd y\int_0^y\dd x\\[2pt]
 04&(x,y,w)=(s_{21},s_{31},s_{41})&\int_0^1\dd x\int_0^1\dd y\int_0^1\dd w
\end{array}
\end{equation}
The simplex measure remains part of $J_r$. The endpoint and contact faces are
retained when performing integration by parts.

\clearpage
\subsection{3e-3g-01 and 3e-3g-02}
Reflection of the ordered domain, together with $a\leftrightarrow c$, gives
$I_{01}=I_{02}$. The common exact expression is
\begin{equation}
 I_{01}=I_{02}=\frac{G_2(\eps)}{32\eps^4}
 \left[2\Gamma(1+\eps)\Gamma(1-\eps)
 -\frac{\Gamma(1-\eps)^2}{\Gamma(2-2\eps)}\right],
 \qquad G_2=e^{-2\gamma_E\eps}\Gamma(1-2\eps).
\end{equation}
There is no remaining $z$ dependence. The seven Laurent coefficients are
\begin{equation}
\begin{aligned}
@COEFF01@
\end{aligned}
\end{equation}
The coefficients $t_k$ on this page apply to both 01 and 02. The expansion
through $\eps^2$ uses the Gamma combinations through $\eps^6$.

\clearpage
\subsection{3e-3g-03}
The numerator is a path derivative,
$P_{03}=(1+z)(\partial_y-\partial_x)S_{03}/4$.
The boundary $x=y$ contributes and is included in the exact result
\begin{equation}
 I_{03}=\frac{G_2(\eps)}{8\eps^2\cosh(\eps\ell)}
 \left[\mathcal V_M(0,\eps;z)-
 \frac{\Gamma(1-\eps)^2}{\Gamma(2-2\eps)}\mathcal V_M(\eps,\eps;z)\right].
 \label{eq:tgv03-exact}
\end{equation}
Here $\mathcal V_M$ is the one-dimensional moment source defined in the
function dictionary. Its zeroth-order terms cancel. Thus
$t_{-4}=t_{-3}=t_{-2}=0$, and
\begin{equation}
\begin{aligned}
t_{-1}&=-\frac{M_{0,1}}4-\frac{M_{1,1}}8,\\[3pt]
t_0&=\left(\frac{\zeta_2}{8}-\frac12\right)M_{0,1}
 -\frac{M_{1,1}}4-\frac{M_{2,1}}{16},\\[5pt]
t_1&=-\frac{M_{3,1}+M_{1,3}}{48}-\frac{M_{2,1}}8
 -\frac{M_{0,1}+4\zeta_2+4}{8}M_{1,1}\\
&\quad-\frac5{24}M_{0,1}^2+
 \left(\frac{\zeta_3}4-\zeta_2-1\right)M_{0,1},\\[5pt]
t_2&=-\frac{M_{4,1}}{192}-\frac{M_{2,3}}{96}
 -\frac{M_{3,1}+M_{1,3}}{24}
 -\left(\frac{M_{0,1}}{16}+\frac{\zeta_2}4+\frac14\right)M_{2,1}\\
&\quad-\left(\frac{M_{0,1}}4+\zeta_2+\frac{\zeta_3}{12}+1\right)M_{1,1}
 +\left(\frac{5\zeta_2}{48}-\frac5{12}\right)M_{0,1}^2\\
&\quad+\left(\frac{59\zeta_2^2}{80}-2\zeta_2-\frac{\zeta_3}{6}-2\right)M_{0,1}.
\end{aligned}
\end{equation}
The identities $\ell^2=-2M_{0,1}-6\zeta_2$ and $M_{0,3}=-M_{0,1}^2$
have been used only to shorten the coefficients. All functions are evaluated
at the same positive real $z$.

\clearpage
\subsection{3e-3g-04}
Path integration by parts splits the original five-dimensional integral
into six endpoint faces and a nonzero body remainder. Endpoint subtraction
is performed before expanding the integrable coefficient kernels.
The Gamma and crossing factors in \eqref{eq:tgv-normalization} are then
convolved with the resulting $J$ coefficients.

For this representative $t_{-4}=t_{-3}=0$, and the complete result is
\begin{equation}
\begin{aligned}
t_{-2}&=\frac{\zeta_2}{8},\qquad
t_{-1}=\frac{-M_{1,1}+\zeta_3}{8},\\[4pt]
t_0&=\frac{M_{0,1}^2}{48}-\frac{M_{2,1}}8+\frac9{10}\zeta_2^2,\\[6pt]
t_1&=-\frac{5M_{3,1}}{96}-\frac{M_{1,3}}{96}
 -\frac{M_{0,1}M_{1,1}}{16}-\frac{3\zeta_2M_{1,1}}{16}\\
&\quad+\frac{3\zeta_3M_{0,1}}4+\frac58\zeta_5
 -\frac{37}{24}\zeta_2\zeta_3-\frac{N_{1,2}}{16},\\[6pt]
t_2&=\frac{223M_{0,1}^3}{11520}-\frac{13\zeta_2M_{0,1}^2}{1920}
 +\frac{169\zeta_2^2M_{0,1}}{320}\\
&\quad-\frac{5M_{0,1}M_{2,1}}{64}+\frac{M_{1,1}^2}{32}
 -\frac{39M_{4,1}}{1024}+\frac{25M_{2,3}}{1536}-\frac{N_{2,2}}{32}\\
&\quad-\frac{21\zeta_2M_{2,1}}{64}+\frac53\zeta_3M_{1,1}
 -\frac{19}{24}\zeta_3^2+\frac{7381}{1680}\zeta_2^3
 +\frac18\tau^+_{D,4}(z).
\end{aligned}
\end{equation}
The additional weight-six function $\tau^+_{D,4}$ is explicitly defined and
given as a finite GPL expression in the function dictionary. It was already
defined in the pure-exchange family; it is not a name assigned to an unknown
remainder. The coefficient at $\eps^1$ belongs to the $M/N$ language.
At $\eps^2$, the stated rational-coefficient, fixed-weight $M/N$ polynomial
space requires this additional direction. This finite-order comparison does
not establish an all-order minimal function basis.

The two positive-order coefficients above have been rechecked by exact
substitution into their full stored Hlog expressions. The original
multidimensional integrations and historical numerical controls are reused.
'''.replace('@COEFF01@',aligned_coeffs(I01,4))
(OUT/'tgv-results.tex').write_text(tex,encoding='utf-8')

# The complete tau expression is a word table, with all equal words combined.
byword={}
for term in tau_data['hlog_terms']:
    word=tuple(term['word']);cp=sp.Rational(term['coefficient'])
    for n,power in zip((2,3,5),term['zeta_powers']):cp*=Z[n]**power
    byword[word]=byword.get(word,0)+cp
items=list(sorted(byword.items()))
lines=[]
half=(len(items)+1)//2
for i in range(half):
    entries=[]
    for idx in (i,i+half):
        if idx<len(items):
            word,cp=items[idx]
            label=r'\varnothing' if not word else '('+','.join(map(str,word))+')'
            entries.extend([f'${label}$',f'${latex(cp)}$'])
        else:entries.extend(['',''])
    lines.append(' & '.join(entries)+r'\\')
dictionary=r'''% Generated from the existing exact tauDq4 Hlog data.
\subsection{The additional function at weight six}
\label{app:tau-d4}
With $D_z(t,s)=(1-t)(1-s)+zts$, define
\begin{equation}
 \begin{aligned}
 X_+&=\log\frac{D_z(t,s)\,ts}{\sqrt z},&
 X_-&=\log\frac{D_z(t,s)(1-t)(1-s)}{\sqrt z},\\
 \tau^+_{D,4}(z)&=-\frac{1+z}{48}
 \int_0^1\!\int_0^1\frac{X_+^4+X_-^4}{D_z(t,s)}\,\dd t\,\dd s.
 \end{aligned}
\end{equation}
This convergent moment definition fixes the normalization and branch for
$z>0$. For direct symbolic evaluation, the same function has the finite
hyperlogarithmic representation
\begin{equation}
 \tau^+_{D,4}(z)=\sum_{w}c_w G(w;z),\qquad
 G(a_1,\ldots,a_n;z)=\int_0^z\frac{\dd u}{u-a_1}
 G(a_2,\ldots,a_n;u).
\end{equation}
Here $G(\varnothing;z)=1$ and
$G(\underbrace{0,\ldots,0}_{n};z)=\log^n z/n!$.
Every letter is $0$ or $-1$; no contour crosses a nonzero letter for $z>0$.
The following table is the complete list of nonzero $c_w$.
\begin{center}
\small\setlength{\tabcolsep}{5pt}\renewcommand{\arraystretch}{1.24}
\begin{tabular}{@{}llll@{}}
\toprule
$w$ & $c_w$ & $w$ & $c_w$\\
\midrule
@TAUROWS@
\bottomrule
\end{tabular}
\end{center}
'''.replace('@TAUROWS@','\n'.join(lines))
(OUT/'tgv-functions.tex').write_text(dictionary,encoding='utf-8')

for rel in ['benchmarks/tgv/3e-3g/all-epsilon-status.md','benchmarks/tgv/3e-3g/high-order-04.md',
            'workflow/checks/tgv04-high-order/compact.py']:
    sources.add(ROOT/rel)
sources.add(Path(__file__).resolve())
result={'format':'two-loop-report-coefficients/1','family':'TGV','domain':'z>0',
        'normalization':'Current full I, Euler factor once, C4 rotations included',
        'through_power':2,'remainder':'O(eps^3)',
        'integrals':{name:{str(q):{'expression':str(v),'latex':latex(v)} for q,v in cs.items()}
                     for name,cs in data.items()},'formal':'not_run'}
(HERE/'tgv-data.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
validation={'passed':True,'checks':checks,'fresh_integration':False,'formal':'not_run',
 'tau_hlog_terms':len(tau),'tau_unique_words':len(byword),
 'source_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in sorted(sources)},
 'outputs':[str((OUT/f).relative_to(ROOT)).replace('\\','/') for f in ('tgv-results.tex','tgv-functions.tex')]}
(HERE/'tgv-checks.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'passed':True,'checks':len(checks),'tau_unique_words':len(byword),'integrals':list(data)}))
