"""Render the current exact three-loop result as an algebraically checked article."""
from pathlib import Path
import argparse,hashlib,json,re
import sympy as s

L,rho,CA,CF,T=s.symbols('L rho CA CF T');M01=s.Symbol('M_0_1')
Z={k:s.Symbol('zeta_'+str(k)) for k in (2,3,4,5,6)}
NAMES={'AA':r'C_FC_A^2','Af':r'C_FC_AT','Ff':r'C_F^2T','ff':r'C_FT^2'}
FIELDS=('U','bare','post_counterterm')
ROW=' '+chr(92)*2+'\n'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(text):return s.sympify(text,locals={'L':L,'rho':rho,'CA':CA,'CF':CF,'T':T})
def normal(expr):return s.expand(expr.subs(M01,-(rho**2+s.pi**2)/2))
def zero(expr):return normal(expr)==0
def zeta_form(expr):
    expr=s.expand(expr).subs({s.pi**6:216*Z[2]**3,s.pi**4:36*Z[2]**2,s.pi**2:6*Z[2],s.zeta(3):Z[3],s.zeta(5):Z[5]})
    return s.expand(expr)
def unzeta(expr):return s.expand(expr.subs({Z[k]:s.zeta(k) for k in Z}))
def latex(expr):
    x=s.latex(s.expand(expr));x=re.sub(r'M_{(\d+)}_(\d+)',r'M_{\1,\2}',x)
    x=re.sub(r'N_{(\d+)}_(\d+)',r'N_{\1,\2}',x)
    # SymPy's default printer renders M_0_1 as M_{0 1}; use exact symbol names.
    for symbol in expr.free_symbols:
        m=re.fullmatch(r'([MN])_(\d+)_(\d+)',str(symbol))
        if m:x=x.replace(s.latex(symbol),m[1]+'_{'+m[2]+','+m[3]+'}')
    return x
def signed_lines(expr,max_terms=5):
    terms=s.Add.make_args(s.expand(expr));terms=sorted(terms,key=s.default_sort_key)
    chunks=[]
    for i in range(0,len(terms),max_terms):
        part=sum(terms[i:i+max_terms],s.S.Zero);text=latex(part)
        if i and not text.startswith('-'):text='+ '+text
        chunks.append(text)
    return (ROW+'&'+r'\quad ').join(chunks)

def generate(root):
    source=root/'benchmarks/current/results/3loop.json';raw=read(source);channels={c:{key:({int(p):parse(e) for p,e in value.items()} if isinstance(value,dict) else parse(value)) for key,value in row.items()} for c,row in raw['channels'].items()}
    low={n:read(root/f'benchmarks/current/results/{n}loop.json')['coefficients'] for n in (1,2)}
    b0=(11*CA-4*T)/12;b1=(34*CA**2-20*CA*T-12*CF*T)/48
    colors={'AA':CA**2,'Af':CA*T,'Ff':CF*T,'ff':T**2};checks=[];rendered=[]
    def check(label,expr):
        assert zero(expr),(label,s.expand(expr));checks.append({'name':label,'exact_zero':True})
    def formula(pointer,expr):
        canonical=zeta_form(expr);check('zeta conversion '+pointer,unzeta(canonical)-expr)
        body=latex(canonical);rendered.append({'json_pointer':pointer,'exact_expression':str(canonical),'latex':body});return canonical
    text=r'''# Three-loop Wilson-loop exponent: four-color analytic results

The three-loop exponent is resolved into four color channels. This article gives their complete Laurent coefficients through the finite term, the renormalized finite functions, and the exact connection between the scale-stripped sum and the physical result. The coefficients are written in the established $M,N$ functions and zeta values.

**Read:** [finite functions](#finite-functions) · [all-laurent coefficients](#all-laurent-coefficients) · [function definitions](#function-definitions) · [exact data](../../benchmarks/current/results/3loop.json) · [reproduce](../../docs/reproduction.md).

## Quantity and normalization

Set $D=4-2\epsilon$, $\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi$ and $T=T_Fn_f$, with $T_F=1/2$. The renormalized connected exponent is

$$\log W_R=\sum_{n\ge1}\lambda^nF_n,\qquad
F_3=C_F\left[C_A^2F_{3,AA}+C_ATF_{3,Af}+C_FTF_{3,Ff}+T^2F_{3,ff}\right].$$

The same color decomposition applies to $U_3$, $B_3$ and $P_3$. An entry labelled $AA$, for example, is the coefficient after removing $C_FC_A^2$. These are coefficients of the connected logarithm.

For the null quadrilateral, set $s=2v_1\!\cdot v_2$, $t=2v_1\!\cdot v_4$ and $z=s/t$. Use $\rho=\log z$ and $L=(L_s+L_t)/2$, so $L_s=L+\rho/2$ and $L_t=L-\rho/2$. The original source domain, signs and Euler factors follow the [single physical convention](../../workflow/current-conventions.md). The displayed real-function representation uses $z>0$.

The calculation restores the external scale, adds coupling counterterms and separates the pole and finite parts:

$$S_n=e^{n\epsilon(L+\rho/2)}+e^{n\epsilon(L-\rho/2)},\qquad
B_n=S_nU_n,\qquad P_3=B_3+CT_3=(\log Z)_3+F_3+O(\epsilon),$$
$$CT_3=-\frac{2b_0}{\epsilon}B_2+
\left(\frac{b_0^2}{\epsilon^2}-\frac{b_1}{2\epsilon}\right)B_1,
\qquad b_0=\frac{11C_A-4T}{12},\quad
b_1=\frac{34C_A^2-20C_AT-12C_FT}{48}.$$

Each lower-loop $B_n$ keeps its own $S_n$. The complete lower-loop input series are [one loop](../../benchmarks/current/results/1loop.json) and [two loops](../../benchmarks/current/results/2loop.json).

## Function definitions

Write $\zeta_n=\zeta(n)$ and suppress the argument $z$ on $M_{k,j}$ and $N_{k,j}$ in the formulas. Their raw-moment definitions are

$$r=\frac{z}{1+z},\qquad W_\beta=\log[\beta(1-\beta)],\qquad
\ell_\beta=\rho+\log\frac{1-\beta}{\beta},$$
$$M_{k,j}(z)=\int_0^1\frac{W_\beta^k\ell_\beta^j}{\beta-r}\,d\beta,
\qquad N_{k,j}(z)=M_{k,j+1}(z)-\rho M_{k,j}(z),\qquad j\ge1.$$

The numerator vanishes at $\beta=r$ and defines the removable value of the integrand. Powers are unnormalized: the definition includes no factorials. In particular,

$$M_{0,1}=-\frac{\rho^2+\pi^2}{2},\qquad M_{0,3}=-M_{0,1}^2,\qquad \pi^2=6\zeta_2.$$

The stored common-$z$ HLOG representation uses $G(a,\mathbf w;z)=\int_0^zG(\mathbf w;t)\,dt/(t-a)$, $G(\varnothing;z)=1$, with letters $0,-1$ and the standard logarithmic endpoint convention $G(0;z)=\log z$. The current four-color totals close on $M,N,\zeta_n$ and the physical scale $L$.

## Finite functions

Let $f_{3,c}(z)=F_{3,c}(0,z)$. The expressions below separate the finite shape $f_{3,c}$ from its full scale dependence.

'''
    for c,a in channels.items():
        text+=f'### Finite {c}: ${NAMES[c]}$\n\n'
        f=formula(f'/channels/{c}/f_at_L0',a['f_at_L0']);F=formula(f'/channels/{c}/F',a['F'])
        check(c+' finite slice',a['F'].subs(L,0)-a['f_at_L0'])
        poly=s.Poly(s.expand(F-f),L);parts=[]
        for degree in range(poly.degree(),0,-1):
            value=poly.coeff_monomial(L**degree)
            if value:
                if len(s.Add.make_args(value))==1:
                    part=latex(value*L**degree);parts.append(part if part.startswith('-') else '+ '+part)
                else:parts.append(('+ L' if degree==1 else f'+ L^{{{degree}}}')+f'\\left({latex(value)}\\right)')
        text+='$$\n\\begin{aligned}\nF_{3,'+c+'}(L,z)={}&f_{3,'+c+'}(z)'+ROW+'&'+(ROW+'&').join(parts)+' .\n\\end{aligned}\n$$\n\n'
        text+='$$\n\\begin{aligned}\nf_{3,'+c+'}(z)={}&'+signed_lines(f)+' .\n\\end{aligned}\n$$\n\n'
        check(c+' finite constant is P0',a['F']-a['post_counterterm'][0])
    text+=r'''## All-Laurent coefficients

For each color channel $c\in\{AA,Af,Ff,ff\}$ define

$$U_{3,c}=\sum_{p=-6}^{0}u_{c,p}\epsilon^p+O(\epsilon),\quad
B_{3,c}=\sum_{p=-6}^{0}b_{c,p}\epsilon^p+O(\epsilon),\quad
P_{3,c}=\sum_{p=-6}^{-1}q_{c,p}\epsilon^p+F_{3,c}+O(\epsilon).$$

Every coefficient, including zero entries, is displayed below. The counterterm coefficient is $[\epsilon^p]CT_{3,c}=q_{c,p}-b_{c,p}$ for $p<0$ and $F_{3,c}-b_{c,0}$ for $p=0$; its complete evaluated values are also stored in the exact data.

'''
    for c,a in channels.items():
        text+=f'### Laurent {c}: ${NAMES[c]}$\n\n'
        for field,label,symbol,powers in [('U','Scale-stripped sum','u',range(-6,1)),('bare','Scale-restored bare exponent','b',range(-6,1)),('post_counterterm','Pole coefficients after coupling renormalization','q',range(-6,0))]:
            text+='**'+label+'.**\n\n'
            for p in powers:
                e=formula(f'/channels/{c}/{field}/{p}',a[field][p]);text+='$$\n\\begin{aligned}\n'+f'{symbol}_{{{c},{p}}}='+r'{}&'+signed_lines(e)+' .\n\\end{aligned}\n$$\n\n'
        text+=f'The zero-order coefficient is $q_{{{c},0}}=F_{{3,{c}}}$, given above.\n\n'
        for p in range(-6,1):
            check(f'{c} scale restoration {p}',a['bare'][p]-sum((3**n/s.factorial(n))*((L+rho/2)**n+(L-rho/2)**n)*a['U'][p-n] for n in range(p+7)))
            check(f'{c} counterterm assembly {p}',a['post_counterterm'][p]-a['bare'][p]-a['counterterm'][p])
            ct=-2*b0*parse(low[2]['bare'].get(str(p+1),'0'))+b0**2*parse(low[1]['bare'].get(str(p+2),'0'))-b1/2*parse(low[1]['bare'].get(str(p+1),'0'))
            color=s.Poly(s.expand(ct),CA,CF,T).coeff_monomial(colors[c]);check(f'{c} lower-loop counterterm {p}',a['counterterm'][p]-color)
            if p<0:check(f'{c} pole/logZ {p}',a['post_counterterm'][p]-a['logZ'][p])
    text+=r'''## Renormalization-group coefficients

The cusp and noncusp coefficients have the same four-color decomposition, $k_3=C_F\sum_c\mathcal C_ck_{3,c}$ and $h_3=C_F\sum_c\mathcal C_ch_{3,c}$. In the convention above,

$$\partial_LF_3-2b_0F_2-b_1F_1=-2k_3L-h_3.$$

'''
    totalF=sum(colors[c]*a['F'] for c,a in channels.items());totalk=sum(colors[c]*a['k3'] for c,a in channels.items());totalh=sum(colors[c]*a['h3'] for c,a in channels.items())
    check('complete four-color finite RG',s.diff(totalF,L)-2*b0*parse(low[2]['F'])-b1*parse(low[1]['F'])+2*totalk*L+totalh)
    for c,a in channels.items():
        k=formula(f'/channels/{c}/k3',a['k3']);h=formula(f'/channels/{c}/h3',a['h3'])
        text+=f'**{c}.**\n\n$$k_{{3,{c}}}={latex(k)},\\qquad h_{{3,{c}}}={latex(h)}.$$\n\n'
    text+=r'''## Exact data and calculation record

The [three-loop JSON](../../benchmarks/current/results/3loop.json) stores these formulas under `channels.AA`, `channels.Af`, `channels.Ff` and `channels.ff`: `U`, `bare`, `counterterm`, `post_counterterm`, `logZ`, `F`, `f_at_L0`, `k3` and `h3`. `native_common_z_HLOG_channels` supplies the corresponding HLOG construction. The [current source catalogue](../../benchmarks/current/catalog.json) connects all163 source orbits to their inputs and selected integral objects.

The [source/joint ledger](../2026-11-02-current/source-ledger.json) and [matter ledger](../2026-11-02-current/matter-ledger.json) record their consumption. The accepted integral derivations, exact function dictionary, source assembly and RG checks form the calculation record. Read the [integral derivations](../../integrals/README.md), [calculation validation](../../validation/README.md) and [formula-presentation check](../../validation/results-presentation-20261103/README.md) for the corresponding steps.

This article is generated from the current exact data by [the formula renderer](../../tools/render_four_color_results.py). Its check reproduces every displayed coefficient, the finite slices, scale restoration, lower-loop counterterms, pole separation and the four-color finite RG relation.

```text
python -B tools/render_four_color_results.py --check
python -m workflow.current_calculation result 3
```

Earlier releases and dated research records are indexed in the [historical archive](../../docs/archive/README.md).
'''
    return text,{'status':'passed','source_sha256':sha(source),'source':'benchmarks/current/results/3loop.json','lower_loop_sha256':{str(n):sha(root/f'benchmarks/current/results/{n}loop.json') for n in (1,2)},'displayed_exact_expressions':len(rendered),'checks':checks,'rendered_formulas':rendered,'scope':'Exact presentation and algebraic assembly checks of the current accepted result.'}
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--check',action='store_true');a=p.parse_args();root=a.root.resolve()
    text,report=generate(root);article=root/'results/current/four-color-results.md'
    if a.check:assert article.read_text(encoding='utf-8')==text,'Article differs from the current exact renderer'
    else:
        article.parent.mkdir(parents=True,exist_ok=True);article.write_text(text,encoding='utf-8')
        report['article_sha256']=sha(article);out=root/'validation/results-presentation-20261103/formula-check.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print({'status':'passed','displayed_expressions':report['displayed_exact_expressions'],'exact_checks':len(report['checks']),'article':str(article)})
if __name__=='__main__':main()
