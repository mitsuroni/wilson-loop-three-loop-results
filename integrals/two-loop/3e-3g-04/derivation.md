# 3e-3g-04: all six faces, the regularized recoil and the complete coefficient assembly
We use $e=\epsilon$, $D=4-2e$, fixed real $z>0$ and real logarithms, and set
\[
\lambda=1-2e,\quad G_2=e^{-2\gamma_Ee}\Gamma(1-2e),\quad
[d\alpha]_3=da\,db\,dc\,\delta(1-a-b-c).
\]
The parameter root is the complete prepared $C/J$ representative, rather than a newly rederived momentum tensor. Its upstream tree three-gluon contraction and Fourier integration give
\[
\mathcal A=CJ,\quad C=\frac{C_AC_Fg_0^4}{64\pi^{4-2e}}
e^{-2\gamma_Ee}\Gamma(2-2e)(-2v_{14}+i0)^{2e},\quad
J=\int_{\mathcal D}d^3s\int[d\alpha]_3(abc)^{-e}P S^{-2+2e}.
\]
The input binds the original parameter order and the full color representative. There is exactly one Euler compensation in $C$. With the four-element cyclic orbit, the definition of the normalized scalar is
\[
I=\frac{\lambda G_2}{2}\,
\frac{J(z,e)+z^{2e}J(1/z,e)}{1+z^{2e}}.\tag{N}
\]
Thus neither the cyclic factor nor a second Gamma factor is to be inserted later. The dimensionless positive-real calculation below strips the common physical scale and its $i0$; it does not discard the $i0$ in $C$. All changes of order and IBP are first made on a common absolutely convergent regulator strip, $1/2<\Re e<1$, or on a shared auxiliary-regulated version of a split endpoint combination, and then meromorphically continued near zero.


The three attachments are $(2,3,4)$. Use $x=s_{21},y=s_{31},w=s_{41}$, each in $[0,1]$, in that original outer-to-inner order, and the simplex $a+b+c=1$. The complete prepared root is
\[
S=aby(1-x)+acD_z(x,w)+zbcw(1-y),\qquad D_z(x,w)=(1-x)(1-w)+zxw,
\]
\[
4P=ab[-(1+z)(1-x)+zy]+ac[1+3z-(1+z)x+z(1+z)w]+zbc[1-y-(1+z)w].
\]
The parameter root, rather than a result table, is the starting point. A direct polynomial identity isolates the integrable boundary and body contributions:
\[
P=\frac{-zS_x-(1+z)S_y+S_w}{4}+P_R,
\qquad P_R=\frac{1+z}{2}c[a(1-x+zw)-zbw].\tag{1}
\]
No term in the original numerator is discarded. Define $F_{x,b}$ as the integral of $(abc)^{-e}S^{-1+2e}|_{x=b}$ over the simplex and the two remaining unit path intervals, and similarly for $y,w$. Let $\mathcal O$ be the original five-dimensional integral with $P_R$ in place of $P$. IBP gives
\[
J=\frac{z(F_{x,1}-F_{x,0})+(1+z)(F_{y,1}-F_{y,0})-(F_{w,1}-F_{w,0})}{4\lambda}+\mathcal O.\tag{2}
\]
The endpoint polynomials are explicit:

| Face | Restricted denominator |
|---|---|
| $x=0$ | $ac(1-w)+aby+zbcw(1-y)$ |
| $x=1$ | $zcw[a+b(1-y)]$ |
| $y=0$ | $acD_z(x,w)+zbcw$ |
| $y=1$ | $a[b(1-x)+cD_z(x,w)]$ |
| $w=0$ | $a(1-x)(c+by)$ |
| $w=1$ | $aby(1-x)+zacx+zbc(1-y)$ |

Use $\mathcal C=F_{w,0}$, $\mathcal H(e,z)=F_{x,0}$ and $\mathcal D(e,z)=e(F_{y,0}-F_{y,1})$. Parameter exchange and path reflection give $F_{x,1}=z^{2e-1}\mathcal C$ and $F_{w,1}=z^{2e-1}\mathcal H(e,1/z)$, hence
\[
\boxed{J=\frac{(1+z^{2e})\mathcal C-(1+z)\mathcal D/e-z\mathcal H(e,z)-z^{2e-1}\mathcal H(e,1/z)}{4(1-2e)}+\mathcal O(e,z).}\tag{3}
\]
These are named **specific integrals**, not unknown functions introduced to hide the finite term. The following steps evaluate exactly the jets required from them.

On the cusp face take $(a,b,c)=(1-u,u(1-v),uv)$, Jacobian $u$. The scale and $x$ integrals give $1/e$ and $1/(2e)$, while the remaining path primitive gives $[1-v^{2e}]/[2e(1-v)]$. The complete Beta difference is
\[
\mathcal C=\frac1{4e^3}\int_0^1v^{-e}(1-v)^{-1-e}(1-v^{2e})dv
=\frac{\Gamma(1+e)\Gamma(1-e)-\Gamma(1-e)^2/\Gamma(1-2e)}{4e^4}.\tag{4}
\]
Its leading powers are $\zeta_2/(2e^2)+\zeta_3/(2e)+2\zeta_2^2/5$, because the fourth- and third-order apparent poles cancel in the full numerator.

For $\mathcal D$, the two face reductions should be combined before expanding. On $y=0$ take $(a,b,c)=(u(1-v),uv,1-u)$, integrate $u$ and then $x$; on $y=1$ take $(a,b,c)=(1-u,u(1-v),uv)$, integrate $u$ and then $w$. In both cases the scale integral is $1/e$. The remaining unit variable is $\xi=w$ or $x$, respectively. Define
\[
\eta=\frac1{1+z},\quad Q=\eta(1-v)-\xi(\eta-v),\quad U=\log(Q/\eta),\quad V=\log(z\xi),\quad W_0=\log(1-\xi),\quad h=\log[v(1-v)].
\]
Then the actual reduced integral is
\[
\mathcal D=\int_0^1d\xi\int_0^1dv\,
\frac{\eta e^{-eh}}{\eta-\xi}
\left[\frac{e^{2eU}-e^{2eV}}{2e(1-v)}-\frac{e^{2eW_0}-e^{2eU}}{2ev}\right].\tag{5}
\]
$Q>0$ inside the square. At $\xi=\eta$, $U=V=W_0$, so the apparent pole is removable. Let
\[
E_n(A,B)=\sum_{j=0}^n\frac{(-h)^j2^{n-j}}{j!(n-j+1)!}(A^{n-j+1}-B^{n-j+1}).
\]
Every required coefficient is now an explicit integrable two-variable logarithmic kernel:
\[
D_n=[e^n]\mathcal D=\int_0^1d\xi\int_0^1dv\,
\frac\eta{\eta-\xi}\left[\frac{E_n(U,V)}{1-v}-\frac{E_n(W_0,U)}v\right],\quad n=0,1,2,3.\tag{6}
\]
Their actual HyperInt integration order is $v,\xi$; the endpoint differences stay together throughout. The two lowest evaluated outputs are
\[
D_0=-\frac{\zeta_2\ell+\ell^3/6}{1+z},
\]
\[
D_1=\frac{17\zeta_2^2/10-3\zeta_2G_{0,0}-5G_{0,0,0,0}+2\zeta_3G_0+2\zeta_2G_{-1,0}+2G_{-1,0,0,0}}{1+z}.
\]
The complete outputs $D_2,D_3$, with their endpoint constants, are printed in [the finite HLOG-jet appendix](appendices/hlog-jets.md).

For $\mathcal H$, take
\[
(a,b,c)=\frac{(r,s,1)}T,\ T=1+r+s,\quad y=\frac u{1+u},\ w=\frac v{1+v},
\quad F_H=r(1+u)+rsu(1+v)+zsv.
\]
All four variables range over $(0,\infty)$. The simplex Jacobian is $T^{-3}$, and the path Jacobians are $(1+u)^{-2}(1+v)^{-2}$. Including them gives
\[
\mathcal H=\int_0^\infty dr\,du\,ds\,dv\,
(rs)^{-e}T^{-1-e}(1+u)^{-1-2e}(1+v)^{-1-2e}F_H^{-1+2e}
=\int f_{H,0}e^{eL_H},
\]
\[
f_{H,0}=\frac1{T(1+u)(1+v)F_H},\quad
L_H=-\log r-\log s-\log T-2\log(1+u)-2\log(1+v)+2\log F_H.
\]
Thus $H_n=\int f_{H,0}L_H^n/n!$ for $n=0,1,2$. Its genuine first primitive, with $A=1+u+su(1+v)$, $B=zsv$, $C=1+s$, is
\[
\int_0^\infty\frac{dr}{(r+C)(Ar+B)}=\frac{\log[B/(AC)]}{B-AC}.
\]
The coincident-root limit $B=AC$ is finite. Continuing the actual order $r,u,s,v$ gives
\[
H_0=\frac2z[G_{-1,-1,0,0}+\zeta_3G_{-1}-\zeta_2G_{-1,0}+3\zeta_2G_{-1,-1}],
\]
and the full $H_1,H_2$ in the appendix. These follow from $f_{H,0}L_H$ and $f_{H,0}L_H^2/2$, respectively, not by differentiating a known final $I$.

For the recoil add $x=t/(1+t)$ in the same simplex chart. The fifth path Jacobian is $(1+t)^{-2}$. Define
\[
X=(1+t)(1+u)(1+v),\quad
F=rsu(1+v)+r(1+u)(1+ztv)+zsv(1+t),
\]
\[
b_e=(rs)^{-e}T^{-1-e}X^{-1-2e}F^{-2+2e},\qquad
R=\frac{1+z}{2}(1+u)[r(1+v)+z(r-s)v(1+t)].
\]
The original recoil density including all Jacobians is $Rb_e$. Its required subtraction is described in [the Euler regularization proof](appendices/euler-regularization.md). Specifically, for $E=r\partial_r-u\partial_u$,
\[
R_{\rm reg}=R-\frac{ER+R E\log b_e}{e}=A/e+B,
\]
\[
d_0=-r/T+u/(1+u)-2EF/F,\quad d_1=-1-r/T+2u/(1+u)+2EF/F,
\quad A=-ER-Rd_0,\ B=R(1-d_1).
\]
There is no ambiguity about the kernels hidden in $A,B$: $R,F,E,d_0,d_1$ are explicit rational polynomials/operators above. Set
\[
L_b=-\log r-\log s-\log T-2\log X+2\log F,\qquad b_e=b_0e^{eL_b}.
\]
The actual necessary coefficient integrals are
\[
O_{-1}=\int b_0A,\quad O_0=\int b_0(B+AL_b),\quad
O_1=\int b_0(BL_b+AL_b^2/2),\quad
O_2=\int b_0(BL_b^2/2+AL_b^3/6).\tag{7}
\]
All five integrations run over $(0,\infty)$ in the historical actual order $r,u,s,t,v$. The coefficient kernels were divergence-checked before integration. The evaluated pole is
\[
O_{-1}=-\tfrac12G_{-1,0,0}-\tfrac32\zeta_2G_{-1}-\tfrac12\zeta_3+\tfrac12\zeta_2\ell.
\]
$O_0,O_1,O_2$ are printed completely in the appendix. In the original Maple workflow, `convert(...,HlogRegInf)` first expresses logs as regularized hyperlogarithms, `hyperInt` does the successive rational-logarithmic integrations, and `fibrationBasis(...,[z])` places the answer on the same $z$ branch. The divergence and abort checks are enabled. A typical simple-pole primitive is $G(a,\mathbf w;q)$ for $G(\mathbf w;q)/(q-a)$; repeated poles use IBP, products use shuffle, and endpoint periods are kept in the output. The included explicit first $r$ primitive and full output jets specify what those CAS stages actually evaluate. The historical Maple calculation evaluated and assembled these jets. The new Python replay instead assigns the accepted closed $I_{-2},I_{-1},I_0$ formulas and reads `compact.json` for the already assembled normalized $I_1,I_2$, then compares their complete historical HLOG coordinates. It does not live-convolve $D/H/O$; equations (8)-(9) below expose the historical mathematical assembly.

For clarity about the order budget: $[e]I$ requires $D_2,H_1,O_1$ and $[e^2]I$ requires $D_3,H_2,O_2$. It would be wrong to stop $\mathcal D$ at $D_2$ because its prefactor contains $1/e$.

To show the highest-order assembly explicitly, let $\mathscr B_q$ denote the coefficients of the bracket in (3), and write $H_n^\iota=H_n(1/z)$. Besides
\[
\mathscr B_{-2}=\zeta_2,\quad
\mathscr B_{-1}=\zeta_2\ell+\zeta_3-(1+z)D_0,
\]
\[
\mathscr B_0=\zeta_2\ell^2+\zeta_3\ell+4\zeta_2^2/5-(1+z)D_1-zH_0-z^{-1}H_0^\iota,
\]
we need
\[
\begin{aligned}
\mathscr B_1={}&\tfrac23\zeta_2\ell^3+\zeta_3\ell^2+\tfrac45\zeta_2^2\ell+3\zeta_5-\zeta_2\zeta_3
-(1+z)D_2-zH_1-z^{-1}(H_1^\iota+2\ell H_0^\iota),\\
\mathscr B_2={}&\tfrac13\zeta_2\ell^4+\tfrac23\zeta_3\ell^3+\tfrac45\zeta_2^2\ell^2
+(3\zeta_5-\zeta_2\zeta_3)\ell+\tfrac{11}{14}\zeta_2^3-\zeta_3^2\\
&-(1+z)D_3-zH_2-z^{-1}(H_2^\iota+2\ell H_1^\iota+2\ell^2H_0^\iota).
\end{aligned}\tag{8}
\]
Then, for $-2\le q\le2$,
\[
j_q=[e^q]J=\tfrac14\sum_{k=0}^{q+2}2^k\mathscr B_{q-k}+O_q,\quad O_{-2}=0,
\]
so $j_2=(\mathscr B_2+2\mathscr B_1+4\mathscr B_0+8\mathscr B_{-1}+16\mathscr B_{-2})/4+O_2$.
The exact reflection $(a,b,c;x,y,w;z)\mapsto(c,b,a;1-w,1-y,1-x;1/z)$ gives $S\mapsto S/z$, $P\mapsto P/z^2$ and $J(1/z)=z^{-2e}J(z)$. Thus
\[
I=\frac{(1-2e)G_2}{1+z^{2e}}J=\sum n_ke^k\sum j_qe^q,
\qquad a_2=\tfrac12j_2+n_1j_1+n_2j_0+n_3j_{-1}+n_4j_{-2}.\tag{9}
\]
For reproducible exact evaluation, $n_k=[e^k](1-2e)\exp(\sum_{m=2}^4 2^m\zeta_me^m/m)/(1+\exp(2e\ell))$ through $k=4$. In particular $n_1=-1-\ell/2$, $n_2=\ell+\zeta_2$. This finite rule includes all pole times positive-order prefactors. The historical assembly scripts implement precisely (8)-(9), then convert the HLOG words and the separately evaluated moment dictionary to the same basis and compare every coordinate.

Use the evaluated $M,N$ convention in [the function appendix](../3e-3g-03/appendices/evaluated-moments.md). At order two the previously solved pure-exchange family also contributes
\[
\tau^+_{D,4}=-\frac{1+z}{48}\int_0^1dt\int_0^1ds\,
\frac{X_+^4+X_-^4}{D_z(t,s)},\quad
X_+=\log\frac{D_z(t,s)ts}{\sqrt z},\ X_-=\log\frac{D_z(t,s)(1-t)(1-s)}{\sqrt z}.
\]
Its adopted HLOG evaluation is the already evaluated double-endpoint moment in [the original tau dictionary](../../../benchmarks/2ne/theory/generation/indexed-extra-sources-20261016/tau-e-dictionary.md), on $z>0$ with the same HLOG convention. Its appearance below is a change of function coordinates for the fully evaluated answer, not an unevaluated recoil. The entire required order is already given by the finite HLOG jets and equations (8)-(9), independently of how the final compact functions are named.

## Complete Laurent coefficients


Write $I=\sum_{q=-4}^2a_q(z)e^q+O(e^3)$. The complete table is:

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=\frac{\pi^{2}}{48}.$$

$$a_{-1}(z)=- \frac{M_{1,1}}{8} + \frac{\zeta\left(3\right)}{8}.$$

$$a_{0}(z)=\frac{M_{0,1}^{2}}{48} - \frac{M_{2,1}}{8} + \frac{\pi^{4}}{40}.$$

$$a_{1}(z)=- \frac{M_{0,1} M_{1,1}}{16} + \frac{3 M_{0,1} \zeta\left(3\right)}{4} - \frac{\pi^{2} M_{1,1}}{32} - \frac{M_{1,3}}{96} - \frac{5 M_{3,1}}{96} - \frac{N_{1,2}}{16} - \frac{37 \pi^{2} \zeta\left(3\right)}{144} + \frac{5 \zeta\left(5\right)}{8}.$$

$$a_{2}(z)=\frac{223 M_{0,1}^{3}}{11520} - \frac{13 \pi^{2} M_{0,1}^{2}}{11520} - \frac{5 M_{0,1} M_{2,1}}{64} + \frac{169 \pi^{4} M_{0,1}}{11520} + \frac{M_{1,1}^{2}}{32} + \frac{5 M_{1,1} \zeta\left(3\right)}{3} - \frac{7 \pi^{2} M_{2,1}}{128} + \frac{25 M_{2,3}}{1536} - \frac{39 M_{4,1}}{1024} - \frac{N_{2,2}}{32} + \frac{\tau^+_{D,4}}{8} - \frac{19 \zeta^{2}\left(3\right)}{24} + \frac{7381 \pi^{6}}{362880}.$$

