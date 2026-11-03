# 3e-3g-03: ordered contact faces become evaluated moments
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


Let $(t,x,y)=(s_{21},s_{41},s_{42})$ with outer-to-inner order $t,y,x$, $0<x<y<1$, and define
\[
D_z(t,u)=(1-t)(1-u)+ztu,\quad S=a[bD_z(t,x)+cD_z(t,y)],
\quad P=\frac{1+z}{4}a(b-c)[1-(1+z)t].
\]
The complete numerator identity $P=(1+z)(\partial_y-\partial_x)S/4$ implies the following actual triangular IBP, using $\partial S^{-\lambda}=-\lambda S^{-2+2e}\partial S$:
\[
J=-\frac{1+z}{4\lambda}\int dt\int[d\alpha]_3(abc)^{-e}
\left[\int_0^1dx\,S(t,x,1)^{-\lambda}+\int_0^1dy\,S(t,0,y)^{-\lambda}
-2\int_0^1dq\,S(t,q,q)^{-\lambda}\right].\tag{1}
\]
The two copies of the diagonal have separate orientations: the lower $y$ limit of $\partial_y$ and the upper $x$ limit of $-\partial_x$ both give a negative contact flux. A square-domain IBP would miss them.

Set $(a,b,c)=(1-u,uv,u(1-v))$, with Jacobian $u$, so $S=u(1-u)D_z(t,vx+(1-v)y)$. For the power in (1), the scale integral is $\int_0^1u^0(1-u)^{e-1}du=1/e$. Push the other two faces to the affine coordinate $q=vx+(1-v)y$: $y=1$ gives $q\in(1-v,1)$ with Jacobian $1/v$, and $x=0$ gives $q\in(0,1-v)$ with Jacobian $1/(1-v)$. After $v\mapsto1-v$ in the first expression, their combined weight, including the diagonal, is
\[
F_e(q)=\int_0^q v^{-e}(1-v)^{-1-e}dv+
\int_q^1v^{-1-e}(1-v)^{-e}dv-2B_e
=\frac{[q(1-q)]^{-e}-B_e}{e}.\tag{2}
\]
For a proof of the last equality, differentiate: the difference of the first two densities is $(2q-1)[q(1-q)]^{-1-e}$, which is also the derivative of the right side. Its integral over $q\in(0,1)$ is zero because $\int[q(1-q)]^{-e}=B_e$; the two face integrals have total average $2B_e$, fixing the constant. Auxiliary regulators are removed in the combined expression, never separately on divergent faces.

The five-variable original integral is now the explicit two-variable kernel
\[
J=-\frac{1+z}{4\lambda e^2}\int_0^1dt\int_0^1dq\,
D_z(t,q)^{2e-1}\{[q(1-q)]^{-e}-B_e\}.\tag{3}
\]
The remaining path primitive is elementary:
\[
\int_0^1D_z(t,q)^{2e-1}dt=
\frac{(zq)^{2e}-(1-q)^{2e}}{2e[(1+z)q-1]}.
\]
At $(1+z)q=1$ it is the continuous divided difference. With $\beta=1-q$, use the moment definitions in [the evaluated-moment appendix](appendices/evaluated-moments.md) to rewrite it exactly as
\[
-(1+z)z^{-e}\int_0^1D_z(t,1-\beta)^{2e-1}dt
=\frac{e^{eW}\sinh(eL)}{e(\beta-r)}.
\]
Consequently the weighted and unweighted terms in (3) are $\mathcal V(0,e)$ and $\mathcal V(e,e)=V(2e)$, respectively. The actual domain reflection gives $J(1/z)=z^{-2e}J(z)$, so (N) yields
\[
\boxed{I=\frac{G_2}{8e^2\cosh(e\ell)}[\mathcal V(0,e)-B_e\mathcal V(e,e)]}.\tag{4}
\]

There is an independent closed evaluation of the first source. Set $x=(1-\beta)/\beta$ in its defining integral. Decompose $(1+z)/[(1+x)(1-zx)]=1/(1+x)+z/(1-zx)$. The combined endpoint-difference Mellin integral uses
\[
\mathcal I(a)=\operatorname{PV}\int_0^\infty x^a\left[\frac1{1+x}+\frac z{1-zx}\right]dx
=\frac{\pi[z^{-a}\cos(\pi a)-1]}{\sin(\pi a)},
\]
obtained by the Beta reflection formulas for the two rational kernels, with the same temporary prescription for both. Their original sinh combination vanishes at $zx=1$, so its final value is an ordinary integral. Thus
\[
\mathcal V(0,e)=\frac{z^e\mathcal I(e)-z^{-e}\mathcal I(-e)}{2e}
=\frac{\pi[\cos(\pi e)-\cosh(e\ell)]}{e\sin(\pi e)}.
\]
This fixes both the function and its period constants, rather than merely its derivative.
For an explicit reproducible coefficient rule define $\mathcal E_n(l_1,\ldots,l_n)=[x^n]\exp(\sum_{k=1}^n l_kx^k)$ by
\[
\mathcal E_0=1,\qquad \mathcal E_n=\frac1n\sum_{k=1}^n k l_k\mathcal E_{n-k}.
\]
The only logarithmic series needed here are
\[
G_2=\exp\sum_{n\ge2}\frac{2^n\zeta_n}{n}e^n,\quad
F=\Gamma(1+e)\Gamma(1-e)=\exp\sum_{\substack{n\ge2\\n\text{ even}}}\frac{2\zeta_n}{n}e^n,
\]
\[
B_e=B(1-e,1-e)=\frac1{1-2e}\exp\sum_{n\ge2}\frac{(2-2^n)\zeta_n}{n}e^n.
\]
The finite recurrence, followed by finite convolution, specifies every displayed coefficient without numerical fitting.


For a finite all-order assembly let $b_i=[e^i]B_e$, $p_i=[e^i](G_2/\cosh(e\ell))$, and
\[
v_n=\sum_{j=0}^{\lfloor n/2\rfloor}\frac{M_{n-2j,2j+1}}{(n-2j)!(2j+1)!},\quad
u_n=\begin{cases}M_{0,n+1}/(n+1)!,&n\text{ even},\\0,&n\text{ odd},\end{cases}
\quad d_n=u_n-\sum_{i=0}^n b_i v_{n-i}.
\]
Then $d_0=0$ and, with negative-index coefficients zero,
\[
a_q=\frac18\sum_{i=0}^{q+2}p_i d_{q+2-i}\quad(-2\le q\le2),\qquad a_{-4}=a_{-3}=0.
\]
For the most demanding displayed order, $a_2=(d_4+p_2d_2+p_3d_1)/8$, since $p_0=1,p_1=d_0=0$. The jets through degree four are therefore sufficient, and the $M_{0,5}$ terms cancel in $d_4$ before final moment identities. This is the actual finite coefficient convolution, including pole times positive powers.

## Complete Laurent coefficients


Here $I=\sum_{q=-4}^2a_q(z)e^q+O(e^3)$ at fixed positive real $z$. Every order, including zeros, is explicit.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=- \frac{M_{0,1}}{4} - \frac{M_{1,1}}{8}.$$

$$a_{0}(z)=- \frac{M_{0,1}}{2} + \frac{\pi^{2} M_{0,1}}{48} - \frac{M_{1,1}}{4} - \frac{M_{2,1}}{16}.$$

$$a_{1}(z)=- \frac{5 M_{0,1}^{2}}{24} - \frac{M_{0,1} M_{1,1}}{8} - \frac{\pi^{2} M_{0,1}}{6} - M_{0,1} + \frac{M_{0,1} \zeta\left(3\right)}{4} - \frac{\pi^{2} M_{1,1}}{12} - \frac{M_{1,1}}{2} - \frac{M_{1,3}}{48} - \frac{M_{2,1}}{8} - \frac{M_{3,1}}{48}.$$

$$a_{2}(z)=- \frac{5 M_{0,1}^{2}}{12} + \frac{5 \pi^{2} M_{0,1}^{2}}{288} - \frac{M_{0,1} M_{1,1}}{4} - \frac{M_{0,1} M_{2,1}}{16} - \frac{\pi^{2} M_{0,1}}{3} - 2 M_{0,1} - \frac{M_{0,1} \zeta\left(3\right)}{6} + \frac{59 \pi^{4} M_{0,1}}{2880} - \frac{\pi^{2} M_{1,1}}{6} - M_{1,1} - \frac{M_{1,1} \zeta\left(3\right)}{12} - \frac{M_{1,3}}{24} - \frac{\pi^{2} M_{2,1}}{24} - \frac{M_{2,1}}{4} - \frac{M_{2,3}}{96} - \frac{M_{3,1}}{24} - \frac{M_{4,1}}{192}.$$

