# 3e-3g-01: separate the original ordered kernel

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

The attachments are $(3,4,4)$. Put $(t,x,y)=(s_{31},s_{41},s_{42})$, with original outer-to-inner path order $t,y,x$ and $0<x<y<1$. The simplex root is
\[
S=az(1-t)(bx+cy),\qquad P=az^2(1-t)(b-c)/4.
\]
Use $a=1-u$, $b=uv$, $c=u(1-v)$ and $x=yq$. The simplex and ordered-path Jacobians are respectively $u$ and $y$, and all five new variables range over $(0,1)$. Including both Jacobians, not merely substituting $S$, gives
\[
J=\frac{z^{2e}}4\int_0^1dt(1-t)^{-1+2e}
\int_0^1du(1-u)^{-1+e}\int_0^1dy\,y^{-1+2e}
\int_0^1dv[v(1-v)]^{-e}(2v-1)
\int_0^1dq[1-v+vq]^{-2+2e}.
\]
The first three integrations give $1/(2e)$, $1/e$, $1/(2e)$, and the actual primitive in $q$ is $[1-v+vq]^{-1+2e}/[v(-1+2e)]$. Consequently
\[
\frac J{z^{2e}}=\frac1{16e^3}\int_0^1dv\,v^{-e-1}(1-v)^{-e}(2v-1)
\frac{1-(1-v)^{-1+2e}}{-1+2e}.
\]
The endpoint difference is kept intact. To evaluate it, a common auxiliary endpoint regulator may first be added to the split Beta terms. Their combined meromorphic limit is
\[
K=\frac{2B(1-e,1-e)-B(-e,1-e)-2B(1-e,e)+B(-e,e)}{2e-1}.
\]
Gamma recursion gives $B(-e,1-e)=-(1-2e)B_e/e$, $B(1-e,e)=F/e$ and the corresponding combined-limit value $B(-e,e)=0$. The last equality is a limit of the Beta ratio, not a claim that an unregulated divergent integral equals zero. Therefore
\[
K=\frac{2F-B_e}{e(1-2e)},\qquad
J=\frac{z^{2e}(2F-B_e)}{16(1-2e)e^4}.
\]
For $J=z^{2e}j$, the two terms in the numerator of (N) are $z^{2e}j$ and $j$, respectively: explicitly $J(z)+z^{2e}J(1/z)=(1+z^{2e})j$. Thus the actual result is
\[
\boxed{I=\frac{G_2(2F-B_e)}{32e^4}}.
\]
The displayed intermediate crossing statement fixes why no $z$ remains. Apparent singularities at $e=1/2$ in the separated formulas are treated by the full meromorphic limit.

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

Because the pole is $e^{-4}$, both Gamma sources must be retained through degree six for the $e^2$ answer. Precisely, $a_q=[e^{q+4}]G_2(2F-B_e)/32$, $-4\le q\le2$. For example $2F-B_e=1-2e+(-4+3\zeta_2)e^2+\cdots$; multiplying $G_2=1+2\zeta_2e^2+\cdots$ yields $a_{-2}=(-4+5\zeta_2)/32$.

## Complete Laurent coefficients


Write $I=\sum_{q=-4}^{2}a_q(z)\epsilon^q+O(\epsilon^3)$ at fixed positive real $z$. Every coefficient, including the vanishing ones, is displayed below.

$$a_{-4}(z)=\frac{1}{32}.$$

$$a_{-3}(z)=- \frac{1}{16}.$$

$$a_{-2}(z)=- \frac{1}{8} + \frac{5 \pi^{2}}{192}.$$

$$a_{-1}(z)=- \frac{1}{4} - \frac{\pi^{2}}{96} + \frac{7 \zeta\left(3\right)}{48}.$$

$$a_{0}(z)=- \frac{1}{2} - \frac{\pi^{2}}{48} - \frac{\zeta\left(3\right)}{24} + \frac{119 \pi^{4}}{11520}.$$

$$a_{1}(z)=-1 - \frac{\pi^{2}}{24} - \frac{7 \pi^{4}}{5760} - \frac{\zeta\left(3\right)}{12} + \frac{31 \zeta\left(5\right)}{80} + \frac{23 \pi^{2} \zeta\left(3\right)}{288}.$$

$$a_{2}(z)=-2 - \frac{\pi^{2}}{12} - \frac{7 \pi^{4}}{2880} - \frac{\zeta\left(3\right)}{6} - \frac{\pi^{2} \zeta\left(3\right)}{144} - \frac{\zeta\left(5\right)}{40} + \frac{31 \zeta^{2}\left(3\right)}{144} + \frac{569 \pi^{6}}{161280}.$$

