# 3e-3g-02: the actual ordered-domain map to the evaluated representative

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

Here $(x,y,t)=(s_{31},s_{32},s_{41})$, $0<x<y<1$, and $t\in(0,1)$. The original ordered integral has
\[
S_{02}=czt[b(1-y)+a(1-x)],\qquad P_{02}=cz^2t(b-a)/4.
\]
To use the already evaluated integral, make the full map
\[
(t,x,y;a,b,c)_{01}=(1-t,1-y,1-x;c,b,a)_{02}.
\]
The triangle $0<x<y<1$ maps to itself, all absolute Jacobians equal one, and $(abc)^{-e}$ and the simplex constraint are unchanged. Direct pullback gives
\[
S_{01}=az(1-t)(bx+cy)\longmapsto czt[b(1-y)+a(1-x)]=S_{02},
\]
\[
P_{01}=az^2(1-t)(b-c)/4\longmapsto cz^2t(b-a)/4=P_{02}.
\]
The numerator sign, denominator, path bounds and simplex measure therefore agree together. This proves equality first in the common convergence domain and then by meromorphic continuation. Both representatives also have the same $C$ and cyclic multiplicity. Hence $J_{02}=J_{01}$ and $I_{02}=I_{01}$; equality is not inferred from a picture or a final table.

The exact adopted representative lemma is the integral in [01](../3e-3g-01/derivation.md): for the displayed $S_{01},P_{01}$ on its displayed ordered domain,
\[
J_{01}=\frac{z^{2e}[2\Gamma(1+e)\Gamma(1-e)-B(1-e,1-e)]}{16(1-2e)e^4}.
\]
Its proof explicitly performs the simplex substitution $(a,b,c)=(1-u,uv,u(1-v))$, the path substitution $x=yq$, three scale integrals, the $q$ primitive and the combined Beta endpoint limit. Applying that proved lemma with the map above, then (N), gives
\[
I_{02}=\frac{G_2}{32e^4}[2\Gamma(1+e)\Gamma(1-e)-B_e].
\]
This adoption is an exact change-of-variables identity for this graph, not an assumption of equality of Laurent coefficients.

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

As in 01, $a_q=[e^{q+4}]G_2(2F-B_e)/32$. Degree six is necessary and sufficient for $q=2$; the same recurrence performs the explicit expansion and no new quadrature is required for 02.

## Complete Laurent coefficients


Write $I=\sum_{q=-4}^{2}a_q(z)\epsilon^q+O(\epsilon^3)$ at fixed positive real $z$. Every coefficient, including the vanishing ones, is displayed below.

$$a_{-4}(z)=\frac{1}{32}.$$

$$a_{-3}(z)=- \frac{1}{16}.$$

$$a_{-2}(z)=- \frac{1}{8} + \frac{5 \pi^{2}}{192}.$$

$$a_{-1}(z)=- \frac{1}{4} - \frac{\pi^{2}}{96} + \frac{7 \zeta\left(3\right)}{48}.$$

$$a_{0}(z)=- \frac{1}{2} - \frac{\pi^{2}}{48} - \frac{\zeta\left(3\right)}{24} + \frac{119 \pi^{4}}{11520}.$$

$$a_{1}(z)=-1 - \frac{\pi^{2}}{24} - \frac{7 \pi^{4}}{5760} - \frac{\zeta\left(3\right)}{12} + \frac{31 \zeta\left(5\right)}{80} + \frac{23 \pi^{2} \zeta\left(3\right)}{288}.$$

$$a_{2}(z)=-2 - \frac{\pi^{2}}{12} - \frac{7 \pi^{4}}{2880} - \frac{\zeta\left(3\right)}{6} - \frac{\pi^{2} \zeta\left(3\right)}{144} - \frac{\zeta\left(5\right)}{40} + \frac{31 \zeta^{2}\left(3\right)}{144} + \frac{569 \pi^{6}}{161280}.$$

