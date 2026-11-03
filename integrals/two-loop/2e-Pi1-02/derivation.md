# 2e-Pi1-02: keep the four corners and integrate the opposite-edge source

We take $D=4-2e$, $z>0$, $Z=e^{-\gamma_Ee}\Gamma(1-e)$, and adopt the full one-loop transverse operator
\[
\Pi^{ab}_{\mu\nu}(k)=i\delta^{ab}\mathcal A_1(e)(k^2g_{\mu\nu}-k_\mu k_\nu)(-k^2-i0)^{-e},
\]
\[
\mathcal A_1=\frac{g^2\Gamma(e)\Gamma(1-e)^2}{2(4\pi)^{2-e}\Gamma(2-2e)}\mathcal B_1,
\quad\mathcal B_1=\frac{(10-6e)C_A-8(1-e)T_Fn_f}{3-2e}.
\]
The inserted momentum root is $\int ds\,dt\,(igv_i^\mu t^a)D_{0,\mu\rho}^{ab}(k)\Pi^{bc,\rho\sigma}(k)D_{0,\sigma\nu}^{cd}(k)(igv_j^\nu t^d)e^{-ik\cdot X}$ integrated over $d^Dk/(2\pi)^D$, where $D_{0,\mu\nu}^{ab}=-i\delta^{ab}g_{\mu\nu}/(k^2+i0)$ and $X=x_i(s)-x_j(t)$. The $k^2g$ and $kk$ pieces are both retained. This calculation adopts the inner one-loop operator; it does not recompute its quark/gluon/ghost subdiagrams.

For clarity, the bridge to the parameter root uses
\[
\mathcal F_a(X)=\int\frac{d^Dk}{(2\pi)^D}e^{-ikX}(-k^2-i0)^{-a}
=i\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}(-X^2+i0)^{a-D/2}.
\]
Contracting $k^2v_i\cdot v_j-(v_i\cdot k)(v_j\cdot k)$ uses $-\partial_\mu\partial_\nu\mathcal F_{2+e}$; define the Fourier coefficients by $\mathcal F_{1+e}=ih_1Q^{-1+2e}$ and $\mathcal F_{2+e}=ih_2Q^{2e}$. Their ratio is $h_2/h_1=-1/[8e(1+e)]$. For $Q=-X^2+i0$,
$\partial_\mu\partial_\nu Q^{2e}=-4eg_{\mu\nu}Q^{-1+2e}+8e(2e-1)X_\mu X_\nu Q^{-2+2e}$.
Set $r=v_{14}$, $f=X^2/(2r)$, $d=v_i\cdot v_j/r$, $u=v_i\cdot X/r$, $w=v_j\cdot X/r$. Removing the common scale $-2r+i0$ gives the complete parameter root
\[
J=\int ds\,dt\,\mathcal K,\qquad
\mathcal K=\frac{(1+2e)df+(1-2e)uw}{2(1+e)}(f+i0)^{-2+2e}.
\]
Gamma recursion, $\Gamma(e)/\Gamma(1+e)=1/e$ and $\Gamma(1-2e)/\Gamma(2-2e)=1/(1-2e)$, gives the prepared physical prefactor
\[
C=\left(\frac{\alpha_s}{\pi}\right)^2 C_F\mathcal B_1(-2\pi v_{14}+i0)^{2e}R_1,
\qquad R_1=-\frac{Z^2}{16e(1-2e)}.
\]
It contains exactly one two-loop Euler compensation. With cyclic multiplicity $m$, the normalized scalar is
\[
I=\frac m2R_1\frac{J(z)+z^{2e}J(1/z)}{1+z^{2e}}.\tag{P}
\]
This is the current native physical r002 sign: assembly uses $+\mathcal B_1$, with no additional minus sign. Both color coefficients of $\mathcal B_1$ must subsequently be convolved as rational functions of $e$. The scalar $I$ itself has no color or dimensionful scale. The common $i0$ is retained in $C$; the dimensionless evaluation below is on the positive-real branch.

For edges 2 and 4 let $(s,t)=(s_{21},s_{41})$, outer-to-inner in that order, both on $[0,1]$. Direct geometry gives
\[
f=(1-s)(1-t)+zst,\quad d=-(1+z),\quad u=-1+(1+z)t,\quad w=1-(1+z)s.
\]
Substitution into the full tensor root gives
\[
\mathcal K=\frac{-2(1+z)f+(1-2e)z}{2(1+e)}f^{2e-2}.
\]
The crucial exact identity is
\[
\mathcal K=-(1+z)f^{2e-1}+\frac{\partial_s\partial_tf^{2e}}{4e(1+e)}.\tag{1}
\]
To check it, use $f_{st}=1+z$ and $f_sf_t=(1+z)f-z$ in
$\partial_s\partial_tf^{2e}=2ef^{2e-1}f_{st}+2e(2e-1)f^{2e-2}f_sf_t$.

First take $\Re e>0$; at a vanishing corner $f$ is a positive linear radial function to leading order, and the tensor term is locally $r^{2\Re e-2}$ with two-dimensional area $r\,dr$. It is therefore integrable. The complete corner bracket is
\[
[f^{2e}]_{s=0}^{s=1}\bigg|_{t=0}^{t=1}
=f(1,1)^{2e}-f(0,1)^{2e}-f(1,0)^{2e}+f(0,0)^{2e}
=z^{2e}-0-0+1.
\]
The zero values at $(0,1)$ and $(1,0)$ are taken in this initial regulator domain; they are not defined by substituting $e=0$ first. Thus
\[
J_O=-(1+z)\int_0^1ds\int_0^1dt\,f^{2e-1}
+\frac{1+z^{2e}}{4e(1+e)}.\tag{2}
\]
The first term also has an explicit path primitive:
\[
\int_0^1dt\,f^{2e-1}=
\frac{(zs)^{2e}-(1-s)^{2e}}{2e[(1+z)s-1]}.
\]
With $\beta=1-s$, its divided difference becomes precisely $z^e\mathcal V(e,e)=z^e V(2e)$ as defined and evaluated in [the moment lemma](../3e-3g-03/appendices/evaluated-moments.md). This is the actual bridge from the opposite-edge square to that lemma, with its positive-real branch and removable interior point unchanged. Equation (2) is therefore
\[
J_O=z^e V(2e)+\frac{1+z^{2e}}{4e(1+e)}.
\]
Its domain reflection $(s,t)\mapsto(1-s,1-t)$ gives $J_O(1/z)=z^{-2e}J_O(z)$. With $m=2$ in (P),
\[
\boxed{I_O=-\frac{Z^2}{16e(1-2e)}\left[\frac{V(2e)}{\cosh(e\ell)}+\frac1{2e(1+e)}\right]}.\tag{3}
\]
The second term is an evaluated corner contribution and changes the double pole; it is not an optional subtraction.

For the finite assembly let $v_n=[e^n]V(2e)$, so
\[
v_0=M_{0,1},\ v_1=M_{1,1},\ v_2=\frac{M_{2,1}}2+\frac{M_{0,3}}6,
\quad v_3=\frac{M_{3,1}+M_{1,3}}6.
\]
Set $P(e)=Z^2/[(1-2e)\cosh(e\ell)]$ and $T(e)=Z^2/[(1-2e)(1+e)]$. The full explicit coefficient rule is
\[
a_q=-\frac1{16}\sum_{i=0}^{q+1}P_i v_{q+1-i}-\frac{T_{q+2}}{32},\qquad -2\le q\le2,
\]
where a sum with negative upper limit is zero. $a_{-4}=a_{-3}=0$. Here $Z^2=\exp(\sum_{n\ge2}2\zeta_ne^n/n)$, and each series is evaluated by the finite exponential recurrence given for the adjacent graph. At the highest order, $a_2=-(P_0v_3+P_1v_2+P_2v_1+P_3v_0)/16-T_4/32$. Thus the moment source needs degree three and the corner source degree four; they have different pole budgets. Low-index moment identities from the adopted lemma produce the compact table below.

## Complete Laurent coefficients


Here $I=\sum_{q=-4}^2a_q(z)e^q+O(e^3)$ at fixed positive real $z$. Every order, including zeros, is explicit.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=- \frac{1}{32}.$$

$$a_{-1}(z)=- \frac{M_{0,1}}{16} - \frac{1}{32}.$$

$$a_{0}(z)=- \frac{M_{0,1}}{8} - \frac{M_{1,1}}{16} - \frac{3}{32} - \frac{\pi^{2}}{192}.$$

$$a_{1}(z)=- \frac{5 M_{0,1}^{2}}{96} - \frac{\pi^{2} M_{0,1}}{24} - \frac{M_{0,1}}{4} - \frac{M_{1,1}}{8} - \frac{M_{2,1}}{32} - \frac{5}{32} - \frac{\pi^{2}}{192} - \frac{\zeta\left(3\right)}{48}.$$

$$a_{2}(z)=- \frac{5 M_{0,1}^{2}}{48} - \frac{M_{0,1} M_{1,1}}{16} - \frac{\pi^{2} M_{0,1}}{12} - \frac{M_{0,1}}{2} - \frac{M_{0,1} \zeta\left(3\right)}{24} - \frac{\pi^{2} M_{1,1}}{24} - \frac{M_{1,1}}{4} - \frac{M_{1,3}}{96} - \frac{M_{2,1}}{16} - \frac{M_{3,1}}{96} - \frac{11}{32} - \frac{\pi^{2}}{64} - \frac{7 \pi^{4}}{11520} - \frac{\zeta\left(3\right)}{48}.$$

