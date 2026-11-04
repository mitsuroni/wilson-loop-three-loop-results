# Serial Pi1 and single Pi2: tensor source to the finite coefficient

The four entries have two endpoint geometries but different inserted operators.
Use the [common normalization](common-normalization.md), \(e=\epsilon\),
\(q=-k^2-i0\), and the transverse tensor \(T^{\mu\nu}=k^2g^{\mu\nu}-k^\mu k^\nu\).
The original contribution is the Wilson trace of
\[
\mathcal A_{ij}=\int ds\,dt\,(igv_i^{\mu}t^a)(igv_j^{\nu}t^b)
\int\frac{d^Dk}{(2\pi)^D}e^{-ikX_{ij}}\,
\begin{cases}(D_0\Pi_1D_0\Pi_1D_0)^{ab}_{\mu\nu},&\text{chain},\\
(D_0\Pi_2D_0)^{ab}_{\mu\nu},&\text{single Pi2},\end{cases}\tag{1}
\]
where \(D_0=-i\delta g/(k^2+i0)\). There is one exchanged momentum; the
internal one- or two-loop operator has already been integrated. Dyson
coefficients are +1. Labelling identical serial insertions does not produce
another \(1/2!\).

## Contract the full tensor before doing endpoint integrals

The accepted one-loop operator, including its scalar bubble evaluation, is
\[
\Pi_1^{ab,\mu\nu}=i\delta^{ab}\mathcal A_1(e)T^{\mu\nu}q^{-e},\qquad
\mathcal A_1=\frac{g^2}{2(4\pi)^{2-e}}
\frac{\Gamma(e)\Gamma(1-e)^2}{\Gamma(2-2e)}\mathcal B_1,
\]
\[
\mathcal B_1=aC_A-bT_f,\quad T_f=T_Fn_f,\qquad
a=\frac{10-6e}{3-2e},\quad b=\frac{8(1-e)}{3-2e}.\tag{2}
\]
Because \(TgT=k^2T\) and \(k_\mu T^{\mu\nu}=0\), serial contraction
keeps one full transverse tensor. The colour deltas leave \(C_F\), so its
external colour operator is \(C_F\mathcal B_1^2\). It is the scalar operator
that is squared, not a coordinate integral or a completed single-Pi1 I.

The Pi2 scalar is a distinct accepted full-D bare 1PI source. In this same
insertion convention
\[
\Pi_2=i\delta\mathcal A_2Tq^{-2e},\quad
\mathcal A_2=-\frac{g^4}{(4\pi)^D}
\left[\frac{\eta(e)}{1-2e}\right]^2P_2(e),\quad
\eta=\frac{\Gamma(1+e)\Gamma(1-e)^2}{\Gamma(1-2e)}.\tag{3}
\]
Its precise \(P_2\) is evaluated in (15) below. This source adoption does not
claim that this reader calculation recomputes all internal two-loop self-energy
diagrams. The old full-D literature-to-native convention is preserved in the
[original Pi2 theory](runtime/snapshots/benchmarks/self-energy/pi2/theory.md).

Direct contraction of (1), retaining both Lorentz terms, gives
\[
\mathcal A_{ij}=-ig^2C_F\mathcal A_X
\int ds\,dt\int\frac{d^Dk}{(2\pi)^D}e^{-ikX}
\left[(v_i\cdot v_j)q^{-1-2e}+(v_i\cdot k)(v_j\cdot k)q^{-2-2e}\right],\tag{4}
\]
with \(\mathcal A_X=\mathcal A_1^2\) for the chain and \(\mathcal A_X=\mathcal A_2\)
for Pi2. The signs follow from the explicit source factors and common
Feynman boundary values, rather than a later sign chosen from a result table.
The accepted Fourier identity is
\[
\mathcal F_\alpha(X)=\int\frac{d^Dk}{(2\pi)^D}e^{-ikX}q^{-\alpha}
=i\frac{\Gamma(D/2-\alpha)}{4^\alpha\pi^{D/2}\Gamma(\alpha)}
Q^{\alpha-D/2},\qquad Q=-X^2+i0.
\]
It is interpreted by dimensional analytic continuation, not by asserting a
common real convergence strip for the separately split Fourier integrals.
Write \(\mathcal F_{1+2e}=ih_1Q^{-1+3e}\),
\(\mathcal F_{2+2e}=ih_2Q^{3e}\). Gamma recursion gives
\[
h_1=\frac{\Gamma(1-3e)}{4^{1+2e}\pi^{2-e}\Gamma(1+2e)},\qquad
\frac{h_2}{h_1}=-\frac1{12e(1+2e)}.
\]
Multiplication by \(k_\mu k_\nu\) corresponds to \(-\partial_\mu\partial_\nu\), and
\[
\partial_\mu\partial_\nu Q^{3e}
=-6e g_{\mu\nu}Q^{-1+3e}+12e(3e-1)X_\mu X_\nu Q^{-2+3e}.
\]
Substitution in (4) produces the complete coordinate tensor
\[
\mathcal A_{ij}=g^2C_F\mathcal A_Xh_1\int ds\,dt
\left[\frac{1+4e}{2(1+2e)}(v_i\cdot v_j)Q^{-1+3e}
-\frac{1-3e}{1+2e}(v_i\cdot X)(v_j\cdot X)Q^{-2+3e}\right].\tag{5}
\]
No longitudinal/contact contribution is removed. Its D-dimensional trace
coefficient is \(D-1\), using \(D=4-2e\).

## Reduce the actual two original domains

The adjacent source is ordered \(s_{31}\) then \(s_{41}\); set
\(x=s_{31},y=s_{41}\). Its source separation is
\(X=-[(1-x)v_3+yv_4]\). The opposite source is ordered \(s_{21}\) then
\(s_{41}\), with \(x=s_{21},y=s_{41}\) and
\(X=(1-y)v_4+v_1+xv_2\). Both parameters lie in \([0,1]\), with Jacobian1.
The null-quadrilateral Gram matrix, divided by \(v_{14}\), is
\[
\begin{pmatrix}0&z&-(1+z)&1\\z&0&1&-(1+z)\\-(1+z)&1&0&z\\1&-(1+z)&z&0\end{pmatrix}.
\]
Thus contracting (5) with these actual separations gives exactly
\[
J=\int_0^1dx\int_0^1dy\frac{N}{2(1+2e)}f^{-2+3e},\qquad
\begin{array}{ll}
f_A=z(1-x)y,&N_A=(2+e)zf_A,\\
f_O=D=(1-x)(1-y)+zxy,&N_O=-(2+e)(1+z)D+(1-3e)z.
\end{array}\tag{6}
\]
More explicitly, set \(r_0=v_{14}\), \(f=X^2/(2r_0)\),
\(d=v_i\cdot v_j/r_0\), \(u=v_i\cdot X/r_0\),
\(w=v_j\cdot X/r_0\). Factoring the common scale and the source's scalar
coefficient from (5) leaves
\[
N=(1+4e)df+(1-3e)uw.
\]
For the adjacent source \((d,u,w)=(z,-zy,-z(1-x))\), so \(uw=zf\).
For the opposite source \((d,u,w)=(-(1+z),-1+(1+z)y,1-(1+z)x)\),
so \(uw=z-(1+z)f\). Substituting these two identities gives precisely the
two numerators in (6); no extra transverse correction is inserted afterwards.
The same geometric J applies to the chain and Pi2; their extracted source
coefficients differ. For \(0<\Re e<1/3\), the adjacent integrals use
\(\int_0^1u^{3e-1}du=1/(3e)\), so
\[
J_A=\frac{(2+e)z^{3e}}{18e^2(1+2e)}.\tag{7}
\]
For the opposite geometry direct differentiation, using the full N, gives
\[
\frac{N_OD^{-2+3e}}{2(1+2e)}
=-(1+z)D^{-1+3e}+\frac{\partial_x\partial_yD^{3e}}{6e(1+2e)}.\tag{8}
\]
The four signed corner values are
\[
D(1,1)^{3e}-D(1,0)^{3e}-D(0,1)^{3e}+D(0,0)^{3e}
=z^{3e}-0-0+1.
\]
At each zero corner \(D\) is comparable to the sum of two nonnegative
distances; \(D^{-2+3e}\) is locally integrable for \(\Re e>0\).
Applying the fundamental theorem first with cutoffs and then dominated
convergence justifies the complete mixed-derivative boundary. Zero-corner
values are taken on that strip before continuation, not set to zero at e=0.
With the evaluated mother H below, (8) becomes
\[
J_O=z^{3e/2}H(3e,z)+\frac{1+z^{3e}}{6e(1+2e)}.\tag{9}
\]

## The retained opposite integral is actually evaluated

Define
\[
H(\delta,z)=-(1+z)z^{-\delta/2}\int_{[0,1]^2}
D_z(x,y)^{-1+\delta}dxdy.
\]
The corner bound makes H holomorphic for a small disk around zero.
Integrating y first gives the removable divided difference
\[
\int_0^1[D_z(x,y)]^{\delta-1}dy
=\frac{(zx)^\delta-(1-x)^\delta}{\delta[(1+z)x-1]}.
\]
It stays whole across \((1+z)x=1\). With \(\beta=1-x\),
\(r=z/(1+z)\), \(W=\log[\beta(1-\beta)]\),
\(L=\rho+\log[(1-\beta)/\beta]\), this proves
\[
H=\int_0^1\frac{e^{\delta W/2}}{\beta-r}
\frac{\sinh(\delta L/2)}{\delta/2}\,d\beta
=M_{0,1}+\frac\delta2M_{1,1}
+\delta^2\left(\frac18M_{2,1}+\frac1{24}M_{0,3}\right)+O(\delta^3).\tag{10}
\]
Logarithmic derivatives are dominated by a slightly stronger corner power,
which justifies coefficient extraction. To evaluate the mother as a standard
function, Mellin--Barnes on the two terms of D and both actual Beta integrals
gives, initially \(0<\delta<1\), \(\delta-1<\Re w<0\),
\[
H=-\frac{(1+z)z^{-\delta/2}}{\Gamma(1-\delta)\Gamma(1+\delta)^2}
\frac1{2\pi i}\int dw\,
\Gamma(-w)\Gamma(\delta-w)^2\Gamma(1-\delta+w)\Gamma(1+w)^2z^w.
\]
This is the defining Meijer G integral, hence
\[
H=-\frac{(1+z)z^{-\delta/2}}{\Gamma(1-\delta)\Gamma(1+\delta)^2}
G_{3,3}^{3,3}\!\left(z\;\middle|\;\begin{matrix}\delta,0,0\\0,\delta,\delta\end{matrix}\right).\tag{11}
\]
This Meijer G is distinct from the finite HLOG notation. The solved moments
needed in (10), adopted in their established positive-real branch, are the
following complete finite HLOG expressions. Put \(G_{a_1\cdots a_n}=G(a_1,\ldots,a_n;z)\):
\[
M_{0,1}=-3\zeta_2-G_{00},\qquad
M_{1,1}=-G_{000}+2\zeta_3+6\zeta_2G_{-1}-3\zeta_2G_0+2G_{-1,0,0},
\]
\[
\begin{aligned}
M_{2,1}={}&12\zeta_2G_{-1,0}+4G_{0,-1,0,0}-24\zeta_2G_{-1,-1}
-8G_{-1,-1,0,0}-\frac{33}5\zeta_2^2-4\zeta_2G_{00}\\
&-8\zeta_3G_{-1}+4\zeta_3G_0-2G_{0000}+4G_{-1,0,0,0}
+12\zeta_2G_{0,-1},\\
M_{0,3}={}&-9\zeta_2^2-6\zeta_2G_{00}-6G_{0000}=-M_{0,1}^2.
\end{aligned}\tag{12}
\]
There is no remaining unknown subintegral in (10). The exact original four
text outputs are under
`runtime/snapshots/workflow/checks/single-opposite/runs/20260916T090538777022Z/`;
the current Python parser adopts these evaluated words rather than rerunning
their historical integration. Their kernels, branch and concrete substitution
are specified here; agreement with a final I table is a separate check.

## Match C and assemble all required geometric coefficients

Define the actual Gamma ratio
\[
U=\frac{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)}
{\Gamma(1-2e)^2\Gamma(1+2e)},\qquad P=\mathcal G^3U.
\]
Removing \(C_F\mathcal B_1^2\) from the chain source gives
\(R_{\rm chain}=-P/[128e^2(1-2e)^2]\); removing
\(C_F\mathcal B_2\), \(\mathcal B_2=4e^2P_2\), from Pi2 gives its positive.
These signs follow from the bound native C rational components.
Substitute (7)/(9) into (N1) with m=4/2 to obtain the complete chain I:
\[
I_A^-=-\frac{P(2+e)}{1152e^4(1-2e)^2(1+2e)},\qquad
I_O^-=-\frac{P}{384e^3(1-2e)^2(1+2e)}
-\frac{PH(3e,z)}{128e^2(1-2e)^2\cosh(3e\rho/2)}.\tag{13}
\]
For Pi2, \(I_A^+=-I_A^-\), \(I_O^+=-I_O^-\).
This shared geometric I is the same unique physical normalization; it
excludes the complete external operator \(\phi_c\). Geometric equality of J
and opposite signs of I do not identify the two inserted operators.

For an explicit finite rule put \(v(e)=(1-2e)^{-2}\),
\(r(e)=v(e)/(1+2e)\), and set
\[
\alpha_n=[e^n]P(2+e)r,\quad \beta_n=[e^n]Pr,\quad
\chi_n=[e^n]\frac{PvH(3e,z)}{\cosh(3e\rho/2)}.
\]
All are finite sums of the Gamma recurrence (N3), the elementary rational
jets, and the four explicit HLOG expressions (12); no undefined integral
has been introduced. Negative coefficient indices mean zero. With
\(\sigma=-1\) for the chain and +1 for Pi2,
\[
c_k^A=\sigma\frac{\alpha_{k+4}}{1152},\qquad
c_k^O=\sigma\left(\frac{\beta_{k+3}}{384}+\frac{\chi_{k+2}}{128}\right),\qquad -6\le k\le0.\tag{14}
\]
This also specifies every leading zero. For c0 the adjacent expression needs
P/r through4; the opposite boundary needs through3 and its bulk through2.
The latter uses \(H(3e,z)\) only through2, with coefficients
\(M_{01},3M_{11}/2,9(M_{21}/8+M_{03}/24)\).

## Full Pi2 operator, its cancellation and colour finite terms

The accepted bare Feynman-gauge full-D operator is
\[
P_2=C_A^2p_{AA}+C_AT_fp_{Af}+C_FT_fp_{Ff},\qquad
p_c=A_c+B_c/U,\qquad\phi_c=4e^2p_c,\tag{15}
\]
with the actual six rational functions
\[
\begin{aligned}
A_{AA}&=\frac{2e^2-4e-1}{2e^3(2e-3)},&
B_{AA}&=\frac{(2e-1)(6e^5-e^4-44e^3+79e^2-30e-2)}
{2e^3(e+1)(2e-3)(3e-2)(3e-1)},\\
A_{Af}&=\frac{2(e-1)(2e^2-3e+2)}{e^3(2e-3)},&
B_{Af}&=-\frac{2(e-1)(2e-1)(2e^4-11e^3-10e^2+17e-4)}
{e^3(e+1)(2e-3)(3e-2)(3e-1)},\\
A_{Ff}&=-\frac{4(e-1)(2e^2-e+2)}{e^3(2e-3)},&
B_{Ff}&=\frac{8(e-1)(2e-1)^2(e^2-2e+2)}
{e^3(2e-3)(3e-2)(3e-1)}.
\end{aligned}\tag{16}
\]
Thus the physical monomials are \(C_A^2C_F,C_AC_FT_f,C_F^2T_f\).
The last is a genuine non-first-coordinate contribution and is retained;
there is no \(T_f^2\) Pi2 term. In contrast (2) squared gives the chain
multipliers \(a^2,-2ab,b^2\) of \(C_A^2C_F,C_AC_FT_f,C_FT_f^2\).

To calculate rather than name the needed Pi2 jet, Gamma expansion gives
\[
\log U=\sum_{n\ge2}\frac{2(-1)^n+1+3^n-2\,2^n-(-2)^n}{n}\zeta_ne^n
=6\zeta_3e^3+9\zeta_4e^4+42\zeta_5e^5+O(e^6).
\]
Use (N3) for \(u_n=[e^n]U\) and \(\bar u_n=[e^n]U^{-1}\), through5.
For \(a_{c,j}=[e^j]A_c\), \(b_{c,j}=[e^j]B_c\), both beginning at -3,
the exact finite-sum jet is
\[
\phi_{c,k}=4\left(a_{c,k-2}+\sum_{j=-3}^{k-2}b_{c,j}\bar u_{k-2-j}\right),
\quad
\psi_{c,k}=4\left(b_{c,k-2}+\sum_{j=-3}^{k-2}a_{c,j}u_{k-2-j}\right),
\quad \psi=U\phi.\tag{17}
\]
The apparent k=-1 terms cancel exactly. The surviving leading terms are
\(\phi_{AA,0}=-25/3\), \(\phi_{Af,0}=20/3\),
\(\phi_{Ff,0}=0\), \(\phi_{Ff,1}=8\).
Adjacent colour finite terms are
\(\sum_{k=0}^4\phi_{c,k}c^+_{-k}\); opposite ones need k through3.
For the chain replace \(\phi_c\) by \(a^2,-2ab,b^2\) and \(c^+\) by
\(c^-\). This is the actual finite-term rule for all four entries.

In (17), the adjacent k=4 term includes \(b_{c,-3}\bar u_5\); cutting U
at order4 would change a finite physical answer. A second organization gives
the same answer directly as \(\psi_c\) times the geometry with U removed.
The actual implementation checks these two full Laurent convolutions for
every Pi2 colour and coefficient before comparing with historical results.
The fully evaluated seven scalar/geometric coefficients and coloured finite
terms are printed in each entry's own article.

The supported [chains_pi2.py](runtime/scripts/chains_pi2.py) reproduces (6)–(17)
using native C/J, exact finite jets and the declared evaluated M words. It
does not read the historical complete-I or old operator-expansion table while
generating. No new internal two-loop integration or formal proof is claimed;
the accepted full-D operator and previously evaluated M dictionary retain
their original evidence scope.
