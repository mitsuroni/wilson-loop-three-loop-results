# 2e-Pi1+Pi1-02: Two serial one-loop insertions, opposite edges

Work with the actual prepared source at $D=4-2e$, $e=\epsilon$, fixed real $z>0$ and $\rho=\log z$.
The target is the complete normalized scalar I through $e^0$, with every position $-6,\ldots,0$ explicit.
The inherited Fourier/scale prescription is retained. The [input](input.json) binds the graph, native full-colour source and ordered C/J, while the [shared derivation](../../families/three-loop-insertions/chains-pi2.md) gives the exact physical and mathematical lemmas used below.

## This momentum operator and this ordered kernel

The source is the normalized Wilson trace of
\[
\mathcal A_{ij}=\int_0^1ds\int_0^1dt\,(igv_i^\mu t^a)(igv_j^\nu t^b)
\int\frac{d^Dk}{(2\pi)^D}e^{-ikX_{ij}}
(D_0\Pi_1D_0\Pi_1D_0)^{ab}_{\mu\nu}.
\]
Both transverse one-loop operators use the same exchanged momentum. The identity $TgT=k^2T$ leaves one transverse tensor; the external scalar is $C_F\mathcal B_1(e)^2$. The Dyson coefficient is1, with no extra factorial. It is not the square of a coordinate integral. Its full insertion scalar is given by equation(2) of the shared derivation.

The complete source tensor gives two coordinate terms, proportional to
$(v_i\cdot v_j)Q^{-1+3e}$ and $(v_i\cdot X)(v_j\cdot X)Q^{-2+3e}$, with respective coefficients $(1+4e)/[2(1+2e)]$ and $-(1-3e)/(1+2e)$. The shared Fourier derivation obtains these coefficients from $h_2/h_1=-1/[12e(1+2e)]$ before taking any endpoint or $D=4$ limit.

Here $(i,j)=(2,4)$, $x=s_{21}$, $y=s_{41}$ and $X=(1-y)v_4+v_1+xv_2$. The source order is $s_{21}$ outermost, then $s_{41}$; both intervals are $[0,1]$ and the Jacobian is1. The actual complete kernel is
\[
D=(1-x)(1-y)+zxy,\qquad N=-(2+e)(1+z)D+(1-3e)z,
\qquad J_O=\int_0^1dx\int_0^1dy\,\frac{ND^{-2+3e}}{2(1+2e)}.
\]
Using the mixed derivative of the *complete* distance power gives
\[
\frac{ND^{-2+3e}}{2(1+2e)}=-(1+z)D^{-1+3e}
+\frac{\partial_x\partial_yD^{3e}}{6e(1+2e)}.
\]
For $0<\Re e<1/3$ the original zero-corner model is $(u+v)^{-2+3e}$ and is absolutely integrable. Cutoff integration and dominated convergence therefore retain all four oriented corner terms:
\[
D(1,1)^{3e}-D(1,0)^{3e}-D(0,1)^{3e}+D(0,0)^{3e}
=z^{3e}-0-0+1.
\]
Consequently
\[
J_O=z^{3e/2}H(3e,z)+\frac{1+z^{3e}}{6e(1+2e)}.\tag{1}
\]
The boundary term is essential: it supplies the scalar's leading pole. The centered opposite integral H has been evaluated, not left as an unknown mother. Equation(11) of the shared derivation gives its Meijer-G value; equation(12) prints every HLOG word and endpoint period of the four required moments. Its exact needed jet is
\[
H(3e,z)=M_{0,1}+\frac{3e}{2}M_{1,1}
+9e^2\left(\frac18M_{2,1}+\frac1{24}M_{0,3}\right)+O(e^3).
\]
These are the raw unnormalised M moments; their complete defined functions and positive-real branch are specified in the shared derivation. No old complete-I table substitutes for (1).

## Match this source C and keep its complete finite budget

Write $P=\mathcal G^3U$, with $\mathcal G=e^{-\gamma_Ee}\Gamma(1-e)$ and
\[
U=\frac{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)}{\Gamma(1-2e)^2\Gamma(1+2e)}.
\]
The residual native C after extracting $C_F\mathcal B_1^2$ is
\[
R=-\frac{P}{128e^2(1-2e)^2}.
\]
This sign is in the actual native source; it is not a correction chosen from a historical result. The cyclic multiplicity is $m=2$. Crossing/orbit extraction gives
\[
I=\frac{m}2R\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.
\]
Exactly one prepared $e^{-3\gamma_Ee}$ compensation has been consumed by $\mathcal G^3$. The complete scalar is therefore
\[I=-\frac{P}{384e^3(1-2e)^2(1+2e)}-\frac{PH(3e,z)}{128e^2(1-2e)^2\cosh(3e\rho/2)}.\tag{2}\]

To make coefficient extraction finite and explicit, let
\[
r=(1-2e)^{-2}(1+2e)^{-1},\quad\beta_n=[e^n]Pr,\qquad
\chi_n=[e^n]\frac{P(1-2e)^{-2}H(3e,z)}{\cosh(3e\rho/2)}.
\]
Gamma expansion gives $P_0=1$, $P_1=0$, $P_2=3\zeta_2/2$, $P_3=7\zeta_3$. The complete coefficients needed for the scalar finite term are therefore
\[
\beta_3=16+3\zeta_2+7\zeta_3,
\]
\[
\chi_2=\left(12+\frac32\zeta_2-\frac98\rho^2\right)M_{0,1}
+6M_{1,1}+\frac98M_{2,1}+\frac38M_{0,3}.
\]
These follow from the displayed solved H jet, the rational coefficients of $(1-2e)^{-2}$ and the second-order sech term. Lower coefficients use their same finite Cauchy convolution; negative indices mean zero. None is an uncomputed subintegral.
For this entry,
\[c_k=-1\left(\frac{\beta_{k+3}}{384}+\frac{\chi_{k+2}}{128}\right),\quad
c_0=-1\left(\frac{\beta_3}{384}+\frac{\chi_2}{128}\right).\]
The boundary product is required through3 and the bulk through2. This is why the pole times positive H orders cannot be discarded.

## All seven scalar/geometric coefficients

Write $I(e,z)=\sum_{k=-6}^{0}c_k(z)e^k+O(e)$ at fixed $z>0$. The explicit coefficients obtained from (2) are:

\[c_{-6}(z)=0.\]

\[c_{-5}(z)=0.\]

\[c_{-4}(z)=0.\]

\[c_{-3}(z)=- \frac{1}{384}.\]

\[c_{-2}(z)=- \frac{M_{0,1}}{128} - \frac{1}{192}.\]

\[c_{-1}(z)=- \frac{M_{0,1}}{32} - \frac{3 M_{1,1}}{256} - \frac{1}{48} - \frac{\pi^{2}}{1536}.\]

\[c_{0}(z)=\frac{9 M_{0,1} \rho^{2}}{1024} - \frac{3 M_{0,1}}{32} - \frac{\pi^{2} M_{0,1}}{512} - \frac{3 M_{0,3}}{1024} - \frac{3 M_{1,1}}{64} - \frac{9 M_{2,1}}{1024} - \frac{1}{24} - \frac{7 \zeta\left(3\right)}{384} - \frac{\pi^{2}}{768}.\]

## Complete external colours and the actual finite convolution

For $a=(10-6e)/(3-2e)$, $b=8(1-e)/(3-2e)$, the three channel multipliers of this very scalar are $\phi_{AA}=a^2$, $\phi_{Af}=-2ab$, $\phi_{ff}=b^2$. They multiply $C_A^2C_F$, $C_AC_FT_f$, $C_FT_f^2$, respectively. There is no second Wilson colour factor from the intermediate propagator.
The complete finite coefficient for any listed channel is the explicit sum
\[F_c(0)=-1\left[\frac1{384}\sum_{j=0}^3\phi_{c,j}\beta_{3-j}
+\frac1{128}\sum_{j=0}^2\phi_{c,j}\chi_{2-j}\right].\tag{3}\]
The evaluated results of this full convolution are:

**AA**, monomial `CA^2*CF`:
\[F_{AA}(0)=\frac{25 M_{0,1} \rho^{2}}{256} - \frac{947 M_{0,1}}{864} - \frac{25 \pi^{2} M_{0,1}}{1152} - \frac{25 M_{0,3}}{768} - \frac{155 M_{1,1}}{288} - \frac{25 M_{2,1}}{256} - \frac{5845}{11664} - \frac{175 \zeta\left(3\right)}{864} - \frac{5 \pi^{2}}{324}.\]

**Af**, monomial `CA*CF*TF*nf`:
\[F_{Af}(0)=- \frac{5 M_{0,1} \rho^{2}}{32} + \frac{5 \pi^{2} M_{0,1}}{144} + \frac{161 M_{0,1}}{108} + \frac{5 M_{0,3}}{96} + \frac{7 M_{1,1}}{9} + \frac{5 M_{2,1}}{32} + \frac{13 \pi^{2}}{648} + \frac{35 \zeta\left(3\right)}{108} + \frac{899}{1458}.\]

**ff**, monomial `CF*(TF*nf)^2`:
\[F_{ff}(0)=\frac{M_{0,1} \rho^{2}}{16} - \frac{M_{0,1}}{2} - \frac{\pi^{2} M_{0,1}}{72} - \frac{M_{0,3}}{48} - \frac{5 M_{1,1}}{18} - \frac{M_{2,1}}{16} - \frac{133}{729} - \frac{7 \zeta\left(3\right)}{54} - \frac{\pi^{2}}{162}.\]

The complete per-power coloured series is in [result.json](result.json). Equation(3), rather than just $\phi_c(0)c_0$, explains every coloured finite expression. M/HLOG coordinate substitution is exact; no numerical fit fixes these coefficients.

The scalar relation to the other operator's same geometry is $I_{\Pi_2}=-I_{\rm chain}$, because their actual residual C signs differ and their original J agrees. Their external operators and full coloured graph contributions remain different. No raw/unprojected source trace is set to zero by that relation.

See [reproduction](reproduce.md) and [validation](validation.md) for actual implementation and adopted-input scopes. The statement retains the source's positive-real-z branch and local regulator continuation; it adds no arbitrary complex-kinematics or direct $z=-1$ theorem.
