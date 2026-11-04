# 2e-Pi1+Pi1-01: Two serial one-loop insertions, adjacent edges

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

Here $(i,j)=(3,4)$, $x=s_{31}$, $y=s_{41}$ and $X=-[(1-x)v_3+yv_4]$. The source order is $s_{31}$ outermost, then $s_{41}$; both intervals are $[0,1]$ and the Jacobian is1. Its actual reduced kernel is
\[
f=z(1-x)y,\qquad N=(2+e)zf,\qquad
J_A=\int_0^1dx\int_0^1dy\,\frac{Nf^{-2+3e}}{2(1+2e)}.
\]
The two ordered endpoint integrations use the actual powers $(1-x)^{3e-1}$ and $y^{3e-1}$. On $0<\Re e<1/3$ their primitives give
\[
J_A=\frac{(2+e)z^{3e}}{2(1+2e)}
\left[\frac1{3e}\right]\left[\frac1{3e}\right]
=\frac{(2+e)z^{3e}}{18e^2(1+2e)}.\tag{1}
\]
The endpoint values at zero are taken in this strip; local continuation about $e=0$ follows after these definite integrals. No residual endpoint integral remains.

## Match this source C and keep its complete finite budget

Write $P=\mathcal G^3U$, with $\mathcal G=e^{-\gamma_Ee}\Gamma(1-e)$ and
\[
U=\frac{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)}{\Gamma(1-2e)^2\Gamma(1+2e)}.
\]
The residual native C after extracting $C_F\mathcal B_1^2$ is
\[
R=-\frac{P}{128e^2(1-2e)^2}.
\]
This sign is in the actual native source; it is not a correction chosen from a historical result. The cyclic multiplicity is $m=4$. Crossing/orbit extraction gives
\[
I=\frac{m}2R\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.
\]
Exactly one prepared $e^{-3\gamma_Ee}$ compensation has been consumed by $\mathcal G^3$. The complete scalar is therefore
\[I=-\frac{P(2+e)}{1152e^4(1-2e)^2(1+2e)}.\tag{2}\]

To make coefficient extraction finite and explicit, let
\[
r=(1-2e)^{-2}(1+2e)^{-1},\qquad\alpha_n=[e^n]P(2+e)r.
\]
The exact logarithmic Gamma recurrence gives
\[
P_0=1,\quad P_1=0,\quad P_2=3\zeta_2/2,\quad P_3=7\zeta_3,
\quad P_4=39\zeta_4/4+9\zeta_2^2/8.
\]
The elementary rational factor $(2+e)r$ has coefficients $2,5,18,40,112$ through order4. Thus the actual highest needed coefficient is the explicit sum
\[
\alpha_4=112+18P_2+5P_3+2P_4.
\]
The other coefficients use the same finite Cauchy convolution, with negative indices zero. This adjacent geometry needs no opposite-kernel M evaluation.
For this entry, $c_k=-1\,\alpha_{k+4}/1152$ and specifically
\[c_0=-1\,\alpha_4/1152.\]
The regular scalar product must therefore be retained through order4 before division by $e^4$.

## All seven scalar/geometric coefficients

Write $I(e,z)=\sum_{k=-6}^{0}c_k(z)e^k+O(e)$ at fixed $z>0$. The explicit coefficients obtained from (2) are:

\[c_{-6}(z)=0.\]

\[c_{-5}(z)=0.\]

\[c_{-4}(z)=- \frac{1}{576}.\]

\[c_{-3}(z)=- \frac{5}{1152}.\]

\[c_{-2}(z)=- \frac{1}{64} - \frac{\pi^{2}}{2304}.\]

\[c_{-1}(z)=- \frac{5}{144} - \frac{7 \zeta\left(3\right)}{576} - \frac{5 \pi^{2}}{4608}.\]

\[c_{0}(z)=- \frac{7}{72} - \frac{\pi^{2}}{256} - \frac{35 \zeta\left(3\right)}{1152} - \frac{67 \pi^{4}}{276480}.\]

## Complete external colours and the actual finite convolution

For $a=(10-6e)/(3-2e)$, $b=8(1-e)/(3-2e)$, the three channel multipliers of this very scalar are $\phi_{AA}=a^2$, $\phi_{Af}=-2ab$, $\phi_{ff}=b^2$. They multiply $C_A^2C_F$, $C_AC_FT_f$, $C_FT_f^2$, respectively. There is no second Wilson colour factor from the intermediate propagator.
The complete finite coefficient for any listed channel is the explicit sum
\[F_c(0)=-1\frac1{1152}\sum_{j=0}^4\phi_{c,j}\alpha_{4-j}.\tag{3}\]
The evaluated results of this full convolution are:

**AA**, monomial `CA^2*CF`:
\[F_{AA}(0)=- \frac{120923}{104976} - \frac{707 \pi^{2}}{15552} - \frac{2765 \zeta\left(3\right)}{7776} - \frac{335 \pi^{4}}{124416}.\]

**Af**, monomial `CA*CF*TF*nf`:
\[F_{Af}(0)=\frac{67 \pi^{4}}{15552} + \frac{469 \zeta\left(3\right)}{972} + \frac{61 \pi^{2}}{972} + \frac{19597}{13122}.\]

**ff**, monomial `CF*(TF*nf)^2`:
\[F_{ff}(0)=- \frac{3179}{6561} - \frac{7 \pi^{2}}{324} - \frac{77 \zeta\left(3\right)}{486} - \frac{67 \pi^{4}}{38880}.\]

The complete per-power coloured series is in [result.json](result.json). Equation(3), rather than just $\phi_c(0)c_0$, explains every coloured finite expression. M/HLOG coordinate substitution is exact; no numerical fit fixes these coefficients.

The scalar relation to the other operator's same geometry is $I_{\Pi_2}=-I_{\rm chain}$, because their actual residual C signs differ and their original J agrees. Their external operators and full coloured graph contributions remain different. No raw/unprojected source trace is set to zero by that relation.

See [reproduction](reproduce.md) and [validation](validation.md) for actual implementation and adopted-input scopes. The statement retains the source's positive-real-z branch and local regulator continuation; it adds no arbitrary complex-kinematics or direct $z=-1$ theorem.
