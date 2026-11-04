# 2e-Pi2-01: One two-loop 1PI insertion, adjacent edges

Work with the actual prepared source at $D=4-2e$, $e=\epsilon$, fixed real $z>0$ and $\rho=\log z$.
The target is the complete normalized scalar I before its external Pi2 colour operator through $e^0$, with every position $-6,\ldots,0$ explicit.
The inherited Fourier/scale prescription is retained. The [input](input.json) binds the graph, native full-colour source and ordered C/J, while the [shared derivation](../../families/three-loop-insertions/chains-pi2.md) gives the exact physical and mathematical lemmas used below.

## This momentum operator and this ordered kernel

The source is the normalized Wilson trace of
\[
\mathcal A_{ij}=\int_0^1ds\int_0^1dt\,(igv_i^\mu t^a)(igv_j^\nu t^b)
\int\frac{d^Dk}{(2\pi)^D}e^{-ikX_{ij}}
(D_0\Pi_2D_0)^{ab}_{\mu\nu}.
\]
There is exactly one bare two-loop1PI operator, not two one-loop insertions. Its full-D scalar is $P_2=C_A^2p_{AA}+C_AT_fp_{Af}+C_FT_fp_{Ff}$, with $p_c=A_c+B_c/U$. All six rational A/B functions and the Gamma U are explicitly given in equations(15)–(17) of the shared derivation. This internal operator is an accepted source input; the two internal loops are not re-integrated here. The $C_F^2T_f$ contribution remains part of this full graph.

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
The residual native C after extracting $C_F\mathcal B_2$, $\mathcal B_2=4e^2P_2$ is
\[
R=+\frac{P}{128e^2(1-2e)^2}.
\]
This sign is in the actual native source; it is not a correction chosen from a historical result. The cyclic multiplicity is $m=4$. Crossing/orbit extraction gives
\[
I=\frac{m}2R\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.
\]
Exactly one prepared $e^{-3\gamma_Ee}$ compensation has been consumed by $\mathcal G^3$. The complete scalar is therefore
\[I=+\frac{P(2+e)}{1152e^4(1-2e)^2(1+2e)}.\tag{2}\]

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
For this entry, $c_k=+1\,\alpha_{k+4}/1152$ and specifically
\[c_0=+1\,\alpha_4/1152.\]
The regular scalar product must therefore be retained through order4 before division by $e^4$.

## All seven scalar/geometric coefficients

Write $I(e,z)=\sum_{k=-6}^{0}c_k(z)e^k+O(e)$ at fixed $z>0$. The explicit coefficients obtained from (2) are:

\[c_{-6}(z)=0.\]

\[c_{-5}(z)=0.\]

\[c_{-4}(z)=\frac{1}{576}.\]

\[c_{-3}(z)=\frac{5}{1152}.\]

\[c_{-2}(z)=\frac{\pi^{2}}{2304} + \frac{1}{64}.\]

\[c_{-1}(z)=\frac{5 \pi^{2}}{4608} + \frac{7 \zeta\left(3\right)}{576} + \frac{5}{144}.\]

\[c_{0}(z)=\frac{67 \pi^{4}}{276480} + \frac{35 \zeta\left(3\right)}{1152} + \frac{\pi^{2}}{256} + \frac{7}{72}.\]

## Complete external colours and the actual finite convolution

Use the actual three source multipliers $\phi_c=4e^2(A_c+B_c/U)$, $c=AA,Af,Ff$, from the six explicit rational functions in the shared derivation. The corresponding monomials are $C_A^2C_F$, $C_AC_FT_f$ and $C_F^2T_f$. The absent $T_f^2$ is a property of this 1PI source, not a reason to substitute the serial operator. Its geometric I excludes $p_c$; three operator coefficients multiply the same I.

For $a_{c,j}=[e^j]A_c$, $b_{c,j}=[e^j]B_c$, $\bar u_n=[e^n]U^{-1}$,
\[
\phi_{c,k}=4\left(a_{c,k-2}+\sum_{j=-3}^{k-2}b_{c,j}\bar u_{k-2-j}\right).
\]
The apparent $e^{-1}$ term cancels; the actual constant multipliers are $-25/3$, $20/3$, $0$ for AA/Af/Ff and the first Ff term is $8e$. All three complete colour series, including these structural leading zeros, are retained.
For the finite adjacent term $\phi_{c,4}$ contains $4b_{c,-3}\bar u_5$. The inverse Gamma ratio is therefore required through order5 before the artificial internal cubic poles cancel. A regular jet cut at order4 before that cancellation would be insufficient. The implementation uses both $\phi I$ and $(U\phi)(I/U)$ as exact cross-checks.
The complete finite coefficient for any listed channel is the explicit sum
\[F_c(0)=+1\frac1{1152}\sum_{j=0}^4\phi_{c,j}\alpha_{4-j}.\tag{3}\]
The evaluated results of this full convolution are:

**AA**, monomial `CA^2*CF`:
\[F_{AA}(0)=- \frac{1187401}{559872} - \frac{2941 \pi^{2}}{62208} + \frac{\pi^{2} \zeta\left(3\right)}{576} + \frac{71 \pi^{4}}{276480} + \frac{7 \zeta\left(5\right)}{144} + \frac{605 \zeta\left(3\right)}{5184}.\]

**Af**, monomial `CA*CF*TF*nf`:
\[F_{Af}(0)=\frac{223 \pi^{4}}{207360} + \frac{\pi^{2} \zeta\left(3\right)}{72} + \frac{199 \pi^{2}}{7776} + \frac{7 \zeta\left(5\right)}{18} + \frac{73 \zeta\left(3\right)}{216} + \frac{97655}{139968}.\]

**Ff**, monomial `CF^2*TF*nf`:
\[F_{Ff}(0)=- \frac{557 \zeta\left(3\right)}{648} - \frac{7 \zeta\left(5\right)}{9} - \frac{\pi^{2} \zeta\left(3\right)}{36} - \frac{\pi^{4}}{324} + \frac{23 \pi^{2}}{864} + \frac{7993}{3888}.\]

The complete per-power coloured series is in [result.json](result.json). Equation(3), rather than just $\phi_c(0)c_0$, explains every coloured finite expression. M/HLOG coordinate substitution is exact; no numerical fit fixes these coefficients.

The scalar relation to the other operator's same geometry is $I_{\Pi_2}=-I_{\rm chain}$, because their actual residual C signs differ and their original J agrees. Their external operators and full coloured graph contributions remain different. No raw/unprojected source trace is set to zero by that relation.

See [reproduction](reproduce.md) and [validation](validation.md) for actual implementation and adopted-input scopes. The statement retains the source's positive-real-z branch and local regulator continuation; it adds no arbitrary complex-kinematics or direct $z=-1$ theorem.
