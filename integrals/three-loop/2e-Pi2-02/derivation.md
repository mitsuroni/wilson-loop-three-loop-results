# 2e-Pi2-02: One two-loop 1PI insertion, opposite edges

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
The residual native C after extracting $C_F\mathcal B_2$, $\mathcal B_2=4e^2P_2$ is
\[
R=+\frac{P}{128e^2(1-2e)^2}.
\]
This sign is in the actual native source; it is not a correction chosen from a historical result. The cyclic multiplicity is $m=2$. Crossing/orbit extraction gives
\[
I=\frac{m}2R\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.
\]
Exactly one prepared $e^{-3\gamma_Ee}$ compensation has been consumed by $\mathcal G^3$. The complete scalar is therefore
\[I=+\frac{P}{384e^3(1-2e)^2(1+2e)}+\frac{PH(3e,z)}{128e^2(1-2e)^2\cosh(3e\rho/2)}.\tag{2}\]

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
\[c_k=+1\left(\frac{\beta_{k+3}}{384}+\frac{\chi_{k+2}}{128}\right),\quad
c_0=+1\left(\frac{\beta_3}{384}+\frac{\chi_2}{128}\right).\]
The boundary product is required through3 and the bulk through2. This is why the pole times positive H orders cannot be discarded.

## All seven scalar/geometric coefficients

Write $I(e,z)=\sum_{k=-6}^{0}c_k(z)e^k+O(e)$ at fixed $z>0$. The explicit coefficients obtained from (2) are:

\[c_{-6}(z)=0.\]

\[c_{-5}(z)=0.\]

\[c_{-4}(z)=0.\]

\[c_{-3}(z)=\frac{1}{384}.\]

\[c_{-2}(z)=\frac{M_{0,1}}{128} + \frac{1}{192}.\]

\[c_{-1}(z)=\frac{M_{0,1}}{32} + \frac{3 M_{1,1}}{256} + \frac{\pi^{2}}{1536} + \frac{1}{48}.\]

\[c_{0}(z)=- \frac{9 M_{0,1} \rho^{2}}{1024} + \frac{\pi^{2} M_{0,1}}{512} + \frac{3 M_{0,1}}{32} + \frac{3 M_{0,3}}{1024} + \frac{3 M_{1,1}}{64} + \frac{9 M_{2,1}}{1024} + \frac{\pi^{2}}{768} + \frac{7 \zeta\left(3\right)}{384} + \frac{1}{24}.\]

## Complete external colours and the actual finite convolution

Use the actual three source multipliers $\phi_c=4e^2(A_c+B_c/U)$, $c=AA,Af,Ff$, from the six explicit rational functions in the shared derivation. The corresponding monomials are $C_A^2C_F$, $C_AC_FT_f$ and $C_F^2T_f$. The absent $T_f^2$ is a property of this 1PI source, not a reason to substitute the serial operator. Its geometric I excludes $p_c$; three operator coefficients multiply the same I.

For $a_{c,j}=[e^j]A_c$, $b_{c,j}=[e^j]B_c$, $\bar u_n=[e^n]U^{-1}$,
\[
\phi_{c,k}=4\left(a_{c,k-2}+\sum_{j=-3}^{k-2}b_{c,j}\bar u_{k-2-j}\right).
\]
The apparent $e^{-1}$ term cancels; the actual constant multipliers are $-25/3$, $20/3$, $0$ for AA/Af/Ff and the first Ff term is $8e$. All three complete colour series, including these structural leading zeros, are retained.
For the opposite finite term the conservative implementation computes U and its inverse through5; through4 is sufficient for this geometry's highest required multiplier $\phi_{c,3}$. It checks the equivalent route $(U\phi)(I/U)$ as well as $\phi I$ after the intermediate pole cancellation.
The complete finite coefficient for any listed channel is the explicit sum
\[F_c(0)=+1\left[\frac1{384}\sum_{j=0}^3\phi_{c,j}\beta_{3-j}
+\frac1{128}\sum_{j=0}^2\phi_{c,j}\chi_{2-j}\right].\tag{3}\]
The evaluated results of this full convolution are:

**AA**, monomial `CA^2*CF`:
\[F_{AA}(0)=\frac{75 M_{0,1} \rho^{2}}{1024} - \frac{14311 M_{0,1}}{13824} - \frac{25 \pi^{2} M_{0,1}}{1536} + \frac{M_{0,1} \zeta\left(3\right)}{32} - \frac{25 M_{0,3}}{1024} - \frac{583 M_{1,1}}{1536} - \frac{75 M_{2,1}}{1024} - \frac{177067}{248832} - \frac{283 \pi^{2}}{27648} + \frac{\pi^{4}}{5760} + \frac{55 \zeta\left(3\right)}{1152}.\]

**Af**, monomial `CA*CF*TF*nf`:
\[F_{Af}(0)=- \frac{15 M_{0,1} \rho^{2}}{256} + \frac{5 \pi^{2} M_{0,1}}{384} + \frac{M_{0,1} \zeta\left(3\right)}{4} + \frac{1961 M_{0,1}}{3456} + \frac{5 M_{0,3}}{256} + \frac{101 M_{1,1}}{384} + \frac{15 M_{2,1}}{256} + \frac{\zeta\left(3\right)}{32} + \frac{41 \pi^{2}}{6912} + \frac{\pi^{4}}{720} + \frac{14597}{62208}.\]

**Ff**, monomial `CF^2*TF*nf`:
\[F_{Ff}(0)=- \frac{M_{0,1} \zeta\left(3\right)}{2} + \frac{55 M_{0,1}}{96} + \frac{3 M_{1,1}}{32} - \frac{\pi^{4}}{360} - \frac{7 \zeta\left(3\right)}{36} + \frac{\pi^{2}}{192} + \frac{1195}{1728}.\]

The complete per-power coloured series is in [result.json](result.json). Equation(3), rather than just $\phi_c(0)c_0$, explains every coloured finite expression. M/HLOG coordinate substitution is exact; no numerical fit fixes these coefficients.

The scalar relation to the other operator's same geometry is $I_{\Pi_2}=-I_{\rm chain}$, because their actual residual C signs differ and their original J agrees. Their external operators and full coloured graph contributions remain different. No raw/unprojected source trace is set to zero by that relation.

See [reproduction](reproduce.md) and [validation](validation.md) for actual implementation and adopted-input scopes. The statement retains the source's positive-real-z branch and local regulator continuation; it adds no arbitrary complex-kinematics or direct $z=-1$ theorem.
