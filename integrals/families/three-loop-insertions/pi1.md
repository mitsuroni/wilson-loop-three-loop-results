# One transverse Pi1 insertion and one bare exchange

The eleven sources have four ordered Wilson endpoints. Their two independent exchange momenta carry, respectively, the accepted full operator $D_0\Pi_1D_0$ and a bare propagator. The bubble adds the third loop. The [single project convention](../../../workflow/current-conventions.md) determines the physical sign, colour projection and Euler compensation before integration. Work at $D=4-2e$, $e=\epsilon$, fixed $z>0$, real $\rho=\log z$, initially $0<\Re e<1/8$, followed by local meromorphic continuation.

## The transverse root and the normalization

If the insertion connects endpoints $p,r$ and the bare line connects $q,u$, the representative source is
\[
\mathscr A_i=-(ig_0)^4\int_{\mathcal D_i}d^4s\,
\prod_{j=1}^4v_{e_j}^{\mu_j}\frac{\operatorname{tr}(t^{a_1}t^{a_2}t^{a_3}t^{a_4})}{N_c}
\int\frac{d^D\ell}{(2\pi)^D}[D_0\Pi_1D_0]_{\mu_p\mu_r}^{a_pa_r}(\ell)e^{-i\ell\Phi_A}
\int\frac{d^Dk}{(2\pi)^D}[D_0]_{\mu_q\mu_u}^{a_qa_u}(k)e^{-ik\Phi_B}.
\]
The source's phase/import convention fixes the displayed overall sign; the internal transverse operator is an adopted physical root, not recomputed here. Each individual article identifies its insertion/bare endpoints and actual distance polynomials. The regulated Fourier transform gives precisely
\[
J_i=\int_{\mathcal D_i}d^4s\;
\frac{N_i}{2(1+e)}F_i^{-2+2e}d_{0,i}F_{0,i}^{-1+e}.
\tag{1}
\]
For the two actual insertion endpoints $s,t$, direct differentiation proves
\[
\frac{N}{2(1+e)}F^{-2+2e}
=dF^{-1+2e}+\frac{\partial_s\partial_tF^{2e}}{4e(1+e)}.
\tag{2}
\]
For an adjacent insertion $N=2dF$, so the complete density becomes $dF^{-1+2e}/(1+e)$. For an opposite insertion $F=D_z(s,t)=(1-s)(1-t)+zst$, $d=-(1+z)$ and $N=-2(1+z)D_z+(1-2e)z$. Keeping the derivative in (2) is essential; its rectangle faces are calculated below.

The raw colour is $C_F(C_F-C_A/2)B_1(e)$. The accepted connected-web projection removes the disconnected $C_F^2$ component and retains $a_1=-1/2$ in the $C_AC_F$ coordinate. This does not assert that the full raw trace vanishes in other coordinates. With
\[
B_1(e)=\frac{(10-6e)C_A-8(1-e)T_Fn_f}{3-2e},\quad
Z=e^{-\gamma_Ee}\Gamma(1-e),\quad R_3=\frac{Z^3}{32e(1-2e)},
\]
the complete connected scalar, before the regular $B_1$ multiplier, is
\[
\boxed{I_i=\frac{m_i a_1}{2}R_3\,
\frac{J_i(z)+z^{3e}J_i(1/z)}{1+z^{3e}}.}\tag{3}
\]
Here $m=4$ for01–06/09/10 and $m=2$ for07/08/11. The current native $C$ has the positive rational Fourier factor $1/[32e(1-2e)]$; its projected colour supplies the negative scalar prefactor. No historical family sign is appended. Euler compensation is consumed exactly once in $Z^3$.

On the starting strip adjacent powers have positive endpoint exponents. For a complete opposite tensor $N$ remains bounded and generally tends to the nonzero value $(1-2e)z$ at a zero corner. Its absolute density is therefore bounded by a constant times $(x+y)^{-2+2\Re e}$; the two-dimensional radial measure gives $r^{-1+2\Re e}dr$, integrable strictly for $\Re e>0$. Endpoint variables belonging to the two different lines are distinct. A full-cube product bound therefore justifies the particular regroupings used here on that strip. It does not prove holomorphy of the full Pi1 source at zero and does not authorize arbitrary permutation of unrelated prepared roots. The subsequent evaluated formulas provide the local meromorphic continuation.

## The evaluated weighted opposite kernel

Set $\Lambda=1+z$ and
\[
B_{\nu,\lambda}=\int_0^1ds\int_0^1dt\;t^\lambda D_z(s,t)^{\nu-1},\qquad
H(\nu,\lambda;z)=-\Lambda z^{-\nu/2}B_{\nu,\lambda}.
\tag{4}
\]
First integrate the variable $s$, which appears linearly:
\[
\int_0^1D_z(s,t)^{\nu-1}ds=
\frac{(zt)^\nu-(1-t)^\nu}{\nu[(1+z)t-1]}.
\]
The divided difference is removable at $t=1/(1+z)$; its numerator must remain combined. Substitute $\beta=1-t$, $r=z/(1+z)$, and define
\[
W=\log[\beta(1-\beta)],\quad L_\beta=\rho+\log\frac{1-\beta}{\beta},\quad
M_{k,j}=\int_0^1\frac{W^kL_\beta^j}{\beta-r}d\beta\quad(j\ge1).
\]
These are raw powers, without factorial normalization. Since $L_\beta(r)=0$, no principal-value prescription is introduced. For $\lambda=c\nu$, put $A_c=(c+1)W/2+c(L_\beta-\rho)/2$. The actual integrated difference gives the finite Taylor rule
\[
\boxed{[\nu^q]H(\nu,c\nu;z)=\frac1{(q+1)!}
\int_0^1\frac{(A_c+L_\beta/2)^{q+1}-(A_c-L_\beta/2)^{q+1}}{\beta-r}d\beta.}\tag{5}
\]
Expand this finite polynomial and replace $W^kL_\beta^j$ by the evaluated $M_{k,j}$. The [explicit moment appendix](appendices/pi1-moments.md) supplies every needed HLOG value and its precise adoption scope, including all endpoint periods. Thus (5) gives actual solved functions, not a newly named remaining integral. For example
\[
H(\nu,0)=M_{0,1}+\tfrac\nu2M_{1,1}
+\nu^2(M_{2,1}/8+M_{0,3}/24)+O(\nu^3).
\]
Crossing sends $\rho\mapsto-\rho$, $M_{k,j}\mapsto(-1)^{j+1}M_{k,j}$. Keep both channels until the physical crossing average is taken.

For an all-parameter evaluation of (4), Mellin–Barnes separation of the two terms in $D_z$ and the two actual Euler integrals give
\[
B_{\nu,\lambda}=\frac1{2\pi i}\int dw\,z^w
\frac{\Gamma(-w)\Gamma(1-\nu+w)\Gamma(1+w)\Gamma(1+\lambda+w)\Gamma(\nu-w)^2}
{\Gamma(1-\nu)\Gamma(1+\nu)\Gamma(1+\nu+\lambda)}.
\]
For the present nonnegative slopes choose small positive real $\nu,\lambda$ and $\nu-1<\Re w<0$; the Gamma poles are separated and both Euler integrals converge. Matching the defining standard kernel yields
\[
B_{\nu,\lambda}=
\frac{G_{3,3}^{3,3}\!\left(z\mid {\nu,0,-\lambda\atop0,\nu,\nu}\right)}
{\Gamma(1-\nu)\Gamma(1+\nu)\Gamma(1+\nu+\lambda)}.
\]
This Meijer G is distinct from the iterated HLOG $G$ of the appendix. The small-parameter integral is holomorphic: corner power bounds with a margin absorb every needed logarithmic derivative. Taylor extraction in (5) is therefore legitimate.

## Adjacent sources: the actual Beta differences

For01 put $x=1-s_{31}$, $y=1-s_{32}$, $a=s_{41}$, $b=s_{42}$. The domain is $0<y<x<1$, $0<a<b<1$, Jacobian1. Its complete kernel separates as
\[
J_{01}=\frac{z^{3e}}{1+e}
\left[\int_{0<y<x<1}x^{2e-1}y^{e-1}dxdy\right]
\left[\int_{0<a<b<1}a^{2e-1}b^{e-1}dadb\right].
\]
The first bracket is $1/(3e^2)$ and the second $1/(6e^2)$. For02 the insertion and bare exponents exchange; the two brackets become $1/(6e^2)$ and $1/(3e^2)$. Hence both actual roots evaluate to $z^{3e}/[18e^4(1+e)]$.

For05 take $(a,u,b,t)=(s_{21},s_{32},s_{31},s_{41})$ with $b<u$. Integrating $a,t,b$ gives
\[
J_{05}=\frac{z^e}{2e^3(1+e)}\int_0^1u^{2e-1}[1-(1-u)^e]du
=\frac{z^e(1-Q)}{4e^4(1+e)},\quad
Q=\frac{\Gamma(1+e)\Gamma(1+2e)}{\Gamma(1+3e)}.
\]
For06 the roles exchange: the remaining integral is $\int_0^1u^{e-1}[1-(1-u)^{2e}]du=(1-Q)/e$, while the preceding factor is $z^{2e}/[4e^3(1+e)]$. The different bare roots therefore have $z^e$ and $z^{2e}$ phases. Their complete crossing average, with their own equal $m,a_1$, is the same. The cancellation $1-Q=O(e^2)$ follows from these Beta differences, not a final coefficient fit.

## Opposite insertion: the P3 derivative boundary

For03 put $x=s_{21}$, $t=s_{41}$. The bare variables are $s_{31}\in[0,1]$ and $s_{42}\in[t,1]$; their actual integrations produce $z^e(1-t^e)/e^2$. The metric part is $-\Lambda(B_{2e,0}-B_{2e,e})$. For the derivative part, first integrating $x$ gives
\[
\int_0^1\partial_x\partial_tD_z^{2e}dx
=\partial_t[(zt)^{2e}-(1-t)^{2e}].
\]
The remaining full boundary mass is
\[
\int_0^1(1-t^e)\partial_t[(zt)^{2e}-(1-t)^{2e}]dt
=\frac{z^{2e}}3+1-Q.
\]
Consequently the fully reduced root is
\[
J_{03}=\frac{z^e}{e^2}\left[-\Lambda(B_{2e,0}-B_{2e,e})
+\frac{1+z^{2e}/3-Q}{4e(1+e)}\right].\tag{6}
\]
The10 map, in the04-style ordering of four symbols, is
\[
(s_{21},s_{41},s_{42},s_{31})_{03}
=(1-s_{41},1-s_{22},1-s_{21},1-s_{31})_{10}
\]
with the03 integrand evaluated at $1/z$. It preserves the insertion/bare labels, bounds and absolute Jacobian. Each line contributes its own power, giving $J_{10}(z)=z^{3e}J_{03}(1/z)$ and hence $I_{10}=I_{03}$.

## Adjacent insertion and bare opposite: P4

For04 retain the bare opposite endpoints $x=s_{21},t=s_{41}$. Integrating the full adjacent insertion over $s_{31}\in[0,1]$, $s_{42}\in[t,1]$ gives $z^{2e}(1-t^{2e})/[4e^2(1+e)]$. Thus
\[
J_{04}=-\frac{\Lambda z^{2e}}{4e^2(1+e)}(B_{e,0}-B_{e,2e}).\tag{7}
\]
For09 use the same displayed complement map as above with03 replaced by04 and10 by09. Its own adjacent lower primitive is $[1-(1-s_{22})^{2e}]/(2e)$; the free endpoint supplies $1/(2e)$. The result is $J_{09}=z^{3e}J_{04}(1/z)$, so the physical I values agree without identifying the bare roots.

## The coupled P5 rectangle and all eight corners

For07/08 set $s_1=s_{21}$, $s_2=s_{22}$, $t_1=s_{41}$, $t_2=s_{42}$, with $0<s_1<s_2<1$, $0<t_1<t_2<1$. Define
\[
T(a,b;z)=\Lambda^2z^{-(a+b)/2}\int_{\Omega_2}
D_z(s_1,t_1)^{a-1}D_z(s_2,t_2)^{b-1}d^4s.
\tag{8}
\]
The metric terms of the two sources are respectively $z^{3e/2}T(2e,e)$ and $z^{3e/2}T(e,2e)$. The derivative rectangle of07 is $[0,s_2]\times[0,t_2]$; that of08 is $[s_1,1]\times[t_1,1]$. Their signed corner masses are exactly

| Rectangle | $+F(\mathrm{upper},\mathrm{upper})^{2e}$ | $-F(\mathrm{lower},\mathrm{upper})^{2e}$ | $-F(\mathrm{upper},\mathrm{lower})^{2e}$ | $+F(\mathrm{lower},\mathrm{lower})^{2e}$ |
|---|---|---|---|---|
|07| $D_z(s_2,t_2)^{2e}$ | $(1-t_2)^{2e}$ | $(1-s_2)^{2e}$ | $1$ |
|08| $z^{2e}$ | $(zs_1)^{2e}$ | $(zt_1)^{2e}$ | $D_z(s_1,t_1)^{2e}$ |

Each mass is multiplied by $-\Lambda D_z^{e-1}/[4e(1+e)]$ on the remaining full square. Put $H_0^{(a)}=H(ae,0;z)$, $H_2=H(e,2e;z)$ and $\bar H_2=H(e,2e;1/z)$. Directly evaluating the eight square integrals gives
\[
B_{07}+B_{08}=\frac{z^{3e/2}}{2e(1+e)}A,
\quad A=H_0^{(3)}+\cosh(e\rho)(H_0^{(1)}-H_2-\bar H_2)
-\sinh(e\rho)(H_2-\bar H_2).\tag{9}
\]
For example the two $D^{2e}$ terms give $2z^{3e/2}H_0^{(3)}$ before the common $1/[4e(1+e)]$; the constants give $(1+z^{2e})z^{e/2}H_0^{(1)}$, and the four remaining powers give $-2z^{e/2}\bar H_2-2z^{5e/2}H_2$. This displays every finite boundary source.

Complementing both ordered chains reverses the small and large endpoints and proves $J_{08}=z^{3e}J_{07}(1/z)$. Equation(3) therefore uses $J_{07}+J_{08}$, giving
\[
I_{07}=I_{08}=-\frac{Z^3}{128e(1-2e)\cosh(3e\rho/2)}
\left[T(2e,e)+T(e,2e)+\frac{A}{2e(1+e)}\right].\tag{10}
\]
Only the first jet of the coupled body is needed. Holomorphy of(8) follows by the product of the two independent opposite-corner bounds on the full cube. Thus
\[
T(2e,e)+T(e,2e)=2\rho_0+3e\rho_1+O(e^2),\quad
\rho_q=[e^q]T(e,e),\quad \rho_0=M_{0,1}^2/6.\tag{11}
\]
This first-derivative identity is not an equality at generic regulator. The complete thirteen-term HLOG value of $\rho_1$ is printed in the [appendix](appendices/pi1-moments.md), with its actual kernel and adopted evaluation. It is not the logarithm $\rho$ and not an unknown remaining period.

## Separate opposite channels: P6

For11 reflect only the insertion endpoints $s_{11},s_{31}$ to turn $\widetilde D_z=z(1-s_{11})(1-s_{31})+s_{11}s_{31}$ into $D_z$. Its square has four derivative corner values $z^{2e},0,0,1$ on the starting strip. The integrated complete insertion and bare exchange are therefore
\[
z^eH_0^{(2)}+\frac{1+z^{2e}}{4e(1+e)},\qquad z^{e/2}H_0^{(1)}.
\]
The domain really is a product of these two squares. The complete paired root is
\[
\frac{J_{11}(z)+z^{3e}J_{11}(1/z)}{1+z^{3e}}
=\frac{H_0^{(1)}}{\cosh(3e\rho/2)}
\left[H_0^{(2)}+\frac{\cosh(e\rho)}{2e(1+e)}\right].\tag{12}
\]
The second term supplies the double pole; factorization does not remove the transverse boundary.

## Finite coefficient and complete colour convolution

Expand the actual paired root as $K_i(e,z)=\sum k_{i,r}e^r$. Its exact expression follows from(6)–(12), or the adjacent evaluations. All needed coefficients of $H$ and $T$ are now explicit evaluated functions. Let
\[
Z^3=\sum g_ne^n,\quad g_0=1,\quad
g_n=\frac3n\sum_{j=2}^n\zeta_jg_{n-j},\quad
\frac{Z^3}{1-2e}=\sum a_ne^n,\quad a_n=\sum_{j=0}^n2^{n-j}g_j.
\]
With $A_i=m_i a_1/64$, the complete coefficient rule and its finite specialization are
\[
c_{i,p}=A_i\sum_{n\ge0}a_nk_{i,p+1-n},\qquad
c_{i,0}=A_i\sum_{n\ge0}a_nk_{i,1-n}.
\tag{13}
\]
The sums are finite, using the actual lower bound of $K_i$. This incorporates the source $1/e$ and every necessary positive-order subintegral. A safe budget is:

|Route|actual I pole|inputs needed for finite I|
|---|---:|---|
|P1|5|Gamma/rational through5|
|P2|3|$Q$ through5, remaining regular factors through5 conservatively|
|P3|4|$H$ through3; $Q$ through4; Gamma/rational through4|
|P4|2|$H$ through3; remaining regular factors through2|
|P5|1|$H$ through2 and mixed $T$ sum through1|
|P6|2|$H_0$ through2; remaining regular factors through2|

For the two connected colour channels write $\phi_A=(10-6e)/(3-2e)$ and $\phi_f=-8(1-e)/(3-2e)$. Their coefficients are $\phi_{A,0}=10/3$, $\phi_{f,0}=-8/3$, and for $n\ge1$,
\[
\phi_{A,n}=\tfrac29(\tfrac23)^{n-1},\qquad
\phi_{f,n}=\tfrac89(\tfrac23)^{n-1}.
\]
The actual finite colour coefficient is $\sum_{j=0}^{p_i}\phi_{X,j}c_{i,-j}$ for scalar pole order $p_i$; general order $p$ is $\sum_{j\ge0}\phi_{X,j}c_{i,p-j}$. This is why $B_1(0)c_{i,0}$ is insufficient. The other two connected channels are zero here. Each individual article prints all seven scalar coefficients and both full coloured finite coefficients; result.json retains the complete coloured series.

The actual Python program expands these six already derived route formulas, adopts the evaluated $M/\rho_q$ data, then performs(13) and the regular colour convolution. It does not itself independently derive the six endpoint/face reductions or reintegrate the archived functions. The current reading chain exposes those mathematical connections. Original computational evidence is preserved in [the runtime](runtime/README.md), particularly its snapshots/workflow/checks/self-energy-pi1-* producers; they retain old paths and are not invoked as current writers. These research calculations do not supply an arbitrary complex-kinematics or integral-level formal proof.
