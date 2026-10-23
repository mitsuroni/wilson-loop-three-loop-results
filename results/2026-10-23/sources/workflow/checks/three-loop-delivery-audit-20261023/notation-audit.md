# Notation and conversion audit for the three loop delivery

Date: 23 October 2026. Task: WILSON-5.15. Source workspace: `D:/wilson-agent`, read only. This note audits the 18 October delivery and the subsequent source definitions. It records existing formulas and their stated scope; no CAS calculation, new identity proof, source change or result selection was performed.

## Main finding

The earlier package already defines the basic Hlog, M/N, tau, E and Omega conventions correctly. Its dictionary does not yet cover the newer F/Q/R hierarchy, Euler operations, Green operations, their coefficient factorials, or several incompatible uses of the same symbol. A richer package should attach an explicit symbol dictionary to each displayed expression and publish a conversion table with an order and a domain. It must not describe all 149 three-loop orbit results as belonging to M/N, or treat the newer pure-exchange hierarchy as a reduction of every other family.

The principal distinction is between four types of statement:

| Type | Example | Meaning |
|---|---|---|
| Literal notation alias | `M01` and $M_{0,1}$ in the registered moment convention | The same defined object |
| Exact function identity at fixed weight | $\tau^+_{0,2}=M_{2,1}/8-M_{0,1}^2/24$ | A conversion of an explicitly defined coefficient function; not an all-order source reduction |
| Exact all-order source relation | The Euler lift of $F_\epsilon[t^\epsilon s^\epsilon]$ | Equality of parameter-dependent integrals with the stated domain and continuation |
| Open reduction | Expressing every extended integral using a minimal finite M/N-only set | Not established by a coefficient conversion, numerical agreement or a constructive all-order representation |

## 1. Hlog, GPL and endpoint conventions

The project convention is

\[
H_{a_1,\ldots,a_m}(z)=\operatorname{Hlog}(z,[a_1,\ldots,a_m])
=G(a_1,\ldots,a_m;z),
\]
\[
G(a_1,\ldots,a_m;z)=\int_0^z\frac{dt}{t-a_1}
G(a_2,\ldots,a_m;t),\qquad G(\varnothing;z)=1,
\]
\[
G(0^m;z)=\frac{\log^m z}{m!}.
\]

The **first listed letter is the outermost kernel**. The argument is first in `Hlog(z,word)` and last in the mathematical GPL notation. In particular,

\[
G(-1;z)=\log(1+z),\quad G(1;z)=\log(1-z)\quad(0<z<1),
\qquad G(0,-1;z)=-\operatorname{Li}_2(-z).
\]

The letter $1$ therefore has the opposite kernel sign to the often-used harmonic-polylogarithm kernel $dt/(1-t)$. The origin is tangentially regularized where required. General continuation and interior singularities retain the source path and prescription. These definitions are explicit in [the Gamma1-03 source](/D:/wilson-agent/benchmarks/self-energy/gamma1/workspaces/3e-Gamma1-03/derivation.md:233), [the leg-insertion source](/D:/wilson-agent/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-07/derivation.md:209) and [the weight-nine Hlog reader](/D:/wilson-agent/benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-eps3-hlog.md:6).

Several coefficient formats require a more specific reader:

| Representation | Required interpretation |
|---|---|
| `hlog_terms` with `word`, `coefficient`, `zeta_powers` | Sum the exact rational coefficient times the declared period monomial times Hlog at the declared argument. An empty word is 1. Read the period order from that schema; do not import another format's order. |
| Omega6 `period_powers` | The registered order is $[\zeta_2,\zeta_3,\zeta_5,\log2,A_4,A_6]$. Letters include rational $-1/2$ and $-2$ as exact numbers. |
| Higher-order `GPL_words` with product rows | A word-table index selects a GPL, not a letter. Rows may be products of several GPLs. Preserve argument 1 and the substitutions $Z=-z$, $IZ=-1/z$. |
| Algebraic-Chen or `RegInf` data | Preserve the basepoint, algebraic-root embedding, continuation path, endpoint variable, and the logarithm whose constant coefficient is extracted. It is not an ordinary GPL value obtained by dropping endpoint terms. |

For example, the exact higher-order row

```text
[c,[i1,...,ik],p,q,r,s]
```

in [the epsilon-four reader](/D:/wilson-agent/benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-eps4-gpl.md:5) means

\[
c\prod_jG(\mathrm{GPL\_words}[i_j];1)
\log^p z\,\log^q\!\frac{z}{(1+z)^2}\,\zeta_2^r\zeta_3^s,
\]

with zero-based indices. This is not already a canonical Hlog-at-$z$ expression. Likewise, the source defines $\zeta(5,3)$ with descending summation indices, whereas the Omega source identifies its alternating periods with HyperInt's small-index-first `zeta[1,-3]` and `zeta[1,-5]`. A multiple-zeta label needs its own summation convention.

The graph-25 source is a concrete endpoint warning: its documented initial substitution of `RegInf` by a GPL at 1 was rejected by local comparisons. The corrected construction retains a regularized endpoint and original logarithm. See [the basepoint correction](/D:/wilson-agent/benchmarks/double-tgv/4e-3g+3g/workspaces/4e-3g+3g-25/derivation.md:897) and [the common-endpoint definition](/D:/wilson-agent/benchmarks/double-tgv/4e-3g+3g/workspaces/4e-3g+3g-25/derivation.md:671).

## 2. M and N indices, factorials and generators

Use $z>0$ for the primary real definition and write

\[
\ell=\log z,\quad r_z=\frac{z}{1+z},\quad
W_\beta=\log[\beta(1-\beta)],\quad
R_\beta=\log\frac{1-\beta}{\beta},\quad L_\beta=\ell+R_\beta.
\]
\[
M_{i,j}(z)=\int_0^1\frac{W_\beta^iL_\beta^j}{\beta-r_z}\,d\beta,
\qquad
N_{i,j}(z)=\int_0^1\frac{W_\beta^iR_\beta L_\beta^j}{\beta-r_z}\,d\beta,
\qquad i\ge0,\ j\ge1.
\]

These are **raw logarithmic moments, without factorial denominators**. `M01`, `M_0_1` and $M_{0,1}$ are aliases only when the file adopts this dictionary. The first index counts $W_\beta$, the second counts $L_\beta$; no index shift has been found between the old and current definitions. The point $\beta=r_z$ is removable because $L_\beta(r_z)=0$. The current base-language note explicitly disallows introducing $N_{i,0}$ merely to absorb a bare logarithm.

The exact relations include

\[
N_{i,j}=M_{i,j+1}-\ell M_{i,j},\qquad
M_{0,1}=-\frac{\ell^2+\pi^2}{2},\qquad
M_{0,3}=-M_{0,1}^2.
\]

The common generating source carries the factorials:

\[
\mathscr U(a,b;z)=\int_0^1
\frac{e^{aW_\beta}(e^{bL_\beta}-1)}{\beta-r_z}\,d\beta
=\sum_{i\ge0,j\ge1}\frac{a^ib^j}{i!j!}M_{i,j}(z).
\]

Consequently $\partial_a^i\partial_b^j\mathscr U|_{0,0}=M_{i,j}$, while coefficient extraction gives $[a^ib^j]\mathscr U=M_{i,j}/(i!j!)$. The two-source combinations are

\[
V_M(a,b;z)=\int_0^1\frac{e^{aW_\beta}}{\beta-r_z}
\frac{\sinh(bL_\beta)}b\,d\beta
=\sum_{i,k\ge0}\frac{a^ib^{2k}}{i!(2k+1)!}M_{i,2k+1},
\]
\[
V_N(a,b;z)=\int_0^1\frac{e^{aW_\beta}R_\beta}{\beta-r_z}
\frac{\cosh(bL_\beta)-1}{b^2}\,d\beta
=\sum_{i,k\ge0}\frac{a^ib^{2k}}{i!(2k+2)!}N_{i,2k+2}.
\]

The old one-parameter functions are exact slices:

\[
V_M^{\rm old}(\epsilon,z)=H_0(\epsilon,z)=V_M(\epsilon/2,\epsilon/2;z),
\qquad V_N^{\rm old}(\epsilon,z)=V_N(\epsilon,\epsilon;z).
\]

These are all-order identities near the source parameter origin with continuation as stated in the source. Replacing epsilon in a generator does not rescale any external Gamma, Euler or crossing factor. Sources: [the original source dictionary](/D:/wilson-agent/benchmarks/2ne/families/single-opposite/three-loop.md:43), [the two-source catalogue](/D:/wilson-agent/benchmarks/renormalization/three-loop/slide-data/generating-structures.md:7), and [the current base language](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/base-language/derivation.md:117).

## 3. Explicit registered conversions

### M/N to Hlog in the current Af finite bridge

The fixed-moment expansion files below give actual Hlog-at-$z$ formulas with letters $\{0,-1\}$ and explicit zeta values. They are not just definitions of M/N. The [current matter check](/D:/wilson-agent/workflow/checks/matter-rg-20261017/check.py:284) loads all six files and forms $N_{1,2}=M_{1,3}-\log z\,M_{1,2}$. The [saved Af compact back-substitution](/D:/wilson-agent/workflow/checks/w3-ca-nf/compact.mpl:16) tests its finite and shape M/N expressions against the GPL expressions after applying the relevant moment substitutions. This audit reads those existing artifacts; it does not rerun the check.

| Moment | Explicit M-to-Hlog expansion |
|---|---|
| $M_{0,1}$ | [M_0_1.txt](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt): $-3\zeta_2-G(0,0;z)$ |
| $M_{1,1}$ | [M_1_1.txt](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt) |
| $M_{2,1}$ | [M_2_1.txt](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt) |
| $M_{3,1}$ | [M_3_1.txt](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_3_1.txt) |
| $M_{1,2}$ | [M_1_2.txt](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_2.txt) |
| $M_{1,3}$ | [M_1_3.txt](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_3.txt) |

Each file is a raw moment equality with no extra factorial. The [generation source](/D:/wilson-agent/workflow/checks/single-opposite/runs/20260916T090538777022Z/compute_moments.mpl:28) uses the complementary integration coordinate $x=1-\beta$ and $a=1/(1+z)$; that change of variables preserves the M definition in Section 2. These fixed-index expansions and $N_{i,j}=M_{i,j+1}-\log z\,M_{i,j}$ give a concrete reading route from an M/N formula to Hlog. They do not establish that an arbitrary Hlog expression, extended block or diagram can be reduced back to M/N.

### Other registered conversions and their scope

The table below contains identities, not numerical fits. Each finite-weight row converts the specified coefficient function only. The complete historical conversion list contains the registered tau/eta, Sigma and double-kernel rho/chi names; it is preferable to a rule guessing a definition from a symbol's spelling.

| Registered symbol or expression | Exact conversion | Scope and source |
|---|---|---|
| `tau0q0` | $\tau^+_{0,0}=M_{0,1}$ | Fixed coefficient function; [conversion table](/D:/wilson-agent/benchmarks/2ne/theory/weight-six/conversions.md:133) |
| `tau0q1` | $\tau^+_{0,1}=M_{1,1}/2$ | Same table, line 141 |
| `tau0q2` | $\tau^+_{0,2}=M_{2,1}/8-M_{0,1}^2/24$ | Same table, line 149 |
| `tau0q3` | $\tau^+_{0,3}=(M_{1,3}+M_{3,1})/48$ | Same table, line 157 |
| `tau0q4` | $\tau^+_{0,4}=M_{0,1}^3/1440-\zeta_2M_{0,1}^2/240+M_{2,3}/192+M_{4,1}/384$ | Same table, line 165 |
| `tauDq2` | $\tau^+_{D,2}=-M_{0,1}^2/8+9M_{2,1}/8$ | Same table, line 374 |
| `tauDq3` | $\tau^+_{D,3}=-4\zeta_3M_{0,1}+2\zeta_2M_{1,1}+M_{1,3}/48+9M_{3,1}/16-N_{1,2}/2$ | Same table, line 381 |
| $N_{0,2}$ | $M_{0,1}^2/3+4\zeta_2M_{0,1}$ | Same table, line 45 |
| $\rho_1=[\epsilon]Q$ in the four-endpoint comparison | $M_{0,1}M_{1,1}/2-N_{1,2}+4\zeta_2M_{1,1}-8\zeta_3M_{0,1}$ | [Block summary](/D:/wilson-agent/benchmarks/renormalization/three-loop/slide-data/block-summary.md:59); finite coefficient, not arbitrary unequal-parameter ordered T |
| Old weight-six $E(z)$ | $\widehat Q^{\ell-}_{010}(z)$ | Exact equality of the two integral definitions; [weight-six blocks, Eq. 10](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/weight-six-blocks/derivation.md:135) |
| $M_{0,1}$ | $-\widehat F^+_{000}$ | [Weight-six blocks, Eq. 6](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/weight-six-blocks/derivation.md:89) |
| $M_{1,1}$ | $-2\widehat F^+_{001}$ | Same equation |
| $M_{2,1}$ | $-8\widehat F^+_{002}+7(\widehat F^+_{000})^2/3-6\zeta_2\widehat F^+_{000}$ | Same equation |
| Old $\tau^+_{D,4}$ in the new single-kernel moments | Explicit ten-term polynomial in $\widehat F^+_{000},\widehat F^+_{001},\widehat F^+_{002},\widehat F^+_{004},\widehat F^+_{103},\widehat F^+_{112},\widehat F^+_{211}$ and periods | [Complete exact formula](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/weight-six-blocks/derivation.md:311); weight six only |
| New Euler $\mathcal E_4$ versus old $\tau^+_{D,4}$ | $\mathcal E_4\equiv-\tau^+_{D,4}\pmod{\mathcal A_{MN}}$ | [Euler readout](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/simple-kernels/all-orders.md:256); **not** equality without the M/N correction |

For the normalized unweighted ordered double kernel, the established all-order generator identity is

\[
Q=\frac12V_M(\epsilon/2,\epsilon/2;z)^2-2V_N(\epsilon,\epsilon;z)
+\frac{2\chi(\epsilon)}\epsilon V_M(\epsilon,\epsilon;z),
\quad \chi(\epsilon)=\gamma_E+\psi(1+2\epsilon).
\]

This identity generates the rho1 conversion above. It does not replace the unequal-parameter source $T(2\epsilon,\epsilon)+T(\epsilon,2\epsilon)$ by the diagonal $T(\epsilon,\epsilon)$.

## 4. Definitions of the retained extended coefficient blocks

Set $D_z(t,s)=(1-t)(1-s)+zts$. The old weight-six function is

\[
\tau^+_{D,4}=[\epsilon^4]\frac{H_D(\epsilon,z)+H_D(\epsilon,1/z)}2,
\]
\[
H_D(\epsilon,z)=-(1+z)z^{-\epsilon/2}
\int_0^1\!\int_0^1t^\epsilon s^\epsilon D_z(t,s)^{\epsilon-1}\,dt\,ds.
\]

The phase $z^{-\epsilon/2}$ and the coefficient-extraction factorial are part of the definition. With $\Omega_2=\{0<t_1<t_2<1,0<s_1<s_2<1\}$,

\[
E(z)=\frac{\ell(1+z)^2}{2}\int_{\Omega_2}
\frac{\log(s_2/(1-s_1))}{D_z(t_1,s_1)D_z(t_2,s_2)}\,d^2t\,d^2s.
\]

Both are crossing even. The derivative block is $\Theta\tau^+_{D,4}=z\,d\tau^+_{D,4}/dz$, with any coefficient multiplying it retained. Definitions are in [the unified weight-six basis, section 2](/D:/wilson-agent/benchmarks/2ne/theory/unified-weight-six.md:31). The newer exact single-kernel conversion of tau does not eliminate E or prove all-order M/N closure.

The T4 function has the exact definition

\[
T_6(z)=\int_0^z\frac{\Phi_5(t)}{t+1/2}\,dt,\quad T_6(0)=0,
\qquad \Omega_6(z)=T_6(z)+T_6(1/z).
\]

The full weight-five source $\Phi_5$, the evaluated word tables and the two alternating periods are specified in [the T4 basis definition](/D:/wilson-agent/benchmarks/tgv/families/t4-opposite/basis-functions.md:5). The alphabet is $\{0,-1,-1/2,-2\}$. In particular,

\[
A_4=\sum_{n\ge1}\frac{(-1)^nH_{n-1}}{n^3},\qquad
A_6=\sum_{n\ge1}\frac{(-1)^nH_{n-1}}{n^5}.
\]

These are fixed periods, not free boundary constants. The documented new-direction claim is in the stated finite-weight coordinate space, not a global function-minimality theorem.

## 5. Current F, Q and R hierarchy

### Bare kernels and normalized graph blocks

The general single kernel is

\[
\widetilde F_\lambda(\alpha,\beta;\gamma,\delta;z)
=\int_0^1t^\alpha(1-t)^\beta s^\gamma(1-s)^\delta
D_z(t,s)^{\lambda-1}\,dt\,ds.
\]

The parameters are actual exponents. In particular,

\[
F_\lambda(A,B;z)=\widetilde F_\lambda(A,0;B,0;z),\qquad
F^{(c)}_{ab}=\widetilde F_{c\epsilon}(a\epsilon,0;b\epsilon,0;z).
\]

$c$ is a kernel-exponent multiplier, not the loop count. The more compact notation $F_\lambda[w]$ puts the full weight $w(t,s)$ inside the integral. None of these bare kernels includes $g_\epsilon=-(1+z)z^{-\epsilon/2}$. Definitions and the full pure-exchange J/I map are in [the unified tables, Eqs. 1-3](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/unified-tables/derivation.md:11).

For $m=2,3$, define the ordered domain $\Omega_m$ by separate ascending t and s chains. Bare Q and R are products of two or three D kernels on that domain; their pairing permutation is essential. In the unified-table notation,

\[
\mathcal Q^\pi_{ab}=\int_{\Omega_2}s_1^{a\epsilon}s_2^{b\epsilon}
\prod_{i=1}^2D_z(t_i,s_{\pi(i)})^{\epsilon-1}\,d^2t\,d^2s,
\]
\[
\mathcal R_\pi=\int_{\Omega_3}\prod_{i=1}^3
D_z(t_i,s_{\pi(i)})^{\epsilon-1}\,d^3t\,d^3s.
\]

For Q, 12/id is same-order pairing and 21/sw is crossed pairing. A superscript c denotes complementary endpoint weights, not complex conjugation. Their reflected endpoint indices swap when the variables are put back into ascending order. See [the unified tables, Eqs. 12-15](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/unified-tables/derivation.md:201).

**Notation changes between papers.** The unified tables use plain $Q^\pi_{ab}=g_\epsilon^2\mathcal Q^\pi_{ab}$ and $R_\pi=g_\epsilon^3\mathcal R_\pi$ for the older normalized kernels. The later triple-R paper instead calls the bare kernels $Q_\pi[w]$ and $R_\pi$, and defines positive graph blocks

\[
\mathcal R_{20}=(1+z)^3z^{-3\epsilon/2}R_{123},\qquad
\mathcal R_{21}=\frac{(1+z)^3z^{-3\epsilon/2}}2(R_{132}+R_{213}).
\]

Its actual integral readout is

\[
I_{6e-20}=\frac{Z^3\mathcal R_{20}}{16c_3},\qquad
I_{6e-21}=I_{6e-22}=\frac{Z^3\mathcal R_{21}}{32c_3},
\quad Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad
c_3=\cosh(3\epsilon\ell/2).
\]

The positive sign is explicitly opposite to earlier negatively centered objects. Thus a label such as R123 must be accompanied by its source definition and prefactor. See [triple R, Eqs. 1-2 and its readout table](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/triple-r/derivation.md:11).

For the newer Q graph blocks, the subscript is a graph number, not the pairing:

\[
\mathcal Q_{10}=\epsilon^{-1}\operatorname{Even}_z
[(1+z)^2z^{-\epsilon/2}Q_{12}[1-s_2^\epsilon]],
\]

with analogous expressions for graphs 11 and 12. Their full definition is in [the negative-one boundary note, Eq. 6](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/minus-one/derivation.md:61). A hatted $\widehat{\mathcal Q}_g$ is the seed before Even projection. R's first derivative requires its **Odd** projection, so substituting the already-even physical Q block loses necessary information; [triple R, Eq. 5](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/triple-r/derivation.md:90) states this explicitly.

### Single-kernel coefficient moments

The newer coefficient moments include factorials in their definition:

\[
\widehat F_{pqr}(z)=\frac{1+z}{p!q!r!}\int_0^1\!\int_0^1
\frac{\log^pt\,\log^qs\,\log^rD_z}{D_z}\,dt\,ds.
\]
\[
\widehat F^+_{pqr}=\frac{\widehat F_{pqr}(z)+\widehat F_{pqr}(1/z)}2,
\qquad
\widehat F^-_{pqr}=\frac{\widehat F_{pqr}(z)-\widehat F_{pqr}(1/z)}{2\ell}.
\]

Both are even under inversion; the superscript minus does **not** make the latter an odd function because it also divides by ell. At z=1 use the removable limit. The denominator need not cancel within a finite polynomial Hlog algebra. See [the weight-six definitions](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/weight-six-blocks/derivation.md:11). These three indices are endpoint/kernel log powers, not M's two indices or a Laurent power.

### Whole-source connections and Euler moments

The single-weight F family has the exact connection

\[
-\lambda(1+z)z^{(a-\lambda)/2}F_\lambda[s^a]
=\mathscr U\left(\frac{a+\lambda}2,\frac{a+\lambda}2;z\right)
-\mathscr U\left(\frac{a+\lambda}2,\frac{a-\lambda}2;z\right).
\]

This is an all-order relation, not a finite-weight fit. The complete difference is taken before a degenerate-parameter limit. Source: [base language, Eq. 7](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/base-language/derivation.md:165).

For the same-direction double-weight slice, define

\[
\mathcal T_\rho[f](z)=\int_0^1\frac{x^\rho\,dx}{(1-x)^{1+\rho}}
\left[\frac{1+z}{1+zx}f(zx)-f(z)\right].
\]

Then the recorded all-order lift is

\[
-(1+z)F_\epsilon[t^\epsilon s^\epsilon]
=\frac{g_1}{\epsilon}\mathscr U(\epsilon,\epsilon;z)
-g_2\mathcal T_\epsilon[\mathscr U(\epsilon,\epsilon;\cdot)](z),
\]
\[
g_1=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
g_2=\frac{\Gamma(1+\epsilon)}{\Gamma(1+2\epsilon)\Gamma(1-\epsilon)}.
\]

This uses an explicit Euler operation on the old source; it does not assert a finite polynomial in untransformed M/N. The coefficient operation is

\[
\mathcal T_k^+[f]=[\rho^k]\operatorname{Even}_z
\{z^{\rho/2}\mathcal T_\rho[f]\},\qquad
\mathbb M^+_{ijq}=\frac{\mathcal T_{q-1}^+[M_{ij}]}{i!j!}.
\]

$\mathcal T_k^+$ already contains $1/k!$, hence the three-index moment includes the corresponding $1/(q-1)!$. A continuous subscript rho and a coefficient index k are different. The new $\mathcal E_n$ is a specified finite sum of these moments, not the old double-kernel E. See [the all-order Euler source, Eqs. 2-8](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/simple-kernels/all-orders.md:24). In the latest notation $\mathcal E_4=\sum_{i+j+q=5}\mathbb M^+_{ijq}$; the index shift does not create an additional function.

### Green operations and the negative-one branch

The latest basepoint is -1, reached on the specified upper or lower analytic continuation and connected to positive z without passing through 0. Set

\[
\Lambda(z)=2\operatorname{Log}(1+z)-\operatorname{Log}z-\log4,
\qquad
\mathcal K_c[f](z)=e^{-c\Lambda(z)}\int_{-1}^z
\frac{dy}{y}e^{c\Lambda(y)}f(y).
\]

The earlier basepoint-1 operator has the same letter but different boundary data. The newer source fixes the required value and derivative boundaries from the original contour; inversion exchanges upper and lower continuation branches. One cannot move Even through a single fixed branch as if it were the earlier real basepoint-1 operation.

The coefficient operators are

\[
G_k[f]=\frac1{k!}\int_{-1}^z\frac{dy}{y}[\Lambda(y)-\Lambda(z)]^kf(y),
\]
\[
L_k[f]=\frac1{k!}\int_{-1}^z\frac{dy}{y}
[\operatorname{Log}z-\operatorname{Log}y][\Lambda(y)-\Lambda(z)]^kf(y).
\]

Thus $\mathcal K_{a\epsilon}=\sum(a\epsilon)^kG_k$ and $\mathcal K_{a\epsilon}^2=\sum(a\epsilon)^kL_k$. These G symbols are operators, not GPLs. The exact parameter relation

\[
\mathcal K_c[\chi\mathcal K_a[f]]
=\frac{\mathcal K_a[f]-\mathcal K_c[f]}{c-a},\qquad
\chi=\frac{z-1}{z+1},
\]

uses the same branch and zero-boundary conditions; its coincident-parameter limit is $-\partial_c\mathcal K_c[f]$. Sources: [negative-one boundary, Eqs. 8-10](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/minus-one/derivation.md:92) and [higher weights, Eqs. 1-5](/D:/wilson-agent/benchmarks/theory/integral-geometry-20261017/adjacent-tower/higher-weights/derivation.md:25).

The recent finite relation

\[
[\epsilon^0](I_{6e-21}+I_{6e-20}/2)=-M_{0,1}^3/1152
\]

must carry its epsilon-zero label. The higher-weight note reports a nonzero correction to the natural all-order extension. Also distinguish it from the source-level combination $\mathcal R_{21}+\mathcal R_{20}/2$, whose **full physical** counterpart has the coefficient $I_{6e-21}+I_{6e-20}/4$ because the two I prefactors differ.

## 6. Same-name conflicts to resolve in the delivery

The following conflicts were found in the audited source set. Use source-qualified identifiers internally and show the local definition with the formula.

| Printed name | Different meanings requiring separation |
|---|---|
| $L$ | Mean renormalization scale log; log z in some old tables; $L_\beta=\ell+R_\beta$ inside M/N; the new integral operator $\mathcal L=L_0$ |
| $r$, $\rho$ | $z/(1+z)$ in M/N; $1/(1+z)$ in the TGV-Pi1 Hlog argument; log z; a Gamma ratio in the unified table; $z/(1+z)^2$ in a high-order GPL reader; an Euler-source parameter |
| $b(\epsilon)$ | Negative signed matter coefficient in the two-loop convention; positive magnitude in the three-loop Af convention. Prefer $B_f=-8(1-\epsilon)/(3-2\epsilon)$ and record the mapping below. Neither is a beta-function coefficient $b_0$ or $b_1$. |
| $\lambda$ | Coupling $\alpha_s/\pi$ in RG; kernel exponent in $F_\lambda$ |
| $G$, $G_k$ | GPL; Euler-compensated Gamma factor in older normalization notes; Green coefficient operator |
| $H_0$ | One-parameter M generator; an Hlog zero-letter function in other contexts |
| $W$ | Wilson loop; beta-log weight; ordered endpoint function $W_\epsilon(s)$ |
| $D$ | Spacetime dimension; kernel polynomial $D_z$; the label in tau_D; differential operator $\mathcal D_\epsilon$ |
| $C$ | Amplitude prefactor; endpoint-weight alias in the single-kernel table; Euler split operator $\mathcal C_\rho$ |
| $U$ | Full M source $\mathscr U$; bare two-endpoint F alias; scalar Gamma factor in self-energy conventions |
| $V$ | Complementary double-weight F alias; normalized unweighted F source; $V_M,V_N$ generators |
| $E$ | Old weight-six double-kernel function; unified-table alias $F_{2\epsilon}[s^\epsilon]$; temporary even source $E_\epsilon$ in triple-R; Euler generator $\mathcal E$ and its $\mathcal E_n$ coefficients |
| $Q$, $\mathcal Q$ | Bare ordered kernel; g-normalized kernel; one-parameter diagonal generator; regular graph block with graph-number subscript |
| $R$, $\mathcal R$ | M/N beta-log ratio; crossing denominator; branch prefactor; bare or g-normalized triple kernel; positive centered graph block; regularized endpoint operation in graph 25 |
| $q$, $r$, $g_2$ | Gamma ratios in one note; normalized ordered kernel q in another. In particular g2 in triple-R is not g2 in the Euler lift. |
| $\chi$ | $\gamma_E+\psi(1+2\epsilon)$ in the old Q formula; $(z-1)/(z+1)$ in the Green hierarchy; registered historical double-kernel coefficient names |
| $\tau$ | Numerous old source coefficients; the specifically retained weight-six tau_D,4+; not every tau label denotes a new function |
| $\mathcal K_c$ | Basepoint-1 Green operator; basepoint-minus-one Green operator with different continuation/boundary data |
| Superscript minus or a bar | An odd projection in some sources; an even reflection difference divided by ell for F-hat-minus; complementary weight or crossed pairing rather than conjugation in Q notes |
| $\Omega_2$, $\Omega_6$ | An ordered integration domain; a named weight-six function |

The Pi1 coefficient has one unambiguous color decomposition if the signed coefficient is named explicitly:

\[
a(\epsilon)=\frac{10-6\epsilon}{3-2\epsilon},\qquad
B_f(\epsilon)=-\frac{8(1-\epsilon)}{3-2\epsilon},\qquad
B_1=a C_A+B_f T,\qquad T=T_Fn_f.
\]

The [two-loop convention](/D:/wilson-agent/benchmarks/renormalization/two-loop/preread/CONVENTIONS.md:29) uses $b_{\mathrm{2L}}=B_f$ and its matter sum is $b_{\mathrm{2L}}I_{\Pi}$. The [three-loop Af derivation](/D:/wilson-agent/benchmarks/renormalization/three-loop/workspaces/CF-CA-TF-nf/derivation.md:26) uses $b_{\mathrm{3L}}=-B_f$ and writes $B_1=aC_A-b_{\mathrm{3L}}T$. The common signed decomposition is independently explicit in the [two-loop source audit](/D:/wilson-agent/benchmarks/renormalization/two-loop/source-audit.md:72); the [matter-RG note](/D:/wilson-agent/workflow/checks/matter-rg-20261017/README.md:24) already names $B_f$ with this sign. These are equivalent conventions, not inconsistent physics. Thus a three-loop weight $-b_{\mathrm{3L}}$ becomes $B_f$, and $-2ab_{\mathrm{3L}}$ becomes $2a B_f$. This substitution does not supply or remove a separate family-specific native/historical source-to-physical sign bridge.

The parent [assembly audit](assembly-audit.md) explicitly defines the positive $b_{\mathrm{3L}}$ and is consistent with that choice; the [two-loop audit](two-loop-workflow.md) retains the negative $b_{\mathrm{2L}}$. Their formulas must not share an unqualified `b` identifier in a generated dictionary. Similarly, the TGV-Pi1 argument $r=1/(1+z)$ is fixed in its [local definition](/D:/wilson-agent/benchmarks/self-energy/tgv-pi1/workspaces/3e-3g+Pi1-07/derivation.md:120), and must not be replaced by the complementary M/N coordinate merely because both are printed as r.

The parent overview and assembly audit were also reviewed read-only for notation claims. Their separation of diagram coefficients, color sums, finite functions and limited validation is sound. The distinction $\mathcal Z_{3,c}=\operatorname{Pole}_{\epsilon}\widetilde w_{3,c}$ must remain explicit: it is the relevant coefficient of **$\log Z_4$**, as stored under `log_Z4_coefficients` in the [selected Af result](/D:/wilson-agent/benchmarks/renormalization/three-loop/results/CF-CA-TF-nf/r001.json:55), in the displayed subtraction convention. The multiplicative $Z_4$ has additional lower-order products. The parent assembly audit has adopted this label clarification.

## 7. Scope of the latest relations

The new F/Q/R tables concern the **47 pure-exchange inputs** and their 39 three-loop members. They do not provide a general conversion for the 149-orbit collection. Within that pure-exchange scope, the base language supplies an all-order M/N coefficient rule for its stated 36-input subset. The remaining branches use specified Euler and Green operations, whose current source and boundary connections are recorded through 23 October. A finite-operation generating chain does not establish a smallest constant-coefficient basis or reduce every operation to a finite list of standard functions.

The weight-six F-hat conversions are exact at the displayed weight. They do not assert that a fixed finite set of F-hat moments suffices at every order. The old E-to-Q-hat identity is a definition-level equality; the Euler-E4-to-tau relation is only a quotient-space statement unless its M/N correction is also included. The high-order source includes actual epsilon and epsilon-squared formulas with explicit operation rules, but those are not automatically a new selected Laurent revision for every diagram.

For the remaining families retain the current selected expression and its actual classification: M/N, an extended block, a complete Hlog expression, a finite GPL product, algebraic-Chen data, or an explicitly unknown coefficient. Gamma1 gluon/quark channels and Pi2 or chain color channels remain separate. Neither a shared symbol nor a common alphabet proves a full compact reduction.

## 8. Minimum record needed to read one diagram unambiguously

Each displayed Laurent expression should have the following attached data, preferably filled from the selected result and checked against its source rather than retyped:

1. Orbit ID, immutable revision, source file and exact expression/JSON pointer; separately identify the color channel.
2. The mathematical object being expanded: raw J, normalized I, a kernel/source block, a bare exponent contribution, or the finite renormalized color coefficient.
3. The convention $I_g(z,\epsilon)=\sum_q\epsilon^qI_{g,q}(z)$ and the actual known q range. A coefficient is not an epsilon derivative; no extra q! is implied.
4. Variable substitutions and the real or complex domain. Preserve z versus r versus z1, the meaning of L/rho, and every external normalization factor.
5. A list of all non-elementary identifiers used by that expression, with one exact definition or a direct definition pointer for each. Local aliases must resolve to that list.
6. Hlog/GPL argument, ordered letters, kernel sign, endpoint prescription and period-index convention; for machine rows, include the row-to-formula rule and period order.
7. Any alternate compact expression, paired with the specific exact conversion and its permitted order, parameter slice and coefficient field. Keep a canonical exact expression available when no compact conversion is established.
8. For the F/Q/R presentation, the bare/centered/physical layer, permutation or graph-ID meaning of each subscript, projection, factorial convention and source-to-I prefactor.
9. For Chen or Green formulas, the basepoint, path, root branch, continuous logarithm and boundary data. These are part of the function definition.

The finite coefficient should be printed explicitly wherever the selected record provides it. A label such as “M/N”, a list of basis names, or a link to a giant JSON file alone is insufficient for the professor to read the result. For large exact GPL/Chen expressions, a compact mathematical finite-sum rule plus the complete indexed data and its definition is an unambiguous alternative to thousands of printed terms.

## Audit limits and handoff

Only this English audit note was written. Main-workspace scientific files were read without modification. The audit checked definitions and their documented conversion scopes; it did not rerun HyperInt, Maple, numerical integrations, Lean, or historical identity checks. The original 18 October package was inspected as a dated snapshot, not promoted to a current completeness statement. Before building a new package, reselect current per-diagram inputs and freeze a new manifest; preserve the older package unchanged.
