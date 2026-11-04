# Two serial one-loop insertions and one two-loop insertion

These four source orbits share two endpoint geometries. Sharing that geometry does not identify their internal operators or their full colours. We evaluate the native ordered J, consume the Gamma/Euler factor already present in C, and only then convolve the complete epsilon-dependent operator coefficients. The active program is [chains_pi2.py](../scripts/chains_pi2.py). It reads no final result table and does not invoke a CAS process.

## What is inserted on the exchanged gluon

Use the common convention `wilson-qcd-action-path-fourier/1.0`, D=4-2epsilon, Fourier factor exp(-ikX), and the already defined 1PI operator

\[
\Pi_L^{ab,\mu\nu}(k)
=-i\delta^{ab}(k^2g^{\mu\nu}-k^\mu k^\nu)\mathcal P_L(k^2),
\qquad D_0^{ab,\mu\nu}=-\frac{i\delta^{ab}g^{\mu\nu}}{k^2+i0}.
\]

The Dyson products have coefficient +1. The chain is D0 Pi1 D0 Pi1 D0, whereas the two-loop source is D0 Pi2 D0. Since T=k^2g-kk satisfies T g T=k^2T, the serial chain retains a single transverse tensor and one exchanged momentum; it is not the product of two independent coordinate integrals. There is no extra factorial for interchanging the two identical Pi1 labels. Their colour deltas leave CF, multiplied by the full B1 squared:

\[
\mathcal B_1(\epsilon)=a(\epsilon)C_A-b(\epsilon)T_Fn_f,
\quad a=\frac{10-6\epsilon}{3-2\epsilon},
\quad b=\frac{8(1-\epsilon)}{3-2\epsilon}.
\]

Thus the chain has AA, Af and ff multipliers a^2, -2ab and b^2, with colour monomials CA^2 CF, CA CF TF nf and CF(TF nf)^2. They are epsilon-dependent multipliers of the same I; replacing them by their epsilon=0 values would lose finite contributions.

Pi2 instead uses the native full-D two-loop 1PI scalar. Its accepted input is the explicitly recorded bare Feynman-gauge literature operator, including the internal gluon, ghost and quark contributions. The program does not newly integrate those internal two-loop graphs. It preserves its three distinct AA, Af and Ff coefficients

\[
p_c(\epsilon)=A_c(\epsilon)+\frac{B_c(\epsilon)}{U(\epsilon)},
\qquad \phi_c(\epsilon)=4\epsilon^2p_c(\epsilon),
\]

where the rational A_c and B_c are the actual native operator fields, not fitted Laurent coefficients. The Ff monomial CF^2 TF nf is retained even though it is outside the first fixed-basis colour coordinate. The complete source has one Pi2 insertion and must not be replaced by B1 squared.

## The original ordered two-endpoint integrals

Both geometries have unit-square bounds and unit Jacobian, with the exact source order retained. The adjacent orbit has s31 outermost and s41 innermost; write x=s31, y=s41. The opposite orbit has s21 outermost and s41 innermost; there write x=s21, y=s41. On the positive real branch z>0,

\[
J=\int_0^1dx\int_0^1dy\,
\frac{N(x,y;\epsilon)}{2(1+2\epsilon)}f(x,y)^{-2+3\epsilon}.
\]

The original distance and full tensor numerator are

\[
f_A=z(1-x)y,\qquad N_A=(2+\epsilon)zf_A,
\]
\[
D=f_O=(1-x)(1-y)+zxy,\quad
N_O=-(2+\epsilon)(1+z)D+(1-3\epsilon)z.
\]

The evaluator reconstructs each of these polynomials and the complete original integrand from native source fields. The two operators have the same J. In the adjacent case the prescribed iterated integrals are immediately Euler integrals:

\[
J_A=\frac{(2+\epsilon)z^{3\epsilon}}{2(1+2\epsilon)}
\left(\frac1{3\epsilon}\right)^2
=\frac{(2+\epsilon)z^{3\epsilon}}{18\epsilon^2(1+2\epsilon)}.
\]

The endpoint limits are justified for Re epsilon>0. For the opposite case use the exact identity

\[
N_O D^{-2+3\epsilon}/[2(1+2\epsilon)]
=-(1+z)D^{-1+3\epsilon}
+\frac{\partial_x\partial_y D^{3\epsilon}}
{6\epsilon(1+2\epsilon)}.
\]

This retains the complete numerator, not just its transverse scalar part. The four source corners give D(0,0)=1, D(1,1)=z, and D(0,1)=D(1,0)=0. At Re epsilon>0 the zero-corner values vanish before continuation, so the complete mixed-derivative boundary is 1+z^(3epsilon). Define the already evaluated centered mother

\[
H_0(\delta,z)=-(1+z)z^{-\delta/2}
\int_0^1dx\int_0^1dy\,D^{-1+\delta}.
\]

Then the source reduction is

\[
J_O=z^{3\epsilon/2}H_0(3\epsilon,z)
+\frac{1+z^{3\epsilon}}{6\epsilon(1+2\epsilon)}.
\]

An opposite corner is locally comparable to a sum of its two positive distances. The derivative term is therefore locally absolutely integrable for Re epsilon>0, as is the centered D^(-1+3epsilon) term. This supplies the relevant integral/boundary connection. It does not assert arbitrary Fubini changes of any upstream momentum integral. The native branch contract initially uses 0<Re epsilon<1/3, followed by local meromorphic continuation around zero at fixed positive z. A new physical complex continuation is not established here.

## Consume C and obtain the physical I

Put rho=log z and

\[
\mathcal G=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad
U=\frac{\Gamma(1+\epsilon)^2\Gamma(1-\epsilon)\Gamma(1-3\epsilon)}
{\Gamma(1-2\epsilon)^2\Gamma(1+2\epsilon)},\quad P=\mathcal G^3U.
\]

The native C already contains its single exp(-3 gamma_E epsilon) compensation. The appearance of Gcal^3 consumes that factor once. The source C rational coefficient gives

\[
R_{\rm chain}=-\frac{P}{128\epsilon^2(1-2\epsilon)^2},
\qquad R_{\Pi_2}=+\frac{P}{128\epsilon^2(1-2\epsilon)^2}
\]

after the complete B1 squared or B2=4epsilon^2 P2 is extracted. These signs are read from native C; there is no runtime sign repair or historical-version selection. With a1=1 and source orbit multiplicity m=4 for adjacent or m=2 for opposite,

\[
I=\frac m2R\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.
\]

The chain therefore has the complete expressions

\[
I_{\rm chain,A}=-\frac{P(2+\epsilon)}
{1152\epsilon^4(1-2\epsilon)^2(1+2\epsilon)},
\]
\[
I_{\rm chain,O}=-\frac{P}{384\epsilon^3(1-2\epsilon)^2(1+2\epsilon)}
-\frac{P H_0(3\epsilon,z)}
{128\epsilon^2(1-2\epsilon)^2\cosh(3\epsilon\rho/2)}.
\]

For the same endpoint geometry the current Pi2 I is the negative of the chain I, as follows directly from its C and unchanged J. The historical schema name `I_geom` denotes this Pi2 physical I before multiplying phi_c; it does not introduce a second normalization. Its complete coloured coefficients are I_c=phi_c I_geom. In particular equality of geometric kernels does not identify Pi2 with two serial Pi1 operators.

## Explicit coefficient budgets and complete HLOG representation

The source-bound, previously evaluated M-moment data give

\[
H_0(3\epsilon,z)=M_{0,1}
+\frac{3\epsilon}{2}M_{1,1}
+9\epsilon^2\left(\frac18M_{2,1}+\frac1{24}M_{0,3}\right)
+O(\epsilon^3).
\]

The symbols are the original unnormalised moments. Their complete 0,-1 HLOG words are loaded from the four unchanged M text outputs; they are not new names for unevaluated subintegrals. Rho is the shape logarithm, distinct from the physical average scale L. Products are expanded by the existing exact rational shuffle algebra, including zeta periods. M03 may equivalently be reduced using M03=-M01^2; no new function family is needed.

| Required object | Adjacent | Opposite |
| --- | --- | --- |
| Physical geometric leading pole | epsilon^-4 | epsilon^-3, including its boundary |
| Regular Gamma/rational jet sufficient for c0 | Order4 | Order3 (order4 is also computed) |
| Centered H0 jet | Not required | Order2, because its coefficient has epsilon^-2 |
| External a^2,-2ab,b^2 or phi_c jet | Order4 | Order3 (order4 is also computed) |
| Pi2 U and inverse U before internal cancellation | Order5 | Same conservative order5 |

The last line matters. A_c and B_c individually contain cubic poles. In 4epsilon^2 B_c/U, the fifth-order term of inverse U can contribute at order4. Expanding only to the eventual regular order would miss a finite piece. The program first combines A+B, explicitly proves the artificial negative phi coefficient is zero, and only then forms the regular phi jet. It independently expands U phi=4epsilon^2(AU+B) and verifies its convolution with the geometry lacking U against the full phi times I_geom route in every colour and Laurent position.

Gamma jets are computed from exact logarithmic Gamma series and finite recurrence. All seven requested positions -6,...,0 are exported, including the established leading zeros, with O(epsilon) remainder at fixed z. The full endpoint reductions, source-derived signs, operator multipliers, consumed moment hashes and order budgets are stored with each result. None of the final old I tables or the old Pi2 operator-expansion table is an input to generation.

## Scope of this reproduction

The new work repeats source identities, endpoint/IBP boundary reduction, exact Gamma expansion and all epsilon-dependent colour convolutions. The H0/M integration is explicitly reused through its actual evaluated HLOG data, and the internal full-D Pi2 scalar is an adopted literature input. The work does not claim a second independent calculation of that two-loop self-energy or new HyperInt integrations. Post-generation table comparison, independently scoped source checks and formal proof status are separate; formal=not_run throughout.
