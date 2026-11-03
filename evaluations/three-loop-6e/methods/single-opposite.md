# One opposite propagator: fifteen three-loop integrals

This route covers 6e-04/05/06/13/14/15/16/26/27/28/29/30/31/36/39. The opposite kernel couples two endpoints, while the four adjacent endpoints can be integrated along the actual ordered chains. The purpose of the reduction is to keep that coupling in a small set of two-variable moments.

Let r count adjacent q=z factors. Mark the opposite endpoints by O and use A/B for adjacent endpoint densities as in [the chain method](empty.md). The normalized chain integrations give

\[
J=-(1+z)z^{r\epsilon}\epsilon^{-4}
\int_0^1dt\,ds\,D_z(t,s)^{\epsilon-1}W(t,s;\epsilon).
\tag{4}
\]

The factor epsilon^-4 accounts for the four integrated adjacent endpoints; W includes their ordering factorials and any independent Beta/Gamma word. For example, the original chains for 6e-04 give W=(1-s^epsilon)^2/4. The weights for 6e-36 and 6e-39 are respectively (1-t^epsilon)(1-s^epsilon) and [1-(1-t)^epsilon][1-(1-s)^epsilon]. Each weight follows from its own bounds.

## Retaining ordered insertions correctly

Write a=x^epsilon, b=(1-x)^epsilon and B_A(x)=epsilon^2 F_BA(0,x). Lower-interval integrations include

\[
\widehat F_A=a,\quad\widehat F_B=1-b,\quad
\widehat F_{AA}=a^2/2,\quad\widehat F_{BB}=(1-b)^2/2,
\]
\[
\widehat F_{BA}=B_A,\qquad\widehat F_{AB}=a(1-b)-B_A.
\]

Here a hat includes one epsilon for each adjacent endpoint. Path splitting and reflection give

\[
B_A(x)+B_A(1-x)=a+b-ab-R_2.
\tag{5}
\]

The BA insertion is eliminated only when the actual physical combination permits averaging the entire weight with its reflection. In the present corpus those cases have equal crossing phase weights. Equation (5) then removes B_A at the I level. It would be incorrect to claim the unsymmetrized J had been reduced in this way. The evaluator records this scope explicitly.

## The centered moments and their evaluated representation

Define the common moment

\[
H_{ab;cd}(\epsilon,z)=-(1+z)z^{-\epsilon/2}
\int_0^1dt\,ds\,D_z(t,s)^{\epsilon-1}
t^{a\epsilon}(1-t)^{b\epsilon}s^{c\epsilon}(1-s)^{d\epsilon}.
\tag{6}
\]

For fixed positive z it is holomorphic near epsilon=0. A small disk suffices: at a zero corner all endpoint powers can be bounded by a slightly stronger but integrable radial power, and logarithms from finite derivatives can be absorbed by another arbitrarily small power. The six moment types required after the actual weight reduction are

| Name | (a,b;c,d) |
| --- | --- |
| H0 | (0,0;0,0) |
| H1 | (0,0;1,0) |
| H2 | (0,0;2,0) |
| HW | (0,0;1,1) |
| HD | (1,0;1,0) |
| HX | (1,0;0,1) |

Their full-epsilon value can also be written as a standard Meijer G:

\[
H_{ab;cd}=-\frac{(1+z)z^{-\epsilon/2}}
{\Gamma(1-\epsilon)\Gamma(1+(a+b+1)\epsilon)\Gamma(1+(c+d+1)\epsilon)}
G^{3,3}_{3,3}\!\left(z\;\middle|\begin{matrix}
\epsilon,-a\epsilon,-c\epsilon\\0,(b+1)\epsilon,(d+1)\epsilon
\end{matrix}\right).
\tag{7}
\]

Mellin--Barnes separation of the two terms of D followed by the two Euler integrals gives the six Gamma factors of this G. This standard-function identity is part of the existing analytic derivation. In the replay, the finite Taylor coefficients are read again from the original HyperInt text and parsed; neither evaluation of the Meijer G formula nor a new Mellin--Barnes/CAS integration is claimed.

The reflection law H_ab;cd(1/z)=H_ba;dc(z), together with exchange of the t and s integrations, supplies the raw/reflected channels. In particular H0, HW and HX are crossing even.

## Assembly, budget and compact functions

Expand the actual weight as W=sum_j b_j(epsilon) t^(a_j epsilon)(1-t)^(b'_j epsilon) s^(c_j epsilon)(1-s)^(d_j epsilon), after any legal I-only reflection average. Then (4) becomes

\[
J=\epsilon^{-4}z^{(r+1/2)\epsilon}\sum_j b_jH_j.
\]

With delta=r-1, the physical result is

\[
I=\frac{\kappa\mathcal G^3}{\epsilon^4\cosh(3\epsilon\rho/2)}
\frac{e^{\delta\epsilon\rho}\sum_j b_jH_j(z)
+e^{-\delta\epsilon\rho}\sum_j b_jH_j(1/z)}2.
\tag{8}
\]

The regular product in (8) must be known through epsilon^4 to determine c0. This includes the moment jets q=0,...,4, the Gamma factors, any Beta coefficients b_j and both phase weights. Computing only H_j at zero would determine only the potential leading pole. Lower powers epsilon^-6 and epsilon^-5 vanish by this bound; weight cancellations can remove further poles.

The familiar coefficient functions are ordinary Taylor coefficients, with no extra q!:

\[
\tau_{m,q}=[\epsilon^q]H_m^+\quad(m=0,1,2),\qquad
\eta_{i,q}=\rho[\epsilon^q]H_i^-\quad(i=1,2,D).
\]

HW/HD/HX use the corresponding even coefficient functions. Here H2 equals twice the older ordered-AA mother, so that factor must be retained. Every tau/eta used in the finite result has an already evaluated HLOG value; these names are compact coordinates for solved functions. Conversion to M/N and zeta is exact substitution in those coordinates. Some finite coefficients, including 36/39, require the existing extended tau/eta representation within the tested compact spaces. Failure of a restricted M/N fit does not invalidate a complete HLOG answer or prove non-representability in every possible M/N expression.

Historical raw moment outputs are reused only with matched definition, endpoint weights, order, branch, script and successful execution. The fresh part of P2 is recognition, raw parsing and finite physical convolution for every orbit. Comparison with the saved final table is performed afterwards.

The preserved [historical basis derivation](../evidence/historical/benchmarks/2ne/families/single-opposite/basis.md) fixes the six moment conventions, and [the evaluated coefficient definitions](../evidence/historical/benchmarks/2ne/families/single-opposite/basis-functions.md) give their complete HLOG expansions.
