# Two opposite propagators in one channel: six three-loop integrals

The relevant orbits are 6e-10/11/12/23/24/25. Their two opposite propagators share two ordered supporting edges. The remaining adjacent propagator can be integrated exactly, but its position on the supporting edge leaves a gap weight. Keeping that weight is the distinction from a product of independent one-opposite moments.

Let

\[
\Omega_2=\{0<t_1<t_2<1,\quad0<v<w<1\}.
\]

For id pairing use D_z(t1,v)D_z(t2,w); for swap use D_z(t1,w)D_z(t2,v). Define

\[
V_\pi[W]=(1+z)^2z^{-\epsilon}
\int_{\Omega_2}D_1^{\epsilon-1}D_2^{\epsilon-1}W(v,w;\epsilon)
\,dt_1dt_2dv\,dw.
\tag{9}
\]

Two actual gap weights are

\[
W_A=1-w^\epsilon,\qquad W_B=w^\epsilon-v^\epsilon.
\]

Integrating an OOA edge gives W_A/epsilon; OAO gives W_B/epsilon. The other endpoint of the adjacent propagator gives another 1/epsilon. Thus the canonical roots are

\[
J_{10}=z^{2\epsilon}V_{\mathrm{id}}[W_A]/\epsilon^2,\quad
J_{11}=z^{2\epsilon}V_{\mathrm{id}}[W_B]/\epsilon^2,\quad
J_{12}=z^{2\epsilon}V_{\mathrm{swap}}[W_A]/\epsilon^2.
\tag{10}
\]

Their kappa values are respectively -1/8,-1/16,-1/16. Simultaneous reflection of both supporting edges uses D_z(1-t,1-s)=z D_(1/z)(t,s), exchanges the appropriate endpoint positions and gives

\[
J_{23}(z)=z^{3\epsilon}J_{10}(1/z),\quad
J_{25}(z)=z^{3\epsilon}J_{11}(1/z),\quad
J_{24}(z)=z^{3\epsilon}J_{12}(1/z).
\]

The kappa factors match within each pair. Equation (1) in [the common normalization](common-normalization.md) then proves the complete physical-I equalities 10=23,11=25,12=24. It does not identify the unreflected J functions.

## Why the apparent double pole is reduced

For fixed z>0, each D vanishes linearly at its two isolated corners. To continue (9) around epsilon=0, retain the whole gap difference and use, for |epsilon|<=delta,

\[
\left|\frac{1-w^\epsilon}{\epsilon}\right|
\le |\log w|w^{-\delta},\qquad
\left|\frac{w^\epsilon-v^\epsilon}{\epsilon}\right|
\le |\log v|v^{-\delta}.
\]

With a sufficiently small delta, for example delta<1/4, these weights and two independent stronger corner kernels remain integrable. One may enlarge the ordered region to a cube for this bound. Finite logarithmic derivatives are covered by a slightly stronger power. Hence V[W] has a zero at epsilon=0 and

\[
V_\pi[W]=\epsilon U_1+\epsilon^2U_2+O(\epsilon^3).
\]

The complete three-loop I therefore has at most a simple pole, even though an epsilon^-2 factor appears in (10). The replay retains the order-0 slot and checks its cancellation rather than dropping it by an assumed diagram rule.

## The integrals that were actually evaluated

First integrate t1 from 0 to u=t2. Put s=v,r=w for id and s=w,r=v for swap, and write

\[
D_1=D_z(u,s),\quad D_2=D_z(u,r),\quad B=1-s,\quad d=(1+z)s-1.
\]

The exact slice is

\[
\int_0^uD_z(t,s)^{\epsilon-1}dt
=\frac{D_1^\epsilon-B^\epsilon}{\epsilon d}.
\]

At d=0 the quotient is removable. It is not two separately convergent pole integrals. Define

\[
K_0=(1+z)^2\frac{\log(D_1/B)}{dD_2},\qquad
K_\log=\log D_2+\tfrac12(\log D_1+\log B)-\log z.
\]

If W=epsilon W1+epsilon^2 W2+..., the required subintegrals are

\[
U_1=\int_0^1du\int_0^1dw\int_0^w dv\,K_0W_1,
\]
\[
U_2=\int_0^1du\int_0^1dw\int_0^w dv\,K_0(W_2+W_1K_\log).
\tag{11}
\]

For A, W1=-log w and W2=-(log w)^2/2. For B, W1=log(w/v) and W2=[(log w)^2-(log v)^2]/2. Historical HyperInt integration used a=1/(1+z) and finite integrations v, then u, then w, keeping divergence checks active. Its complete outputs were converted back to the z alphabet {0,-1}.

The current general engine instead emits the directly reduced ordered four-variable gap kernel and its full jet q=0,1,2. This is the same root and epsilon budget, with a different intermediate organization. P2 replays that engine's exact archived request/script/raw output. It does not claim to have launched a new computation of (11).

Writing Uq^+=(Uq(z)+Uq(1/z))/2 and Uq^-=(Uq(z)-Uq(1/z))/2, the final two nonzero candidate slots are

\[
c_{-1}=\kappa U_1^+,\qquad
c_0=\kappa\left(U_2^++\frac{\log z}{2}U_1^-\right).
\tag{12}
\]

Gcal^3 and sech(3 epsilon log z/2) have no linear epsilon term, but the odd crossing coefficient in (12) must be kept. It contributes to the finite result.

The historical compact notation SigmaA/SigmaB/SigmaC represents the three already evaluated finite coefficients. These are not unsolved remainders. Complete HLOG values and exact conversions to the previously defined tau/eta and E directions are preserved in the representation dictionaries. Whether a shorter M/N expression exists is a representation question, distinct from having all words and constants in (12).

The original detailed derivation is retained as [historical method evidence](../evidence/historical/benchmarks/2ne/families/same-channel-two/method.md).
