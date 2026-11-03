# A shared rectangle kernel for the intrinsic C and D sources

2026-10-16. Independent analytic derivation. The normalizations below refer to the bare geometric integrals C and D, without color, Gamma or crossing prefactors. This note derives a reusable one-propagator rectangle before introducing any three-loop moments. It does not prove minimality or reduction to the previously selected P/Q/R slices. Formal status: `not_run`.

## 1. Evaluate a rectangle with one zero corner

Let A=exp(eta), B=exp(-eta), so AB=1, and let

\[
0<\kappa\leq A+B=2\cosh\eta,\qquad
g(x,y)=Ax+By-\kappa xy.
\]

The shape polynomial is positive in the open square. It has a zero corner at (0,0), and at the limiting value kappa=A+B also at (1,1). For Re(epsilon)>-1, define the one-propagator rectangle

\[
\boxed{\mathcal B_\epsilon(\eta,\kappa)
=-\kappa\int_0^1dx\int_0^1dy\,
g(x,y)^{\epsilon-1}.}
\tag{K1}
\]

This integral is evaluated explicitly as follows. First integrate x:

\[
\mathcal B_\epsilon=-\frac\kappa\epsilon\int_0^1
\frac{[A+(B-\kappa)y]^\epsilon-(By)^\epsilon}{A-\kappa y}\,dy.
\tag{K2}
\]

The difference quotient is continuous when its denominator vanishes. Set

\[
u=1-\kappa y/A,\quad q=1-\kappa/A,
\quad r=1-\kappa/B.
\]

Because AB=1, the two bases become kappa^{-1}(1-ru) and kappa^{-1}(1-u). With

\[
J_\epsilon(x)=\int_0^x\frac{(1-v)^\epsilon-1}{v}\,dv,
\]

the remaining primitive is kappa^{-epsilon}[J_epsilon(ru)-J_epsilon(u)]/epsilon. Evaluation from u=1 to u=q gives

\[
\boxed{\mathcal B_\epsilon(\eta,\kappa)
=-\frac{\kappa^{-\epsilon}}\epsilon
\left[J_\epsilon(q)+J_\epsilon(r)-J_\epsilon(qr)-J_\epsilon(1)\right].}
\tag{K3}
\]

Here q,r<1 and

\[
qr=1-\kappa(2\cosh\eta-\kappa)\leq1.
\]

Consequently all four J values use the real branch inherited from K1. The coordinate u can cross zero, but the difference of primitives is regular there. No principal-value prescription is introduced.

K1 also proves analyticity near epsilon=0. If raw parameter derivatives are used,

\[
\mathcal B_n(\eta,\kappa)
=\left.\partial_\epsilon^n\mathcal B_\epsilon\right|_0
=-\kappa\int_0^1dx\int_0^1dy\,
\frac{\log^n g(x,y)}{g(x,y)},
\qquad
\mathcal B_\epsilon=\sum_{n\geq0}\frac{\epsilon^n}{n!}\mathcal B_n.
\tag{K4}
\]

These coefficients need not be left as two-dimensional integrals. Put

\[
j_m(x)=\frac1{m!}\int_0^x\frac{\log^m(1-v)}v\,dv,
\qquad J_\epsilon(x)=\sum_{m\geq1}\epsilon^m j_m(x).
\]

Then the explicit finite coefficient rule is

\[
\boxed{\mathcal B_n
=-n!\sum_{m=1}^{n+1}
\frac{(-\log\kappa)^{n+1-m}}{(n+1-m)!}
[j_m(q)+j_m(r)-j_m(qr)-j_m(1)].}
\tag{K5}
\]

In particular,

\[
\boxed{\mathcal B_0=
\operatorname{Li}_2(q)+\operatorname{Li}_2(r)
-\operatorname{Li}_2(qr)-\zeta_2.}
\tag{K6}
\]

Two elementary checks are

\[
\mathcal B_{\epsilon=1}(\eta,\kappa)=-\kappa,
\qquad
\mathcal B_\epsilon(0,1)=\frac{J_\epsilon(1)}\epsilon
=-\frac{\gamma_E+\psi(1+\epsilon)}\epsilon.
\]

The first statement evaluates the complete generating function at epsilon=1; it is not a statement about the raw derivative B_1.

## 2. The complete one-loop M source is a boundary of this kernel

The x/y symmetry gives

\[
\mathcal B_\epsilon(-\eta,\kappa)=\mathcal B_\epsilon(\eta,\kappa).
\tag{K7}
\]

At eta=-log(z)/2 and kappa=(1+z)/sqrt(z)=2cosh(eta), q=-z, r=-1/z, qr=1. Therefore

\[
\boxed{\mathcal B_\epsilon\!\left(-\frac{\log z}{2},
\frac{1+z}{\sqrt z}\right)=\mathcal V_M(\epsilon,z).}
\tag{K8}
\]

This is the existing bare one-loop generating function, including its normalization, rather than a new definition of it. The rectangle kernel thus extends a one-loop source by a genuine second shape parameter kappa. Its boundary and epsilon derivatives connect directly to the original M moments.

## 3. Map the complementary rectangle of the ordered graph

Use

\[
h(t,s)=(1-t)(1-s)+zts,\quad H=h/\sqrt z,\quad
c=(1+z)/\sqrt z,\quad \alpha=(1+z)^{-1}.
\]

The propagator is k_epsilon(t,s)=-c H(t,s)^{epsilon-1}. For the complementary rectangle of D,

\[
G_\epsilon(t,s)=\int_t^1du\int_0^sdv\,k_\epsilon(u,v),
\]

make the actual affine substitution u=1-(1-t)x, v=sy. Let

\[
a=1-t,\quad b=s,\quad
\eta(t,s)=\tfrac12\log\frac{1-t}{zs},\quad
\kappa(t,s)=(1+z)\sqrt{\frac{(1-t)s}{z}}.
\]

Then H(u,v)=sqrt(ab) g(x,y), including the same Jacobian as K1, and

\[
\boxed{G_\epsilon(t,s)=[(1-t)s]^{\epsilon/2}
\mathcal B_\epsilon(\eta(t,s),\kappa(t,s)).}
\tag{K9}
\]

The scale extracted from the rectangle is [(1-t)s]^{epsilon/2}; it is independent of z. Moreover

\[
2\cosh\eta-\kappa=H(t,s)/\sqrt{(1-t)s}>0
\]

inside the square, so the physical parameters lie in the K1 domain. Four original end-point integrations will have been completed once K3 is inserted in the formulas below.

## 4. The two linear primitives of D

The original ordered-domain reduction is

\[
D=\int_0^1dt\int_0^1ds\,
\left[\int_s^1k_\epsilon(t,v)dv\right]
\left[\int_0^tk_\epsilon(u,s)du\right]G_\epsilon(t,s).
\tag{K10}
\]

Define endpoint bases A_0=sqrt(z)t, B_0=(1-s)/sqrt(z). The two linear primitives are exactly

\[
-\frac{A_0^\epsilon-H^\epsilon}{\epsilon(t-\alpha)},
\qquad
-\frac{H^\epsilon-B_0^\epsilon}{\epsilon(s-\alpha)}.
\tag{K11}
\]

Their endpoint product A_0 B_0=t(1-s) is independent of z. Introduce the local logarithms

\[
L=\tfrac12\log\frac{zt}{1-s},\qquad
V=\log\frac{H}{\sqrt{t(1-s)}},\qquad
U=\tfrac12\log[t(1-t)s(1-s)]+\log H.
\tag{K12}
\]

Using the geometric mean sqrt(A_0 B_0), the product of power differences becomes

\[
(A_0^\epsilon-H^\epsilon)(H^\epsilon-B_0^\epsilon)
=2[t(1-s)]^\epsilon e^{\epsilon V}
[\cosh(\epsilon L)-\cosh(\epsilon V)].
\tag{K13}
\]

Combining K9-K13 gives the complete two-dimensional outer source

\[
\boxed{D=2\int_0^1dt\int_0^1ds\,
\frac{e^{\epsilon U}[\cosh(\epsilon L)-\cosh(\epsilon V)]}
{\epsilon^2(t-\alpha)(s-\alpha)}
\mathcal B_\epsilon(\eta(t,s),\kappa(t,s)).}
\tag{K14}
\]

This is a direct reduction from the graph. No global rho factor is stripped off and reattached under a new name. All four integrated variables are represented by elementary primitives and the explicitly evaluated independent rectangle kernel.

At t=alpha one has H=A_0 and V=L; at s=alpha one has H=B_0 and V=-L. Thus the two apparent interior poles are removable. More specifically,

\[
\cosh(\epsilon L)-\cosh(\epsilon V)
=2\sinh\frac{\epsilon(L+V)}2\sinh\frac{\epsilon(L-V)}2,
\]

where L-V=log(A_0/H) vanishes at t=alpha and L+V=log(H/B_0) vanishes at s=alpha. These are separate divided differences, so their simultaneous intersection is regular too. This regularity should be preserved in numerical implementations.

## 5. Indexed coefficients with constant assembly coefficients

For m,n>=0 and j>=1 define the two-dimensional moments

\[
\mathfrak X_{m,n,j}(z)=\int_0^1dt\int_0^1ds\,
\frac{U^m(L^{2j}-V^{2j})}{(t-\alpha)(s-\alpha)}
\mathcal B_n(\eta(t,s),\kappa(t,s)).
\tag{K15}
\]

The individual moments retain the removable-pole factors, since L^{2j}-V^{2j} is divisible by (L-V)(L+V). Substitution of the evaluated K5 gives a two-dimensional integral of explicitly known one-variable GPL functions. These are not new unresolved four-dimensional integrals.

The complete Taylor/Laurent-regular coefficients are

\[
\boxed{D_q=2\sum_{\substack{m,n\geq0,\ j\geq1\\m+n+2j=q+2}}
\frac{\mathfrak X_{m,n,j}}{m!n!(2j)!}.}
\tag{K16}
\]

Thus

\[
\begin{aligned}
D_0&=\mathfrak X_{0,0,1},\\
D_1&=\mathfrak X_{1,0,1}+\mathfrak X_{0,1,1},\\
D_2&=\tfrac12\mathfrak X_{2,0,1}+\mathfrak X_{1,1,1}
+\tfrac12\mathfrak X_{0,2,1}+\tfrac1{12}\mathfrak X_{0,0,2},\\
D_3&=\tfrac16\mathfrak X_{3,0,1}+\tfrac12\mathfrak X_{2,1,1}
+\tfrac12\mathfrak X_{1,2,1}+\tfrac16\mathfrak X_{0,3,1}
+\tfrac1{12}(\mathfrak X_{1,0,2}+\mathfrak X_{0,1,2}),\\
D_4&=\tfrac1{24}\mathfrak X_{4,0,1}+\tfrac16\mathfrak X_{3,1,1}
+\tfrac14\mathfrak X_{2,2,1}+\tfrac16\mathfrak X_{1,3,1}
+\tfrac1{24}\mathfrak X_{0,4,1}\\
&\quad+\tfrac1{24}\mathfrak X_{2,0,2}+\tfrac1{12}\mathfrak X_{1,1,2}
+\tfrac1{24}\mathfrak X_{0,2,2}+\tfrac1{360}\mathfrak X_{0,0,3}.
\end{aligned}
\tag{K17}
\]

All assembly coefficients are rational constants. There is no external log(rho) convolution. The GPL weight of B_n is n+2; the natural total weight of X_m,n,j is m+n+2j+4, so K16 has weight q+6.

Under z->1/z, t->1-s, s->1-t, H,U,V,kappa remain unchanged, L and eta change sign, and the product of denominators is unchanged. K7 therefore shows that every X_m,n,j is crossing-even individually.

## 6. The same rectangle kernel also gives C

For the reversed pairing, the original domain gives

\[
C=\int_0^1dt\int_0^1ds\,k_\epsilon(t,s)
G_\epsilon(t,s)G_\epsilon(s,t).
\]

Consequently, with exactly the same U,

\[
\boxed{C=-\int_0^1dt\int_0^1ds\,\frac cH e^{\epsilon U}
\mathcal B_\epsilon(\eta(t,s),\kappa(t,s))
\mathcal B_\epsilon(\eta(s,t),\kappa(s,t)).}
\tag{K18}
\]

This independent second graph is important: B_epsilon is a shared single-propagator building block, not a name for the complete D integrand. The corresponding moments with U^m B_n(t,s) B_p(s,t) give C_q by m+n+p=q and the constant factors 1/(m!n!p!).

## 7. What has and has not been obtained

K14-K18 establish a natural, indexed, crossing-adapted two-dimensional source with no external rho powers. Its one-propagator building block is evaluated in J, its coefficient kernels are explicit finite GPL combinations, and its complete-rectangle boundary is the original M generator. This is a structural reformulation of the already derived full-epsilon C/D integrals; it is not a new independent numerical value of any coefficient.

These formulas do not prove that the two-dimensional indexed moments reduce to a finite selection of M,N,P,Q,R, that they are analytically independent of those families, or that two dimensions are necessary. A further comparison with the previously isolated quotient directions is needed before calling any selected moment a minimal new master family.

Original ordered-domain sources: [D](../pure-three-cd-20261016/D/derivation.md), [C](../pure-three-cd-20261016/C/derivation.md). Existing complete M source and two-loop N reduction: [4e-04](../4e04-direct-mn.md).
