# One opposite kernel: integrate ordered endpoints, then the solved HLOG jets

This shared lemma covers the fifteen individually derived sources
6e-04/05/06/13/14/15/16/26/27/28/29/30/31/36/39. A label is not an answer
rule: each article displays its own J, bounds, opposite variables, adjacent
q factors, endpoint words, weight and kappa before applying this lemma.
The source is fixed positive real z, real logs and local epsilon continuation.

## Endpoint primitives retain the actual placement of the opposite point

Let t,s be the two opposite endpoints. Put
A(x)=x^(epsilon-1),B(x)=(1-x)^(epsilon-1) and mark a retained endpoint O.
The other four endpoints belong to adjacent propagators. The concrete root
has numerator -(1+z) times the adjacent q factors, so for r adjacent q=z

$$J=-\frac{(1+z)z^{r\epsilon}}{\epsilon^4}
\int_{[0,1]^2}D_z(t,s)^{\epsilon-1}W(t,s;\epsilon)dt ds,\quad
D_z=(1-t)(1-s)+zts.\tag{S1}$$

The factor epsilon^-4 counts the four actually integrated adjacent endpoints.
W includes their ordering factorials and independent Beta-chain factors.
No Jacobian from an old cube export is inserted. For Re epsilon>0, the
full-cube adjacent densities are integrable and D is comparable to the sum
of the two distances at each of its zero corners. An integrable product
majorant therefore justifies the specific chain regrouping and separations;
all original bounds survive in the per-entry inequalities. This is an
analytic bridge for these roots, not a formal general Fubini receipt.

Write a=x^epsilon,b=(1-x)^epsilon. A hat means the integral multiplied by
one epsilon for each integrated adjacent endpoint. Direct primitives give

$$\widehat F_A(0,x)=a,\quad \widehat F_B(0,x)=1-b,\quad
\widehat F_A(x,1)=1-a,\quad\widehat F_B(x,1)=b.$$

For equal-factor words the ordered-square symmetry gives half their squares:
AA lower=a²/2, AA upper=(1-a)²/2, BB lower=(1-b)²/2, BB upper=b²/2.
For one endpoint on each side of O the factors multiply. Whole unmarked
AB or BA chains contribute R2 or1-R2, by the evaluated Beta calculation in
[adjacent-chains](adjacent-chains.md).

The only non-elementary-looking partial word can also be kept exactly:

$$B_A(x)=\epsilon^2\int_{0<u<v<x}(1-u)^{\epsilon-1}v^{\epsilon-1}du dv
=x^\epsilon-\epsilon\int_0^xv^{\epsilon-1}(1-v)^\epsilon dv.$$

The lower AB word is a(1-b)-B_A; path composition gives upper AB=R2-a+B_A
and upper BA=a+b-ab-R2-B_A. To prove the reflection identity, differentiate
the two incomplete Beta terms at x and1-x. Their epsilon-weighted sum has
derivative d(ab)/dx, and its x=0 value is R2. Hence

$$\boxed{B_A(x)+B_A(1-x)=a+b-ab-R_2.}\tag{S2}$$

For30/31 the actual phase is r-1=0. Reflection of both retained endpoints
therefore permits the *physical paired I* to replace W by(W+RW)/2, removing
B_A via(S2). This does not evaluate or replace the unsymmetrized bare J.
Other entries use their actual weights without this operation.

## Evaluate the retained kernel as a standard function

For nonnegative integers a,b,c,d define the centered kernel

$$H_{ab;cd}(\epsilon,z)=-(1+z)z^{-\epsilon/2}
\int_{[0,1]^2}D_z^{\epsilon-1}
 t^{a\epsilon}(1-t)^{b\epsilon}s^{c\epsilon}(1-s)^{d\epsilon}dt ds.\tag{S3}$$

The definition alone is not a claimed solution. Separate zuv and the other
term in D by Mellin--Barnes, with0<epsilon<1/4 and epsilon-1<Re w<0. The
actual t,s Beta integrals produce

$$\frac1{2\pi i}\int dw\,z^w
\frac{\Gamma(-w)\Gamma(1-\epsilon+w)
\Gamma(1+a\epsilon+w)\Gamma((b+1)\epsilon-w)
\Gamma(1+c\epsilon+w)\Gamma((d+1)\epsilon-w)}
{\Gamma(1-\epsilon)\Gamma(1+(a+b+1)\epsilon)
\Gamma(1+(c+d+1)\epsilon)}.$$

Multiplication by -(1+z)z^(-epsilon/2) gives H. This contour separates both
sets of Gamma poles and the Euler exponents have positive real parts. The
six Gamma factors match the standard MeijerG definition, giving the complete
all-epsilon evaluation

$$H_{ab;cd}=-\frac{(1+z)z^{-\epsilon/2}}
{\Gamma(1-\epsilon)\Gamma(1+(a+b+1)\epsilon)\Gamma(1+(c+d+1)\epsilon)}
G_{3,3}^{3,3}\!\left(z\;\middle|\;\begin{matrix}
\epsilon,-a\epsilon,-c\epsilon\\0,(b+1)\epsilon,(d+1)\epsilon
\end{matrix}\right).\tag{S4}$$

This MeijerG is distinct from the HLOG G used for finite coefficients. The
formula is an evaluated standard function, not an unknown integral renamed H.
Reflection(t,s)->(1-t,1-s) has unit absolute Jacobian and gives
H_ab;cd(1/z)=H_ba;dc(z); exchanging t and s gives H_ab;cd=H_cd;ab.
Consequently all powers in the present paired weights reduce to

| Name | Exact parameters(a,b;c,d) |
|---|---|
| H0 | (0,0;0,0) |
| H1 | (0,0;1,0) |
| H2 | (0,0;2,0) |
| HW | (0,0;1,1) |
| HD | (1,0;1,0) |
| HX | (1,0;0,1) |

H0,HW,HX are crossing even. H2 includes the full s^(2epsilon) weight; it is
twice an ordered-AA mother and must not acquire a second factorial.

## The actual finite Taylor integrals and their HLOG evaluation

H is holomorphic in a sufficiently small epsilon disk at each fixed z>0.
The corner power bound has a margin, and finite log insertions are absorbed
by a slightly stronger integrable power. Set h_ab;cd,q=[epsilon^q]H.
Use the source-specific rational map

$$t=\frac{x}{1+x},\quad s=\frac{y}{1+y},\qquad
 dt ds=\frac{dxdy}{(1+x)^2(1+y)^2},\qquad (x,y)\in(0,\infty)^2.$$

Since D=(1+zxy)/[(1+x)(1+y)], the *actual* coefficient integral is

$$h_{ab;cd,q}=-\frac{1+z}{q!}\int_0^\infty\!dx\int_0^\infty\!dy\,
\frac{\Lambda_{ab;cd}(x,y)^q}{(1+x)(1+y)(1+zxy)},\tag{S5}$$
$$\Lambda=\log(1+zxy)-\ell/2+a\log x+c\log y
-(1+a+b)\log(1+x)-(1+c+d)\log(1+y).$$

This specifies the exact normalized object, Taylor order, domain and kernel
for each adopted jet. It agrees term for term with wbBase and wbLog in the
native compute.mpl. The historical HyperInt operation is

```maple
hyperInt(convert(wbBase*wbLog^q/q!,HlogRegInf),[x,y]);
fibrationBasis(wbAnswer,[z]);
```

Here each unspecified endpoint in[x,y] is0..infinity. The primitive in x has
only the letters0,-1,-1/(zy); evaluation uses its regularized0/infinity
limits. HyperInt's log-polynomial conversion and fibration combine these
limits before y is integrated; raw endpoint pieces need not converge
separately. The remaining y endpoints are also0/infinity. Finishing the
fibration puts the solution into z letters0,-1 with the stored period table.
Divergence checks and abort-on-divergence were enabled, and nonlinear
polynomials were not ignored. Reflected outputs substitute z=1/zz into the
*evaluated* raw expression, then refibrate in zz before renaming it z.

As a concrete primitive at q=0, before the second integration,

$$\int_0^\infty\frac{dx}{(1+x)(1+zxy)}
=\frac{\log(zy)}{zy-1}.$$

The zy=1 point is removable. The remaining real integral gives
h_i0=-(ell²+pi²)/2=M01 for all six types. At higher q the logarithm Lambda
and1/q! in(S5) fix all endpoint constants; a symbol-only integration or a
numeric fit would not suffice.

The [evaluated jet appendix](appendices/single-opposite-jets.md) prints every
h_iq for the six types and q=0,...,4 as complete HLOG/zeta expressions,
transcribed from the native solved raw outputs. Thus no coefficient is left
as an unevaluated(S5) integral. The exact raw/reflected files, successful logs,
request definitions and generated scripts remain byte-bound in the runtime;
the current replay rechecks and parses them. It does not rerun HyperInt.
This is explicit adoption of already evaluated subintegrals.

## Substitution into each physical I and finite coefficient rule

Expand the particular W as a finite polynomial in
AT=t^epsilon,BT=(1-t)^epsilon,AS=s^epsilon,BS=(1-s)^epsilon.
For a monomial AT^a BT^b AS^c BS^d replace its integral by(S3), preserving
its coefficient, including any R2. Let T[W] be that finite linear combination
of the six evaluated H functions and their reflections. From(S1),

$$J=\epsilon^{-4}z^{(r+1/2)\epsilon}\mathcal T[W](\epsilon,z).$$

In the [single normalization](common-normalization.md), set delta=r-1. The
complete physical result is

$$\boxed{I=\frac{\kappa Z^3}{\epsilon^4\cosh(3\epsilon\ell/2)}
\frac{e^{\delta\epsilon\ell}\mathcal T[W](z)
+e^{-\delta\epsilon\ell}\mathcal T[W](1/z)}2.}\tag{S6}$$

Each article states its concrete W and the paired polynomial, making this
substitution explicit. If W=sum b_alpha(epsilon)AT^a BT^b AS^c BS^d and
b_alpha=sum b_alpha,j epsilon^j, define

$$d_m=\frac12\sum_\alpha\sum_{j+u+v=m}
b_{\alpha,j}\frac{(\delta\ell)^u}{u!}
[h_{\alpha,v}(z)+(-1)^u h_{\alpha,v}(1/z)].\tag{S7}$$

The same expression applies after an allowed I-only BA average, with that
averaged polynomial. With p_n from the normalization recurrence,

$$c_{-6}=c_{-5}=0,\qquad
c_k=\kappa\sum_{n=0}^{k+4}p_n d_{k+4-n}\quad(-4\le k\le0),$$
$$\boxed{c_0=\kappa[d_4+p_2d_2+p_3d_1+p_4d_0].}\tag{S8}$$

Thus q<=4 jets, R2 through4, phase exponential through4 and prefactor through4
are a sufficient finite budget. Exact cancellations may improve it but are
not assumed to eliminate a source before it is computed. All seven requested
powers and their zeroes follow from these finite sums, not from fitting an old
I table. Their explicit compact results appear in each article.

For comparison with established named functions, tau_mq is the ordinary
coefficient of Hm^+ for m=0,1,2; tau_D,q is the coefficient of HD^+;
eta_iq=ell times the coefficient of Hi^- for i=1,2,D. No extra q! is used.
In particular tau_D4=[h_HD,4(z)+h_HD,4(1/z)]/2, an explicitly solved function
in the appendix. Its appearance in36/39 is not an unknown remainder. M/N
compactification is an exact dictionary substitution in complete HLOG/period
coordinates; it is separate from the integration. The current dictionary
contains no general proof that all such functions reduce to M/N, and a failed
restricted fit is not a proof of impossibility.
