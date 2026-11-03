# Endpoint conversion and analytic constants for4e-04

This appendix proves equations(7)--(9) of the [main derivation](../derivation.md).
It does not introduce a different Q. All notation V,Y,Yc,V0,E,J_a,F,Psi is
as defined there. Write sigma=z/(1+z)², r=z/(1+z), ell=log z and theta=z d/dz.
The proof starts in a common absolutely convergent parameter domain and then
continues to a small epsilon neighborhood. The answer is restricted to z>0.

## A source differential identity from an actual total derivative

For general small a,b set

$$\mathscr U(a,b;z)=\int_0^1
\frac{e^{aW}(e^{bL}-1)}{\beta-r}d\beta,\qquad
C(a,b)=\frac{\Gamma(1+a-b)\Gamma(1+a+b)}{\Gamma(1+2a)}.$$

The numerator vanishes at beta=r. To derive its equation, first take r outside
the integration segment and introduce

$$\mathcal C_{s,t}(r)=\int_0^1\frac{x^s(1-x)^t}{x-r}dx.$$

Integrate the total derivative of x(1-x)x^s(1-x)^t/(x-r); the endpoints
vanish initially for Re s,Re t>0. Collecting terms gives

$$r(1-r)\mathcal C'_{s,t}
=[s-(s+t)r]\mathcal C_{s,t}
-\frac{\Gamma(1+s)\Gamma(1+t)}{\Gamma(1+s+t)}.$$

The two Cauchy transforms forming U have equal internal residues. Take the
mean of the two boundary values only as an auxiliary calculation: their
residues cancel in U, so the final object is precisely the ordinary removable
integral, with no new principal-value prescription on the physical source.
Using theta r=r(1-r) now gives

$$\boxed{\theta(\sigma^{-a}\mathscr U(a,b;z))
=\sigma^{-a}[C(a,0)-C(a,b)z^b].}\tag{A1}$$

For A=C(epsilon,0)sigma^(-epsilon), D=sigma^(-epsilon)V0 and
B=sigma^(-epsilon)Y, the b=0 subtraction and b derivative in(A1) imply

$$\theta D=1-A,\qquad\theta B=-\ell A.\tag{A2}$$

The ordered Beta source F of the main text satisfies its coefficientwise
Euler equation. The Gauss reflection combination is

$$G_\epsilon(z)+G_\epsilon(1/z)-1=A,\qquad
G_\epsilon(z)={}_2F_1(1,-\epsilon;1+\epsilon;-z).$$

This follows from Euler's Beta integral for the Gauss function, or its
standard connection formula, with the negative real arguments retained.
Consequently theta²Psi(epsilon,epsilon)=A-1. D is crossing odd and Psi is
crossing even; evaluating thetaPsi at z=1 sets its odd constant to zero, so

$$\theta\Psi(\epsilon,\epsilon)=-D.$$

Equations(A2) then show that
Psi-(B-ell D)-ell²/2 is independent of z. The boundary calculation below fixes
this constant to pi²/2. Since M01=-(ell²+pi²)/2, this proves equation(7).

For the other conversion differentiate(A1) with respect to b at
b=epsilon and b=-epsilon. At those points

$$C(\epsilon,\pm\epsilon)=1,\quad
\partial_bC(\epsilon,\epsilon)=-d,\quad
\partial_bC(\epsilon,-\epsilon)=d,\quad d=J_{2\epsilon}(1).$$

It follows that

$$\theta\{\sigma^{-\epsilon}[2Y_c+2\epsilon dV(2\epsilon)]\}
=-\ell\sigma^{-\epsilon}(z^\epsilon+z^{-\epsilon}).\tag{A3}$$

Direct differentiation of the series defining F and the endpoint J gives
the identical derivative for
Psi(2epsilon,0)-ell[J_{2epsilon}(-z)-J_{2epsilon}(-1/z)]+2M01.
Their difference is a z-independent constant. The next subsection computes
it as zero, establishing equation(8).

## Determine both constants without reading a Laurent table

Initially0<epsilon<1/4 and let u=log z tend to minus infinity. This ensures
both positive exponents epsilon and2epsilon lie between0 and1. There are two
uniform asymptotic tools. The first comes from the Euler integral and Gauss
connection formula:

$$\mathcal C_{s,t}(r)=B(s,t+1)
{}_2F_1(1,-s-t;1-s;r)-\pi\cot(\pi s)r^s(1-r)^t.\tag{B1}$$

It is initially stated for noninteger0<Re s<1, then continued with removable
poles canceled. In particular the s->0 limit must combine the Beta and cot
terms before expanding. The second is the Mellin transform

$$\int_0^\infty y^{v-1}{}_2F_1(1,-a;1+b;-y)dy
=\frac{\Gamma(v)\Gamma(1-v)\Gamma(-a-v)\Gamma(1+b)}
{\Gamma(-a)\Gamma(1+b-v)}.\tag{B2}$$

Two Euler/Beta integrals prove(B2) for Re a<0,Re b>0 and
0<Re v<min(1,-Re a), after which one continues parameters. F is the double
Euler primitive of the Gauss function minus1. The v=0 triple pole controls
its quadratic log and its constant; continuation to0<a<1 retains the moving
power pole as well. Computing those residues gives

$$F(a,b;-y)=\frac{\Gamma(1+a)\Gamma(1+b)}{\Gamma(1+a+b)}
\frac{y^a}{a^2}-\frac12\log^2y+c_{ab}\log y+d_{ab}+o(1),$$
$$c_{ab}=\psi(1+b)-\psi(-a),\quad
d_{ab}=-\frac12[c_{ab}^2+\psi_1(-a)-\psi_1(1+b)+2\zeta_2].\tag{B3}$$

The arguments -y stay on the negative real axis, away from the standard
hypergeometric cut[1,infinity). The small-parameter differentiated remainders
also tend to zero in this domain.

For the first constant set P=pi cot(pi epsilon),
E2=pi² csc²(pi epsilon). Formula(B1) and its b=0 parameter derivative give

$$D=A/\epsilon+u-P+o(1),\qquad
B=A(u/\epsilon+1/\epsilon^2)-E2+o(1).$$

Using(B3), F(epsilon,epsilon;1)=[psi1(1+epsilon)-zeta2]/2 and the polygamma
reflection identity gives

$$\Psi(\epsilon,\epsilon)=A/\epsilon^2-u^2/2+Pu-E2+\pi^2/2+o(1).$$

Hence Psi-(B-uD)-u²/2 tends to pi²/2. Its exact derivative was zero, so this
is the required global constant, not a finite-order numerical adjustment.

For the second constant set a=2epsilon,d=J_a(1),P=pi cot(pi a),
E2=pi² csc²(pi a),T1=T_a(1), and

$$d_T=-\tfrac12[(d-P)^2+E2-\psi_1(1+a)+\zeta_2].$$

The left member of equation(8), after division by sigma^epsilon, has from(B3)
the expansion

$$(1+1/z)^a(u/a+1/a^2)-u^2/2+d_T-2T1-\pi^2+o(1).$$

The right member, computed from the b=+/-epsilon derivatives of(B1), has the
same power term and the same quadratic log. Its constant is

$$[d^2-\psi_1(1+a)-5\zeta_2]/2-E2+dP.$$

In the b=epsilon endpoint the apparent s=0 pole is canceled before taking
r->0. The necessary local expansion is

$$B(s,1+a-s)=1/s+d+s[d^2+\zeta_2+\psi_1(1+a)]/2+O(s^2),$$

multiplied by exp(-su) and combined with -pi cot(pi s). In particular the
5zeta2 above cannot be dropped. Finally

$$T1=[\psi_1(1+a)-\zeta_2-d^2]/2,\qquad E2-P^2=\pi^2$$

makes the two constants exactly equal. This fixes the constant of(A3) to0.
The ordinary removable integrals analytically continue these identities from
the starting real epsilon interval to the local complex epsilon neighborhood.

## Eliminate the last temporary source by an explicit linear substitution

Expand the two powers in U(epsilon,+/-epsilon):

$$e^{\epsilon W}e^{\epsilon L}=[z(1-\beta)^2]^\epsilon,\qquad
e^{\epsilon W}e^{-\epsilon L}=[\beta^2/z]^\epsilon.$$

Subtract and add sigma^epsilon before integrating. For the first power take
t=(beta-r)/(1-r), with d beta=(1-r)dt and bounds[-z,1]. For the second take
t=1-beta/r, with d beta=-r dt and bounds[1,-1/z]. Both transformed numerators
are complete removable differences. The endpoint primitives are exactly

$$\mathscr U(\epsilon,\epsilon)
=\sigma^\epsilon[d-J_{2\epsilon}(-z)]-V_0,$$
$$\mathscr U(\epsilon,-\epsilon)
=\sigma^\epsilon[J_{2\epsilon}(-1/z)-d]-V_0.$$

Their sum is E by its definition in the main text. Adding cancels d and
proves equation(9). Substitution back into the same Q produces VN; no
auxiliary unevaluated integral is retained in the final M/N mother formula.

The mathematical source of this appendix is the existing local
[direct M/N derivation](../../../../benchmarks/2ne/theory/generation/4e04-direct-mn.md),
including its analytic-constant calculation. This new text makes that adopted
lemma and its actual source bridge readable here. The source's historical
regression and current runtime replay remain separate evidence categories.
