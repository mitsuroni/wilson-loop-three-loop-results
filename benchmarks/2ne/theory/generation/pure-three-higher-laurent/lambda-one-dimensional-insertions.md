# Complete Lambda: one existing Horn parent and one-dimensional J insertions

2026-10-16. The all-epsilon identity below is derived from the complete original sector before comparing high-order coefficients. Lambda can be represented by the same Horn parent used for W, together with explicit one-dimensional insertions. This does not establish a minimal physical master basis or M/N closure. `formal=not_run`.

## 1. Preserve the complete original sector

Use the [complete Lambda single series](../pure-three-20261015/lambda-single-series.md), with `e=epsilon`:

\[
J_e(x)=\int_0^x\frac{(1-t)^e-1}{t}\,dt,
\quad \Delta J=J_e(Cuv)-J_e(Au)-J_e(Bv)+J_e(D).
\]
\[
\mathcal S_e(A,B,C;D)=\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
\Delta J\,[(1-Cuv)^e-(1-Au)^e][(1-Cuv)^e-(1-Bv)^e].
\tag{C1}
\]

The three oriented physical sectors have `(A,B,C,D)=(1,1,-z,-1/z),(-1/z,-1/z,-1/z,-1/z),(1,-1/z,1,-1/z)` with coefficients `1,1,-2`. Their sum times `-rho^(3e/2)/e^3` gives L; adding its reciprocal-z partner gives Lambda. All four corner terms remain in C1.

Reuse the existing fixed W parent:

\[
\mathscr F_e(C;A,B)=\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
(1-Au)^e(1-Bv)^eJ_e(Cuv).
\tag{C2}
\]

This is the equal-exponent slice of the [three-parameter W Horn source](../pure-three-20261015/w-source-and-coefficients.md), defined independently of any current remainder.

## 2. Eliminate endpoint moments before summing

Let `b_n(c)=(-c)_n/n!`, `h_n(c)=[x^n](1-x)^c J_e(x)`, and write

\[
B_A=\int_0^1u^{n-1}(1-Au)^e du,\quad
H_A=\int_0^1u^{n-1}(1-Au)^eJ_e(Au)du,
\quad J_A=J_e(A),
\]

Use the corresponding notation for B. The old complete moment formula contains a moment which integrates directly:

\[
Q_A=\frac{J_A}{n}-\frac{B_A}{n}+\frac1{n^2}.
\]

Substituting this expression and retaining all corners gives the complete coefficient at nonzero C index:

\[
\begin{aligned}
s_n={}&\frac{h_n(2e)}{n^2}-\frac{h_n(e)}n(B_A+B_B)\\
&+b_n(2e)\left[\frac{B_A+B_B+J_D-J_A-J_B}{n^2}-\frac2{n^3}\right]\\
&+b_n(e)\left[\frac{H_A+H_B-B_AB_B+(J_B-J_D)B_A+(J_A-J_D)B_B}{n}
+\frac{B_A+B_B}{n^2}\right].
\end{aligned}
\tag{C3}
\]

The only product of two endpoint moments is `-b_n(e) B_A B_B/n`. Summing it gives **minus the existing parent C2**. All other terms require a single endpoint coordinate and therefore assemble into one-dimensional integrals. This follows from the original moment expansion. The [actual Wolfbook/MCP check](../../../../../workflow/checks/pure-three-higher-laurent-20261016/lambda-moment-wolfbook.json) verifies C3 and rejects the wrong parent sign.

## 3. Complete one-dimensional correction

Define primitives of the known J and two endpoint insertion sources:

\[
K_e(x)=\int_0^x\frac{J_e(t)}t\,dt,\qquad
L_e(x)=\int_0^x\frac{K_e(t)}t\,dt,
\]
\[
V_e(C)=\int_0^C\frac{dx}{x}\int_0^x\frac{(1-t)^{2e}J_e(t)}t\,dt,
\quad
Y_e(A,C)=\int_0^1\frac{(1-Au)^eJ_e(Cu)}u\,du,
\]
\[
Z_e(A,C)=\int_0^1\frac{(1-Au)^e}{u}
\left[J_e(Au)J_e(Cu)-\frac12J_e(Cu)^2+K_{2e}(Cu)\right]du.
\tag{C4}
\]

`J` is the explicitly known `-e*x*3F2(1,1,1-e;2,2;x)`; K and L are its further Mellin primitives. None of the source coefficients below is left as an unevaluated integral.

Then the complete sector is

\[
\boxed{\begin{aligned}
\mathcal S_e(A,B,C;D)={}&J_AJ_B\left[J_D-\frac{J_A+J_B}{2}\right]-\mathscr F_e(C;A,B)\\
&+Z_e(A,C)+Z_e(B,C)\\
&+(J_B-J_D)Y_e(A,C)+(J_A-J_D)Y_e(B,C)\\
&+V_e(C)+(J_D-J_A-J_B)K_{2e}(C)-2L_{2e}(C).
\end{aligned}}
\tag{C5}
\]

The first line includes the entire C-zero mode. One useful simplification in deriving C4 is

\[
\int_0^x\frac{(1-t)^eJ_e(t)}t\,dt=\frac{J_e(x)^2}{2}+K_e(x).
\]

It cancels an apparent separate K_e contribution to Z. Every term of C5 has now been specified by an already solved Horn parent or a fixed one-variable source. Individual terms can start below e cubed; they must be recombined before dividing by e cubed in the physical expression.

The coefficient proof first holds for small `|A|,|B|,|C|`, with a common analytic neighborhood of `e=0`. At the physical arguments all power bases are positive in the open interval. Endpoints with A, B or C equal to 1 are integrable for sufficiently small `|Re(e)|`, for example below 1/8. J starts linearly at zero; it removes all apparent lower endpoint denominators. These bounds also control finite parameter derivatives. They extend C5 along the original real branch to the physical sectors, without treating a Taylor series outside its disk as convergent.

## 4. The new natural source is a product of two J insertions

For independent a,b,c define

\[
\mathcal B_{X;Y,Z}(a,b,c)=\int_0^1\frac{(1-Xu)^aJ_b(Yu)J_c(Zu)}u\,du.
\tag{C6}
\]

In physical sectors X,Y,Z take values 1 or -z, together with their reciprocal partners. The two J arguments can be unordered because their product is symmetric, leaving six argument patterns. A zero elementary exponent is allowed.

An additional integration by parts uses `u*dJ_e(Au)/du=(1-Au)^e-1` and `u*dK_2e(Cu)/du=J_2e(Cu)`. It gives

\[
\int_0^1\frac{(1-Au)^eK_{2e}(Cu)}u\,du
=J_e(A)K_{2e}(C)-\mathcal B_{1;A,C}(0,e,2e)+L_{2e}(C).
\]

The subscript 1 on the middle term is immaterial because its first exponent is zero. Consequently

\[
\boxed{Z_e(A,C)=\mathcal B_{A;A,C}(e,e,e)
-\frac12\mathcal B_{A;C,C}(e,e,e)-\mathcal B_{1;A,C}(0,e,2e)
+J_e(A)K_{2e}(C)+L_{2e}(C).}
\tag{C7}
\]

Thus the new nonlinear one-dimensional insertion in C5 is a specific slice of a three-parameter family. Its source position, boundary and parameter values are fixed before coefficient fitting.

Its every Taylor coefficient is a finite GPL expression:

\[
[a^ib^jc^k]\mathcal B_{X;Y,Z}
=G\!\left(0,\operatorname{sh}\big[(1/X)^i,\,(0,(1/Y)^j),\,(0,(1/Z)^k)\big];1\right),
\quad i\ge0,\ j,k\ge1.
\tag{C8}
\]

Here `sh` is the shuffle sum with multiplicities; the notation means prepend 0 to each shuffled word. Its weight is `i+j+k+3`. Repeated-letter GPLs already include the Taylor factorials. The same argument gives fixed finite rules for Y,V,K,L. At weight at most nine, all six patterns give 158 coefficients after the exact interchange symmetry of equal J arguments. These coefficients are generated from C8, with no target data.

## 5. Why ordinary double-Beta convolution was tested first

The complete-Beta slice of C2 has an ordinary `5F4` value. Releasing its two lower endpoint indices naturally gives the independent double-Euler source

\[
\mathcal E_2(a,b,c,d,e;z)=\int_0^1\!\int_0^1
u^{b-1}(1-u)^{d-b}v^{c-1}(1-v)^{e-c}[(1+zuv)^{-a}-1]du\,dv.
\]
\[
{}_3F_2(a,b,c;1+d,1+e;-z)-1
=bc\,g(b,d)g(c,e)\mathcal E_2,\quad
g(b,d)=\frac{\Gamma(1+d)}{\Gamma(1+b)\Gamma(1+d-b)}.
\tag{C9}
\]

All five parameters are independent. Its third-order ODE and origin boundary generate every jet:

\[
[\theta(\theta+d)(\theta+e)+z(\theta+a)(\theta+b)(\theta+c)]F=0,\quad F(0)=1.
\]

The new calculation generates all 462 raw jets through weight nine, reduced by upper S3 and lower S2 symmetry to 97, plus 120 Mellin jets. Twelve independent Pochhammer/Mellin comparisons and four original Beta moments pass. Adding their full homogeneous product algebra leaves the previous weight-nine even rank at 373; the two physical targets still increase it to 375. [Exact declared comparison](../../../../../workflow/checks/pure-three-higher-laurent-20261016/beta-convolution-runs/20261016T034627077157Z/comparison.json).

This is a scoped negative result for ordinary complete-Beta convolution. It neither tests all independently weighted ordered two-propagator kernels nor excludes one-dimensional sources such as C6. The latter have two independently inserted J functions, a different coupling from merely releasing ordinary 3F2 parameters.

## 6. Fixed single-index generators for the six J patterns

The one-dimensional integrals C6 also admit explicit coefficient sums with a fixed rational recurrence. This avoids introducing an unspecified new integral at each epsilon order. Define

\[
\beta_n(r)=\frac{(-r)_n}{n!},\quad j_n(r)=\frac{\beta_n(r)}n\ (n\ge1),\quad j_0=0,
\quad h_n(r,s)=[x^n](1-x)^rJ_s(x),\ h_0=0.
\]
\[
\beta_n(r)=\frac{n-1-r}{n}\beta_{n-1}(r),\qquad
n h_n(r,s)+(r-n+1)h_{n-1}(r,s)=\beta_n(r+s+1)-\beta_n(r+1).
\tag{C10}
\]

The coefficient of a product of two J functions needs no convolution at evaluation time:

\[
d_n(b,c)=[x^n]J_b(x)J_c(x)
=\frac{h_n(b,c)+h_n(c,b)-j_n(b)-j_n(c)}n.
\tag{C11}
\]

For the remaining elementary power set `g_n(a;b,c)=[x^n](1-x)^a J_b(x)J_c(x)` and

\[
r_n=h_n(a+b,c)+h_n(a+c,b)-h_n(a,c)-h_n(a,b).
\]
\[
n g_n+(a-n+1)g_{n-1}=r_n-r_{n-1},\qquad g_0=0.
\tag{C12}
\]

These identities follow by differentiating the specified products, so their initial data are fixed. The required endpoint moments also have explicit initial data and transport:

\[
H_n(a,b)=\int_0^1u^{n-1}(1-u)^aJ_b(u)du,
\quad H_1=\frac{\psi(a+2)-\psi(a+b+2)}{a+1},
\]
\[
(n+a+1)H_{n+1}=nH_n+B(n,a+b+2)-B(n,a+2).
\tag{C13}
\]
\[
Q_n(b)=\frac{J_b(1)}n-\frac{B(n,b+1)}n+\frac1{n^2},
\quad J_b(1)=-\gamma_E-\psi(1+b),
\]
\[
P_n(b,c)=\frac{J_b(1)J_c(1)-H_n(b,c)-H_n(c,b)+Q_n(b)+Q_n(c)}n.
\tag{C14}
\]

Here P is the already evaluated moment of `J_b(u)J_c(u)`. Beta values obey their usual rational n transport. Every coefficient in C10–C14 is therefore a finite Gamma/digamma evaluation or a fixed rational recurrence with explicit initial data.

Write `B_stu` for C6 with `X_0=1,X_1=-z`. The six complete sources become

| Pattern | Explicit value, initially for abs(z)<1 |
|---|---|
| 000 | `sum(n>=1) g_n(a;b,c)/n` (independent of z) |
| 001 | `sum(n>=1) j_n(c) H_n(a,b) (-z)^n` |
| 011 | `sum(n>=1) d_n(b,c) B(n,a+1) (-z)^n` |
| 100 | `B_000(0,b,c) + sum(n>=1) beta_n(a) P_n(b,c) (-z)^n` |
| 101 | `sum(n>=1) h_n(a,c) Q_n(b) (-z)^n` |
| 111 | `sum(n>=1) g_n(a;b,c) (-z)^n/n` |

The period in the 100 row is the explicit convergent unit series in the 000 row, with a=0. It is not a free boundary constant. The derivation uses an ordinary power expansion of the original elementary/J factors, followed by endpoint moments. The zero lower endpoints fix the boundaries of the rows with a moving J.

For general positive z, the original sources are holomorphic in `Re(z)>-1` near the parameter origin, with the same integrable bounds as C6. Set `w=z/(z+2)`. If `sum c_n(-z)^n` denotes any row, its negative-axis continuation is

\[
c_0+\sum_{k\ge1}w^k\sum_{n=1}^k(-2)^n\binom{k-1}{n-1}c_n.
\tag{C15}
\]

The map sends the half-plane to the unit disk, hence specifies a convergent series for every positive z. Unit-period sums can converge slowly near zero exponents, and the binomial transform can suffer numerical cancellation; no uniform performance or precision bound is asserted.

## 7. Validation and open structural question

Three independent symbolic polynomial checks at integer e=1,2,3 preserve arbitrary A,B,C,D and verify C5 directly from the original two-variable integral. Wrong Horn signs and omitted zero modes fail. Two noninteger cases e=0.37 and -0.13 compare the complete original sector with C5; their observed differences are below 4.25e-61 at 60 digits. These numbers are empirical discrepancies, not certified error bounds.

The fresh [all-z HyperInt comparison](../../../../../workflow/checks/pure-three-higher-laurent-20261016/lambda-compression-fibration-N7/status.json) checks each physical sector independently, retaining common endpoint limits. All 21 differences at source degrees one through seven vanish exactly, including the degree required for physical epsilon four. Comparing unsimplified GPL word dictionaries alone was insufficient: nonzero word differences at fixed physical letters vanished after the full fibration transformation. That preliminary word test is retained and is not reported as an exact failed integral identity.

All 158 C6 moments through weight nine and both crossing components have now been computed; 16 independent original-integral controls pass. The single-series rules also pass six noninteger-source comparisons, exact rational convolution checks, and two negative-axis Cayley tests at z=2. [Actual generator checks](../../../../../workflow/checks/pure-three-higher-laurent-20261016/bilinear-j-single-series-verification.json).

Adding the full weighted bilinear-J family closes the two previously missing weight-nine physical/source directions. After removing all old three-kernel and E/tau generators, the same declared algebra still reconstructs all physical and bare orders zero through three exactly. Without the elementary endpoint power, the tested subspace remains rank 381 to 383; the full family is 399 to 399. [Natural-source relations](natural-source-reconstruction.md) distinguish these finite statements from the all-epsilon identity C5.

A fixed all-epsilon replacement of W's physical Horn slice by bilinear J sources remains unproved. The full freely decorated ordered two-propagator Q family has not been exhaustively compared. No conclusion about the minimal number of master functions or minimal integration dimension is made here.

The subsequent [shared physical generator](physical-z-transport.md) connects the complete Lambda parent sum directly to the V_n square of the existing W moment system; its reciprocal is the U_n square. It gives fixed one-dimensional parameter curves for every remainder term, the complete normalized I20/I21 formulas, and a first-order physical-z Fuchsian realization with explicit boundary data. Thus a shared all-epsilon representation of the two physical integrals is now available without assuming a BJ-only identity. The finite BJ reconstructions and the stronger unresolved minimality questions retain their scopes above.
