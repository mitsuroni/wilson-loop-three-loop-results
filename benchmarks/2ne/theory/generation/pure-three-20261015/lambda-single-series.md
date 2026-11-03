# The complete three-kernel Lambda source as one explicit series

2026-10-15. This note evaluates every coefficient of a single-index representation of the complete source `Lambda`, including all four corner terms and all three sign sectors. The remaining sum is a specified analytic function with rational coefficient transport and explicit initial data. It is not asserted to be a finite combination of ordinary univariate hypergeometric functions, a minimal master basis, or a newly verified complete physical three-loop result.

The starting point is the full source (E3) in [the existing three-kernel derivation](../three-kernel-eps2.md). Existing results and their normalization are unchanged. All calculations below precede comparison with any stored Hlog coefficients. Exact and noninteger controls are in [lambda_check.py](../../../../../workflow/checks/pure-three-generators-20261015/lambda_check.py); the actual report is [lambda-verification-r001.json](../../../../../workflow/checks/pure-three-generators-20261015/lambda-verification-r001.json). Formal proof status is `not_run`.

## 1. Preserve the complete physical source

Write `e=epsilon`, `rho=z/(1+z)^2`, and

$$
J_e(x)=\int_0^x\frac{(1-t)^e-1}{t}\,dt
=-ex\,{}_3F_2(1,1,1-e;2,2;x).
$$

The actual source is

$$
\mathcal L(e,z)=-\frac{\rho^{3e/2}}{e^3}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}
\left[J_e(-pq/z)-J_e(p/z)-J_e(q/z)+J_e(-1/z)\right]
$$
$$
\hspace{9mm}\times
\left[(1+pq/z)^e-(1-p/z)^e\right]
\left[(1+pq/z)^e-(1-q/z)^e\right],
\qquad \Lambda(e,z)=\mathcal L(e,z)+\mathcal L(e,1/z).
\tag{L1}
$$

The apparent factors `1/p` and `1/q` are removable in the complete integrand. The two power differences and the corner difference start at order `e`, so the entire numerator starts at `e^3`. Terms must be recombined before interpreting any Laurent pole.

Introduce a sector function for independent arguments:

$$
\mathcal S_e(A,B,C;d)=\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
\left[J_e(Cuv)-J_e(Au)-J_e(Bv)+J_e(d)\right]
$$
$$
\hspace{8mm}\times
\left[(1-Cuv)^e-(1-Au)^e\right]
\left[(1-Cuv)^e-(1-Bv)^e\right].
\tag{L2}
$$

This integral is the derivation of the function, not its final unevaluated representation. Equations (L7)-(L11) below give its value without coefficient integrals. Splitting `p,q` at zero and retaining the orientations gives the complete source:

$$
\boxed{\mathcal L(e,z)=-\frac{\rho^{3e/2}}{e^3}
\left[\mathcal S_e(1,1,-z;-1/z)
+\mathcal S_e(-1/z,-1/z,-1/z;-1/z)
-2\mathcal S_e(1,-1/z,1;-1/z)\right].}
\tag{L3}
$$

The minus sign of the mixed sector follows from `dp/p` on `[-1,0]`. The interchange symmetry of (L2) combines its two mixed sectors. Thus the four corners and every region of the original square remain in (L3).

## 2. Three coefficient sequences suffice for the product-variable dependence

Define the explicit coefficients

$$
b_n(c)=\frac{(-c)_n}{n!},\qquad
j_0=0,\quad j_n=\frac{b_n(e)}n\quad(n\ge1),
$$
$$
h_n(c)=\sum_{r=1}^n\frac{b_r(e)}r b_{n-r}(c),\qquad h_0(c)=0.
\tag{L4}
$$

They are the Taylor coefficients of `(1-x)^c`, `J_e(x)` and `(1-x)^c J_e(x)`. In particular all sums in (L4) are finite; no integral or unknown period is hidden in a coefficient.

For evaluation the convolution is unnecessary. The first-order identity

$$
\left[(1-x)\theta_x+cx\right]\left[(1-x)^cJ_e(x)\right]
=(1-x)^{c+e+1}-(1-x)^{c+1}
$$

gives the rational transport rules

$$
b_n(c)=\frac{n-1-c}{n}b_{n-1}(c),\qquad
nh_n(c)+(c-n+1)h_{n-1}(c)
=b_n(c+e+1)-b_n(c+1).
\tag{L5}
$$

Only `c=e` and `c=2e` are required, driven by the five binomial sequences with exponents `e,2e,e+1,2e+1,3e+1`.

## 3. Every endpoint moment has already been evaluated

For positive integer `n`, temporarily write the endpoint moments as

$$
B_n^c(A)=\int_0^1u^{n-1}(1-Au)^cdu
=\frac1n{}_2F_1(-c,n;n+1;A),
$$
$$
Q_n(A)=\int_0^1u^{n-1}J_e(Au)du
=\frac{J_e(A)}n-\frac{B_n^e(A)}n+\frac1{n^2},
$$
$$
H_n(A)=\int_0^1u^{n-1}(1-Au)^eJ_e(Au)du.
\tag{L6}
$$

The last integral also has a finite standard-function value. Direct differentiation gives

$$
\frac{d}{dv}
\frac{v^{a+1}J_e(1-v)+J_{a+1}(1-v)-J_{a+e+1}(1-v)}{a+1}
=v^aJ_e(1-v).
$$

Substitute `v=1-Au` and expand the integer power `(1-v)^(n-1)`. The upper endpoint `v=1` contributes zero, and hence

$$
\boxed{H_n(A)=-A^{-n}\sum_{m=0}^{n-1}(-1)^m\binom{n-1}{m}
\frac{(1-A)^{e+m+1}J_e(A)+J_{e+m+1}(A)-J_{2e+m+1}(A)}{e+m+1}.}
\tag{L7}
$$

This is a finite combination of specified `3F2` functions for every integer `n`. At `A=1` take the common endpoint limit; for `Re(e)>0` all powers displayed here vanish. At `A=0` use the removable limit. The physical arguments are nonzero.

An alternative exact evaluation uses rational transport, keeping the boundary terms:

$$
A(n+c+1)B_{n+1}^c(A)=nB_n^c(A)-(1-A)^{c+1},
\qquad
B_1^c(A)=\frac{1-(1-A)^{c+1}}{A(c+1)},
$$
$$
\boxed{A(n+e+1)H_{n+1}(A)=nH_n(A)
+B_n^{2e+1}(A)-B_n^{e+1}(A)
-(1-A)^{e+1}J_e(A),}
$$
$$
H_1(A)=-\frac{(1-A)^{e+1}J_e(A)+J_{e+1}(A)-J_{2e+1}(A)}{A(e+1)}.
\tag{L8}
$$

Both recurrences follow by integrating derivatives of `u^n(1-Au)^(c+1)` or `u^n(1-Au)^(e+1)J_e(Au)`. They preserve the actual upper endpoint. Upward recurrence can be numerically ill-conditioned for small `|A|`; the finite formula or a series may then be preferable. Exact coefficient transport is unaffected by this numerical issue.

## 4. The complete single-index sum

Abbreviate

$$
J_A=J_e(A),\quad J_B=J_e(B),\quad J_d=J_e(d),\quad
B_A=B_n^e(A),\quad B_B=B_n^e(B),
$$

and let `Q_A,Q_B,H_A,H_B` denote the corresponding order-`n` quantities in (L6)-(L8). The complete coefficient is

$$
\begin{aligned}
s_n(A,B;d)={}&\frac{h_n(2e)}{n^2}
-\frac{h_n(e)}n(B_A+B_B)+j_n B_AB_B\\
&-\frac{b_n(2e)}n(Q_A+Q_B)
+\frac{b_n(e)}n(H_A+H_B)\\
&+b_n(e)(Q_AB_B+B_AQ_B)\\
&+J_d\left[\frac{b_n(2e)}{n^2}-\frac{b_n(e)}n(B_A+B_B)\right].
\end{aligned}
\tag{L9}
$$

The sum of all zero-product-index terms can be evaluated without introducing another source. The identity

$$
\int_0^A\frac{(1-x)^eJ_e(x)}x dx
=\frac12J_e(A)^2+\int_0^A\frac{J_e(x)}x dx
$$

cancels the two primitive terms in the complete zero mode. Its result is

$$
s_0=J_AJ_B\left[J_d-\frac{J_A+J_B}{2}\right].
\tag{L10}
$$

Therefore the full sector is

$$
\boxed{\mathcal S_e(A,B,C;d)=
J_AJ_B\left[J_d-\frac{J_A+J_B}{2}\right]
+\sum_{n=1}^{\infty} C^n s_n(A,B;d).}
\tag{L11}
$$

Every coefficient in this one sum is given by finite standard-function evaluations or the fixed recurrences (L5),(L8). No remaining `u,v` coefficient integral occurs. The factorization into endpoint moments follows by expanding only the product variable `Cuv` in the four terms

$$
(1-Cuv)^{2e}-(1-Au)^e(1-Cuv)^e
-(1-Bv)^e(1-Cuv)^e+(1-Au)^e(1-Bv)^e,
$$

and multiplying each by all four corner terms in (L2). The apparently divergent zero-index terms are combined before integration. At nonzero `n`, the factors `u^(n-1),v^(n-1)` make each moment separately integrable.

## 5. Domain, boundary data, and what this establishes

One may first take `Re(e)>0`, real `A,B<1`, and `|C|<1`, with the original branch of each positive-base power. Endpoint values `A=1` or `B=1` follow by their integrable limits. The coefficient formula analytically continues the endpoint arguments independently. Its `C` continuation is along the real interval from zero to the physical value; the original integrand fixes this branch. Positive `C=1` is reached by its integrable endpoint limit, while negative values cross no physical cut. A bare Taylor sum outside its convergence disk is not the value of (L11).

The physical sector arguments in (L3) require negative values `C=-z,-1/z` and the endpoint `C=1`. Negative-axis continuation can be given as an explicit convergent series rather than left as an unspecified operation. Set

$$
w=\frac{C}{C-2},\qquad C=-\frac{2w}{1-w},\qquad
\widehat s_k=\sum_{n=1}^k(-2)^n\binom{k-1}{n-1}s_n.
$$
$$
\boxed{\mathcal S_e(A,B,C;d)=s_0+\sum_{k\ge1}w^k\widehat s_k.}
\tag{L12}
$$

The complete sector integral is holomorphic in `Re(C)<1` in the initial convergence domain: the real part of `1-Cuv` stays positive, and its powers and the chosen `J` branch are holomorphic with a locally uniform integrable bound. The Cayley map above sends the unit disk to this half-plane. Consequently (L12) converges for every negative real `C`, since then `0<w<1`. Its coefficients follow from expanding `C^n=(-2)^n w^n(1-w)^(-n)`, so no new boundary data are needed. The independent endpoint arguments `A,B` retain their values; they are not replaced by `w`.

At physical `C=1`, use the original series (L11). For `Re(e)>0`, the relevant binomial sequences and their `J` convolutions are absolutely summable, while the finite endpoint moments have an `O(1/n)` bound; this yields absolute convergence of (L9)-(L11) at that endpoint. The complete `e=0` value and coefficients follow by the same meromorphic continuation and cancellations as the original source. A termwise Laurent expansion must still retain the full `e^3` onset.

Thus all positive-real `z` in (L3) have an explicit convergent series prescription in the initial `e` domain. Numerical efficiency is separate: the binomial transform can have severe cancellation, and `w` approaches one for a large negative `C`. No uniform fast-summation or low-precision stability claim is made. The actual two-sector check at `z=2,e=0.37` is recorded in [lambda-cayley-r001.json](../../../../../workflow/checks/pure-three-generators-20261015/lambda-cayley-r001.json).

Products of the finitely many sequences in (L5),(L8), together with rational factors in `n`, admit a finite rational recurrence realization by tensor products. Consequently the complete `C` series has a finite differential realization with explicit initial data in the stated function field. The recurrences and their initial values are explicit here; the minimal scalar differential operator and the pullback to the physical one-variable `z` slice have not been computed in this note. A matrix dimension of an unreduced realization would not be a count of independent physical master integrals.

The physical graphs are recovered only after combining this source with the independently treated `W` and the known one-kernel `H`:

$$
I_{20}=\frac{P}{16}\left(\mathcal W+\Lambda-\frac{H^3}{2}\right),\qquad
I_{21}=I_{22}=\frac{P}{32}\left(\frac{H^3}{3}-\mathcal W-\frac\Lambda2\right),
\quad P=\frac{[e^{-\gamma_Ee}\Gamma(1-e)]^3}{\cosh(3e\log z/2)}.
$$

This work establishes a full-source single-series representation and exact finite transport for `Lambda`. It does not replace the existing separately verified epsilon-squared coefficients or by itself close the remaining physical compact-generation milestone.

## 6. Actual validation scope

The main script checks the exact mixed endpoint primitive, both transport identities, twelve symbolic convolution recurrences, and complete integer-`e` sector polynomials. It includes all three physical sector configurations at `z=1/2,1,2`, and a wrong zero-mode control. Independent noninteger endpoint quadratures check the mixed moments, and a complete interior sector is compared with the single series. These are analytic/algebraic and numerical controls, not a formal proof.

At `e=1/2,z=1`, the full original `p,q` source (all three oriented sectors, not an isolated subblock) was compared independently with 2000/4000 single-series terms. Resolving the original square-root endpoints by `u=1-x^2,v=1-y^2` gives

$$
\mathcal L_{\rm series}=-3.40123973062680275835\ldots,\qquad
\mathcal L_{\rm original}=-3.40123973062743250641\ldots.
$$

Their observed difference is about `6.30e-13`; the 2000/4000 change is `6.50e-12`, while the 48/96-node original-kernel quadrature change is `1.08e-20`. These differences are numerical convergence evidence, not rigorous remainder bounds. The actual record is [lambda-physical-r001.json](../../../../../workflow/checks/pure-three-generators-20261015/lambda-physical-r001.json). The initial direct-coordinate quadrature was insufficient at its requested tolerance; its unsuccessful report is preserved separately and was not accepted as a passing check.
