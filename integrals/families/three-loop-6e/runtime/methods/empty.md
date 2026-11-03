# Adjacent chains: ten three-loop integrals

This route covers 6e-01/02/03/07/08/09/17/18/19/35. With no opposite propagator, each adjacent factor splits into two one-variable endpoint powers. After checking the actual bounds, read every endpoint chain from its smallest to largest parameter and assign A(x)=x^(epsilon-1) or B(x)=(1-x)^(epsilon-1). The problem becomes a product of short ordered words.

For a word w=c1...cm define

\[
F_w(\epsilon)=\int_{0<x_1<\cdots<x_m<1}
\prod_{j=1}^{m}f_{c_j}(x_j)\,d^mx,\qquad F_\varnothing=1.
\]

If r adjacent propagators carry q=z, the source integral is

\[
J(\epsilon,z)=z^{r\epsilon}\prod_{\mathrm{chains}}F_w(\epsilon).
\tag{2}
\]

For Re epsilon>0, the product of unrestricted endpoint integrals bounds the ordered integral by a finite constant times epsilon^(-6). This supplies the absolute-integrability bridge from the original iterated order to (2), with unit Jacobian. No historical cube Jacobian is multiplied into the ordered root.

## Evaluating the words

Put

\[
R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
R_{12}=\frac{\Gamma(1+\epsilon)\Gamma(1+2\epsilon)}{\Gamma(1+3\epsilon)}.
\]

The required words all have length at most three:

| Word | Exact integral |
| --- | --- |
| A or B | 1/epsilon |
| AA or BB | 1/(2 epsilon^2) |
| AAA or BBB | 1/(6 epsilon^3) |
| AB | R2/epsilon^2 |
| BA | (1-R2)/epsilon^2 |
| AAB or ABB | R12/(2 epsilon^3) |
| ABA or BAB | (R2-R12)/epsilon^3 |
| BAA or BBA | (1-2R2+R12)/(2 epsilon^3) |

For example, integrating the smaller A first gives

\[
F_{AB}=\epsilon^{-1}\int_0^1x^\epsilon(1-x)^{\epsilon-1}dx
=R_2/\epsilon^2.
\]

Splitting the square into its two orderings gives F_A F_B=F_AB+F_BA. The same procedure with two smaller A endpoints gives F_AAB=R12/(2 epsilon^3). Reflection xj -> 1-x_(m+1-j) reverses the word and exchanges A/B. It supplies the paired entries in the table without changing the ordered domain.

Write K=product Fw and ell=log z. Substitution of (2) into [the common I definition](common-normalization.md) yields

\[
I=\kappa\mathcal G^3K
\frac{\cosh[(r-3/2)\epsilon\rho]}{\cosh(3\epsilon\rho/2)}.
\tag{3}
\]

For 6e-01 the two words are AAA and BBB, K=1/(36 epsilon^6), r=3 and kappa=-1/8. Equation (3) gives I=-Gcal^3/(288 epsilon^6), reproducing the completed pilot.

For all ten diagrams, put c_j=cosh(j epsilon ell/2). The complete formulas are

\[
I_{01}=-\frac{\mathcal G^3}{288\epsilon^6},\quad
I_{02}=I_{03}=\frac12I_{01},
\]
\[
I_{07}=I_{17}=-\frac{\mathcal G^3c_1}{32\epsilon^6c_3}(R_2-R_{12}),
\]
\[
I_{08}=I_{18}=-\frac{\mathcal G^3c_1}{32\epsilon^6c_3}(1-2R_2+R_{12}),\quad
I_{09}=I_{19}=\frac12I_{08},
\]
\[
I_{35}=-\frac{\mathcal G^3c_1}{16\epsilon^6c_3}(1-R_2)^2.
\]

These are physical-I equalities with their normalization included. They do not replace each diagram's original input.

## The finite-order budget

The conservative budget is the regular scalar jet through epsilon^6. Use

\[
\log\mathcal G=\sum_{k\ge2}\frac{\zeta(k)}k\epsilon^k,\quad
\log R_2=\sum_{k\ge2}\frac{(-1)^k(2-2^k)\zeta(k)}k\epsilon^k,
\]
\[
\log R_{12}=\sum_{k\ge2}\frac{(-1)^k(1+2^k-3^k)\zeta(k)}k\epsilon^k.
\]

For F=exp(sum a_k epsilon^k)=sum f_m epsilon^m, take f0=1 and f_m=(sum_(k=1)^m k a_k f_(m-k))/m. Convolve these jets with the elementary cosh ratio. Apparent epsilon^-6 factors need not survive: R2-R12 starts at epsilon^2, the next combination at epsilon^3, and (1-R2)^2 at epsilon^4. The output records every resulting zero explicitly.

The only kinematic dependence at this order is a polynomial in ell^2, which can be replaced by -2 M01-6 zeta2. This route requires no new function family or CAS integration. The batch actually repeats its chain evaluation and scalar recurrence; the 6e-01 pilot remains an existing independent packet and is not counted as a newly discovered integral.

The preserved [historical chain method](../evidence/historical/benchmarks/2ne/families/empty/method.md) contains the original word derivations and their evaluation checks.
