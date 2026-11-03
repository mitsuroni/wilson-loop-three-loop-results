# Evaluate the actual adjacent endpoint chains

This lemma applies only after every source propagator has been checked to be
q(1-u)v with q=1 or z and distinct endpoints. Its numerator is the product of
those q, so the original density becomes

$$z^{r\epsilon}\prod_i f_{c_i}(x_i),\quad
f_A(x)=x^{\epsilon-1},\quad f_B(x)=(1-x)^{\epsilon-1},$$

where r counts q=z. The bounds of the ten adjacent-only three-loop entries
are disjoint chains, explicitly listed in their individual articles. Each
chain is read from its low endpoint to its high endpoint. Define

$$F_{c_1\ldots c_m}=\int_{0<x_1<\cdots<x_m<1}
\prod_i f_{c_i}(x_i)d^mx,\qquad F_\varnothing=1.$$

For Re epsilon>0 the modulus on the full6-cube is bounded by the product of
six integrable adjacent endpoint densities. Its integral is proportional to
(Re epsilon)^(-6). Extending the ordered-domain indicator by zero therefore
proves the particular Tonelli/Fubini regrouping into edge chains; it is not
an unproved permission to permute unrelated prepared roots. Grouping has
Jacobian1 and retains every chain inequality. The same compact-subdomain
bound proves holomorphic dependence on epsilon before local continuation.

## Compute the necessary chains, including the mixed ones

Equal factors are symmetric on the full square/cube. Dividing by their
permutation count gives

$$F_A=F_B=1/\epsilon,\quad
F_{AA}=F_{BB}=1/(2\epsilon^2),\quad
F_{AAA}=F_{BBB}=1/(6\epsilon^3).$$

The two mixed Gamma ratios are

$$R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
R_{12}=\frac{\Gamma(1+\epsilon)\Gamma(1+2\epsilon)}{\Gamma(1+3\epsilon)}.$$

Integrating the smallest A explicitly gives

$$F_{AB}=\frac1\epsilon\int_0^1x^\epsilon(1-x)^{\epsilon-1}dx
=\frac{R_2}{\epsilon^2}.$$

The square is the disjoint union of the two ordered triangles, so
F_A F_B=F_AB+F_BA; hence F_BA=(1-R2)/epsilon². Similarly

$$F_{AAB}=\frac{B(1+2\epsilon,\epsilon)}{2\epsilon^2}
=\frac{R_{12}}{2\epsilon^3}.$$

For BAB fix the middle A position x. The B endpoint integrals on its two
sides are [1-(1-x)^epsilon]/epsilon and(1-x)^epsilon/epsilon. Thus

$$F_{BAB}=\frac1{\epsilon^2}\int_0^1x^{\epsilon-1}
[1-(1-x)^\epsilon](1-x)^\epsilon dx
=\frac{R_2-R_{12}}{\epsilon^3}.$$

For BBA first integrate the two lower B endpoints:

$$F_{BBA}=\frac1{2\epsilon^2}\int_0^1x^{\epsilon-1}
[1-(1-x)^\epsilon]^2dx
=\frac{1-2R_2+R_{12}}{2\epsilon^3}.$$

The map x_i->1-x_(m+1-i) has unit absolute Jacobian, preserves the ordered
domain and reverses the word while interchanging A/B. Consequently
F_ABB=F_AAB, F_ABA=F_BAB and F_BAA=F_BBA. These statements evaluate every
word actually used by the ten entries. No length-four mixed chain or
arbitrary-loop closure is inferred.

## Restore crossing and derive every coefficient

With the individual source's actual chains,

$$J=z^{r\epsilon}K(\epsilon),\qquad K=\prod_{\rm edges}F_w.$$

Substitution into the [single normalization](common-normalization.md) yields

$$I=\kappa Z^3K(\epsilon)
\frac{\cosh[(r-3/2)\epsilon\ell]}{\cosh(3\epsilon\ell/2)}.\tag{A1}$$

The per-entry word products show explicitly why some integrals share K but
have different kappa. They remain different orbit contributions.

Every relevant K is epsilon^(-6) times a polynomial in R2,R12. Their logs are

$$\log R_2=\sum_{k\ge2}\frac{(-1)^k(2-2^k)\zeta_k}{k}\epsilon^k,$$
$$\log R_{12}=\sum_{k\ge2}\frac{(-1)^k(1+2^k-3^k)\zeta_k}{k}\epsilon^k.$$

For E=exp(sum l_k epsilon^k), E0=1 and
E_n=sum_{k=1}^n k l_k E_(n-k)/n. This is an explicit finite calculation of
R2,R12 through6; |epsilon|<1/3 initially suffices for their log series.
Set T(epsilon)=epsilon^6 K and let
b_q be the ordinary coefficient of cosh[(r-3/2)epsilon ell]T. Then

$$c_k=\kappa\sum_{j=0}^{k+6}p_j b_{k+6-j},\quad-6\le k\le0.$$

In particular c0=kappa[b6+p2 b4+p3 b3+p4 b2+p5 b1+p6 b0]. This includes all
pole times positive-order prefactor products. The explicit source polynomial
may have b0,b1 or further zeros; their actual cancellations determine the
leading pole. The tables in the individual articles are the finite expansion
of this formula, followed only by exact ell²-to-M01 replacement.

The native Gamma/word computation uses the copied exchange.py recognition
and scalar expansion. Independent original-order integer-epsilon integration
checks J without evaluating Gamma(1-epsilon) at its poles. The current
portable replay reexecutes recognition and Gamma convolution; it does not
relabel an old complete-I reader as an integration program.
