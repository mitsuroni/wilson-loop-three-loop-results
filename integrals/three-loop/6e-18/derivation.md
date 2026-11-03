# 6e-18: original ordered J to complete Laurent coefficients

Work at fixed real z>0, ell=log z and D=4-2epsilon, with the inherited source
i0. The accepted root is this representative's original ordered prepared C/J.
The target is the complete physical I=sum_{k=-6}^0 c_k epsilon^k+O(epsilon).
This graph has m=4, a1=1/2 and kappa=-1/8.
The [normalization derivation](../../families/three-loop-6e/common-normalization.md)
shows where each factor enters exactly once.

## The actual source and its bounds

$$J=z\,\int_{0}^{1}ds_{22}\int_{0}^{1}ds_{33}\int_{0}^{1}ds_{41}\int_{0}^{s_{22}}ds_{21}\int_{0}^{s_{33}}ds_{32}\int_{0}^{s_{32}}ds_{31}\;f_{1}^{\epsilon-1}\,f_{2}^{\epsilon-1}\,f_{3}^{\epsilon-1}.$$

$$f_{1}=s_{32} \left(1 - s_{21}\right).$$

$$f_{2}=s_{33} \left(1 - s_{22}\right).$$

$$f_{3}=s_{41} z \left(1 - s_{31}\right).$$

The source outer-to-inner order is `s22,s33,s41,s21,s32,s31`. No cube Jacobian is inserted. Its bounds are equivalent to the following disjoint ordered chains, low to high.

| Variables, low to high | Word | Evaluated integral |
|---|---|---|
| s21 < s22 | BB | $\frac{1}{2 \epsilon^{2}}$ |
| s31 < s32 < s33 | BAA | $\frac{R_{12} - 2 R_2 + 1}{2 \epsilon^{3}}$ |
| s41 | A | $\frac{1}{\epsilon}$ |

A(x)=x^(epsilon-1), B(x)=(1-x)^(epsilon-1), and O retains an opposite endpoint.
Every letter follows from the displayed propagators. Regrouping has Jacobian1.
For Re epsilon>0 the full-cube product bound in the shared lemma is integrable,
so these particular chain separations preserve the original iterated root;
the original bounds are not discarded. Local epsilon continuation follows
only after the integrations have been evaluated.

## Multiply the evaluated chains, then apply crossing

The [chain lemma](../../families/three-loop-6e/adjacent-chains.md) proves the
specific endpoint/Beta formulas in the table. For this graph their actual
product is

$$K=F_{BB}\,F_{BAA}\,F_{A}=- \frac{- R_{12} + 2 R_2 - 1}{4 \epsilon^{6}},\qquad J=z^{1\epsilon}K.$$

Thus r=1 counts the actual adjacent q=z factors. Inserting both scale images into the single physical normalization gives

$$I=- \frac{1}{8}Z^3\left(- \frac{- R_{12} + 2 R_2 - 1}{4 \epsilon^{6}}\right)\frac{\cosh[(- \frac{1}{2})\epsilon\ell]}{\cosh(3\epsilon\ell/2)}.$$

## Expand the regular part deeply enough for c0

$$T=\epsilon^6K=\frac{R_{12}}{4} - \frac{R_2}{2} + \frac{1}{4}.$$

Generate this polynomial's R2,R12 series through epsilon6 from their
explicit Gamma logarithms in the chain lemma. For an exponential series use
E0=1,E_n=sum_{j=1}^n j l_j E_(n-j)/n. If t_j is the ordinary coefficient of
the displayed T, its source-specific crossing numerator produces

$$b_q=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{[(r-3/2)\ell]^{2j}}{(2j)!}t_{q-2j}.$$

Let p_n be the explicit Z³/cosh(3epsilon ell/2) coefficients from the
normalization recurrence. Every requested Laurent coefficient is

$$c_k=\kappa\sum_{j=0}^{k+6}p_j b_{k+6-j}\quad(-6\le k\le0),$$
$$c_0=\kappa[b_6+p_2b_4+p_3b_3+p_4b_2+p_5b_1+p_6b_0].$$

Use this article's numeric r and kappa above. The sufficient regular budget
is order6, including pole times positive-order prefactor products. Actual
source zeros lower the pole only after the displayed Gamma polynomial is
expanded. No parameter subintegral remains. Only after convolution replace
ell²=-2M01-pi² to obtain the compact table below.

## All seven Laurent coefficients of the complete physical I

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=-\frac{1}{16} \zeta_{3}.$$

$$c_{-2}(z)=\frac{13}{160} \zeta_{2} \zeta_{2}.$$

$$c_{-1}(z)=-\frac{1}{8} M_{0,1}(z) \zeta_{3}-\frac{7}{32} \zeta_{2} \zeta_{3}-\frac{15}{16} \zeta_{5}.$$

$$c_{0}(z)=-\frac{1}{2} \zeta_{3} \zeta_{3}+\frac{13}{80} \zeta_{2} \zeta_{2} M_{0,1}(z)+\frac{435}{448} \zeta_{2} \zeta_{2} \zeta_{2}.$$

The remainder is O(epsilon) at fixed z>0. This branch and local continuation do not provide an arbitrary complex-z or z=-1 boundary theorem. See [reproduction](reproduce.md) and [validation](validation.md) for the adopted-input and executed-check boundaries.
