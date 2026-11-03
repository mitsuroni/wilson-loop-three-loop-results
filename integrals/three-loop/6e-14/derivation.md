# 6e-14: original ordered J to complete Laurent coefficients

Work at fixed real z>0, ell=log z and D=4-2epsilon, with the inherited source
i0. The accepted root is this representative's original ordered prepared C/J.
The target is the complete physical I=sum_{k=-6}^0 c_k epsilon^k+O(epsilon).
This graph has m=4, a1=1/2 and kappa=-1/8.
The [normalization derivation](../../families/three-loop-6e/common-normalization.md)
shows where each factor enters exactly once.

## The actual source and its bounds

$$J=z \left(- z - 1\right)\,\int_{0}^{1}ds_{32}\int_{0}^{1}ds_{22}\int_{0}^{1}ds_{42}\int_{0}^{s_{32}}ds_{31}\int_{0}^{s_{22}}ds_{21}\int_{0}^{s_{42}}ds_{41}\;f_{1}^{\epsilon-1}\,f_{2}^{\epsilon-1}\,f_{3}^{\epsilon-1}.$$

$$f_{1}=s_{32} \left(1 - s_{21}\right).$$

$$f_{2}=s_{22} s_{41} z + \left(1 - s_{22}\right) \left(1 - s_{41}\right).$$

$$f_{3}=s_{42} z \left(1 - s_{31}\right).$$

The source outer-to-inner order is `s32,s22,s42,s31,s21,s41`. No cube Jacobian is inserted. Its bounds are equivalent to the following disjoint ordered chains, low to high.

For this normalized-factor table, t=s22,s=s41. Define the endpoint powers

$$A_t=t^\epsilon,\quad B_t=(1-t)^\epsilon,\quad A_s=s^\epsilon,\quad B_s=(1-s)^\epsilon.$$

These normalized powers differ from the integration densities A(x),B(x) below.

| Variables, low to high | Word | Normalized endpoint factor |
|---|---|---|
| s21 < s22 | BO | $1 - B_t$ |
| s31 < s32 | BA | $1 - R_2$ |
| s41 < s42 | OA | $1 - A_s$ |

A(x)=x^(epsilon-1), B(x)=(1-x)^(epsilon-1), and O retains an opposite endpoint.
Every letter follows from the displayed propagators. Regrouping has Jacobian1.
For Re epsilon>0 the full-cube product bound in the shared lemma is integrable,
so these particular chain separations preserve the original iterated root;
the original bounds are not discarded. Local epsilon continuation follows
only after the integrations have been evaluated.

## Integrate the four adjacent endpoints while retaining their placement

Retain t=s22,s=s41. Use the normalized endpoint powers A_t,B_t,A_s,B_s defined above (AT,BT,AS,BS in the native program). In this table, one epsilon multiplies each integrated adjacent endpoint.

The actual lower/upper primitives are epsilon*integral_0^x A=x^epsilon,
epsilon*integral_x^1 A=1-x^epsilon, epsilon*integral_0^x B=1-(1-x)^epsilon
and epsilon*integral_x^1 B=(1-x)^epsilon. Equal repeated endpoints give half
the square of these factors by ordered-square symmetry. Whole AB and BA
edges contribute R2 and1-R2, respectively, after their epsilon² normalization.
Applying exactly these rules at the retained points in the table gives

$$W=- \left(A_s - 1\right) \left(B_t - 1\right) \left(R_2 - 1\right).$$

$$J=-\frac{(1+z)z^{1\epsilon}}{\epsilon^4}\int_0^1dt\int_0^1ds\,D_z(t,s)^{\epsilon-1}\left(- A_s B_t R_2 + A_s B_t + A_s R_2 - A_s + B_t R_2 - B_t - R_2 + 1\right).$$

Here D_z=(1-t)(1-s)+zts. The four epsilon denominators came from the
actual adjacent endpoints. The [shared derivation](../../families/three-loop-6e/single-opposite.md)
proves these partial-word formulas, the convergent regrouping and the
complete evaluation of the remaining opposite kernel.

## Substitute the evaluated opposite kernels and keep the actual phases

The remaining phase is delta=r-1=0. The actual physical paired weight is

$$\mathscr P=- \frac{A_s B_t R_2}{2} + \frac{A_s B_t}{2} + \frac{A_s R_2}{2} - \frac{A_s}{2} - \frac{A_t B_s R_2}{2} + \frac{A_t B_s}{2} + \frac{A_t R_2}{2} - \frac{A_t}{2} + \frac{B_s R_2}{2} - \frac{B_s}{2} + \frac{B_t R_2}{2} - \frac{B_t}{2} - R_2 + 1.$$

For each monomial AT^a BT^b AS^c BS^d use the solved H_ab;cd defined and evaluated in the shared lemma. Its exact reflection and t<->s laws reduce this particular paired polynomial to

$$\mathcal F=- H_{0}(\epsilon,z) R_2 + H_{0}(\epsilon,z) + H_{1}(\epsilon,z) R_2 - H_{1}(\epsilon,z) + H_{1}(\epsilon,1/z) R_2 - H_{1}(\epsilon,1/z) - H_{X}(\epsilon,z) R_2 + H_{X}(\epsilon,z).$$

$$I=- \frac{1}{8}\frac{Z^3}{\epsilon^4\cosh(3\epsilon\ell/2)}\mathcal F.$$

This entry requires H0, H1, HX. Every ordinary coefficient h_iq for q0..4 is printed as a
complete HLOG expression in the [jet appendix](../../families/three-loop-6e/appendices/single-opposite-jets.md).
The actual integrated kernel is wbBase*Lambda_ab;cd^q/q!, obtained by
 t=x/(1+x),s=y/(1+y) with its explicit Jacobian; integration is in[x,y],
0..infinity. The shared derivation gives this full density, the endpoint
primitive, regularized boundary evaluation, MB/MeijerG identity and branch.
The current replay matches and parses those historical raw outputs; no new
Maple integration or old complete-I input is implied.

## Assemble c0 with every required prefactor term

Let d_q be the ordinary coefficient of the explicit solved mathcal F above.
A coefficient b(epsilon), a displayed phase exp(delta epsilon ell), and an
H jet contribute the finite sum

$$d_q=\sum_{j+u+v=q}b_j\frac{(\delta\ell)^u}{u!}h_v(z),$$

with every reflected term and its opposite phase retained exactly as in
mathcal F. For30/31 use the preceding I-only averaged polynomial. No
unevaluated integral is being named d_q: all h_v are the explicit solved
appendix functions, and b_j comes from the displayed polynomial and R2.
Then the complete physical coefficients are

$$c_{-6}=c_{-5}=0,\qquad
c_k=\kappa\sum_{j=0}^{k+4}p_j d_{k+4-j}\quad(-4\le k\le0),$$
$$c_0=\kappa[d_4+p_2d_2+p_3d_1+p_4d_0].$$

A sufficient budget is all needed H jets through4, R2 through4 if present,
phases through4 and the physical prefactor through4. Further zeroes are
calculated from these actual source weights. Compact M/N conversion is an
exact whole-HLOG/period dictionary substitution after integration and
convolution, not a numerical coefficient fit.

## All seven Laurent coefficients of the complete physical I

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=0.$$

$$c_{-2}(z)=0.$$

$$c_{-1}(z)=0.$$

$$c_{0}(z)=-\frac{1}{32} \zeta_{2} M_{2,1}(z)-\frac{1}{96} \zeta_{2} M_{0,1}(z) M_{0,1}(z)-\frac{1}{16} \zeta_{2} \zeta_{2} M_{0,1}(z).$$

The remainder is O(epsilon) at fixed z>0. This branch and local continuation do not provide an arbitrary complex-z or z=-1 boundary theorem. See [reproduction](reproduce.md) and [validation](validation.md) for the adopted-input and executed-check boundaries.
