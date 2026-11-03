# One Pi1 insertion and one bare exchange: eleven integrals

The original source has four ordered Wilson-line insertions, one exchange carrying a one-loop vacuum polarization, and one bare exchange. The internal bubble adds the third loop. The eleven labelled sources share six geometric reductions, conventionally called P1 through P6. Those six labels refer to routes within this branch; they are not a count of independent master integrals or the name of the overall publication batch.

The active evaluator is [pi1.py](../scripts/pi1.py). It reads the current source C/J and reconstructs the finite-order result using the existing evaluated M moments and one ordered double-kernel moment. It does not read a final integral coefficient table. The preserved original computation scripts and moment files have their own byte bindings in the packet. This is a fresh reduction/assembly with explicitly reused integrations, rather than eleven new independent internal-bubble computations.

## Retaining the transverse insertion

Use the current convention `wilson-qcd-action-path-fourier/1.0`, with D=4-2epsilon and z>0. The one-loop transverse operator is part of the initial momentum object D0 Pi1 D0, including its phase and internal closed-loop signs. Its regular colour factor is

\[
B_1(\epsilon)=\frac{(10-6\epsilon)C_A-8(1-\epsilon)T_Fn_f}{3-2\epsilon}.
\]

The quark, gluon and ghost contributions underlying this adopted operator are not all independently reintegrated in the present curation. The source phase ledger fixes the sign before parameter integration. It is not chosen by comparing a result with an old table.

For each actual source, let F be the insertion-line distance polynomial, d its tangent contraction, and F0,d0 those of the bare line. The source supplies the full numerator N and the original dependent bounds. The parameter kernel is

\[
J_i=\int_{\mathcal D_i}d^4s\;
\frac{N_i}{2(1+\epsilon)}F_i^{-2+2\epsilon}
d_{0,i}F_{0,i}^{-1+\epsilon}.
\tag{1}
\]

The full transverse numerator obeys the exact identity

\[
\frac{N}{2(1+\epsilon)}F^{-2+2\epsilon}
=dF^{-1+2\epsilon}
+\frac{\partial_s\partial_t F^{2\epsilon}}{4\epsilon(1+\epsilon)}.
\tag{2}
\]

Here s,t are the actual two insertion endpoints. The evaluator verifies (1), (2), the stored metric-plus-derivative decomposition and every source order before assembly. The mixed derivative must be integrated over the prescribed ordered domain: many of its endpoint faces have nonzero values. It cannot simply be dropped.

On the open domain, all relevant F,F0 are positive. For small positive Re epsilon, an adjacent insertion has the integrable endpoint powers 2epsilon-1, and an opposite insertion near a zero corner has an integrable radial power after the two-dimensional measure is included. The bare line is likewise integrable. Bounding each ordered domain by the full endpoint cube gives absolute integrability for the regroupings used below. The original outer-to-inner iterated order remains part of each root. Fixed-positive-z results are then continued meromorphically to epsilon=0; no new physical complex-z/i0 continuation is claimed.

## Native C, connected colour and the scalar I

The raw Wilson trace is C_F(C_F-C_A/2). The accepted connected-web projection removes the disconnected C_F^2 component and retains a1=-1/2 in the C_A C_F coordinate. This is a connected colour reconstruction, not a large-N limit. The full raw trace is retained in the source data.

The current native C has the positive Fourier rational factor 1/[32epsilon(1-2epsilon)]. Its Euler compensation is already present once. Put

\[
\mathcal G=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\qquad
R_3=\frac{\mathcal G^3}{32\epsilon(1-2\epsilon)},\qquad
\rho=\log z.
\]

The scalar I, before the regular B1 colour factor is restored, is

\[
I_i=\frac{m_i a_1}{2}R_3
\frac{J_i(z)+z^{3\epsilon}J_i(1/z)}{1+z^{3\epsilon}},\qquad a_1=-\frac12.
\tag{3}
\]

The source determines m=4 for P1-P4 and m=2 for P5-P6. Formula (3) is evaluated directly with the current C. The active program has no legacy-version or family-sign switch. The native unprojected input and the accepted ordered projected root are checked separately; their J objects agree, while their colour components have distinct declared roles.

## The common evaluated moments

Define F_z(s,t)=(1-s)(1-t)+zst and the weighted Euler integral

\[
B_{\nu,\lambda}(z)=\int_0^1\!ds\,dt\;t^\lambda F_z(s,t)^{\nu-1},
\qquad
\mathcal H(\nu,\lambda;z)=-(1+z)z^{-\nu/2}B_{\nu,\lambda}(z).
\tag{4}
\]

The historical Euler integration also identifies B as a Meijer G with upper parameters (0,-lambda,nu), lower parameters (0,nu,nu) and denominator Gamma(1-nu)Gamma(1+nu)Gamma(1+nu+lambda). Its finite-order evaluation is most conveniently connected to the original no-factorial M moments.

For beta in (0,1), put

\[
r=\frac{z}{1+z},\quad W=\log[\beta(1-\beta)],\quad
\Lambda=\rho+\log\frac{1-\beta}{\beta},\qquad
M_{k,j}=\int_0^1\frac{W^k\Lambda^j}{\beta-r}\,d\beta.
\]

The numerator vanishes at beta=r for j>=1. Thus this definition uses the removable value there, rather than an additional principal-value prescription. With lambda=c nu, write A_c=(c+1)W/2+c(Lambda-rho)/2. The integrated endpoint difference gives

\[
\mathcal H(\nu,c\nu;z)
=\sum_{q\ge0}\frac{\nu^q}{(q+1)!}
\int_0^1\frac{(A_c+\Lambda/2)^{q+1}-(A_c-\Lambda/2)^{q+1}}
{\beta-r}\,d\beta.
\tag{5}
\]

Every coefficient in (5) is a finite polynomial in the original M functions and rho. For example H0(nu)=H(nu,0) starts as M01+nu M11/2+nu^2(M21/8+M03/24). Crossing sends rho to -rho and M_kj to (-1)^(j+1)M_kj. The code computes these finite polynomials, then substitutes the existing complete HLOG evaluations of the M moments. Products use exact rational shuffle arithmetic. This reuse is a named mathematical dependency, not a new integration of those moments.

## The six source reductions

| Route | Sources | Ordered-source reduction |
| --- | --- | --- |
| P1 | 01,02 | Two same-edge endpoint chains |
| P2 | 05,06 | An intervening ordered endpoint, giving a Beta difference |
| P3 | 03,10 | Opposite insertion with weighted H and all derivative faces |
| P4 | 04,09 | Adjacent insertion plus a bare opposite line |
| P5 | 07,08 | Two coupled ordered opposite kernels |
| P6 | 11 | Opposite kernels on separate channels, with independent endpoints |

All IDs in this table have the prefix 4e-Pi1-. Equalities below are complete physical-I statements under (3). Each entry nevertheless keeps its original J, labelled source and own normalization. Reflection acts on the complete insertion-plus-bare kernel and its bounds.

For P1, the two ordered endpoint factors multiply to 1/(18epsilon^4), giving

\[
J_{01}=J_{02}=\frac{z^{3\epsilon}}{18\epsilon^4(1+\epsilon)},\qquad
I_{01}=I_{02}=-\frac{\mathcal G^3}{576\epsilon^5(1-2\epsilon)(1+\epsilon)}.
\tag{6}
\]

For P2 define Q=Gamma(1+epsilon)Gamma(1+2epsilon)/Gamma(1+3epsilon). The two raw kernels are z^epsilon(1-Q)/[4epsilon^4(1+epsilon)] and z^(2epsilon)(1-Q)/[4epsilon^4(1+epsilon)]. Their crossing quotient supplies the ratio cosh(epsilon rho/2)/cosh(3epsilon rho/2), and

\[
I_{05}=I_{06}=-\frac{\mathcal G^3(1-Q)\cosh(\epsilon\rho/2)}
{128\epsilon^5(1-2\epsilon)(1+\epsilon)\cosh(3\epsilon\rho/2)}.
\tag{7}
\]

The cancellation 1-Q=O(epsilon^2) comes from the actual endpoint difference.

For the remaining routes use C=cosh(3epsilon rho/2). Write H0^(a)=H(aepsilon,0;z), and set Hhalf=H(2epsilon,epsilon;z) for P3, with barred values denoting the actual reciprocal-z branch. The P3 crossing quotient before multiplication by (m a1/2)R3 is

\[
\frac{\cosh(\epsilon\rho/2)[H_0^{(2)}-(H_{\rm half}+\bar H_{\rm half})/2]
-\sinh(\epsilon\rho/2)(H_{\rm half}-\bar H_{\rm half})/2}{\epsilon^2 C}
+\frac{1+3(1-Q)\cosh(\epsilon\rho/2)/C}{12\epsilon^3(1+\epsilon)}.
\tag{8}
\]

The second term is the complete derivative-boundary contribution. J10(z)=z^(3epsilon)J03(1/z), hence I10=I03; the raw J functions are not identified.

For P4 put H2=H(epsilon,2epsilon;z). Its crossing quotient is

\[
\frac{\cosh(\epsilon\rho)[H_0^{(1)}-(H_2+\bar H_2)/2]
-\sinh(\epsilon\rho)(H_2-\bar H_2)/2}
{4\epsilon^2(1+\epsilon)C}.
\tag{9}
\]

The reflection J09(z)=z^(3epsilon)J04(1/z) yields I09=I04. Equations (8),(9) are evaluated with (5), retaining both even and odd crossing terms before coefficient extraction.

For P5 the ordering couples the two opposite kernels. The original scalar body is

\[
T(a,b;z)=(1+z)^2z^{-(a+b)/2}
\int_{0<s_1<s_2<1\atop0<t_1<t_2<1}
F_z(s_1,t_1)^{a-1}F_z(s_2,t_2)^{b-1}\,d^2s\,d^2t.
\tag{10}
\]

The two mixed bodies are T(2epsilon,epsilon) and T(epsilon,2epsilon). The eight nonzero rectangle-corner contributions reduce to

\[
A=H_0^{(3)}+\cosh(\epsilon\rho)(H_0^{(1)}-H_2-\bar H_2)
-\sinh(\epsilon\rho)(H_2-\bar H_2).
\]

Thus the crossing quotient for P5 is [T(2epsilon,epsilon)+T(epsilon,2epsilon)+A/(2epsilon(1+epsilon))]/(2C). Its required mixed-body jet is

\[
T(2\epsilon,\epsilon)+T(\epsilon,2\epsilon)
=2\rho_0+3\epsilon\rho_1+O(\epsilon^2),\quad
\rho_q=[\epsilon^q]T(\epsilon,\epsilon),\quad
\rho_0=\frac{M_{0,1}^2}{6}.
\tag{11}
\]

The existing rho1 evaluation is a complete weight-five HLOG in the supplied structural dictionary. It is not left as an unevaluated parameter integral. This specific original double-kernel integration is reused; the current program computes its coefficient in the new physical assembly. I07=I08 follows from the complete ordered reflection. The full scalar body is not replaced by a product of two unrelated H functions.

The indexed rho_q in (11) is the historical ordered-body coefficient, stored under the literal dictionary IDs `rho_id_0_0_0` and `rho_id_0_0_1`; its displayed family name is tau^[2,id]_(0,0;q). It is distinct from the unindexed shape logarithm rho=log(z). These source IDs do not define another physical scale or a new unevaluated function.

For P6 the original endpoints do separate, but the insertion still has a boundary term. Its crossing quotient is

\[
\frac{H_0^{(1)}}{C}
\left[H_0^{(2)}+\frac{\cosh(\epsilon\rho)}{2\epsilon(1+\epsilon)}\right].
\tag{12}
\]

Multiplying (12) by (m a1/2)R3 with m=2 gives the complete source-connected I11. Dropping the second term would miss its double pole and finite contribution.

## Finite-order budget and colour reconstruction

The original Gamma logarithm is log(Gcal)=sum_(k>=2)zeta(k)epsilon^k/k. The active code uses the exact exponential Taylor recurrence through order5, which safely covers every route. The required subintegral orders before physical convolution are:

| Route | Actual I pole bound | Required regular inputs for c0 |
| --- | ---: | --- |
| P1 | 5 | Gamma/rational jet through5 |
| P2 | 3 | Beta difference Q through5; conservative Gamma/rational jet through5 |
| P3 | 4 | Weighted H through3, Q through4 and Gamma through4 |
| P4 | 2 | Weighted H through3 and Gamma through2 |
| P5 | 1 | H through2 and the mixed T sum through1; Gamma constant suffices |
| P6 | 2 | H0 through2 and Gamma through2 |

The extra explicit poles in (8)-(11) require these orders even when the final physical pole cancels. All seven requested positions epsilon^-6 through0, including established zeros, are exported. No positive epsilon coefficient of I is needed for multiplication by the regular B1 factor.

Restoring the external operator gives the two connected colour contributions

\[
C_A^2C_F\,\frac{10-6\epsilon}{3-2\epsilon}I_i,
\qquad
C_AC_FT_Fn_f\,\frac{-8(1-\epsilon)}{3-2\epsilon}I_i.
\tag{13}
\]

The a1 and orbit factors are already in I_i. They are not inserted again in (13). The other two connected three-loop colour channels are explicitly zero for this branch, while the unprojected raw trace remains available separately. For a pole bound p, the rational operator in (13) is expanded through epsilon^p and convolved with the full geometric Laurent series. Keeping only B1(0) would omit finite contributions.

The result data provides each geometric coefficient and each coloured coefficient both as a complete 0,-1 HLOG/zeta sum and as its finite moment expression. The logarithm rho in those expressions is log(z), distinct from the physical average scale L. The original M/N conversion N_kj=M_k,j+1-rho M_kj remains available for later compact displays.

The actual batch checks and post-generation old-result comparison are recorded separately by the enclosing packet. Shared source checks, function reuse and a successful CAS/Python operation are not integral-level formal proofs. Formal state remains not_run throughout.
