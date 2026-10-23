# Conventions and the moment dictionary

This delivery treats the three matter colors of the fundamental lightlike quadrilateral Wilson-loop **exponent** in massless QCD. All three-loop coefficients have the coupling and the displayed color factor removed. Set

\[
\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi,\qquad T=T_Fn_f,\qquad
\log W=\sum_{n\geq1}\lambda^n w_n,
\]
\[
w_3=C_FC_A^2w_{3,AA}+C_FC_ATw_{3,Af}
+C_F^2Tw_{3,Ff}+C_FT^2w_{3,ff}.
\]

The three delivered colors are Af, Ff and ff. The AA direct diagram sum remains incomplete. Its additional total-reconstruction constant does not enter any delivered matter color. A complete matter reporting chain therefore does not require completion of AA.

## Scale and expansion variables

Use the initial real region \(s,t<0\), \(z=s/t>0\), \(D=4-2\epsilon\), and the common source continuation prescription:

\[
\mu_x^2=\frac{e^{2\gamma_E}}4\mu_{\overline{\mathrm{MS}}}^2,\quad
X_s=\mu_x^2(-s+i0),\quad X_t=\mu_x^2(-t+i0),
\]
\[
\rho=\log z,\quad L_s=\log X_s,\quad L_t=\log X_t,\quad
L=\frac{L_s+L_t}{2},\quad L_s-L_t=\rho,
\]
\[
S_n=X_s^{n\epsilon}+X_t^{n\epsilon}
=2e^{n\epsilon L}\cosh(n\epsilon\rho/2).
\]

The report always uses \(L\) for the physical scale logarithm and \(\rho\) for the shape logarithm. Some native integral files use `L` for \(\log z\); their unmodified expressions and local dictionary are retained in the diagram ledger. They must be translated before combining expressions. Likewise, `r` can mean either \(1/(1+z)\) or \(z/(1+z)\) in different sources. The source argument is never silently replaced by \(z\).

All Laurent coefficients are ordinary coefficients, with no factorial. The delivered target range is \(\epsilon^{-6},\ldots,\epsilon^0\), including explicit zeros, with remainder \(O(\epsilon)\). The required lower-order budget is \(\widehat w_1\) through \(\epsilon^2\) and \(\widehat w_2\) through \(\epsilon^1\). Each lower order retains its own scale factor.

## Which quantity is being reported?

For each color, define

\[
U_{3,c}=\sum_{p=-6}^0u_{c,p}(z)\epsilon^p+O(\epsilon),\quad
\widehat w_{3,c}=S_3U_{3,c},\quad
\widetilde w_{3,c}=\widehat w_{3,c}+\mathrm{CT}_{3,c}.
\]
\[
\mathcal Z_{3,c}=\operatorname{Pole}_{\epsilon}\widetilde w_{3,c},\quad
F_{3,c}(L,\rho)=[\epsilon^0]\widetilde w_{3,c},\quad
f_{3,c}(\rho)=F_{3,c}(0,\rho).
\]

Here \(\mathcal Z_{3,c}\) is the coefficient of **\(\log Z_4\)** for \(W_B=Z_4W_R\); it is not the coefficient of \(Z_4\). Subtracting it leaves \(F_{3,c}+O(\epsilon)\). The phrase “constant term” in a Laurent expansion means \([\epsilon^0]\) of a specified stage. It is generally a function of \(L,z\), and differs between the bare and renormalized stages. The number at \(L=0,z=1\) is a further specialization, included as a numerical sample, not the definition of the full finite term.

The evaluated native integrals already contain their documented orbit/web factors, crossing convention and exactly one Euler compensation. The current assembly is already an exponent assembly. No extra \(W_1W_2\) subtraction or second multiplication by prepared \(C\) is applied.

## Hlog and GPL

\[
\operatorname{Hlog}(z,[a_1,\ldots,a_m])=G(a_1,\ldots,a_m;z),\quad
G(a_1,\ldots,a_m;z)=\int_0^z\frac{dt}{t-a_1}G(a_2,\ldots,a_m;t),
\]
\[
G(\varnothing;z)=1,\quad G(0^m;z)=\log^m(z)/m!,\quad
G(-1;z)=\log(1+z).
\]

The first letter is the outermost kernel. The tangential endpoint and branch prescription are inherited from the source. In the machine formulas, `H_z_m1_0` denotes \(G(-1,0;z)\); `H_r_...` keeps the native argument \(r\). The portable conversion dictionary lists the actual expansions used, so the Hlog-to-moment relationship is an executable identity at the required weights.

## Fixed M and N moments

Let

\[
r_z=\frac z{1+z},\quad W_\beta=\log[\beta(1-\beta)],\quad
R_\beta=\log\frac{1-\beta}{\beta},\quad L_\beta=\rho+R_\beta.
\]
\[
M_{i,j}(z)=\int_0^1\frac{W_\beta^iL_\beta^j}{\beta-r_z}\,d\beta,\quad
N_{i,j}(z)=\int_0^1\frac{W_\beta^iR_\beta L_\beta^j}{\beta-r_z}\,d\beta,
\quad i\geq0,\ j\geq1.
\]

These are raw moments, without factorial denominators. The apparent pole at \(\beta=r_z\) is removable. The aliases `M01`, `M_0_1`, and \(M_{0,1}\) denote the same defined function in the unified tables. Useful exact identities are

\[
N_{i,j}=M_{i,j+1}-\rho M_{i,j},\qquad
M_{0,1}=-\frac{\rho^2+\pi^2}{2},\qquad M_{0,3}=-M_{0,1}^2.
\]

Changing \(\beta\) to \(1-\beta\) gives

\[
M_{i,j}(1/z)=(-1)^{j+1}M_{i,j}(z),\qquad
N_{i,j}(1/z)=(-1)^jN_{i,j}(z).
\]

Thus all moments retained in the channel totals are even under \(z\leftrightarrow1/z\), although intermediate moments such as \(M_{1,2}\) are odd. The explicit \(N_{1,2}\) bridge is \(M_{1,3}-\rho M_{1,2}\).

## Relation to tau and E notation

The required channel totals close in \(M,N,L,\rho\) and zeta values. No independent \(\tau^+\) or \(E\) function occurs through \(\epsilon^0\) in these three totals or in their consumed counterterms. This statement is restricted to these outputs, not every three-loop diagram or higher epsilon order.

For comparison with the existing block notation, the registered fixed-coefficient identities are

\[
\tau^+_{0,0}=M_{0,1},\quad \tau^+_{0,1}=\frac12M_{1,1},\quad
\tau^+_{0,2}=\frac18M_{2,1}-\frac1{24}M_{0,1}^2,
\]
\[
\tau^+_{D,2}=\frac98M_{2,1}-\frac18M_{0,1}^2,\quad
\tau^+_{D,3}=-4\zeta_3M_{0,1}+2\zeta_2M_{1,1}
+\frac{M_{1,3}}{48}+\frac{9M_{3,1}}{16}-\frac{N_{1,2}}2.
\]

These identities do not identify the old weight-six \(E\) with the newer Euler-lift symbols. The more extensive [notation audit](sources/workflow/checks/three-loop-delivery-audit-20261023/notation-audit.md) records those distinct definitions. The two-loop vertex's independent weight-six tau block first enters an epsilon order beyond the counterterm budget consumed here.

## Counterterms and channel dependence

Write \(\widehat w_1=C_Fv_1\) and \(\widehat w_2=C_F(C_Av_{2,A}+Tv_{2,f})\), with all scale dependence retained. Then

\[
b_0=\frac{11C_A-4T}{12},\quad
b_1=\frac{34C_A^2-20C_AT-12C_FT}{48},
\]
\[
\mathrm{CT}_3=-\frac{2b_0}{\epsilon}\widehat w_2+
\left(\frac{b_0^2}{\epsilon^2}-\frac{b_1}{2\epsilon}\right)\widehat w_1.
\]
After removing each three-loop color monomial,

\[
\mathrm{CT}_{Af}=-\frac{11v_{2,f}}{6\epsilon}
+\frac{2v_{2,A}}{3\epsilon}-\frac{11v_1}{18\epsilon^2}+\frac{5v_1}{24\epsilon},
\]
\[
\mathrm{CT}_{Ff}=\frac{v_1}{8\epsilon},\qquad
\mathrm{CT}_{ff}=\frac{2v_{2,f}}{3\epsilon}+\frac{v_1}{9\epsilon^2}.
\]

Thus each three-loop matter channel is autonomous with respect to the missing three-loop AA data, while retaining its specified lower-order dependencies. The physical beta coefficients are shared constants. Color projection commutes with this explicitly polynomial counterterm calculation.

## Native diagram weights

Define the signed insertion coefficient

\[
a(\epsilon)=\frac{10-6\epsilon}{3-2\epsilon},\quad
B_f(\epsilon)=-\frac{8(1-\epsilon)}{3-2\epsilon}.
\]

Af has 11 Pi1 diagrams weighted by \(B_f\), two already color-resolved chains, two Pi2 Af entries, 12 native TGV-leg insertions weighted by \(B_f\), and four native Gamma1 quark entries. A scalar chain would instead receive \(2aB_f\); its delivered color-resolved coefficient already includes this factor. The ff and Ff channels use two chain and two Pi2 entries respectively. There are 31 unique orbit IDs and 35 channel occurrences.

The historical leg/Gamma minus-sign bridge has already been incorporated in the current native inputs. Applying it again would be an error. The ledger records the actual JSON field and external multiplier for each contribution. Multiplication by \(B_f(\epsilon)\) is a Laurent convolution, including all terms needed through epsilon zero.
