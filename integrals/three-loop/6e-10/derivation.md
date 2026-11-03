# 6e-10: from its actual ordered root to every Laurent coefficient

This graph has two ordered same-channel opposite propagators and one adjacent propagator. The adopted initial root is its original prepared $C/J$ on the specified iterated domain. Work at $D=4-2e$, $e=\epsilon$, fixed real $z>0$, and real logarithms $\ell=\log z$. The physical source inherits $(-2\pi v_{14}+i0)^{3e}$; stripping that scale for the positive-real parameter calculation does not erase its source prescription. Upstream graph rules, Fourier propagation and the prepared color projection are accepted root facts, not a newly executed FeynCalc calculation.

The original three momentum factors have these actual source separations:

| Momentum | Endpoint Lorentz labels | Source separation/phase argument |
| --- | --- | --- |
| l1 | mu1, mu4 | `-((1-s41)*v4 + v1 + s21*v2)` |
| l2 | mu2, mu5 | `-((1-s42)*v4 + v1 + s22*v2)` |
| l3 | mu3, mu6 | `(1-s31)*v3 + s43*v4` |


Write the root, with its full numerator retained, as
\[
J(e,z)=\int_{\Omega}^{\mathrm{ordered}}d^6s\,N(z)\prod_{k=1}^3f_k(\mathbf s,z)^{e-1}.
\]
The actual polynomials are
\[
\begin{aligned}
N(z)&=z \left(z + 1\right)^{2},\\
f_1&=s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right),\\
f_2&=s_{22} s_{42} z + \left(1 - s_{22}\right) \left(1 - s_{42}\right),\\
f_3&=s_{43} z \left(1 - s_{31}\right).
\end{aligned}
\]

The original outer-to-inner order is $s_{22}\to s_{31}\to s_{43}\to s_{21}\to s_{42}\to s_{41}$. Its dependent bounds are part of the root:

| Variable | Lower | Upper |
|---|---|---|
| $s_{22}$ | $0$ | $1$ |
| $s_{31}$ | $0$ | $1$ |
| $s_{43}$ | $0$ | $1$ |
| $s_{21}$ | $0$ | $s_{22}$ |
| $s_{42}$ | $0$ | $s_{43}$ |
| $s_{41}$ | $0$ | $s_{42}$ |

The [input](input.json) binds this actual root and its full-color upstream source. No old cube Jacobian is inserted into this ordered integral. In the fixed three-loop basis its selected coordinate is $a_1=1/2$ and its cyclic multiplicity is $m=4$. Thus $\kappa=-1/8$, with the exact physical normalization
\[
C=\frac{a_1}{(-2)^3}\mathcal N S_{14}Z^3,\quad Z=e^{-\gamma_Ee}\Gamma(1-e),\quad
I=\kappa Z^3\frac{J(e,z)+z^{3e}J(e,1/z)}{1+z^{3e}},\qquad \kappa=-ma_1/16.\tag{1}
\]
Here $\mathcal N=C_A^2C_F(\alpha_s/\pi)^3$ and $S_{14}=(-2\pi v_{14}+i0)^{3e}$. There is exactly one Euler compensation in $Z^3$ and no additional color/orbit factor after (1). The [common derivation](../../families/three-loop-6e/common-normalization.md) explains the matching to the orbit sum.

The graph-specific adjacent integrations are now performed, first on $0<\Re e<1/8$ where the absolute-convergence/Fubini bound in [the shared method](../../families/three-loop-6e/same-channel-two.md) applies. Each primitive keeps its actual endpoints:

| Eliminated variable | Density | Interval | Integral |
|---|---|---|---|
| $s_{31}$ | $\left(1 - s_{31}\right)^{e - 1}$ | $(0,1)$ | $\frac{1}{e}$ |
| $s_{43}$ | $s_{43}^{e - 1}$ | $(s_{42},1)$ | $\frac{1 - s_{42}^{e}}{e}$ |

The adjacent source scale is $q=z$. Combining its original numerator and propagator power gives exactly $q^e=z^{e}$. The remaining endpoints retain their original relative order; every regrouping here has the just-stated Fubini bridge. Define the actual **uncentered** regular kernel
\[
\mathcal M(e,z)=\int_{\Omega_{\rm rem}}\frac{W(e)\exp(e\Delta)}{F_1F_2}\,d^4s,
\quad \Delta=\log F_1+\log F_2,\qquad M_j=[e^j]\mathcal M.
\]
Here one-index $M_j$ is this regular-kernel coefficient and is distinct from the two-index established moment $M_{k,j}$. Explicitly,
\[
F_1=s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right),\quad F_2=s_{22} s_{42} z + \left(1 - s_{22}\right) \left(1 - s_{42}\right),\quad W(e)=1 - s_{42}^{e}.
\]

The remaining bounds are $0<s_{21}<s_{22}$, $0<s_{22}<1$, $0<s_{41}<s_{42}$, $0<s_{42}<1$. The actual archived inner-to-outer integration order is $s_{41}$, $s_{42}$, $s_{21}$, $s_{22}$.

The complete weight vanishes at zero. Its required nonzero Taylor data are
\[
W(e)=eW_1+e^2W_2+O(e^3),
\]
\[
W_1=- \log{\left(s_{42} \right)},\qquad W_2=- \frac{\log{\left(s_{42} \right)}^{2}}{2}.
\]
Therefore the actual subintegrals that must be solved are
\[
M_0=0,\qquad M_1=\int_{\Omega_{\rm rem}}\frac{W_1}{F_1F_2}d^4s,\qquad
M_2=\int_{\Omega_{\rm rem}}\frac{W_2+W_1\Delta}{F_1F_2}d^4s.\tag{2}
\]
Both full raw and reflected **evaluated** results are printed in [the finite HLOG appendix](appendices/evaluated-jets.md), with all word coefficients and zeta-period constants. Its computational lemma adopts precisely the kernels (2), the bounds and actual order above, on the fixed positive-real branch. It does not merely assert that an integration script succeeded. Historical HyperInt used ordinary combined endpoint limits, divergence/abort checks, successive rational partial fractions and hyperlogarithmic primitives, then fibration to $z$ and reflected fibration. The shared method supplies the first explicit endpoint primitive and the interpretation of removable interior divided differences.
The original integral is thus exactly
\[
J(e,z)=\frac{(1+z)^2z^{e}}{e^2}\mathcal M(e,z).\tag{3}
\]
In the centered family presentation, $(t_1,t_2)=(s_{21},s_{22})$ and $(v,w)=(s_{41},s_{42})$. The pairing is id, with gap typeA. Its solved centered coefficients are $U_1=(1+z)^2M_1$ and $U_2=(1+z)^2(M_2-\ell M_1)$ because $V=(1+z)^2z^{-e}\mathcal M$. This scale conversion prevents double counting the center phase.

For the actual finite assembly put $A_j=(1+z)^2M_j(z)$ and $\bar A_j=(1+1/z)^2M_j(1/z)$. Substituting (3) into (1), the exact centered phase is
\[
I=\frac{\kappa Z^3}{2e^2\cosh(3e\ell/2)}
\{e^{(r-3/2)e\ell}\sum_jA_je^j+e^{(3/2-r)e\ell}\sum_j\bar A_je^j\},
\]
with $r=1$ for this source. Since $A_0=\bar A_0=0$ and $Z^3/\cosh(3e\ell/2)=1+O(e^2)$, every necessary pole-times-positive-order term is
\[
c_{-6}=\cdots=c_{-2}=0,\quad
c_{-1}=\frac\kappa2(A_1+\bar A_1),\quad
c_0=\frac\kappa2[A_2+\bar A_2+(r-3/2)\ell(A_1-\bar A_1)].\tag{4}
\]
The $M_2$ jets are essential because the original explicit pole is $e^{-2}$; the odd-crossing term in (4) is essential for the finite answer. The complete finite HLOG tables make (4) a finite, reproducible analytic expression before optional function-family compression. No final-reference coefficient is needed to define it.

## Complete normalized Laurent expansion

Write $I(e,z)=\sum_{q=-6}^{0}c_q(z)e^q+O(e)$ at fixed positive real $z$. All coefficients, including every established zero, are:

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=0.$$

$$c_{-2}(z)=0.$$

$$c_{-1}(z)=-\frac{1}{8} M_{0,1}(z) \zeta_{3}-\frac{1}{96} M_{1,3}(z).$$

$$c_{0}(z)=\begin{aligned}&\frac{1}{24} M_{0,1}(z) M_{2,1}(z)-\frac{5}{192} M_{1,1}(z) M_{1,1}(z)+\frac{1}{16} \zeta_{2} M_{2,1}(z)\\
&-\frac{7}{16} \zeta_{3} M_{1,1}(z)+\frac{263}{34560} M_{0,1}(z) M_{0,1}(z) M_{0,1}(z)+\frac{47}{5760} \zeta_{2} M_{0,1}(z) M_{0,1}(z)\\
&+\frac{59}{240} \zeta_{2} \zeta_{2} M_{0,1}(z)-\frac{9}{1024} M_{4,1}(z)-\frac{13}{1536} M_{2,3}(z)\\
&-\frac{1}{64} N_{2,2}(z)+\frac{1}{24} \tau_{D,4}(z)-\frac{1}{12} E(z).\end{aligned}$$

The two-index $M,N$ and the already solved $\tau_{D,4},E$ use the [common exact function definitions](../../families/three-loop-6e/common-normalization.md). No factorial-normalized moment or fresh letter convention is introduced. The finite HLOG coordinates in the mathematical appendix, followed by the explicit physical convolution above, are a complete alternative representation of these coefficients. The compact table is checked by exact full-word/period substitution, not by a fit at sample points. Physical complex continuation and a direct $z=-1$ boundary calculation are separate tasks; `formal=not_run`.
