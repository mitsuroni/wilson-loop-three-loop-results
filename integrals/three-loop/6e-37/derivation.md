# 6e-37: from its actual ordered root to every Laurent coefficient

This graph has two ordered opposite propagators plus one independent reverse-channel propagator. The adopted initial root is its original prepared $C/J$ on the specified iterated domain. Work at $D=4-2e$, $e=\epsilon$, fixed real $z>0$, and real logarithms $\ell=\log z$. The physical source inherits $(-2\pi v_{14}+i0)^{3e}$; stripping that scale for the positive-real parameter calculation does not erase its source prescription. Upstream graph rules, Fourier propagation and the prepared color projection are accepted root facts, not a newly executed FeynCalc calculation.

The original three momentum factors have these actual source separations:

| Momentum | Endpoint Lorentz labels | Source separation/phase argument |
| --- | --- | --- |
| l1 | mu1, mu4 | `-((1-s31)*v3 + v4 + s11*v1)` |
| l2 | mu2, mu5 | `-((1-s41)*v4 + v1 + s21*v2)` |
| l3 | mu3, mu6 | `-((1-s42)*v4 + v1 + s22*v2)` |


Write the root, with its full numerator retained, as
\[
J(e,z)=\int_{\Omega}^{\mathrm{ordered}}d^6s\,N(z)\prod_{k=1}^3f_k(\mathbf s,z)^{e-1}.
\]
The actual polynomials are
\[
\begin{aligned}
N(z)&=- \left(z + 1\right)^{3},\\
f_1&=s_{11} s_{31} + z \left(1 - s_{11}\right) \left(1 - s_{31}\right),\\
f_2&=s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right),\\
f_3&=s_{22} s_{42} z + \left(1 - s_{22}\right) \left(1 - s_{42}\right).
\end{aligned}
\]

The original outer-to-inner order is $s_{11}\to s_{31}\to s_{22}\to s_{42}\to s_{21}\to s_{41}$. Its dependent bounds are part of the root:

| Variable | Lower | Upper |
|---|---|---|
| $s_{11}$ | $0$ | $1$ |
| $s_{31}$ | $0$ | $1$ |
| $s_{22}$ | $0$ | $1$ |
| $s_{42}$ | $0$ | $1$ |
| $s_{21}$ | $0$ | $s_{22}$ |
| $s_{41}$ | $0$ | $s_{42}$ |

The [input](input.json) binds this actual root and its full-color upstream source. No old cube Jacobian is inserted into this ordered integral. In the fixed three-loop basis its selected coordinate is $a_1=1/2$ and its cyclic multiplicity is $m=2$. Thus $\kappa=-1/16$, with the exact physical normalization
\[
C=\frac{a_1}{(-2)^3}\mathcal N S_{14}Z^3,\quad Z=e^{-\gamma_Ee}\Gamma(1-e),\quad
I=\kappa Z^3\frac{J(e,z)+z^{3e}J(e,1/z)}{1+z^{3e}},\qquad \kappa=-ma_1/16.\tag{1}
\]
Here $\mathcal N=C_A^2C_F(\alpha_s/\pi)^3$ and $S_{14}=(-2\pi v_{14}+i0)^{3e}$. There is exactly one Euler compensation in $Z^3$ and no additional color/orbit factor after (1). The [common derivation](../../families/three-loop-6e/common-normalization.md) explains the matching to the orbit sum.

All six endpoints occur in different opposite kernels. At either zero corner a kernel is comparable to $x+y$, and $\int_0^\delta\int_0^\delta(x+y)^{-1-\eta}dxdy$ is finite for $\eta<1$. Enlarging the actual ordered domain to a cube gives the product of three independent such bounds, uniformly on a small disk about $e=0$, including finite logarithmic derivatives by slightly stronger powers. This proves both Fubini for the stated reductions and holomorphy at zero. Therefore no negative Laurent coefficient is hidden in an endpoint finite-part convention:
\[
c_{-6}=\cdots=c_{-1}=0,\qquad
c_0=\frac\kappa2[J_0(z)+J_0(1/z)].\tag{2}
\]
Only the ordinary original ordered six-dimensional integral at $e=0$ is needed. Its full evaluated raw and reflected HLOG expressions, including endpoint periods, are printed in [the finite integral appendix](appendices/evaluated-jets.md). The exact accepted integration lemma specifies the above root, its actual ordinary bounds and the inner-to-outer order
$s_{41}\to s_{21}\to s_{42}\to s_{22}\to s_{31}\to s_{11}$. Historical HyperInt applied divergence/abort checks, integrated the rational/hyperlogarithmic kernels in that order, took the combined ordinary endpoint values, and fibred raw and reflected answers to $z$. The table is $J_0$, not an already normalized $I_0$ and not an assumed centered all-$e$ mother. The new Python batch checks and parses that already integrated original-$J_0$ primitive and repeats (2); it does not execute a new cold CAS integration.

For this actual source the independent reverse-channel square is $(a,b)=(s_{11},s_{31})$, and the ordered forward endpoints are $(t_1,t_2)=(s_{21},s_{22})$, $(s_1,s_2)=(s_{41},s_{42})$, with id pairing. Reflecting $(a,b)$ complements both endpoints with unit absolute Jacobian, converting $\widetilde D_z$ to $D_z$. The domain is a genuine product with the remaining ordered four-variable domain, so [the exact factorization](../../families/three-loop-6e/mixed-opposite.md) gives $J=z^{3e/2}H_0Q_{\rm id}$. This retains the distinct id/swap ordered two-kernel mother; it is not a product of three unrestricted squares.

The evaluated zero-order kernels are $H_0(0,z)=M_{0,1}$, $Q_{\rm id}(0,z)=M_{0,1}^2/6$, and $Q_{\rm swap}(0,z)=M_{0,1}^2/3$. Their actual first path primitives, ordinary endpoint interpretation and adoption scope are in the shared method. Equivalently, the complete original-$J_0$ table printed here simplifies exactly to $M_{0,1}^3/6$. Substitution into(2), with this own $\kappa$, gives $-M_{0,1}^3/96$. Thus the finite term is completely solved even independently of the optional $H_0Q$ interpretation.

The matching finite terms of 37 and 38 do **not** imply equality of the original $J$ or complete physical $I$. At $e=1$, the regular kernels have $Q_{\rm id}=Q_{\rm swap}=(1+z)^2/(4z)$, while complete-$I$ equality would require $Q_{\rm swap}=2Q_{\rm id}$. This counterargument evaluates regular kernels, not Gamma factors at their poles. The actual new independent control also gives the reduced physical difference $1161/266240$ at $e=2,z=1/2$. Only $c_0^{37}=c_0^{38}$ is the accepted equality.

## Complete normalized Laurent expansion

Write $I(e,z)=\sum_{q=-6}^{0}c_q(z)e^q+O(e)$ at fixed positive real $z$. All coefficients, including every established zero, are:

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=0.$$

$$c_{-2}(z)=0.$$

$$c_{-1}(z)=0.$$

$$c_{0}(z)=-\frac{1}{96} M_{0,1}(z) M_{0,1}(z) M_{0,1}(z).$$

The two-index $M,N$ and the already solved $\tau_{D,4},E$ use the [common exact function definitions](../../families/three-loop-6e/common-normalization.md). No factorial-normalized moment or fresh letter convention is introduced. The finite HLOG coordinates in the mathematical appendix, followed by the explicit physical convolution above, are a complete alternative representation of these coefficients. The compact table is checked by exact full-word/period substitution, not by a fit at sample points. Physical complex continuation and a direct $z=-1$ boundary calculation are separate tasks; `formal=not_run`.
