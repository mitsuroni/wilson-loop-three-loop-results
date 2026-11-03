# 6e-20: from its actual ordered root to every Laurent coefficient

This graph has three ordered same-channel opposite propagators. The adopted initial root is its original prepared $C/J$ on the specified iterated domain. Work at $D=4-2e$, $e=\epsilon$, fixed real $z>0$, and real logarithms $\ell=\log z$. The physical source inherits $(-2\pi v_{14}+i0)^{3e}$; stripping that scale for the positive-real parameter calculation does not erase its source prescription. Upstream graph rules, Fourier propagation and the prepared color projection are accepted root facts, not a newly executed FeynCalc calculation.

The original three momentum factors have these actual source separations:

| Momentum | Endpoint Lorentz labels | Source separation/phase argument |
| --- | --- | --- |
| l1 | mu1, mu4 | `-((1-s41)*v4 + v1 + s21*v2)` |
| l2 | mu2, mu5 | `-((1-s42)*v4 + v1 + s22*v2)` |
| l3 | mu3, mu6 | `-((1-s43)*v4 + v1 + s23*v2)` |


Write the root, with its full numerator retained, as
\[
J(e,z)=\int_{\Omega}^{\mathrm{ordered}}d^6s\,N(z)\prod_{k=1}^3f_k(\mathbf s,z)^{e-1}.
\]
The actual polynomials are
\[
\begin{aligned}
N(z)&=- \left(z + 1\right)^{3},\\
f_1&=s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right),\\
f_2&=s_{22} s_{42} z + \left(1 - s_{22}\right) \left(1 - s_{42}\right),\\
f_3&=s_{23} s_{43} z + \left(1 - s_{23}\right) \left(1 - s_{43}\right).
\end{aligned}
\]

The original outer-to-inner order is $s_{23}\to s_{43}\to s_{22}\to s_{42}\to s_{21}\to s_{41}$. Its dependent bounds are part of the root:

| Variable | Lower | Upper |
|---|---|---|
| $s_{23}$ | $0$ | $1$ |
| $s_{43}$ | $0$ | $1$ |
| $s_{22}$ | $0$ | $s_{23}$ |
| $s_{42}$ | $0$ | $s_{43}$ |
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
$s_{41}\to s_{21}\to s_{42}\to s_{22}\to s_{43}\to s_{23}$. Historical HyperInt applied divergence/abort checks, integrated the rational/hyperlogarithmic kernels in that order, took the combined ordinary endpoint values, and fibred raw and reflected answers to $z$. The table is $J_0$, not an already normalized $I_0$ and not an assumed centered all-$e$ mother. The new Python batch checks and parses that already integrated original-$J_0$ primitive and repeats (2); it does not execute a new cold CAS integration.

For the actual geometry, $(t_1,t_2,t_3)=(s_{21},s_{22},s_{23})$ and $(s_1,s_2,s_3)=(s_{41},s_{42},s_{43})$, with pairing $\pi=123$. Define $R_\pi=z^{-3e/2}J_\pi$. Reflection/reversal $t_i\mapsto1-t_{4-i}$, $s_i\mapsto1-s_{4-i}$ preserves the ordered triangles, has absolute Jacobian one, and sends the pairing to123. Each $D$ contributes the correct inverse-$z$ factor, hence $R_{123}(e,1/z)=R_{123}(e,z)$. At zero $R_0=J_0$, so the centered family and direct primitive coincide exactly at the needed order.

The [shared reduction](../../families/three-loop-6e/same-channel-three.md) explicitly integrates $t_1$ and $t_3$, with primitives $L_-/(v-a)$ and $L_+/(v-a)$ on their actual endpoint intervals, to give the complete four-/three-variable logarithmic kernels. Those primitives and removable divided differences explain how the original kernel enters HLOG integration. Its complete evaluated zero-order sources are printed in each own appendix; their exact shuffle-coordinate conversion gives
\[
T_3=2E-\tfrac12M_{0,1}M_{2,1}+\tfrac12M_{1,1}^2+\tfrac{\pi^2}{9}M_{0,3}-\tfrac1{12}M_{0,5},
\quad R_{132}^+(0,z)=M_{0,1}^3/36-T_3.
\]
Here $T_3=J_{20}(0,z)$ is already evaluated, and its physical finite coefficient is $-T_3/16$. The centered123 source is self-crossing.

The $E$ in this expression is the previously solved double-kernel endpoint insertion, with its exact structural definition in the shared method and its established HLOG dictionary in the [common functions](../../families/three-loop-6e/common-normalization.md). It is optional compact notation for a fully printed HLOG answer, not a renamed unsolved part of this graph. Neither the finite relation nor the compact basis is promoted to an all-$e$ relation between the three-kernel mothers.

## Complete normalized Laurent expansion

Write $I(e,z)=\sum_{q=-6}^{0}c_q(z)e^q+O(e)$ at fixed positive real $z$. All coefficients, including every established zero, are:

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=0.$$

$$c_{-2}(z)=0.$$

$$c_{-1}(z)=0.$$

$$c_{0}(z)=\frac{1}{32} M_{0,1}(z) M_{2,1}(z)-\frac{1}{32} M_{1,1}(z) M_{1,1}(z)+\frac{1}{144} M_{0,1}(z) M_{0,1}(z) M_{0,1}(z)-\frac{1}{8} E(z).$$

The two-index $M,N$ and the already solved $\tau_{D,4},E$ use the [common exact function definitions](../../families/three-loop-6e/common-normalization.md). No factorial-normalized moment or fresh letter convention is introduced. The finite HLOG coordinates in the mathematical appendix, followed by the explicit physical convolution above, are a complete alternative representation of these coefficients. The compact table is checked by exact full-word/period substitution, not by a fit at sample points. Physical complex continuation and a direct $z=-1$ boundary calculation are separate tasks; `formal=not_run`.
