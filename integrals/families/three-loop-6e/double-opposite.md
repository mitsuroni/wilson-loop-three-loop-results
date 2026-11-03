# Opposite propagators in two different channels

The 32/33/34 roots each contain $D_z(u,v)=(1-u)(1-v)+zuv$, $\widetilde D_z(u,v)=z(1-u)(1-v)+uv$, and one adjacent propagator. Their own articles write the actual polynomials, bounds and two endpoint primitives. The decisive property is that after those adjacent integrations, the remaining four endpoints form two independent square factors. No ordering inequality remains between the two factors; this fact must be checked for each root rather than inferred from the family name.

For $0<\Re e<1/8$, enlarging to the cube bounds the original adjacent endpoint densities by integrable powers and each opposite kernel by a two-dimensional corner kernel. This proves absolute integrability and the Fubini bridge that permits the stated grouping. The four elementary endpoint positions relative to an opposite endpoint $s$ are

| Placement | Actual adjacent integral |
|---|---|
| AO | $s^e/e$ |
| OA | $(1-s^e)/e$ |
| BO | $[1-(1-s)^e]/e$ |
| OB | $(1-s)^e/e$ |

Both factors $1/e$ remain explicit. Their numerator differences, continued as complete expressions from this initial domain, determine which poles cancel.

Use the **evaluated** centered one-opposite kernels
\[
A(e,z)=H_0=-(1+z)z^{-e/2}\int_{[0,1]^2}D_z(u,v)^{e-1}du\,dv,
\]
\[
B(e,z)=H_1=-(1+z)z^{-e/2}\int_{[0,1]^2}D_z(u,v)^{e-1}v^e du\,dv.
\]
The exact kernel, positive-quadrant Jacobian, actual integration and complete HLOG coefficients through $e^2$ are the adopted lemma in [the one-opposite method](single-opposite.md) and its [finite jet appendix](appendices/single-opposite-jets.md). The first path primitive is $[(zv)^e-(1-v)^e]/[e((1+z)v-1)]$, with its continuous divided-difference value at the apparent interior pole. Thus adoption means a specific integrated kernel with fixed ordinary endpoints and constants, not an unknown symbol $A$ or $B$.

To handle the reverse channel, reflect both endpoints: $\widetilde D_z(1-u,1-v)=D_z(u,v)$, with unit absolute Jacobian. An insertion $(1-v_{\rm old})^e$ becomes $v^e$, and $1-(1-v_{\rm old})^e$ becomes $1-v^e$. This transports both the kernel and the weight. Alternatively $\widetilde D_z=zD_{1/z}$ gives the same centered-scale rule. From the actual endpoints one obtains

| Root | Reverse-channel factor | Forward-channel factor | Complete original $J$ |
|---|---|---|---|
|32|$B$|$A-B$|$z^{2e}B(A-B)/e^2$|
|33|$A-B$|$A-B$|$z^{2e}(A-B)^2/e^2$|
|34|$A-B$|$B$|$z^{2e}(A-B)B/e^2$|

For example in 32, integrating $s_{32}\in(s_{31},1)$ with its B density gives $(1-s_{31})^e/e$, and integrating $s_{42}\in(s_{41},1)$ with its A density gives $(1-s_{41}^e)/e$. Reflection of the first square therefore yields $B$; the second yields $A-B$. In 34 the B endpoint lies before $s_{32}$ and the A endpoint before $s_{42}$, yielding $1-(1-s_{32})^e$ and $s_{42}^e$; the two factors exchange. This proves $J_{32}=J_{34}$ for the complete original integral, and their equal $\kappa=-1/16$ also proves complete physical-$I$ equality. It is stronger than a finite-coefficient alias. For 33 $\kappa=-1/8$.

With $F=B(A-B)$, $K=(A-B)^2$, the [physical normalization](common-normalization.md) gives exactly
\[
I_{32}=I_{34}=-\frac{Z^3}{32e^2\cosh(3e\ell/2)}
[e^{e\ell/2}F(z)+e^{-e\ell/2}F(1/z)],
\]
\[
I_{33}=-\frac{Z^3}{16e^2\cosh(3e\ell/2)}
[e^{e\ell/2}K(z)+e^{-e\ell/2}K(1/z)].
\]
These functions reuse the solved $A,B$ mother kernels. For the finite coefficient only the following fully evaluated jets can enter, with $t_j=[e^j]H_0$:
\[
A=t_0+et_1+e^2t_2+O(e^3),\quad
B^+=t_0+2et_1+e^2(4t_2-\ell^2t_0/24)+O(e^3),
\qquad B^-=-e\ell t_0/6+O(e^2).
\]
$B^\pm=[B(z)\pm B(1/z)]/2$. The $e^2$ odd coefficient of $B$ is not needed: it multiplies an even $A_0$ in $F_2$, so cancels under the symmetric physical assembly; the crossing-odd term that survives comes only from $F_1^-$.

For an explicit highest-order convolution, put $b_1^+=2t_1$, $b_1^-=-\ell t_0/6$, $b_2^+=4t_2-\ell^2t_0/24$. Since $A_0=B_0=t_0$,
\[
F_0=0,\quad F_1^+=-t_0t_1,\quad F_1^-=-t_0b_1^-,
\]
\[
F_2^+=t_0(t_2-b_2^+)+b_1^+(t_1-b_1^+)-(b_1^-)^2.
\]
Thus $c_{-1}^{32}=\kappa F_1^+$ and $c_0^{32}=\kappa(F_2^++\ell F_1^-/2)$, with $\kappa=-1/16$. For 33,
$K_0=K_1=0$, $K_2^+=t_1^2+(b_1^-)^2$ and $c_0^{33}=-K_2^+/8$. These are finite convolutions of actual solved jets, including the otherwise easily missed odd-crossing contribution. Using $\ell^2=-2t_0-6\zeta_2$, they give
\[
c_{-1}^{32}=c_{-1}^{34}=t_0t_1/16,
\]
\[
c_0^{32}=c_0^{34}=3t_0t_2/16+t_1^2/8+7t_0^3/576+7\zeta_2t_0^2/192,
\]
\[
c_0^{33}=-t_1^2/8+t_0^3/144+\zeta_2t_0^2/48.
\]
All other negative-power coefficients vanish. The known conversions $t_0=M_{0,1}$, $t_1=M_{1,1}/2$, $t_2=M_{2,1}/8-M_{0,1}^2/24$ produce the own-entry compact coefficients.

The current runtime uses a different, equally explicit organization: it analytically removes the two adjacent endpoints and adopts the original uncentered four-variable gap jets $M_0,M_1,M_2$, then reassembles the physical prefactors. The complete raw/reflected tables are printed in each own entry; their exact coefficient kernels, domains and actual integration order are stated there. In this route the finite result is obtained by formula(S) in [the two-kernel method](same-channel-two.md), with adjacent $z$ power one. These direct-gap outputs are genuinely subintegrals of $J$, not already normalized $I$ tables. The new Python execution consumes them and repeats the physical convolution; it does not claim a new cold CAS or an independent fresh evaluation of the older $A/B$ route.

The [original research method](../../../benchmarks/2ne/families/double-opposite/method.md) retains the factorization calculation and its historical tests. Real positive $z$, local meromorphic continuation in $e$, and `formal=not_run` are the scope.
