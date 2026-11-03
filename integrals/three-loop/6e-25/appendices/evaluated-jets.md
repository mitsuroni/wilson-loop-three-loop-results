# 6e-25: complete adopted regular-kernel Taylor jets

For $z>0$, write $G_{(a_1,\ldots,a_n)}=G(a_1,\ldots,a_n;z)$, using the first-letter kernel $dt/(t-a_1)$, and $G_{\varnothing}=1$. Each table is the complete finite analytic expression
\[
X(z)=\sum_{\mathbf w\in\mathcal W_X}c_{X,\mathbf w}(z)G_{\mathbf w}(z).
\]
All listed coefficients are exact; unlisted words have zero coefficient. The empty-word row retains the endpoint periods. These are finite HLOG coordinate tables, not definitions of unknown integrals. A zero table means the complete coefficient is zero.
The exact adopted object is $M_q(z)=[e^q]M(e,z)$, the **uncentered** regular four-variable gap kernel in the main derivation. It excludes $(1+z)^2z^{r e}/e^2$, Gamma and physical crossing. The reflected table is the same function evaluated at $1/z$, already converted to words in $z$ on the positive-real branch. No extra scale, crossing or period factor is part of this table.

The actual archived variable order, inner to outer, is $s_{21}$, $s_{41}$, $s_{42}$, $s_{23}$. All dependent bounds and the exact kernels are in the main text. The accepted computational lemma is that successive integration of those stated complete kernels, with ordinary combined endpoints, gives precisely these finite tables. Its condition is fixed real $z>0$ and the local regulator/continuation specified in the main text.

The byte-preserved [request](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/request.json) and [integration script](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/compute.mpl) identify the actual integration rather than a result-only pointer. Historical HyperInt enabled divergence checks, abort-on-divergence and nonlinear-factor rejection, then performed fibration to $z$ and reflected fibration. The new Python batch rechecks and parses these integrated subresults and reassembles physical $I$; it does not launch a new CAS.

## $M_0$

Source: [q0.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/q0.raw.txt).

$$M_0=0.$$

## $M_0^{\iota}$

Source: [q0.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/q0.reflected.txt).

$$M_0^{\iota}=0.$$

## $M_1$

Source: [q1.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/q1.raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{5 \left(2 \zeta_{2} \zeta_{3} + \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0)$ | $- \frac{3}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0,0)$ | $- \frac{1}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $\frac{5 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0)$ | $\frac{2}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $- \frac{\zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $- \frac{3 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{4 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $- \frac{9 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0)$ | $\frac{2}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $\frac{6 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{\zeta_{3}}{\left(z + 1\right)^{2}}$ |

## $M_1^{\iota}$

Source: [q1.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/q1.reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{z^{2} \left(8 \zeta_{2} \zeta_{3} + \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $- \frac{4 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $- \frac{z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{3 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $- \frac{5 z^{2} \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $\frac{9 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $- \frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $- \frac{2 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |

## $M_2$

Source: [q2.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/q2.raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{403 \zeta_{2}^{3} - 60 \zeta_{3}^{2}}{30 \left(z + 1\right)^{2}}$ |
| $(0,-1,0,0)$ | $\frac{18 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1)$ | $\frac{48 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $- \frac{12 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1)$ | $\frac{9 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1)$ | $- \frac{18 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $- \frac{15 \left(2 \zeta_{2} \zeta_{3} + \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $\frac{5 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{97 \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(-1,0)$ | $- \frac{12 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $- \frac{3 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1)$ | $\frac{54 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0)$ | $\frac{47 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{16 \zeta_{2} \zeta_{3} + 11 \zeta_{5}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1)$ | $- \frac{75 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1)$ | $- \frac{39 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0)$ | $\frac{3 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $\frac{18 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $- \frac{25 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0)$ | $- \frac{38 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0)$ | $- \frac{15 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0)$ | $- \frac{8 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{29 \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1)$ | $\frac{27 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0,0)$ | $- \frac{20}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1,0,0)$ | $- \frac{6}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1,0,0)$ | $\frac{16}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1,0,0)$ | $\frac{18}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1,0,0)$ | $- \frac{25}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0,0)$ | $\frac{22}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0,0)$ | $\frac{3}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0,0,0)$ | $- \frac{6}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0,0,0)$ | $- \frac{9}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0,0)$ | $\frac{11}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1,0,0)$ | $\frac{9}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1,0,0)$ | $- \frac{13}{\left(z + 1\right)^{2}}$ |

## $M_2^{\iota}$

Source: [q2.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-381f6aaacf2c411e8787be767b44053f/q2.reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{z^{2} \left(61 \zeta_{2}^{3} + 210 \zeta_{3}^{2}\right)}{105 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0,0)$ | $- \frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1,0,0)$ | $\frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1,0,0)$ | $- \frac{18 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1,0,0)$ | $\frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0,0)$ | $\frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0,0)$ | $\frac{4 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1,0,0)$ | $- \frac{9 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1,0,0)$ | $\frac{5 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $- \frac{18 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $\frac{2 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0)$ | $\frac{16 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0)$ | $\frac{12 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0)$ | $\frac{4 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{51 z^{2} \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1)$ | $- \frac{27 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0)$ | $\frac{5 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1)$ | $\frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $\frac{6 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1)$ | $- \frac{9 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1)$ | $\frac{9 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{3 z^{2} \left(8 \zeta_{2} \zeta_{3} + \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $\frac{8 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $- \frac{8 z^{2} \zeta_{2}^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0)$ | $\frac{15 z^{2} \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $\frac{6 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1)$ | $- \frac{54 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{z^{2} \left(14 \zeta_{2} \zeta_{3} - 3 \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1)$ | $\frac{6 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1)$ | $\frac{15 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0)$ | $\frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
