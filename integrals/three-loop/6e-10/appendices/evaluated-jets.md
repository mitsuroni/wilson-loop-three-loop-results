# 6e-10: complete adopted regular-kernel Taylor jets

For $z>0$, write $G_{(a_1,\ldots,a_n)}=G(a_1,\ldots,a_n;z)$, using the first-letter kernel $dt/(t-a_1)$, and $G_{\varnothing}=1$. Each table is the complete finite analytic expression
\[
X(z)=\sum_{\mathbf w\in\mathcal W_X}c_{X,\mathbf w}(z)G_{\mathbf w}(z).
\]
All listed coefficients are exact; unlisted words have zero coefficient. The empty-word row retains the endpoint periods. These are finite HLOG coordinate tables, not definitions of unknown integrals. A zero table means the complete coefficient is zero.
The exact adopted object is $M_q(z)=[e^q]M(e,z)$, the **uncentered** regular four-variable gap kernel in the main derivation. It excludes $(1+z)^2z^{r e}/e^2$, Gamma and physical crossing. The reflected table is the same function evaluated at $1/z$, already converted to words in $z$ on the positive-real branch. No extra scale, crossing or period factor is part of this table.

The actual archived variable order, inner to outer, is $s_{41}$, $s_{42}$, $s_{21}$, $s_{22}$. All dependent bounds and the exact kernels are in the main text. The accepted computational lemma is that successive integration of those stated complete kernels, with ordinary combined endpoints, gives precisely these finite tables. Its condition is fixed real $z>0$ and the local regulator/continuation specified in the main text.

The byte-preserved [request](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/request.json) and [integration script](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/compute.mpl) identify the actual integration rather than a result-only pointer. Historical HyperInt enabled divergence checks, abort-on-divergence and nonlinear-factor rejection, then performed fibration to $z$ and reflected fibration. The new Python batch rechecks and parses these integrated subresults and reassembles physical $I$; it does not launch a new CAS.

## $M_0$

Source: [q0.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/q0.raw.txt).

$$M_0=0.$$

## $M_0^{\iota}$

Source: [q0.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/q0.reflected.txt).

$$M_0^{\iota}=0.$$

## $M_1$

Source: [q1.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/q1.raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{\zeta_{2} \zeta_{3} + 2 \zeta_{5}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $- \frac{2 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0)$ | $\frac{1}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{3 \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{3 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0)$ | $\frac{1}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $\frac{\zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0)$ | $\frac{1}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $\frac{3 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0)$ | $- \frac{2}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{\zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $\frac{3 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $- \frac{6 \zeta_{2}}{\left(z + 1\right)^{2}}$ |

## $M_1^{\iota}$

Source: [q1.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/q1.reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{5 z^{2} \zeta_{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{2 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{3 z^{2} \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0)$ | $- \frac{21 z^{2} \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0)$ | $\frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $\frac{z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $- \frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0)$ | $\frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $- \frac{z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $- \frac{z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $- \frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $\frac{6 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |

## $M_2$

Source: [q2.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/q2.raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{143 \zeta_{2}^{3} - 294 \zeta_{3}^{2}}{42 \left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $- \frac{11 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0)$ | $- \frac{5 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $\frac{3 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{\zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $- \frac{7 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1)$ | $\frac{6 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0)$ | $\frac{23 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1)$ | $- \frac{21 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0)$ | $\frac{9 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $\frac{12 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0,0,0)$ | $- \frac{2}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1,0,0)$ | $\frac{5}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0,0)$ | $- \frac{9}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1,0,0)$ | $- \frac{11}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1,0,0)$ | $\frac{12}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0,0)$ | $\frac{6}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1,0,0)$ | $\frac{6}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0,0,0)$ | $- \frac{5}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1,0,0)$ | $- \frac{7}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1,0,0)$ | $- \frac{4}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1)$ | $- \frac{12 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{42 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{2 \left(7 \zeta_{2} \zeta_{3} - 9 \zeta_{5}\right)}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1)$ | $\frac{15 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1)$ | $- \frac{15 \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0)$ | $- \frac{2 \left(2 \zeta_{2} \zeta_{3} - 5 \zeta_{5}\right)}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0)$ | $- \frac{8 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1)$ | $\frac{18 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1)$ | $- \frac{33 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0)$ | $- \frac{21 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0)$ | $\frac{5 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0)$ | $\frac{17 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1)$ | $\frac{36 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0,0)$ | $\frac{6}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0,0)$ | $\frac{6}{\left(z + 1\right)^{2}}$ |

## $M_2^{\iota}$

Source: [q2.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-6228fdff75ab4cb38c8547eeb69371e8/q2.reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{z^{2} \left(3391 \zeta_{2}^{3} - 1470 \zeta_{3}^{2}\right)}{210 \left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0,0)$ | $- \frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0,0,0)$ | $\frac{4 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1,0,0)$ | $\frac{7 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1,0,0)$ | $- \frac{12 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1,0,0)$ | $- \frac{6 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0,0,0)$ | $- \frac{5 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1,0,0)$ | $\frac{5 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1,0,0)$ | $\frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1)$ | $\frac{6 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $- \frac{9 z^{2} \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{4 z^{2} \left(8 \zeta_{2} \zeta_{3} - 3 \zeta_{5}\right)}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1)$ | $- \frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1)$ | $- \frac{15 z^{2} \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0)$ | $- \frac{2 z^{2} \left(8 \zeta_{2} \zeta_{3} - 3 \zeta_{5}\right)}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0)$ | $\frac{10 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1)$ | $- \frac{18 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1)$ | $\frac{21 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0)$ | $\frac{15 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0)$ | $\frac{2 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0)$ | $- \frac{9 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $\frac{6 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1)$ | $- \frac{36 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $\frac{7 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0)$ | $- \frac{5 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $- \frac{3 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $- \frac{21 z^{2} \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $\frac{5 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1)$ | $- \frac{6 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0)$ | $\frac{97 z^{2} \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1)$ | $\frac{15 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0)$ | $\frac{2 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $- \frac{12 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
