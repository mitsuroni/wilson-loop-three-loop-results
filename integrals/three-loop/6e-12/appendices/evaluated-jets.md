# 6e-12: complete adopted regular-kernel Taylor jets

For $z>0$, write $G_{(a_1,\ldots,a_n)}=G(a_1,\ldots,a_n;z)$, using the first-letter kernel $dt/(t-a_1)$, and $G_{\varnothing}=1$. Each table is the complete finite analytic expression
\[
X(z)=\sum_{\mathbf w\in\mathcal W_X}c_{X,\mathbf w}(z)G_{\mathbf w}(z).
\]
All listed coefficients are exact; unlisted words have zero coefficient. The empty-word row retains the endpoint periods. These are finite HLOG coordinate tables, not definitions of unknown integrals. A zero table means the complete coefficient is zero.
The exact adopted object is $M_q(z)=[e^q]M(e,z)$, the **uncentered** regular four-variable gap kernel in the main derivation. It excludes $(1+z)^2z^{r e}/e^2$, Gamma and physical crossing. The reflected table is the same function evaluated at $1/z$, already converted to words in $z$ on the positive-real branch. No extra scale, crossing or period factor is part of this table.

The actual archived variable order, inner to outer, is $s_{41}$, $s_{42}$, $s_{21}$, $s_{22}$. All dependent bounds and the exact kernels are in the main text. The accepted computational lemma is that successive integration of those stated complete kernels, with ordinary combined endpoints, gives precisely these finite tables. Its condition is fixed real $z>0$ and the local regulator/continuation specified in the main text.

The byte-preserved [request](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/request.json) and [integration script](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/compute.mpl) identify the actual integration rather than a result-only pointer. Historical HyperInt enabled divergence checks, abort-on-divergence and nonlinear-factor rejection, then performed fibration to $z$ and reflected fibration. The new Python batch rechecks and parses these integrated subresults and reassembles physical $I$; it does not launch a new CAS.

## $M_0$

Source: [q0.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/q0.raw.txt).

$$M_0=0.$$

## $M_0^{\iota}$

Source: [q0.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/q0.reflected.txt).

$$M_0^{\iota}=0.$$

## $M_1$

Source: [q1.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/q1.raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{2 \zeta_{2} \zeta_{3} - 7 \zeta_{5}}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0)$ | $\frac{1}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $- \frac{\zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{\zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{3 \zeta_{2}^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{\zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $\frac{2 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0)$ | $\frac{2}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $\frac{3 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $- \frac{\zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0)$ | $- \frac{1}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $- \frac{3 \zeta_{2}}{\left(z + 1\right)^{2}}$ |

## $M_1^{\iota}$

Source: [q1.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/q1.reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{z^{2} \left(8 \zeta_{2} \zeta_{3} - 5 \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $\frac{2 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{2 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0)$ | $\frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $\frac{3 z^{2} \zeta_{2}^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $- \frac{z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $\frac{2 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0)$ | $\frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $- \frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{z^{2} \zeta_{2}^{2}}{2 \left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $- \frac{z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0,0)$ | $- \frac{2 z^{2}}{\left(z + 1\right)^{2}}$ |

## $M_2$

Source: [q2.raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/q2.raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{319 \zeta_{2}^{3} - 1680 \zeta_{3}^{2}}{210 \left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0,0)$ | $- \frac{3}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0,0,0)$ | $- \frac{10}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1,0,0)$ | $- \frac{5}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0,0)$ | $\frac{3}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1,0,0)$ | $- \frac{3}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1,0,0)$ | $\frac{3}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1,0,0)$ | $\frac{8}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0,0)$ | $\frac{4}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1,0,0)$ | $- \frac{6}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0,0,0)$ | $- \frac{4}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0)$ | $\frac{9 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1)$ | $\frac{24 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $\frac{3 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1)$ | $- \frac{9 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1)$ | $- \frac{18 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $- \frac{8 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(-1,-1)$ | $- \frac{15 \zeta_{2}^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0)$ | $\frac{10 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1)$ | $- \frac{3 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,-1)$ | $- \frac{15 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $- \frac{46 \zeta_{2} \zeta_{3} - 45 \zeta_{5}}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{24 \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $- \frac{7 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0)$ | $- \frac{10 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0)$ | $\frac{3 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $- \frac{6 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1)$ | $\frac{9 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0)$ | $- \frac{\zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{8 \zeta_{2} \zeta_{3} - 23 \zeta_{5}}{2 \left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1)$ | $\frac{9 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0)$ | $- \frac{11 \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0)$ | $\frac{29 \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $\frac{8 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $\frac{3 \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,-1,0,0)$ | $\frac{3}{\left(z + 1\right)^{2}}$ |
| $(-1,0,0,0,0,0)$ | $\frac{14}{\left(z + 1\right)^{2}}$ |

## $M_2^{\iota}$

Source: [q2.reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-moments-9738d43368f44b6084d24135a2bd8665/q2.reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $\frac{z^{2} \left(293 \zeta_{2}^{3} - 840 \zeta_{3}^{2}\right)}{105 \left(z + 1\right)^{2}}$ |
| $(-1,0,0,0)$ | $\frac{z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1)$ | $\frac{3 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1)$ | $- \frac{z^{2} \left(64 \zeta_{2} \zeta_{3} - 39 \zeta_{5}\right)}{2 \left(z + 1\right)^{2}}$ |
| $(0,-1)$ | $\frac{51 z^{2} \zeta_{2}^{2}}{5 \left(z + 1\right)^{2}}$ |
| $(-1,0,0)$ | $- \frac{10 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0)$ | $- \frac{10 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0)$ | $\frac{7 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1)$ | $\frac{6 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1)$ | $- \frac{9 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0)$ | $- \frac{10 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0)$ | $\frac{13 z^{2} \left(\zeta_{2} \zeta_{3} - \zeta_{5}\right)}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0)$ | $\frac{6 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0)$ | $- \frac{49 z^{2} \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(0,0,-1)$ | $- \frac{z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0)$ | $\frac{2 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0)$ | $- \frac{9 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1)$ | $- \frac{3 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0)$ | $- \frac{3 z^{2} \zeta_{3}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1)$ | $\frac{9 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1)$ | $\frac{18 z^{2} \zeta_{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0)$ | $\frac{13 z^{2} \zeta_{2}^{2}}{10 \left(z + 1\right)^{2}}$ |
| $(-1,-1)$ | $- \frac{15 z^{2} \zeta_{2}^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,0,0,0)$ | $\frac{5 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,-1,0,0,0,0)$ | $- \frac{10 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,0,0,0)$ | $- \frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,-1,0,0)$ | $\frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,-1,0,0)$ | $- \frac{3 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,-1,-1,0,0)$ | $- \frac{z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,0,0,0,0,0)$ | $\frac{4 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,0,0,0,0)$ | $\frac{6 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(0,-1,-1,-1,0,0)$ | $\frac{6 z^{2}}{\left(z + 1\right)^{2}}$ |
| $(-1,0,-1,0,0,0)$ | $- \frac{7 z^{2}}{\left(z + 1\right)^{2}}$ |
