# 6e-20: complete adopted original zero-order integral

For $z>0$, write $G_{(a_1,\ldots,a_n)}=G(a_1,\ldots,a_n;z)$, using the first-letter kernel $dt/(t-a_1)$, and $G_{\varnothing}=1$. Each table is the complete finite analytic expression
\[
X(z)=\sum_{\mathbf w\in\mathcal W_X}c_{X,\mathbf w}(z)G_{\mathbf w}(z).
\]
All listed coefficients are exact; unlisted words have zero coefficient. The empty-word row retains the endpoint periods. These are finite HLOG coordinate tables, not definitions of unknown integrals. A zero table means the complete coefficient is zero.
The exact adopted object is $J_0(z)$, the ordinary original ordered six-variable integral at $e=0$ in the main derivation. It includes the original numerator and is **not** an already normalized $I$ or a centered all-e mother. The reflected table is the same function evaluated at $1/z$, already converted to words in $z$ on the positive-real branch. No extra scale, crossing or period factor is part of this table.

The actual archived variable order, inner to outer, is $s_{41}$, $s_{21}$, $s_{42}$, $s_{22}$, $s_{43}$, $s_{23}$. All dependent bounds and the exact kernels are in the main text. The accepted computational lemma is that successive integration of those stated complete kernels, with ordinary combined endpoints, gives precisely these finite tables. Its condition is fixed real $z>0$ and the local regulator/continuation specified in the main text.

The byte-preserved [request](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-constant-14d23b53549f4b6b9547250deaceb13f/request.json) and [integration script](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-constant-14d23b53549f4b6b9547250deaceb13f/compute.mpl) identify the actual integration rather than a result-only pointer. Historical HyperInt enabled divergence checks, abort-on-divergence and nonlinear-factor rejection, then performed fibration to $z$ and reflected fibration. The new Python batch rechecks and parses these integrated subresults and reassembles physical $I$; it does not launch a new CAS.

## $J_0$

Source: [raw.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-constant-14d23b53549f4b6b9547250deaceb13f/raw.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{69 \zeta_{2}^{3} - 20 \zeta_{3}^{2}}{10}$ |
| $(0,0,0,-1,0,0)$ | $-2$ |
| $(0,0,-1,0)$ | $- 6 \zeta_{2}$ |
| $(0,0,0)$ | $- 2 \zeta_{3}$ |
| $(0,0,-1,-1)$ | $12 \zeta_{2}$ |
| $(0,0,-1)$ | $4 \zeta_{3}$ |
| $(0,0,0,-1)$ | $- 6 \zeta_{2}$ |
| $(0)$ | $- 2 \left(3 \zeta_{2} \zeta_{3} + \zeta_{5}\right)$ |
| $(0,0)$ | $- \frac{27 \zeta_{2}^{2}}{10}$ |
| $(0,0,0,0,0,0)$ | $-1$ |
| $(0,0,0,0)$ | $- \zeta_{2}$ |
| $(0,0,-1,0,0,0)$ | $-2$ |
| $(0,0,-1,-1,0,0)$ | $4$ |

## $J_0^{\iota}$

Source: [reflected.txt](../../../families/three-loop-6e/runtime/primitives/objects/cas-runs/ordered-constant-14d23b53549f4b6b9547250deaceb13f/reflected.txt).

| Word $\mathbf w$ | Exact coefficient |
|---|---|
| $\varnothing$ | $- \frac{69 \zeta_{2}^{3} - 20 \zeta_{3}^{2}}{10}$ |
| $(0,0,0,-1,0,0)$ | $-2$ |
| $(0,0,-1,0)$ | $- 6 \zeta_{2}$ |
| $(0,0,0)$ | $- 2 \zeta_{3}$ |
| $(0,0,-1,-1)$ | $12 \zeta_{2}$ |
| $(0,0,-1)$ | $4 \zeta_{3}$ |
| $(0,0,0,-1)$ | $- 6 \zeta_{2}$ |
| $(0)$ | $- 2 \left(3 \zeta_{2} \zeta_{3} + \zeta_{5}\right)$ |
| $(0,0)$ | $- \frac{27 \zeta_{2}^{2}}{10}$ |
| $(0,0,0,0,0,0)$ | $-1$ |
| $(0,0,0,0)$ | $- \zeta_{2}$ |
| $(0,0,-1,0,0,0)$ | $-2$ |
| $(0,0,-1,-1,0,0)$ | $4$ |
