# Validation scope for the three matter channels

The delivery covers the three selected coefficients of the logarithm of the
closed null-quadrilateral Wilson loop: $`C_FC_AT_Fn_f`$,
$`C_F^2T_Fn_f`$, and $`C_F(T_Fn_f)^2`$. The last color contains
$`T_F^2`$. It does not supply a complete fourth, pure-gluon
$`C_FC_A^2`$ finite result or a formal proof.

The new UV audit passes **115 checks, including 6 negative controls**, against
17 byte-bound source files. The exact assertions and residuals are in
`uv-checks.json`; the files, selected JSON fields, and their hashes are in
`uv-source-hashes.json`. The run timestamp is recorded separately from the
dates of the reused scientific sources.

## What is freshly checked

`verify_uv.py` reads the current three-loop selector and checks the selected
file hashes. It reads the explicit `kref` and `href` dictionaries in the
October 19 bridge builder as restricted exact arithmetic syntax. It does
**not** execute that builder or derive the external constants from the
selected three-loop poles. It then checks:

- All three selected $`k_3,h_3`$ values against those external references.
- All six negative Laurent positions, $`\epsilon^{-6}`$ through
  $`\epsilon^{-1}`$, in each channel against an independently predicted
  $`\log Z_4`$ pole template.
- The selected poles, finite functions, and $`L=0`$ shapes against the
  October 19 QCD total, with the source's complete moment-to-GPL maps.
- For the mixed $`C_FC_AT_Fn_f`$ channel, the selected `finite_MN` and
  `shape_MN` against the selected `finite_GPL` and `shape_GPL`, including
  the shuffle identities needed by $`N_{1,2}`$.
- The finite RG equation using the external $`k_3,h_3`$, the accepted
  one- and two-loop finite results, and the full lower-order color algebra.
- The raw-factor coordinate-scale algebra; the moment budget of the
  compact finite functions; and the absence of $`C3_{N4}`$ in all three
  matter finite functions.

The negative controls reject replacing $`h_3=\Gamma_{\Box,3}/2`$ by
$`h_3=\Gamma_{\Box,3}`$, and reject adding a finite constant when
comparing the finite expressions. Separate assertions demonstrate that
the finite RG equation is blind to an added $`L`$-independent even
function such as $`1+\rho^2`$.

The exact arithmetic uses SymPy and registered moment/GPL expressions.
These checks are computational identities, not formal proofs or new
three-loop integrations. The earlier symbolic parser and shuffle helper
are reused, with their source bytes included in the manifest.

## Exact conventions and lower-order dependencies

Set $`T=T_Fn_f`$, $`\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi`$,
$`L=(L_s+L_t)/2`$, and $`\rho=L_s-L_t=\log z`$, on the accepted
$`z>0`$ continuation with the original Feynman-branch convention.
The coordinate scale is

```math

\mu_x^2=\frac{e^{2\gamma_E}}4\mu_{\overline{\mathrm{MS}}}^2,
\qquad
b_0=\frac{11C_A-4T}{12},\qquad
b_1=\frac{34C_A^2-20C_AT-12C_FT}{48}.

```

The source $`I_i`$ already includes its stated web/orbit normalization
and one Euler compensation per loop. A color-stripped, unit-scale graph
sum $`U_{3,c}`$, the scale-restored bare coefficient
$`\widehat w_{3,c}`$, and its coupling-renormalized coefficient
$`\widetilde w_{3,c}`$ are different objects:

```math

S_n=e^{n\epsilon(L+\rho/2)}+e^{n\epsilon(L-\rho/2)},\qquad
\widehat w_{3,c}=S_3 U_{3,c},

```


```math

\widetilde w_3=\widehat w_3-
\frac{2b_0}{\epsilon}\widehat w_2+
\left(\frac{b_0^2}{\epsilon^2}-\frac{b_1}{2\epsilon}\right)
\widehat w_1.

```

The lower $`\widehat w_1,\widehat w_2`$ retain their own $`S_1,S_2`$.
Their counterterm inputs are bare coefficients, not the renormalized
finite functions $`F_1,F_2`$. After extracting a physical color monomial,
the negative powers of $`\widetilde w_{3,c}`$ are the corresponding
coefficient of **$`\log Z_4`$**. Its $`\epsilon^0`$ coefficient is
$`F_{3,c}(L,\rho)`$; $`f_{3,c}(\rho)=F_{3,c}(0,\rho)`$.

The finite counterterm needs $`\widehat w_1`$ through
$`\epsilon^2`$ and $`\widehat w_2`$ through $`\epsilon^1`$.
Each selected three-loop graph is needed through $`\epsilon^0`$.
This is the budget for coefficients of $`\log W`$. It is not the
larger budget needed to reconstruct $`\log W`$ by multiplying ordinary
Wilson-loop coefficients $`W_1,W_2,W_3`$.

More explicitly, write $`w_1=\widehat w_1/C_F`$ and
$`\widehat w_2/C_F=C_Aw_{2,A}+Tw_{2,f}`$. The three counterterms,
with their whole three-loop color monomials removed, are

```math

\mathrm{CT}_{Af}=-\frac{11}{6\epsilon}w_{2,f}
+\frac{2}{3\epsilon}w_{2,A}
-\frac{11}{18\epsilon^2}w_1+\frac{5}{24\epsilon}w_1,

```


```math

\mathrm{CT}_{Ff}=\frac{w_1}{8\epsilon},\qquad
\mathrm{CT}_{ff}=\frac{2w_{2,f}}{3\epsilon}
+\frac{w_1}{9\epsilon^2}.

```

The authoritative inputs are the [one-loop result](../../results/2026-10-23/sources/benchmarks/renormalization/one-loop/result-r001.json),
[bare lower-order budget](../../results/2026-10-23/sources/benchmarks/renormalization/three-loop/lower-order-inputs-r001.json),
[two-loop finite/RG result](../../results/2026-10-23/sources/benchmarks/renormalization/two-loop/result-r002.json),
and [counterterm table](../../results/2026-10-23/sources/benchmarks/renormalization/three-loop/counterterms-r001.json).

## Independent pole prediction and its limits

After removing the common external $`C_F`$, let $`[\cdot]_c`$
extract respectively $`C_AT`$, $`C_FT`$, or $`T^2`$. Define

```math

B=[b_0^2]_c,\quad P=[b_0 k_2]_c,\quad
Q=[b_0 h_2]_c,\quad D=[b_1]_c,

```

where $`k_2,h_2`$ retain their complete lower-order color dependence
and have only the common external $`C_F`$ removed. With external
reference values $`k_{3,c},h_{3,c}`$, the predicted negative Laurent
coefficients are

```math

Z_{3,c}^{[-6]}=Z_{3,c}^{[-5]}=0,\qquad
Z_{3,c}^{[-4]}=-\frac{11}{9}B,

```


```math

Z_{3,c}^{[-3]}=-\frac23 BL+\frac59P+\frac89D,
\qquad
Z_{3,c}^{[-2]}=\frac23(P+D)L+\frac13Q-\frac29k_{3,c},

```


```math

Z_{3,c}^{[-1]}=-\frac{2k_{3,c}L+h_{3,c}}3.

```

Here $`Z_{3,c}^{[p]}`$ denotes a coefficient of $`\log Z_4`$, not
of $`Z_4`$. The finite RG equation checked with the same external
constants is

```math

\partial_LF_{3,c}=[2b_0F_2+b_1F_1]_c-2k_{3,c}L-h_{3,c}.

```

The reference constants come from the preserved [UV bridge builder](../../results/2026-10-23/sources/workflow/checks/n4-w3-normalization-20261019/uv-bridge/build_bridge.py).
The historical literature retrieval is not repeated. The
[bridge record](../../results/2026-10-23/sources/workflow/checks/n4-w3-normalization-20261019/uv-bridge/uv-coordinate-bridge-r002.json)
binds its source documents and the published finite-line relation
$`\Gamma_\Box=4(\Gamma_{\cap}-\Gamma_{\wedge})=4(\gamma_G-2B_\delta)`$.
The three-loop non-cusp comparison is therefore conditional on that
external operator/normalization dictionary. It is independent of
extracting $`h_3`$ from these selected poles, but is not a primitive
three-loop UV integration.

## Why the coordinate-scale bridge matters

The [October 19 derivation](../../results/2026-10-23/sources/benchmarks/renormalization/theory/n4-w3-normalization-20261019/derivation-r002.md)
obtains the scale relation from the raw $`\pi`$, $`2\pi`$, and
$`4\pi`$ factors and a common bare coupling, before comparing any
$`h_n`$:

```math

\mu_\Box^2=\frac{e^{2\gamma_E}}2\mu_{\overline{\mathrm{MS}}}^2
=2\mu_x^2,\qquad
\mu_{\mathrm{FGM,local}}^2=e^{2\gamma_E}\mu_{\overline{\mathrm{MS}}}^2
=2\mu_\Box^2.

```

For $`s=2n\!\cdot\!n^*`$,
$`\mu_\Box^2(n\!\cdot\!n^*)=\mu_x^2s`$. The matching closed-contour
logarithms therefore give $`h_n=\Gamma_{\Box,n}/2`$. Substituting the
local integration scale into the closed-contour logarithm without its
geometric conversion introduces an erroneous $`\log2`$ shift.
The earlier scale discrepancy remains a historical, convention-bound
diagnostic; it is not the final status of the later UV comparison.

The new audit repeats the scale algebra. The original source reading,
primitive one-loop numerical controls, and the finite-line geometry
argument are reused from the October 19 record and are not rerun here.

## Finite-function coverage and remaining independence gap

Through $`\epsilon^0`$, the three compact matter sums need only these
fixed moments, zeta values, and polynomials in $`L,\rho`$:

| Color | Moment symbols in the compact finite expression |
| --- | --- |
| $`C_FC_AT_Fn_f`$ | $`M_{0,1},M_{1,1},M_{2,1},M_{1,3},N_{1,2}`$, including $`M_{0,1}^2`$ |
| $`C_F^2T_Fn_f`$ | $`M_{1,1}`$ |
| $`C_F(T_Fn_f)^2`$ | $`M_{1,1},M_{2,1}`$ |

The source identities are
$`M_{0,1}=-(\rho^2+\pi^2)/2`$ and
$`N_{1,2}=M_{1,3}-\rho M_{1,2}`$. No
$`\tau_{D,4}^{+}`$, $`E`$, or $`\Omega_6`$ is needed by these final
matter sums at this order. The two-loop $`\tau`$ term first enters
the relevant higher-order TGV04 expansion at $`\epsilon^2`$, beyond
the $`\epsilon^1`$ counterterm input budget. This is a statement about
the consumed order and these channel sums, not a claim that every
individual graph or all-order expression has a minimal M/N basis.

The later [QCD total](../../results/2026-10-23/sources/workflow/checks/n4-w3-normalization-20261019/uv-bridge/qcd-total-up-to-C3-r001.json)
has an unknown $`C3_{N4}`$ only in its $`C_FC_A^2`$ finite term.
Each of the three selected matter finite functions is complete in its
accepted computational sense and is explicitly checked to contain no
such unknown constant.

The comparison with that later total is **shared-input consistency**:
its matter terms reuse the same earlier integral information. Exact
M/N-to-GPL agreement and finite RG agreement do not provide an
independent determination of an arbitrary finite shape. In particular,
an independent original-integral validation of the complete
$`C_FC_AT_Fn_f`$ finite shape is not established by this audit.

The earlier [matter-independent checker](../../results/2026-10-23/sources/workflow/checks/matter-independent-20261017/verify.py)
provides a different accepted-$`C/J`$-root reconstruction for the
$`C_F^2T_Fn_f`$ and $`C_F(T_Fn_f)^2`$ channels, with direct scalar
geometry numerical checks. Those results or a fresh replay must retain
their own source bindings and numerical error scope; they are not
performed by `verify_uv.py`. They do not independently derive the
upstream physical operator or establish the mixed-channel finite shape.

The [current matter RG checker](../../results/2026-10-23/sources/workflow/checks/matter-rg-20261017/check.py)
freshly compares the consumed current native coefficients and rebuilds
the scale/counterterm algebra when replayed. For the 31-graph mixed
channel it reuses the historical exact HyperInt $`r`$-to-$`z`$
canonicalization witness. Applying the old extra minus signs again to
the corrected native leg/Gamma values would double-count the sign
conversion. A successful replay is not a fresh evaluation of that
mixed-argument graph sum. Its literal internal date `2026-10-17` is the
checker's historical metadata, not the new replay date.

## Replaying this audit

With Python and SymPy available, run the script with the repository or
the package's `sources` directory as the read-only source root and a
new output directory:

```text
python verify_uv.py --root <source-root> --output-dir <new-output-directory>
```

Existing `uv-checks.json` or `uv-source-hashes.json` files cause an
explicit refusal to overwrite. All consumed source files are hashed
before use and rechecked before publication. The checker changes no
scientific source, selector, stored proof, or historical evidence.
