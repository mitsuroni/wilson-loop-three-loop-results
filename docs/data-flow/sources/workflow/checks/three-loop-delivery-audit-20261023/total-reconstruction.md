# QCD total reconstruction and the remaining finite constant

23 October 2026. This supplement records research completed after the October 18 delivery. It is additional total-level information, while the selected QCD diagram catalogue and its missing coefficients remain unchanged.

## Completed five-family difference

The [five-family calculation](../n4-five-family-20261019/README.md) evaluates
\[
\Delta w_3=w_3^{\rm QCD}-w_3^{\mathcal N=4}
\]
from 31 connected orbit contributions: 11 four-endpoint Pi1, two Pi1 chains, 12 TGV-leg Pi1, four Gamma1 and two Pi2. Shared pure-gluon contributions cancel in this difference. The calculation uses the stated common bare normalization and DRED field content, then applies the QCD coupling counterterms with the original lower-order scales. It does not replace the three-loop calculation by a two-point insertion alone.

The [scale-stripped bare difference](../n4-five-family-20261019/five-family-bare-r001.json) and [post-counterterm difference](../n4-five-family-20261019/renormalized-difference-r001.json) contain four colors through epsilon zero. The three matter colors reproduce the selected QCD results. The AA difference supplies the additional data required to reconstruct the pure QCD total from a correctly matched N4 total.

## Later resolution of the local versus closed-contour scale issue

The [October 19 UV bridge](../../../benchmarks/renormalization/theory/n4-w3-normalization-20261019/derivation-r002.md) distinguishes the coordinate scales obtained from raw pi, 2pi and 4pi factors:
\[
\mu_{\rm closed}^2=\mu_{\rm local}^2/2=2\mu_x^2,
\qquad \mu_{\rm closed}^2(n\cdot n^*)=\mu_x^2 s.
\]
It supplies a research-level interpretation of the earlier log(2) mismatch without fitting the shift from a desired three-loop coefficient. The [bridge data](../n4-w3-normalization-20261019/uv-bridge/uv-coordinate-bridge-r002.json) and [checks](../n4-w3-normalization-20261019/uv-bridge/README.md) record low-order holdouts, 28 exact checks and six original one-loop numerical controls. The four QCD cusp and non-cusp coefficients agree with the independently adopted reference relations in this research convention.

This later result must accompany the [October 17 correction](../../../benchmarks/renormalization/theory/matter-independent-20261017/derivation-r002.md). The older residual is not erased; describing its dictionary as simply unresolved today would omit subsequent work. The adopted three-loop finite-line non-cusp relation is a factorization prediction, not a direct reintegration of every three-loop closed-contour graph or a formal proof.

## Four-color total known up to one constant

The [research total](../n4-w3-normalization-20261019/uv-bridge/qcd-total-up-to-C3-r001.json) has explicit status `complete_poles_and_finite_shape_up_to_one_undetermined_constant_research`, `absolute_F3_complete=false`, and `selector_changed=false`. Its unknown is `C3_N4 in AA finite term only`.

In the documented normalization,
\[
\frac{F_3^{\mathcal N=4}}{C_FC_A^2}
=C_3^{\mathcal N=4}-\frac{11\pi^4}{720}(L^2-\rho^2/4)
-\left(6\zeta_5+\frac{11\pi^2\zeta_3}{36}\right)L,
\]
\[
C_3^{\mathcal N=4}
=\left.[\epsilon^0]\frac{(\log W_{\mathcal N=4})_{\lambda^3}}{C_FC_A^2}\right|_{L=0,\,s/t=1}.
\]
The notation `(log W)_(lambda^3)` means its lambda-cubed coefficient. At that same point,
\[
F_{3,AA}^{\rm QCD}=C_3^{\mathcal N=4}+\Delta F_{3,AA},
\qquad \Delta F_{3,AA}\simeq-0.5259594010446349257640.
\]
The numerical difference is a recorded evaluation, not a value for the unknown total. Shifting the constant changes only the AA epsilon-zero result; RG and reflection do not determine it.

The [latest total-input boundary](../../../benchmarks/renormalization/theory/n4-w3-normalization-20261019/derivation-r012.md) distinguishes known amplitude, insertion-ratio and threshold-soft constants from the same-definition closed-Wilson constant. A complete BDS amplitude or a vanishing remainder does not by itself supply the finite regulator conversion. The threshold-soft total is a concrete additional input, but its finite contour/ordering/regulator map remains to be established; it is not adopted as C3. The intervening scalar-jet work is a possible calculation route, not an accepted substitute total.

## How to present this beside the diagram tables

Give this route a separate supplement with the bare difference, counterterms, four-color pole and finite formulas, source/scale dictionary and exact C3 definition. Explicitly distinguish (1) selected individual diagram results, (2) direct diagram-by-diagram color assembly, and (3) this total reconstruction using additional theory input. Its total-level information does not determine the missing I15/I16/I24 coefficients separately and does not close their proof or calculation tasks.

The two epsilon budgets must also remain separate. The existing connected-exponent/counterterm route needs one loop through epsilon squared and two loops through epsilon one. A direct calculation of `W3-W1W2+W1^3/3` from the ordinary unlogged W coefficients requires W1 through epsilon four and W2 through epsilon two for the finite coefficient under the recorded leading-pole orders. The [current extraction contract](../n4-w3-normalization-20261019/uv-bridge/C3-extraction-contract-r002.json) documents this difference.
