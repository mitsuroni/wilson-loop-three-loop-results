# One- and two-loop scientific reading audit

All 14 pilot articles now connect their actual starting integral to every requested Laurent coefficient. Shared physical roots and evaluated function/CAS lemmas remain explicit adopted inputs. This is a documentation/implementation audit; it is not a second cold integration or a new formal proof.

| Entry | Final reading status | Mathematical chain checked |
|---|---|---|
| [2e-01](../../../../integrals/one-loop/2e-01/derivation.md) | 正文可连续阅读 | Momentum phase and Fourier root → z(1-s31)s41 → two endpoint primitives → orbit/crossing I=-G/e² → Gamma expansion through degree five. |
| [2e-02](../../../../integrals/one-loop/2e-02/derivation.md) | 明确采用已求引理，链已接通 | Momentum geometry → full removable difference quotient → centered moment generator and Mellin-Barnes evaluation → crossing → explicit M coefficients and finite convolution. |
| [4e-01](../../../../integrals/two-loop/4e-01/derivation.md) | 正文可连续阅读 | Two original ordered adjacent chains → symmetric triangle integral → exact J=z^(2e)/(4e⁴) → I=-Z²/(16e⁴) → degree-six Gamma recurrence. |
| [4e-02](../../../../integrals/two-loop/4e-02/derivation.md) | 明确采用已求引理，链已接通 | Ordered upper endpoint gives 1-s^e; its weighted kernel and reflected partner give the exact V(2e) identity → physical crossing → finite d_q/p_q convolution. |
| [4e-03](../../../../integrals/two-loop/4e-03/derivation.md) | 正文可连续阅读 | Independent endpoint integrals leave an actual ordered mixed Beta triangle → R2 Gamma ratio → centered crossing 1/cosh(e log z) → finite expansion through degree six. |
| [4e-04](../../../../integrals/two-loop/4e-04/derivation.md) | 明确采用已求引理，链已接通 | Original four-dimensional ordered root → two linear primitives → antisymmetric zero addition → K Euler derivative/Beta summation with boundary → same-Q substitution → analytically fixed endpoint conversions → evaluated M/N mother → Q2/prefactor convolution. |
| [4e-05](../../../../integrals/two-loop/4e-05/derivation.md) | 明确采用已求引理，链已接通 | Its own lower endpoint gives 1-(1-t)^e; a justified square-variable swap and explicit reflected weighted source prove B2[J05]=B2[J02] → same normalized answer with its distinct original root retained. |
| [4e-06](../../../../integrals/two-loop/4e-06/derivation.md) | 明确采用已求引理，链已接通 | Original disjoint variable pairs → complement map with unit Jacobian → product of evaluated opposite kernels → crossing → finite V²/prefactor convolution. |
| [3e-3g-01](../../../../integrals/two-loop/3e-3g-01/derivation.md) | 正文可连续阅读 | Explicit simplex and ordered-path substitutions with Jacobians → three scale integrals and q primitive → combined regulated Beta difference → crossing → degree-six Gamma/Beta coefficient rule. |
| [3e-3g-02](../../../../integrals/two-loop/3e-3g-02/derivation.md) | 明确采用已求引理，链已接通 | Full orientation-preserving ordered-domain/simplex map → direct S and P pullback → exact 01 lemma → same normalization and degree-six expansion. |
| [3e-3g-03](../../../../integrals/two-loop/3e-3g-03/derivation.md) | 明确采用已求引理，链已接通 | Numerator derivative → triangular boundary/contact IBP → simplex scale and pushed affine-face weights → exact F_e constant → two moment sources → closed Mellin value and finite coefficient convolution. |
| [3e-3g-04](../../../../integrals/two-loop/3e-3g-04/derivation.md) | 明确采用已求引理，链已接通 | Full polynomial numerator → six explicitly restricted faces plus recoil → evaluated cusp/D/H reductions → Euler-regularized recoil → complete adopted HLOG jets → explicit bracket, j2 and a2 convolution → optional M/N/tau compression. |
| [2e-Pi1-01](../../../../integrals/two-loop/2e-Pi1-01/derivation.md) | 正文可连续阅读 | Adopted full transverse self-energy tensor → Fourier/differential parameter kernel → adjacent geometry collapses complete numerator → elementary endpoints → physical current sign and degree-five rational/Gamma expansion. |
| [2e-Pi1-02](../../../../integrals/two-loop/2e-Pi1-02/derivation.md) | 明确采用已求引理，链已接通 | Same full tensor root → mixed-derivative identity → four corner values with regulator strip → actual opposite-kernel primitive → evaluated V(2e) → separate corner and moment pole budgets. |

## Hard examples and actual review findings

For 4e-04 the review checks the actual K derivative, the z=0 constant, the same-Q back-substitution, the two endpoint-conversion constants and the complete Q2-to-I convolution. The corrected M05 constant and N04/M05 conversion make its compact highest-order expression traceable. For TGV04 it checks all six boundary restrictions, the actual D/H/recoil densities and orders, the common-strip Euler cutoff proof, the complete adopted jet tables and the j2-to-a2 prefactor budget. The cutoff proof concerns chart variables r,u; it does not discard the nonzero Wilson-path boundary terms.

The new TGV04 program adopts already assembled normalized I1/I2 and accepted low coefficients, then verifies their representations. Equations (8)–(9) expose the historical mathematical D/H/O assembly; they are not claimed as a new live Python jet convolution. No remaining prose/implementation gap was found after the documented corrections.

The reviewer executed an exact transcription check on every one of the 168 coefficients in the six new HLOG tables; all equal the original raw jet outputs. This checks transcription, not the original integration. See [the check](hlog-transcription-check.json) and [its reproducible script](check_hlog_transcription.py).

## Execution and limits

The integrator's [actual replay](../replay.json) reports 14 integrals and 96 coefficient positions. Its 80 one-loop checks include 27 numerical controls; four negative controls are also reported. The two-loop run reports 149 checks and two negative controls, with 68 separate integral controls. A copied runtime closure also executed successfully and reproduced the result bytes. These batch totals are not independent checks of each original integral.

The source Fourier/QFT and prepared C/J rules, established moment dictionary and historical HyperInt integrations retain their original scope. All articles restrict their conclusion to real z>0 and local meromorphic epsilon continuation. This work establishes no arbitrary physical complex-z branch theorem, z=-1 boundary result, RG sufficiency theorem or formal integral proof. No commit or GitHub publication was performed.

## Standard decision

Retain standard 1.0 unchanged. The trial demonstrates that continuous reasoning is practical for both elementary and hard examples. Precise shared adoption is allowed; unexplained source/result links are not substitutes for graph-specific transformations or endpoint constants. The report evaluates prose, implementation, reproduction, mathematical validation and publication separately.

See the complete [per-entry evidence JSON](audit.json), [before/after baseline findings](baseline-findings.md) and [forward-test assessment](standard-forward-test.md). The JSON identifies concrete mathematical passages and source evidence; headings, page count and file count were not used as scientific acceptance criteria.
