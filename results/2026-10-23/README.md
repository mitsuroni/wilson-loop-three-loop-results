# Three-loop matter channels

**Af, Ff and ff: diagram data, bare sums, renormalization and finite results.** Prepared 23 October 2026, in the same observable and normalization used for the accepted two-loop report.

The three complete matter channels can be carried through the workflow separately from the unfinished three-loop AA channel. Counterterms retain their specified lower-order color dependencies. Each channel has explicit coefficients through epsilon zero at every stage. All three totals can be written using the original M/N moments, logarithms and zeta values; no additional tau or E block is needed at this order.

| Channel | Color factor | Unique orbit contributions | Lower-order dependencies |
|---|---|---:|---|
| [Af](channel-ca-nf.md) | CF CA TF nf | 31 | w1 and both colors of w2 |
| [Ff](channel-cf-nf.md) | CF^2 TF nf | 2 | w1 |
| [ff](channel-nf2.md) | CF (TF nf)^2 | 2 | w1 and the matter color of w2 |

## Reading order

1. [Conventions and function dictionary](conventions.md): observable, scale, Hlog, M/N, tau relations and color weights.
2. [Diagram catalogue](diagrams.md): 31 distinct orbits and 35 channel occurrences, with native and weighted Laurent coefficients.
3. [Af](channel-ca-nf.md), [Ff](channel-cf-nf.md), [ff](channel-nf2.md): scale-stripped U, scale-restored bare exponent, separate CT, post-CT coefficients, log Z4 poles, F and f, k and h.
4. [Consumed lower-order inputs](lower-orders.md): the positive epsilon orders needed for finite counterterms.
5. [Validation and limits](validation-scope.md), [numerical samples](samples.md), and [reproduction instructions](reproduce.md).

## Checks and acceptance boundary

The existing three-color checker was rerun: **452 checks passed**, including 7 negative controls. The compact stage conversion adds **324 exact checks**, including 4 negative controls. The separate UV/finite-source comparison passes its recorded checks in [uv-checks.json](uv-checks.json). The diagram ledger independently checks 217 Af source positions; numerical moment evaluation supplies 18 samples and 24 precision/crossing checks. These suites overlap and should not be added as independent evidence counts.

The data are complete computational results through finite for these three matter channels. The Af diagram sum reuses its bound historical HyperInt argument-conversion witness. The later UV scale bridge supports the cusp/non-cusp comparisons in its research convention. Independent original-input finite reconstruction is available for ff and Ff as [preserved October 17 evidence](sources/workflow/checks/matter-independent-20261017/verification.json), not a newly rerun integral evaluation. Its [raw predictions](sources/workflow/checks/matter-independent-20261017/raw-source-predictions.json) and [handoff](sources/workflow/checks/matter-independent-20261017/HANDOFF-r002.md) are included. A comparable independent Af finite calculation and a same-definition external full-finite three-loop comparison remain open. No formal proof is claimed.

This delivery therefore completes the available assembly, renormalization and reporting chain for the three matter channels, while keeping those further acceptance checks explicit. It does not complete the four-color three-loop result. The unfinished individual AA diagrams and the unknown AA total constant do not appear in these matter formulas.

## What is included for the recipient

This folder includes offline HTML, readable Markdown, exact machine data, the function conversion dictionary, raw source snapshots with SHA-256 hashes, reproduction scripts and check records. It is a new three-channel supplement; the October 18 all-diagram package remains a separate historical delivery.

The finite term is the epsilon-zero coefficient at the named stage. The color-channel decomposition is a separate operation. Both the bare finite coefficient and the renormalized finite function are supplied, so “constant term” cannot be confused with the four-color classification or a number at a special kinematic point.
