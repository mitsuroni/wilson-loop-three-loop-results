# Integral derivations and calculation records

This is the current reading layer under the [integral documentation standard](../workflow/integral-documentation-standard.md). It follows each initial integral through its actual reduction, evaluated subintegrals and Laurent-coefficient assembly. The physical definitions remain those of the [single current convention](../workflow/current-conventions.md).

The first trial covers [two one-loop integrals](one-loop/README.md), through epsilon cubed, and [twelve two-loop integrals](two-loop/README.md), through epsilon squared. Every entry has a derivation, reproduction instructions, validation scope, exact input/result bindings and code/data/check indices. Nontrivial endpoint proofs and large but explicit HLOG jets have mathematical appendices.

The next batch applies the same standard to [all39 three-loop6e integrals](three-loop/README.md), through the finite term. Six shared mathematical routes support39 own-source derivations; their [scientific audit and actual replay](../workflow/checks/integral-documentation-6e-20261103/README.md) are recorded separately. Each entry preserves its own ordered root, normalization and exact relation to any representative.

The following batch adds [15 two-point and exchange insertions](three-loop/README.md#two-point-and-exchange-insertions): eleven Pi1-plus-exchange, two serial Pi1 and two single Pi2 sources. Their [current mathematics, replay and per-entry audit](../workflow/checks/integral-documentation-insertions-20261104/README.md) retain the complete transverse tensors, endpoint terms and external colour operators.

The [trial report](../workflow/checks/integral-documentation-low-loops-20261103/README.md) separates the independent scientific reading audit from the actual [96-coefficient replay](../workflow/checks/integral-documentation-low-loops-20261103/replay.json). Replaying an adopted analytic formula does not establish a new integration of its original kernel. In particular the difficult TGV04 historical integration remains an explicitly adopted computational result.

Shared original programs, primitive outputs and source snapshots live in [low-loop-runtime](families/low-loop-runtime/README.md). Its archived Markdown preserves provenance; use the articles linked here for the current mathematical exposition. Wolfram and Maple files are kept where they actually participated; Python adapters run the supported low-loop batch. There is no requirement to manufacture three versions of an algorithm.

The earlier public packets remain at their original locations. Their packaging/replay acceptance does not imply acceptance under the new continuous-derivation standard; use the new per-entry scientific reviews. The15 insertion articles are included in this release; [public-copy validation](../validation/insertions-20261104/README.md) records the execution scope. Other three-loop families await the same reading audit.
