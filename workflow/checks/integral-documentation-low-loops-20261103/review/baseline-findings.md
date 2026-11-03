# Independent forward review: baseline findings

Reviewer: `/root/lowloop_standard_review`. Date: 2026-11-03.

Scope: the fixed `wilson-integral-derivations` standard 1.0 and its application to all two one-loop and twelve two-loop entries. This review reads mathematical statements and the original derivation/code. It does not perform a new independent integration, alter scientific selections, or grant formal proof status. Runtime replay is owned by the workspace integrator and is assessed separately from mathematical continuity.

The baseline is the frozen public checkout at commit `7e77ddb83d8aa9622c05a553281b262f0e9dc7ca`, under `.wilson/publications/wilson-loop-three-loop-results`. That baseline remains unchanged. The findings below identify missing exposition, not a refutation of its accepted results.

| Entry | Concrete baseline finding | Material that can supply the missing bridge |
|---|---|---|
| 2e-01 | The path phase, distance polynomial, two endpoint primitives, orbit/crossing factors and Gamma expansion are already connected. The new article needs an explicit highest-order budget and precise historical/current check labels. | One-loop source reproducer and existing endpoint/Gamma derivation. |
| 2e-02 | The difference quotient, centered moment generator, Mellin-Barnes/Meijer-G representation and moment compactification are already connected. A finite coefficient rule makes the positive-order endpoint explicit. | Existing one-loop derivation and its moment source. |
| 4e-01 | Elementary ordered endpoint integral is present; complete coefficient extraction needs a displayed finite recurrence and all stated zero terms. | Original ordered integral plus shared Gamma expansion. |
| 4e-02 | `derivation.md:17-18` says the weighted Euler source gives the mother formula without the actual crossing-weight identity. | The author's explicit weighted kernel and two reflected linear primitives. |
| 4e-03 | The Beta-triangle integral is present; the original source domain, exact prefactor and finite expansion budget should be connected locally. | Existing Beta computation and shared Gamma/Beta recurrence. |
| 4e-04 | `derivation.md:33-45` introduces K, asserts its 4F3 evaluation, then jumps to the Q mother formula. The graph-specific derivative, constants and Q back-substitution are absent from the short article. | `benchmarks/2ne/theory/generation/4e04-direct-mn.md`: steps 1–13, Appendix A substitutions, B boundary constants, C K derivative/Beta proof, D linear primitives. This material contains a real derivation and can be exposed rather than replaced by an opaque old-result adoption. |
| 4e-05 | The endpoint primitive is explicit, but reflection to the 4e-02 weighted object does not display the actual change of variables and scale factors. | The explicit weighted-kernel identity also used for 4e-02. |
| 4e-06 | Disjoint parameter factorization and the reflected opposite kernel are explicit. A local moment/prefactor convolution must expose all requested orders. | Existing factorization and the evaluated one-loop moment generator. |
| 3e-3g-01 | `derivation.md:13` jumps from the simplex kernel to a one-dimensional Beta difference. The integrated scales and their Jacobians should be written. | Original `lecture-notes/integrals.tex` and the associated Maple kernel. |
| 3e-3g-02 | `derivation.md:11-18` narrates a reflection without its actual path/simplex pullback. Equality of scalar answers does not by itself prove equality of the original oriented kernels. | Original lecture notes' explicit reflection map and its domain/numerator substitution. |
| 3e-3g-03 | `derivation.md:20` skips the contact-face weight and the remaining kernel's conversion to the two mother functions. | Original lecture notes and boundary calculation. |
| 3e-3g-04 | `derivation.md:30-39` names six required positive-order jets but leaves A, B and the logarithms in the sample recoil density undefined. Merely linking complete scripts does not finish a continuous reader. | Original lecture notes/appendices, regularized coefficient kernels, six raw HLOG outputs and exact prefactor convolution. The Euler/contact manipulations additionally require an actual common convergence strip and boundary argument. |
| 2e-Pi1-01 | The transverse operator, exact elementary kernel and current sign are explicit. Its effective-operator scope and full prefactor convolution should remain visible. | Existing transverse-kernel derivation and source r002. |
| 2e-Pi1-02 | `derivation.md:13-22` says four actual corners produce the rational term without displaying all four values and their signs. | The mixed derivative's endpoint values at D(0,0)=1, D(1,1)=z and the two vanishing corners. |

## Standard forward-test observation

Standard 1.0 correctly distinguishes mathematical exposition, implementation, execution, validation and publication. Its hardest practical question is the boundary between a precisely adopted evaluated lemma and a linked black box. The present hard examples do not justify weakening the standard: actual graph-specific bridge material exists locally and is being written into the new article or a mathematical appendix. A valid adoption must still identify the exact kernel, domain/order, branch, epsilon range, evaluated result, boundary constants and substitution into the present I. Passing a script or matching an old coefficient table is separate evidence.

No header count, formula count or file count is used for the findings above. The final review will reread the new mathematical chain entry by entry and retain any remaining gaps explicitly.
