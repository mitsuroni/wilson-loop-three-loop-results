# Standard 1.0: forward-test assessment

The reviewer applied `workflow/skills/wilson-integral-derivations/SKILL.md` and `references/standard.md` to actual original low-loop articles, research notes and source programs before reviewing the new articles. This is a scientific reading audit, not an independent second evaluation of the original integrals. No new CAS integration, Lean proof, Git commit or publication was performed by the reviewer.

## What the standard usefully required

The one-loop adjacent exchange is a suitable ordinary example: two endpoint primitives already expose the double pole, so a complete article needs the source factor, convergence strip, crossing/orbit normalization and Gamma-series budget. It does not need artificial Maple/Wolfram/Python copies or a long generic catalogue of unused functions.

The difficult two-loop examples reveal why source lists are insufficient. For 4e-04, the old short article's K-to-4F3 assertion and Q-to-M/N jump become a real chain only after including the actual Euler derivative of K, its Beta summation, the endpoint at z=0, the back-substitution into Q and the two analytic conversion constants. Original research material for these steps exists. Thus satisfying the standard requires exposing the research, rather than weakening the requirement to permit an unexplained old-result link.

For 3e-3g-04, a named set of six jets is not their evaluation. The article must identify their regularized densities, actual variable orders and endpoint conventions, give their exact evaluated HLOG values or an explicit finite equivalent, and show which jet/prefactor convolution supplies each target coefficient. Contact/recoil terms and a convergence/boundary argument must remain visible. A successful historical HyperInt run supports this adopted calculation; a new Python representation replay has a different scope.

The weighted 4e-02/05 pair also provides a useful non-CAS test. Displaying the two reflected weighted sources proves their exact crossing identity. Equality of the final coefficient lists alone would not establish that identity for their differently ordered roots.

## Actual findings made during forward review

1. The original short articles omitted particular graph-specific bridges; they were identified by equation and paragraph, as recorded in `baseline-findings.md`. The authors located original detailed research and wrote those bridges, rather than changing the standard.
2. A draft 4e-04 moment identity briefly had the constant `-31*pi^6/126` in M05. The correct fifth Taylor derivative of `pi[cot(pi*a)-exp(ell*a)csc(pi*a)]+ell` also includes the cotangent contribution and is `-pi^6/2`. At z=1 the independent elementary check is `T(a)=-pi*tan(pi*a/2)`. The final article has the corrected constant. The author and reviewer identified this during overlapping work; this is not a change to the frozen accepted result.
3. A draft TGV03 coefficient rule used an undefined Greek nu where the preceding definition used u. This is a continuity defect in an executable mathematical recurrence even though the intended rule is straightforward; it was sent for correction.
4. Initial new `input.json` records identified accepted results but did not bind the original ordered C/J record locally. The integrator added exact graph/native/prepared fields, source hashes, parameter order, bounds and simplex measure. The integrator's `source-binding-check.json` reports all twelve two-loop bindings passed.
5. The new wrappers correctly delegate to a same-loop batch. They are not twelve distinct cold integration programs. Reading the implementation confirms that ordinary entries expand evaluated Gamma/moment mother functions, while TGV04 explicitly adopts integrated HLOG output. This distinction must survive into README, reproduction and validation pages.
6. The easiest exchange article initially repeated a long unused M/N dictionary. The author removed that unnecessary material. The standard's “按需” rule therefore helps both completeness and reading speed; mathematical continuity is not a minimum page count.

## Adoption, implementation and mathematical status

A valid adoption is a stated evaluated lemma with a precise kernel, measure/order, branch, epsilon range, full needed result and a concrete substitution into the present graph. A file path, an old result table or a hash alone is insufficient. Conversely, this task does not require redoing a valid prior integration merely to improve its exposition.

The actual shared runtime was inspected. Its generic engines generate coefficients from mother-function Taylor series and finite convolutions before comparison with old tables. TGV04 is expressly different: it assigns its accepted low coefficients and reads already assembled normalized I1/I2 from `compact.json`, then checks complete HLOG coordinates. It does not perform a new live D/H/O jet convolution. The historical mathematical jet convolution remains exposed in the article. Representation equality shares the moment dictionary. It is useful evidence, but is not an independent integration of the original five-dimensional kernel or a formal proof.

The integrator reports fresh replay runs of 80 one-loop checks (including 27 numerical controls, with four negative controls also reported), 149 two-loop checks with two negative controls, and 68 separate two-loop numerical controls. Those executions are integrator evidence, not reviewer executions. The reviewer separately compared all 168 exact coefficients in the six new HLOG tables with their original raw outputs; this is a transcription check, not an integration. The final per-entry report binds the actual receipt paths. None grants arbitrary complex-z continuation or the z=-1 boundary; the stated pilot domain remains real z>0 with local meromorphic epsilon continuation.

## Decision on the fixed standard

Retain version 1.0 unchanged. It provides a workable distinction between continuous mathematical prose, actual implementation, executed replay, mathematical validation and publication. No concrete case in this pilot warrants lowering its requirement for graph-specific equalities or endpoint constants.

For future application, begin by extracting the actual equation chain and subintegral budget, then write only the needed graph-specific transformations and adopt exact shared lemmas. This is an application recipe, not a relaxed acceptance rule. The full final per-entry findings remain authoritative when any one article has an unresolved gap; a global format or execution pass cannot conceal it.
