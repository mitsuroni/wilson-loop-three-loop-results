# Accepted two-loop workflow and its implications for a three-loop delivery

Information audit, 23 October 2026. Task: WILSON-5.15. Audit worktree base: `34bb1af2727185e38f7c351482975d212ebbcdf3`. Scientific sources were read from `D:/wilson-agent`; source pointers below are relative to that main workspace and include line numbers. No source was edited, no builder or CAS was run, and no historical mathematical check was represented as a new run. This file is the only audit output authored by this reviewer.

The accepted two-loop delivery contains an evaluated integral catalogue **and a separate, completed physical assembly and validation chain**. Its model is therefore more extensive than a directory of diagram coefficients. The current report is the 28-page English r006 edition; the 48-page v5 board guide is a later presentation of the established science with a three-loop progress continuation. Neither presentation revision silently expands the mathematical acceptance scope.

## 1. Authoritative entry points

| Purpose | Source pointer |
| --- | --- |
| Current roadmap and separation of scientific and presentation states | `workflow/master-plan.md:213`, `:258`, `:276`; `workflow/data-flow.md:21`, `:34`, `:44` |
| Source/normalization/assembly narrative | `benchmarks/renormalization/two-loop/derivation.md:5` |
| Current immutable physical result | `benchmarks/renormalization/two-loop/result-r002.json:1`; migration record at `:449` |
| Current English report and section/page map | `benchmarks/renormalization/two-loop/report/README.md:3`, `:9`, `:11` |
| Reader-facing acceptance contract | `benchmarks/renormalization/two-loop/preread/VALIDATION.md:3` |
| Reader-facing conventions and function dictionary | `benchmarks/renormalization/two-loop/preread/CONVENTIONS.md:3` |
| Exact assembly, renormalization, RG and finite-reference formulas | `benchmarks/renormalization/two-loop/report/renormalization.tex:1` |
| External-reference conversion and branch argument | `benchmarks/renormalization/two-loop/reference-comparison.md:3` |
| Recipient data and check-script interfaces | `benchmarks/renormalization/two-loop/preread/README.md:27` |
| Current board-talk scope | `benchmarks/renormalization/two-loop/board-talk-plan/README.md:1` |

The historical source audit and original RG record used legacy Pi1 results with an explicit minus-sign bridge. They are valuable evidence but are not the current assembly instructions. Current r002 Pi1 values already incorporate the physical insertion sign and enter the sum with **plus** sign. The migration preserved all physical bare coefficients, finite functions, contour counterterms and anomalous dimensions. See `workflow/checks/self-energy-signs/README.md:5`, `:7`, `:21`, and `benchmarks/renormalization/two-loop/result-r002.json:449`.

## 2. Objects and diagram coverage

The reporting sequence is

```text
observable and coupling order
  -> allowed attachments, zero/scaleless exclusions and exponent colors
  -> cyclic orbits, multiplicities and ordered endpoints
  -> complete momentum-space amplitude A_r
  -> raw C_r J_r on the original ordered domain
  -> explicit preparation/color/sign conversion
  -> cyclic sum and scale extraction defining normalized I_r
  -> actual evaluation of J_r and Laurent coefficients of I_r
  -> full-color, scale-restored bare exponent
  -> coupling counterterms and contour renormalization
  -> finite exponent and independent physical comparisons
```

This is the explicit chain in `benchmarks/renormalization/two-loop/board-talk-plan/upstream-data-flow-audit.md:7`. The distinction at its final integral step matters: a formula defining I in terms of J is not an evaluation of J.

At order g0^4 the relevant families are six exchange representatives, four tree-level three-gluon-vertex representatives, and two one-loop self-energy representatives. The independent edge-labelled counts are 19 + 16 + 6 = 41 before the C4 quotient, giving 6 + 4 + 2 = 12 representatives. Mirror-related or equal-valued integrals are not deleted. The recorded producer uses the combinatorial `generate_diagrams -> classify_orbits -> definition_exporter` path; the audit did not run the QGRAF executable. See `benchmarks/renormalization/two-loop/board-talk-plan/upstream-data-flow-audit.md:34`, `:44`, `:57`, `:80`.

The exchange exponent color is -CF CA/2. The CF^2 part of the full-color diagrams forms W1^2/2 and is removed when passing to log W. Consequently the twelve normalized entries already assemble the exponent; subtracting W1^2/2 again is an error. See `benchmarks/renormalization/two-loop/board-talk-plan/upstream-data-flow-audit.md:46` and `benchmarks/renormalization/two-loop/integrand-handoff/ASSEMBLY.md:37`.

The original ordered bounds and integration order are retained. C_raw, prepared C, and orbit-normalized I are distinct objects. The report keeps the full amplitude on the left through the momentum/Fourier/parameter transformations, then identifies normalization changes explicitly. The r002 report revision verified 66 phase/Gamma/geometry/kernel connections rather than replacing the source integrals with unexplained normalized formulas. See `workflow/checks/two-loop-report-20261015-r002/README.md:9`.

## 3. Exact conventions and bare assembly

The observable is the fundamental, vacuum-normalized, non-self-intersecting lightlike quadrilateral in massless QCD, in Feynman gauge. The trace is normalized by 1/Nc. With T = TF nf and TF = 1/2,

\[
D=4-2\epsilon,\qquad s=2v_{12}=2v_{34},\qquad t=2v_{14},
\qquad z=s/t>0,\qquad \ell=\log z,
\qquad \lambda=\alpha_s^{\overline{\rm MS}}/\pi.
\]

The initial analytic region is s,t < 0. Scale powers carry the common +i0 prescription, and kernel powers and polylogarithms must continue with it. Positive-real-z checks do not establish arbitrary complex continuation. See `benchmarks/renormalization/two-loop/preread/CONVENTIONS.md:5` and `benchmarks/renormalization/two-loop/integrand-handoff/ASSEMBLY.md:5`.

Every delivered scalar series is

\[
I_r(\epsilon,z)=\sum_{p=-4}^{2}a_{r,p}(z)\epsilon^p+O(\epsilon^3).
\]

These are ordinary coefficients, with no factorial multiplier. The catalogue has 84 scalar positions, including zeros, and 28 further color-weighted Pi1 positions. Equal results retain separate orbit identities. See `benchmarks/renormalization/two-loop/preread/CONVENTIONS.md:15` and `workflow/checks/two-loop-report-20261015/README.md:9`.

Let IE be the sum of six exchange I values, IV the sum of four vertex I values, and IPi the sum of the two current scalar Pi1 values. Then

\[
a(\epsilon)=\frac{10-6\epsilon}{3-2\epsilon},\qquad
b(\epsilon)=-\frac{8(1-\epsilon)}{3-2\epsilon},
\]
\[
\mathcal S_A=I_E+I_V+aI_\Pi,\qquad
\mathcal S_f=bI_\Pi,\qquad
\widehat w_2=C_FS_2(\epsilon)
\left[C_A\mathcal S_A+T\mathcal S_f\right].
\]

The two physical colors are CF CA and CF TF nf. The normalized I already includes its orbit/crossing factors and one Euler compensation. Use either the weighted Pi1 tables or the scalar I times a,b; never both. Prepared C is not another factor to multiply into an already evaluated I. See `benchmarks/renormalization/two-loop/preread/CONVENTIONS.md:23` and `benchmarks/renormalization/two-loop/derivation.md:18`.

The accepted scale dictionary is

\[
\mu_x^2=\frac{e^{2\gamma_E}}4\mu_{\overline{\rm MS}}^2,
\quad X_s=\mu_x^2(-s+i0),\quad X_t=\mu_x^2(-t+i0),
\quad L=\frac{\log X_s+\log X_t}{2},
\]
\[
\log X_s=L+\ell/2,\quad \log X_t=L-\ell/2,
\qquad S_2=X_s^{2\epsilon}+X_t^{2\epsilon}
=2e^{2\epsilon L}\cosh(\epsilon\ell).
\]

L is the overall physical scale logarithm, whereas ell is the shape logarithm. Do not import a local integral's symbol L without translating its meaning. Sources: `benchmarks/renormalization/two-loop/derivation.md:120`, `:159`; `benchmarks/renormalization/two-loop/report/renormalization.tex:9`.

Before counterterms, the dimensionless sums exhibit three explicit cancellations:

\[
[\epsilon^{-4}]\mathcal S_A=-1/16+1/32+1/32=0,
\]
\[
[M_{21}][\epsilon^0]\mathcal S_A=3/32+3/32-1/16-1/8=0,
\]
\[
[M_{01}^2][\epsilon^0]\mathcal S_A
=1/32+1/32-1/48-1/16+1/48=0.
\]

The remaining leading poles are -11/(48 epsilon^3) and 1/(12 epsilon^3) for SA and Sf. The bare finite sum still contains M11; it is not finite before renormalization. The complete two-color sum table is `benchmarks/renormalization/two-loop/derivation.md:87`; the explicit scale-restored bare result is `benchmarks/renormalization/two-loop/derivation.md:135`.

For any summed channel \(\mathcal S=\sum_{p=-3}^\infty s_p\epsilon^p\),

\[
[\epsilon^0](S_2\mathcal S)=2s_0+4Ls_{-1}
+(4L^2+\ell^2)s_{-2}
+\left(\frac83L^3+2L\ell^2\right)s_{-3}.
\]

Thus adding only the finite entries in the diagram table is not the physical finite answer. If one expands individual diagrams before the fourth-order-pole cancellation, their fourth-order scale contributions must also be retained until cancellation. Source: `benchmarks/renormalization/two-loop/report/renormalization.tex:17`.

## 4. Coupling, contour poles, finite result and anomalous dimensions

With \(\beta_0=(11C_A-4T)/3\), the fixed convention is

\[
\lambda_0(4\pi)^\epsilon e^{-\gamma_E\epsilon}
=\mu_{\overline{\rm MS}}^{2\epsilon}\lambda
\left[1-\frac{\beta_0}{4\epsilon}\lambda+O(\lambda^2)\right],
\quad
\mathcal B_2=\widehat w_2-\frac{\beta_0}{4\epsilon}\widehat w_1.
\]

The full one-loop input, with its own scale dependence, is

\[
\frac{\widehat w_1}{C_F}=-\frac2{\epsilon^2}-\frac{2L}{\epsilon}
-L^2+\frac{\ell^2}{4}+\frac{\pi^2}{3}
+\epsilon\left[-\frac{L^3}{3}
+L\left(\frac{\ell^2}{4}+\frac{\pi^2}{3}\right)
-\frac{M_{11}}2-\frac{2\zeta_3}{3}\right]+O(\epsilon^2).
\]

Source: `benchmarks/renormalization/two-loop/report/renormalization.tex:34`. It cannot be replaced by a one-loop coefficient multiplied by S2.

Contour renormalization is \(W_B=Z_4W_R\), with \(Z_{4,n}\) denoting the coefficient of **log Z4**, not of Z4. The strict-pole prescription is

\[
Z_{4,1}=C_F\left(-\frac2{\epsilon^2}-\frac{2L}{\epsilon}\right),
\quad Z_{4,2}=\operatorname{Pole}_\epsilon\mathcal B_2,
\quad F_2=[\epsilon^0]\mathcal B_2.
\]

Extracting one's own poles and subtracting them guarantees finiteness by definition. The independent acceptance test instead predicts the counterterm from reference anomalous dimensions:

\[
Z_{4,2}^{\rm ref}=\frac{3C_F\beta_0}{8\epsilon^3}
+\frac{C_F\beta_0L/4-K_2^{\rm ref}/2}{\epsilon^2}
-\frac{K_2^{\rm ref}L+\gamma_2^{\rm ref}/4}{\epsilon},
\]
\[
\operatorname{Pole}_\epsilon
\left(\widehat w_2-\frac{\beta_0}{4\epsilon}\widehat w_1
-Z_{4,2}^{\rm ref}\right)=0.
\]

All three remaining pole orders in both color channels pass. The poles have no z dependence at fixed L and are at most linear in L. Sources: `benchmarks/renormalization/two-loop/report/renormalization.tex:64`, `:106`, `:166`; `benchmarks/renormalization/two-loop/preread/VALIDATION.md:14`.

The reference anomalous dimensions use

\[
\frac{d\log W_R}{d\log\mu^2}=-2K_F(\lambda)L-\frac12\gamma_4(\lambda),
\quad K_F=C_F\lambda+K_2\lambda^2+\cdots,
\quad\gamma_4=\gamma_2\lambda^2+\cdots,
\]
\[
K_2=C_F\left[C_A\left(\frac{67}{36}-\frac{\pi^2}{12}\right)-\frac59T\right],
\]
\[
\gamma_2=C_F\left[C_A\left(\frac{202}{27}-\frac{11\pi^2}{36}-7\zeta_3\right)
+T\left(-\frac{56}{27}+\frac{\pi^2}{9}\right)\right].
\]

The finite and dimensional-pole RG checks are, respectively,

\[
\partial_LF_2-\frac{\beta_0}{4}F_1=-2K_2L-\frac12\gamma_2,
\qquad F_1=C_F\left[-L^2+\ell^2/4+\pi^2/3\right],
\]
\[
\partial_LZ_{4,2}-2\epsilon Z_{4,2}-\frac{\beta_0}{4}Z_{4,1}
=2K_2L+\frac12\gamma_2.
\]

The classical -epsilon lambda term is essential in the dimensional equation, and coupling running is essential in the finite equation. Source: `benchmarks/renormalization/two-loop/report/renormalization.tex:141`.

The complete accepted F2 is written explicitly, by color, at `benchmarks/renormalization/two-loop/derivation.md:217` and `benchmarks/renormalization/two-loop/report/renormalization.tex:117`; the machine expression is `benchmarks/renormalization/two-loop/result-r002.json:/F2_over_CF`. It is cubic in L, and its only non-elementary shape function through the finite term is M11, with coefficient -CF beta0/8. The M21 and M01^2 cancellations are properties of the assembled finite order, not an all-order simplification theorem. Recovering the Wilson loop itself requires

\[
W_R=1+\lambda F_1+\lambda^2\left(F_2+\frac12F_1^2\right)+O(\lambda^3).
\]

## 5. Independent finite comparison and its limits

Finite RG is blind to \(F_2\mapsto F_2+h(z)\). Consequently agreement of anomalous dimensions does not certify the scale-independent shape function or finite constant. The full finite comparison is an additional acceptance layer, not a decorative reference. See `benchmarks/renormalization/two-loop/report/renormalization.tex:205`.

The independent QCD expression is the authors' `WRen.m` ancillary file for Chen et al., arXiv:2510.07377, Zenodo record 17200473. In its coupling \(a_{\rm anc}=\alpha_s/(4\pi)=\lambda/4\), write \(W_R^{\rm anc}=1+a_{\rm anc}\mathcal W_1+a_{\rm anc}^2\mathcal W_2+\cdots\). Then

\[
F_2^{\rm ref}=\left.\frac{\mathcal W_2-\mathcal W_1^2/2}{16}\right|_{\mathcal D},
\]

where the one-loop comparison fixes

\[
s_a^2/t_a^2=z,\quad
\log(s_a^2\mu_a^2)=L+\ell/2,\quad
\log(t_a^2\mu_a^2)=L-\ell/2.
\]

After matching polylogarithm branches, the complete analytic difference is zero. The original author-file hashes and publisher MD5 values are retained in `workflow/checks/rg-two-loop/reference-source.json`. The comparison does not use the reference finite answer as an integration input. Five polylogarithm identities are checked by differentiation plus boundary constants; the negative-axis logarithm has +i pi and the positive Li3 cut uses the compatible negative-imaginary boundary value. The mapped imaginary terms and CF^2 coefficient vanish. Sources: `benchmarks/renormalization/two-loop/reference-comparison.md:5`, `:21`, `:40`, `:42`; `benchmarks/renormalization/two-loop/report/renormalization.tex:214`.

The cusp reference is arXiv:1510.07803, Eq. (6.12). `ADs.m` provides the independently mapped threshold anomalous dimension; in this convention gamma4 = 2 h2 and the author's gammaW = -2 GammaThr. No extra scale shift is introduced after the one-loop dictionary. See `benchmarks/renormalization/two-loop/reference-comparison.md:44`.

The later N=4-to-QCD reconstruction is a useful additional holdout, distinct from the QCD author-file comparison. It predicts the QCD answer using the corrected N=4 total, the two Pi1 geometries and the one-loop counterterm, without consuming exchange/TGV/QCD-total coefficients until validation. Its exact bridge is

\[
w_{2,B}^{\rm pred}=w_{2,B}^{\mathcal N=4}
+C_FS_2\left(B_1^{\rm QCD}+4C_A\right)(I_{B_A}+I_{B_O}),
\]

with the N=4 DRED operator -4CA including 6+2epsilon real scalars. It passes 28 exact equalities, four deliberately incorrect controls, and three 60-digit formula-evaluation points. It reuses accepted self-energy geometry; it is not independent reintegration of all diagrams. The two-loop success does not justify replacing only two-point self-energy operators at three loops. See `workflow/checks/n4-two-loop-reconstruction-20261018/README.md:3`, `:25`, `:34`, `:58`, `:79`.

## 6. Expansion budgets

| Target | Required inputs | Source |
| --- | --- | --- |
| Two-loop finite exponent F2 | Two-loop scalar integrals through epsilon^0; full one-loop exponent through epsilon^1; Pi1 rational color operator through epsilon^3 because the adjacent scalar starts at epsilon^-3 | `benchmarks/renormalization/two-loop/derivation.md:187`; `benchmarks/renormalization/two-loop/source-audit.md:265`; `workflow/checks/n4-two-loop-reconstruction-20261018/README.md:9` |
| Finite part of the already summed two-loop bare channels | Scale S2 through epsilon^3 after the epsilon^-4 sum cancels | `benchmarks/renormalization/two-loop/report/renormalization.tex:17` |
| Delivered individual two-loop integral catalogue | All twelve scalar I through epsilon^2; 84 scalar plus 28 color-weighted positions, with O(epsilon^3) remainder | `workflow/checks/two-loop-report-20261015/README.md:9` |
| Entire two-loop post-counterterm expression through epsilon^2 | One-loop input through epsilon^3; this extended renormalized expression was **not** the report's physical target | `benchmarks/renormalization/two-loop/report/renormalization.tex:58`; `benchmarks/renormalization/two-loop/preread/CONVENTIONS.md:65` |
| Three-loop finite exponent | One-loop w1 through epsilon^2, two-loop w2 through epsilon^1, three-loop bare w3 through epsilon^0, plus the declared epsilon-dependent operator multipliers and distinct S1/S2/S3 scales | `workflow/data-flow.md:44`; `workflow/project-map.md:49` |

The Pi1 budget is experimentally discriminating, not merely administrative: truncating B1 at epsilon^2 produces a finite error CF(CA+4T)/162 in the N=4 reconstruction. See `workflow/checks/n4-two-loop-reconstruction-20261018/README.md:66`.

The report also distinguishes normalized scalar integrals from a reconstruction of the full bare W2: the square of a one-loop series beginning at epsilon^-2 can require still higher one-loop coefficients. `benchmarks/renormalization/two-loop/integrand-handoff/ASSEMBLY.md:87` gives the specific epsilon^1 target, for which that square requires the one-loop series through epsilon^3. Budgets must therefore name both the target order and the mathematical object.

## 7. What was actually accepted and checked

| Evidence layer | Recorded scope and result | Source |
| --- | --- | --- |
| Source and assembly audit | 382 source/exact checks, including six expected-error controls; independent geometry, C4 counts, colors, ordered-domain shuffle, Dyson sign, scales and tensor endpoint cancellation | `workflow/checks/rg-two-loop-audit/README.md:12` |
| Physical calculation | 334 source/exact checks; 12 selected results, M/Hlog connections, grouped versus per-diagram convolution, full one-loop counterterm, pole and finite RG, correct and incorrect self-energy signs | `workflow/checks/rg-two-loop/README.md:14` |
| Full external comparison | 36 analytic/source checks; coupling/log/object conversion, five polylogarithm identities, zero full finite residual, threshold anomalous dimension | `workflow/checks/rg-two-loop/README.md:15` |
| Numerical controls | 50 controls at 70 working digits: five M11 defining integrals, five reflections, 30 original unreduced author-expression comparisons, ten finite-RG scale differences; maximum normalized discrepancy about 1.05e-41, not an interval certificate | `workflow/checks/rg-two-loop/README.md:16` |
| Current sign migration | 29 unit/regression tests and 855 migration checks; current r002 physical result equals the externally matched legacy result | `workflow/checks/self-energy-signs/README.md:21` |
| English integral report | 427 catalogue/assembly/source/algebra/control checks: 185 exchange, 50 TGV, 92 self-energy, 39 main-text, 61 catalogue checks | `workflow/checks/two-loop-report-20261015/README.md:34` |
| Momentum-to-parameters report revision | 66 checks of phases, Gamma factors, geometry and all twelve prepared kernels; old evidence reused only after source/evidence hash matching | `workflow/checks/two-loop-report-20261015-r002/README.md:13` |
| Renormalization report revision | 52 checks with five negative controls; scale convolution, counterterm, all reference-predicted poles, finite coefficient, dimensional/finite RG, and RG's finite-function blind spot | `workflow/checks/two-loop-report-20261015-r003/README.md:17` |
| Recipient symbolic export | Actual Wolfram Language 14.3/Wolfbook evaluation of the exported portable script: all 24 checks True, no kernel messages | `workflow/checks/two-loop-preread-20261015/README.md:18` |
| Recipient archive | CSV/JSON coefficient equality, orbit metadata, source/evidence hashes, copied bytes, links/includes, manifest, CRC and byte-for-byte ZIP read-back | `workflow/checks/two-loop-preread-20261015/README.md:20` |

These are overlapping evidence collections, not an additive number of independent proofs. The wrong self-energy sign is a particularly useful warning: both signs passed formal cancellation of negative powers in the dimensional derivative. The wrong sign instead fails local cusp structure and the independent reference. See `workflow/checks/rg-two-loop/README.md:19`.

Other recorded negative controls detect repeated Euler normalization, confusing G2 with G^2, omission of the opposite-edge endpoint, repeated orbit multiplicity, repeated subtraction of W1^2/2, incorrect coupling-counterterm sign, omitted classical/finite beta terms, altered cusp coefficients and a scale-independent finite shift. See `workflow/checks/rg-two-loop-audit/README.md:18` and `workflow/checks/two-loop-report-20261015-r003/README.md:17`.

## 8. Scientific acceptance versus presentation acceptance

**Scientific acceptance** binds a fixed source, definition, selected result, normalization and comparison. A complete Laurent expansion, an all-order representation, physical sum agreement and formal proof are separate claims. Sum agreement permits compensating individual-diagram errors and does not validate all positive-epsilon coefficients. The integral catalogue therefore retains its own expansion checks. See `benchmarks/renormalization/two-loop/preread/VALIDATION.md:39`.

**Presentation acceptance** checks transcription, source reuse, rendered content, fonts, clipping, diagrams, links and exact delivered bytes. In the report workflow a build creates a candidate and resets visual review to pending; every final page is rendered and inspected, then publication verifies the reviewed PDF hash and dependencies before updating the stable deliverable. Historical editions remain in `before/`. This is documented at `workflow/checks/two-loop-report-20261015-r003/README.md:23`, `:36` and `workflow/checks/two-loop-report-20261015-r006/README.md:42`.

The board-talk r003 additionally checked the displayed one-parameter VM/VN generating sums: coefficients q=0..12 and six numerical comparisons, while reusing the underlying physical result. The v5 continuation preserves the two-loop source files and equations, shifts unchanged PDF pages, and adds seven three-loop pages. It explicitly adds no complete w3/RG result. Planned 55+15+5 or 40+15+5 minute routes have not been measured by rehearsal. See `workflow/checks/two-loop-board-talk-20261017-r003/README.md:5`, `:9`; `workflow/checks/two-loop-board-talk-20261017-r005/README.md:5`, `:24`.

The recipient supplement contains a PDF, structured JSON and Wolfram data, a twelve-row index, 112 coefficient rows, explicit tau GPL terms, original prepared C/J snapshots, editable TeX, conventions, validation scope, source provenance, a manifest and a portable algebraic check script. Its Wolfram file returns data without installing function definitions; special functions remain symbolic unless a reader supplies an evaluator. See `benchmarks/renormalization/two-loop/preread/README.md:27`, `:54`, `:73`.

## 9. Unresolved matters and non-claims

- Full formal proof remains `not_run`. Packaging, CAS checks, source hashes and task completion do not change that status.
- The report does not independently reintegrate all twelve original multidimensional inputs or prove arbitrary Fourier contour deformations. Internal self-energy species reduction is reused with a stated source bridge.
- General complex-z continuation and all cut-winding paths are outside the accepted positive-real-z comparison.
- The full post-counterterm two-loop answer through epsilon^2 is not claimed, even though each delivered scalar integral is expanded through epsilon^2.
- Finite-order function cancellations do not establish all-order closure or minimality. Vertex 04 already contributes tau_D,4^+/8 at epsilon^2. Its later all-order research does not retroactively change the old finite physical acceptance.
- The earliest source audit's explicit minus sign belongs to legacy Pi1 I; current r002 values use a plus assembly sign. Historical evidence stays bound to its original files.
- External two-loop finite matching establishes its particular coupling, scale, theory and branch dictionary. It is not a blanket license for an arbitrary higher-loop QCD/N=4 finite-scheme replacement.

Sources for these boundaries: `workflow/checks/rg-two-loop/README.md:23`; `workflow/checks/two-loop-report-20261015-r002/README.md:15`; `workflow/checks/two-loop-report-20261015-r003/README.md:25`; `benchmarks/renormalization/two-loop/preread/VALIDATION.md:39`; `workflow/checks/n4-two-loop-reconstruction-20261018/README.md:79`.

## 10. Checklist to carry into the richer three-loop package

1. **Declare the target observable and snapshot.** Identify log W versus W, the scalar-integral normalization, the four physical color structures, selectors and source hashes. Count orbit representatives separately from independent integrals or master functions.
2. **Retain the source-to-integral chain.** For each family include graph/topology, ordered endpoints, momentum/parameter kernel, raw/prepared C/J, normalized I, result revision, coefficient range and exact dependencies. Mark effective insertion operators as such.
3. **Present both individual results and physical color assembly.** Give per-diagram Laurent data, family/color subtotals and each assembled bare channel with exterior scales. Record precisely where multipliers, C4, crossing, color projection and Euler compensation enter.
4. **Show cancellations as explicit identities.** Distinguish cancellations in the dimensionless bare sum, scale-restored expression, coupling-counterterm sum and contour subtraction. Do not call a bare sum finite merely because its highest poles cancel.
5. **Publish the coefficient budget before the calculation.** For finite three-loop renormalization include w1 through epsilon^2, w2 through epsilon^1 and w3 through epsilon^0; separately budget epsilon-dependent insertion operators. Preserve S1, S2 and S3 in their own terms.
6. **Give the full renormalization map.** Show the coupling substitution, the lower-order induced terms, the strict-pole contour factor and finite extraction by physical color. Name log Z rather than an ambiguous Z coefficient.
7. **Separate extraction from independent prediction.** Present the poles/anomalous dimensions extracted from the calculation and the reference-predicted pole template with its own coupling, scale and theory dictionary. A self-subtracted expression is finite by construction.
8. **Report multiple validation scopes.** Include dimensional RG with its classical beta term, finite RG with coupling running, locality, crossing, and any independent original-integral or reference comparisons. Keep numerical precision controls distinct from certified error bounds.
9. **Supply the complete finite function or state its precise remaining gap.** Finite RG does not fix an arbitrary scale-independent shape function. Any available full finite reference or N=4 difference route needs an explicit finite-scheme/scale/branch bridge and a fixed finite constant; missing information must remain visible.
10. **Keep partial channels partial.** A coefficient catalogue can be ready to circulate while a physical channel, finite constant or independent comparison remains unfinished. Use a status matrix separating Laurent completion, physical assembly, compact all-order representation, independent validation and formal proof.
11. **Make the package useful without the workspace.** Include the explanatory PDF, per-diagram catalogue, portable data, original input snapshots, complete special-function definitions/word tables, dependency closure, normalization guide, source manifest and a bounded, portable algebraic checker. Do not silently replace missing coefficients by zero.
12. **Validate delivery independently.** Compare exported coefficients to frozen pointers, check displayed formulas and diagrams, render every PDF page, test offline navigation, verify archive CRC/read-back, and bind the final reviewed files by hash. Record old evidence as reused and new checks as actually run.

The two-loop example supports this checklist as a reporting and acceptance structure. It does not establish that any currently open three-loop channel or finite comparison has already passed those stages; those claims require the companion current-three-loop audit.
