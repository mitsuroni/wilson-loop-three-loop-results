# From diagrams to the renormalized Wilson loop

An illustrated guide to the Wilson project's scientific data flow · 23 October 2026

[中文阅读说明](reading-guide-zh.md) · [Research data](../results/index.md) · [Source register](source-register.md)

Our objective is to calculate the logarithm of a lightlike quadrilateral Wilson loop, including its divergences and finite dependence on the contour. A diagram is the starting specification of a contribution. A Montopar parameter integral, an evaluated scalar integral, a color sum and a renormalized finite function are successive mathematical objects. Keeping their definitions distinct makes the result reproducible and prevents factors from being counted twice.

The diagram below follows the main calculation. Each arrow represents a mathematical transformation whose normalization and source version must be recorded. The numerical and formal checks discussed later accompany these arrows; they do not replace them.

![End-to-end scientific data flow](figures/pipeline.svg)

*Figure 1. Read from left to right within each row, then follow the connecting arrow to the next row. Defining I and evaluating I are separate steps. Lower-loop inputs enter the three-loop counterterms through a second branch.*

## 1. Specify the observable before generating diagrams

We work with a vacuum-normalized Wilson loop in the fundamental representation, on a finite oriented lightlike quadrilateral in massless QCD. The initial analytic region is $`s,t<0`$, with $`z=s/t>0`$, and $`D=4-2\epsilon`$. Path ordering, the Fourier convention and the common $`+i0`$ prescription belong to the definition.

At fixed coupling order, the graph stage specifies Wilson-line attachment points, their edge labels and order, internal connections, vertex orientations, color traces and combinatorial multiplicities. Graphs related by a cyclic relabelling are organized into $`C_4`$ orbits. Reflection relations are recorded separately; equal integral values do not automatically remove different orbit IDs.

The project's `qgraf-pipeline` name should not be read as a claim that the QGRAF executable generated every delivered diagram. The audited two-loop path uses `generate_diagrams -> classify_orbits -> definition_exporter`, a combinatorial generator. It yields 12 orbit representatives: six exchange, four three-gluon-vertex and two self-energy representatives, obtained from 41 allowed edge-labelled configurations. These counts refer to that specified contour calculation, not a theorem about every QCD graph set.

The outputs are structured definitions such as `orbit_definitions.json` and `integrals.json`. They tell Montopar what to calculate. A topology drawing by itself does not fix a mathematical integral: endpoint ordering, phases, numerator and domain are also needed. See the [two-loop upstream audit](sources/benchmarks/renormalization/two-loop/board-talk-plan/upstream-data-flow-audit.md).

## 2. Montopar turns the graph into C times J

For a representative diagram $`r`$, begin with its full momentum-space contribution $`\mathcal A_r`$, obtained from Wilson insertions, propagators, vertices and the ordered color trace. Fourier transforms, internal subgraph integrals and parameterization put it into the form

```math

\mathcal A_r=C_{r,\mathrm{raw}}\,J_{r,\mathrm{raw}},\qquad
J_r=\int_{\mathcal D_r}d\mathbf s\,d\boldsymbol\alpha\;j_r(\mathbf s,\boldsymbol\alpha;\epsilon,z).

```

Here **C** carries factors independent of the remaining integration variables: coupling/color factors, phases, powers of pi, Gamma functions, regulator-dependent rational factors and an extracted kinematic scale, as appropriate to that family. **J** is the remaining parameter integral, including its numerator, denominator powers, integration domain, measure and order. A JSON field called `J.integrand` is only the integrand of J, not its evaluated value.

Different families reach this split by different operations. Exchange diagrams mainly require propagator Fourier transforms. Vertex diagrams also introduce a common internal center and Feynman parameters, integrate that center and retain the vertex numerator. Self-energy diagrams first integrate an internal loop and then carry its tensor operator through the Fourier transformation. A single universal Gamma prefactor would therefore be incorrect.

If integration by parts produces several terms, the output is instead

```math

\mathcal A_r=\sum_t C_tJ_t.

```

The field `combined_C` already includes the overall C. Boundary/contact terms and common regulators must be retained. One must not multiply by the global C a second time or discard a derivative term merely because it looks longitudinal. See the [Montopar contract](sources/benchmarks/theory/montopar/io-contract.md).

![What changes between raw C/J, prepared inputs and I](figures/normalization.svg)

*Figure 2. The parameter integral retains its domain. Color/Euler preparation and orbit normalization act at specified interfaces. Once their factors are included in I, they are not multiplied again during color assembly.*

## 3. Prepare the inputs without changing their meaning

Raw output is not automatically the input convention used by the downstream solver. Preparation records original C, source/code hashes, the adopted sign convention, color transformations and Euler compensation. It preserves the original variables and ordered bounds unless an explicit transformation is separately justified.

The default pure-exchange corpus contains **47 roots** across 2e/4e/6e. Manifest 0.3 uses the pre-hypercube ordered variables and outer-to-inner integration order. Its color operation reconstructs the full fixed-basis coordinates under $`C_A^2=1+2C_AC_F`$, then selects $`a_1`$. This is neither a large-N limit nor a search for a monomial. It is a special upstream projection, not the four-channel physical color decomposition of the complete three-loop answer.

After the Gamma-weight check, prepared C includes exactly one factor $`e^{-n\gamma_E\epsilon}`$. The preserved full-color information is still needed later. Self-energy and vertex families use their own documented native/prepared contracts; all 149 catalogue entries do not pass through the same default 47-root preparation routine.

The original ordered representation and a unit-cube representation are alternatives with different variables, bounds and Jacobians. They are not two contributions to be added. A successful input-schema check also does not prove a change of integration order or an analytic continuation. See [input preparation](sources/preprocessing/README.md).

## 4. Define the orbit-normalized I, then actually evaluate it

The purpose of $`I_r`$ is to package a representative's contribution in a convention suitable for physical assembly. It can include the cyclic orbit sum, selected web color normalization, Euler compensation and removal of the exterior scale. The precise definition belongs to each family/result revision.

The accepted two-loop convention gives a concrete example. Put

```math

\mathcal B[J]=\frac{J(\epsilon,z)+z^{2\epsilon}J(\epsilon,1/z)}{1+z^{2\epsilon}},\qquad
\mathcal G(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon).

```

For an exchange orbit of size $`m_r`$,

```math

I_r=-\frac{m_r}{16}\mathcal G(\epsilon)^2\mathcal B[J_r].

```

Why these factors? A cyclic shift exchanges s and t. Expressed using the original reference scale, the rotated integral carries $`z^{2\epsilon}J(\epsilon,1/z)`$. Summing distinct cyclic images supplies the orbit multiplicity; extracting the symmetric external scale supplies the denominator $`1+z^{2\epsilon}`$. This formula is specific to the stated exchange convention: the vertex and self-energy prefactors differ.

The formula **defines** I. To calculate it, J must still be integrated. The two-loop `4e-01` example closes the entire step:

```math

J_{01}=\frac{z^{2\epsilon}}{4\epsilon^4},\qquad m_{01}=4,
\qquad I_{01}=-\frac{\mathcal G(\epsilon)^2}{16\epsilon^4}.

```

Expanding the last expression gives actual Laurent coefficients. More complicated diagrams require Beta integrals, differential equations, parameter integrations, boundary terms, moment functions or hyperlogarithms. Defining an unevaluated generator is useful, but does not by itself give all coefficients of the original diagram. The [worked normalization derivation](sources/benchmarks/renormalization/two-loop/board-talk-plan/upstream-data-flow-audit.md) supplies the ordered integral and intermediate factors.

## 5. Store evaluated coefficients with their definitions

An evaluated result has the form

```math

I_r(\epsilon,z)=\sum_{p=p_{\min}}^{p_{\max}}a_{r,p}(z)\epsilon^p+O(\epsilon^{p_{\max}+1}).

```

These are ordinary coefficients, without factorials. An omitted coefficient is not automatically zero. Stored zero, zero established by pole-order metadata, a physically vanishing contribution and an uncomputed position must be distinguished.

The result revision fixes the definition, domain, branch, normalization, source bindings and available powers. `current.json` or a family catalogue selects a revision; a larger revision number alone is not an adoption rule. A later numerical or representation update may leave the consumed analytic fields unchanged, but that must be checked at field level.

Hlog, M/N and extended blocks are representations of coefficient functions, not extra color channels. In the delivered matter totals,

```math

N_{i,j}=M_{i,j+1}-\rho M_{i,j},\qquad
M_{0,1}=-\frac{\rho^2+\pi^2}{2},\qquad \rho=\log z.

```

The [function dictionary](../results/conventions.md) defines the integrals and Hlog kernel, and the [conversion data](../../results/2026-10-23/conversion-dictionary.json) supplies explicit maps. The three matter totals through epsilon zero need no independent tau or E block. This does not establish an M/N-only reduction for every diagram or every higher order.

## 6. Assemble the full colors of log W

The target observable is the exponent, with $`\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi`$:

```math

\log W=\lambda w_1+\lambda^2w_2+\lambda^3w_3+\cdots.

```

If starting from the ordinary expansion of W, exponentiation requires $`w_2=W_2-W_1^2/2`$ and $`w_3=W_3-W_1W_2+W_1^3/3`$. Our adopted web/color assembly already forms $`w_n`$, so those products must not be subtracted again.

For $`T=T_Fn_f`$, the three-loop physical color basis is

```math

C_FC_A^2,\quad C_FC_AT,\quad C_F^2T,\quad C_FT^2.

```

For each channel c, first form a scale-stripped sum

```math

U_{3,c}(\epsilon,z)=\sum_r q_{r,c}(\epsilon)I_r(\epsilon,z).

```

The factor $`q_{r,c}`$ must be read from the actual source contract. It can be epsilon dependent, so the product requires a Laurent convolution. A stored color-resolved coefficient may already include it. The internal quark $`C_F`$ contributing to $`C_F^2T`$ is not a reducible external exponentiation product.

The Af channel provides a useful complete trace. Define

```math

a=\frac{10-6\epsilon}{3-2\epsilon},\qquad B_f=-\frac{8(1-\epsilon)}{3-2\epsilon}.

```

Its 31 orbit contributions comprise 11 Pi1 exchange insertions weighted by $`B_f`$, two chains, two Pi2 entries, 12 TGV-leg insertions weighted by $`B_f`$, and four Gamma1 quark entries. A scalar chain would receive $`2aB_f`$; the delivered chain color field already contains this factor. The Pi2 and Gamma entries are also read in their designated color fields. The other matter colors use two chain contributions (ff) and two Pi2 contributions (Ff), giving **35 color occurrences across 31 unique orbit IDs**.

| Example orbit | Meaning of its delivered contribution | Where to inspect it |
|---|---|---|
| `4e-Pi1-01` | Native scalar I, followed by Bf convolution for Af | [Diagram page](../results/diagram-4e-Pi1-01.md) |
| `2e-Pi1+Pi1-02` | Separate already-resolved Af and ff channel coefficients | [Diagram page](../results/diagram-2e-Pi1+Pi1-02.md) |
| `2e-Pi2-02` | Separate Af and Ff operator/color contributions | [Diagram page](../results/diagram-2e-Pi2-02.md) |
| `3e-Gamma1-03` | Native quark-channel coefficient for Af | [Diagram page](../results/diagram-3e-Gamma1-03.md) |

![A Pi2 insertion example](../../results/2026-10-23/assets/diagrams/2e-Pi2-02.svg)

*Figure 3. A Pi2 operator insertion connects two contour edges. The blob is an effective internal subgraph; the picture represents the contour attachment topology. Its Af and Ff contributions retain distinct color fields. It is not an additional multiplicity to apply to the stored integral.*

Historical leg/Gamma expressions used an explicit minus-sign bridge. The current native C/J/I convention already incorporates that sign, so the current external assembly sign is plus. The [ledger](../../results/2026-10-23/diagram-ledger.json) records exact fields and multipliers, preventing a second sign change.

## 7. Restore the scale and construct the counterterms

The symmetric exterior factor is

```math

S_n=e^{n\epsilon(L+\rho/2)}+e^{n\epsilon(L-\rho/2)},\qquad
\widehat w_{3,c}=S_3U_{3,c}.

```

Here $`L=(\log X_s+\log X_t)/2`$, $`X_s=\mu_x^2(-s+i0)`$, $`X_t=\mu_x^2(-t+i0)`$, and $`\mu_x^2=e^{2\gamma_E}\mu_{\overline{\mathrm{MS}}}^2/4`$. The scale logarithm L differs from the shape logarithm rho. Restoring S3 mixes poles into the finite coefficient; summing only the finite entries of individual diagrams would miss these terms.

Coupling renormalization then gives

```math

\widetilde w_3=\widehat w_3-\frac{2b_0}{\epsilon}\widehat w_2
+\left(\frac{b_0^2}{\epsilon^2}-\frac{b_1}{2\epsilon}\right)\widehat w_1,

```


```math

b_0=\frac{11C_A-4T}{12},\qquad b_1=\frac{34C_A^2-20C_AT-12C_FT}{48}.

```

The terms containing lower orders are the counterterm CT3. Their color algebra must be performed before extracting a channel. Each lower order keeps its own S1 or S2; replacing these by S3 changes the result.

![Lower-order dependencies of the three matter channels](figures/dependencies.svg)

*Figure 4. The three matter channels do not consume the missing three-loop AA coefficients. They do consume lower-order data through their own counterterms. “Independent channel” refers to the first fact, not to independence from all lower loops.*

For $`\widehat w_1=C_Fv_1`$ and $`\widehat w_2=C_F(C_Av_{2,A}+Tv_{2,f})`$, the stripped counterterms are

```math

\mathrm{CT}_{Af}=-\frac{11v_{2,f}}{6\epsilon}+\frac{2v_{2,A}}{3\epsilon}
-\frac{11v_1}{18\epsilon^2}+\frac{5v_1}{24\epsilon},

```


```math

\mathrm{CT}_{Ff}=\frac{v_1}{8\epsilon},\qquad
\mathrm{CT}_{ff}=\frac{2v_{2,f}}{3\epsilon}+\frac{v_1}{9\epsilon^2}.

```

Thus the finite three-loop answer needs one loop through epsilon squared, two loops through epsilon one and three loops through epsilon zero. This is the budget for the adopted connected-exponent route. Reconstructing log W directly from unlogged W products has a different, higher budget.

## 8. Separate contour poles from the finite answer

With $`W_B=Z_4W_R`$, minimal contour subtraction gives, channel by channel,

```math

\mathcal Z_{3,c}=\operatorname{Pole}_{\epsilon}\widetilde w_{3,c},\qquad
F_{3,c}(L,\rho)=[\epsilon^0]\widetilde w_{3,c},\qquad
f_{3,c}(\rho)=F_{3,c}(0,\rho).

```

The symbol $`\mathcal Z_{3,c}`$ denotes a coefficient of **log Z4**, not Z4 itself. The bare epsilon-zero coefficient, the renormalized finite function F, the L-independent shape f and the number at L=0,z=1 are four distinct quantities. Each is labelled separately in the delivered data.

In the adopted convention the simple pole obeys

```math

[\epsilon^{-1}]\widetilde w_{3,c}=-\frac{2k_{3,c}L+h_{3,c}}3,

```

where k is the cusp coefficient and $`\gamma_{4,3,c}=2h_{3,c}`$. The finite RG equation is

```math

\partial_LF_{3,c}-[2b_0F_2+b_1F_1]_c+2k_{3,c}L+h_{3,c}=0.

```

The brackets denote physical color extraction; F1 and F2 here carry their full colors. The pole and finite equations provide strong cross-checks connecting data from different stages.

## 9. Distinguish calculation from validation

The diagram equations, symbolic identities, numerical samples, physical comparisons and formal receipts answer different questions. The later UV scale bridge must be read alongside the historical October 17 non-cusp scale caveat: it supplies a research-level matched local/closed-contour scale convention. The new delivery checks k and h against the adopted external references under that convention; the old limitation is not silently erased or reused as current evidence.

| Evidence | What it checks | What it leaves open |
|---|---|---|
| Source hashes and revision selection | Which exact input was consumed | Mathematical correctness |
| Exact M/N-to-Hlog back-substitution | Equality of the stated representations | Independent original-integral evaluation |
| Independent pole prediction and cusp comparison | UV structure and anomalous dimensions | Arbitrary finite functions of rho |
| Finite RG and crossing | Scale dependence and reflection constraints | An added L-independent even function |
| Independent raw-C/J reconstruction | Finite answer from another route using accepted roots | Upstream operator derivation unless separately checked |
| Same-definition external full-finite comparison | The full function after matching conventions | A formal proof of every step |
| Formal receipt | Its exact checked statement and context | Unchecked ancestors or changed revisions |

The three matter channels have complete computational coefficients through finite, and all pass the cusp, pole and finite RG checks. The new package also tests their explicit representation conversions. Independent raw-input finite reconstructions are available for ff and Ff as preserved historical evidence. Af still lacks that comparable independent reconstruction, and a same-definition external complete finite comparison remains open. Adding an even function $`g(\rho)`$ to F would leave its L derivative and the cusp test unchanged; this explains why the remaining check matters.

The direct AA diagram calculation still lacks double-TGV 15/24 finite coefficients and diagram 16 single-pole/finite coefficients. The 15/24 gap is shared by reflection. Diagram 25's further compact-basis reduction is a separate issue, not a missing selected finite result. The newer total-level QCD/N4 route provides information up to an unknown AA constant; it does not solve the missing individual diagrams. See the [total-reconstruction note](sources/workflow/checks/three-loop-delivery-audit-20261023/total-reconstruction.md).

## 10. Where the software and the delivered files sit

The scientific calculation is accompanied by an evidence path. The solver records derivations in Markdown and may submit mathematical steps through the Wilson MCP API. The kernel presents the mathematical task interface; the proof store handles internal proof construction/replay; interfaces control the statement syntax; Lean checks supported formal statements. Existing CAS results and receipts must match the exact context being reused. This architecture does not imply that every published integral has a complete formal proof.

All-order generators, differential systems, higher-weight function research and numerical samples branch from the same integral definitions. They can supply alternative representations or higher coefficients. They update a selected physical result only through an explicit, checked revision; a new research note alone does not change the adopted Laurent table.

| Question | Main project location | Shareable delivery |
|---|---|---|
| Which graph and ordering? | Frozen `sources/` and `orbit_definitions.json` | Source snapshots and diagram schematics |
| What are the raw/prepared C,J? | Family `sources/`, `prepared/`, manifests | Source pointers in `diagram-ledger.json` |
| What has been integrated? | `workspaces/<id>/derivation.md`, `results/<id>/<revision>.json` | Native and weighted per-diagram coefficient maps |
| Which representation is meant? | Family definition and conversion records | `conversion-dictionary.json` and conventions |
| How are the colors assembled? | Renormalization inputs and source-weight ledger | `coefficients.json`, separate U/bare/CT/post/logZ stages |
| What has passed validation? | Fixed-input records under `workflow/checks/` | Explicit validation scope and reproduction records |
| Which result is currently adopted? | `current.json` / family catalogues | A dated snapshot; this repository is not the live selector database |
| What work remains? | Task ledger and current master plan | Stated snapshot limitations |

The 47 default pure-exchange roots, 149 three-loop orbit records in the broader historical catalogue, and 31 unique orbits in the matter delivery count different objects. They are not conflicting counts of the same graph set. Similarly, a “Done” reporting task or passing project-index check is a management result, not an additional integral proof.

To follow an actual result backwards, start from a field of `coefficients.json`, follow its selected result and source ledger, open the referenced per-diagram coefficient field, then follow that result's raw/prepared source and graph definition. To follow it forwards, reverse that sequence while preserving the declared factors and epsilon budgets. The [source register](source-register.md) identifies the documents used for this guide; scientific data in the dated delivery are unchanged.
