# Tree three-gluon sources: momentum, parameters and evaluated functions

The P6 packet connects44 labelled three-loop sources to their complete physical Laurent functions through epsilon^0. They share one three-gluon vertex, an additional bare exchange, and five ordered Wilson endpoints. Eight evaluation families share actual calculation code; every source keeps its own graph, phases, C/J, domains, orders and result. This is a reproducible curation of established primitive evaluations with new exact reassembly, not44 fresh independent integrations.

## From the momentum source to the accepted parameter root

The [physical convention](../snapshots/benchmarks/theory/montopar/shared-sign-conventions.md) fixes the Minkowski metric, Fourier exp(-ipx), Wilson equation U'=igAU, propagators, three-gluon vertex momenta and Wick phases. A source is specified by its Wilson slots, the incidence of four gluon propagators and the cubic vertex, plus the Fourier phase of each loop momentum. Those full components and ordered endpoint bounds are embedded in every source.json, with byte-bound original graph and raw/native records in input.json. They specify the momentum object without substituting a drawing or an integral label for it. Momentum algebra and the upstream Montopar export are accepted input facts in this packet; they are audited, not rerun as a new FeynCalc calculation.

The three propagators meeting at the vertex give the alpha simplex. The stored polynomial S is the quadratic form after its coordinate-centre integration, and P is the contracted derivative numerator. Thus the TGV part is P*(alpha1*alpha2*alpha3)^(-epsilon)/S^(2-2epsilon). The fourth propagator supplies the separate pure-exchange polynomial and tangent contraction. The declared endpoint mapping and Jacobian give the native cube J; J_raw_s retains the corresponding original dependent bounds and explicit outer-to-inner order. Both roots remain present. A polynomial source probe tests this map; it does not establish a finite Laurent coefficient at an integer regulator.

The alpha measure means d alpha1 d alpha2 over alpha1>=0, alpha2>=0, alpha1+alpha2<=1 with alpha3=1-alpha1-alpha2. It is not three independent unit intervals. Subsequent projective changes in the family notes retain their Jacobians and endpoint terms. Each family first works in its stated convergent regulator region and then takes the local meromorphic continuation to epsilon=0.

## C, path ordering and the physical I

Current native C uses the [one current convention](../snapshots/workflow/current-conventions.md), with assembly_sign=1 and exactly one exp(-3 EulerGamma epsilon). The Wilson path order puts later matrices on the left; independent finite-N Fierz contraction of that order gives the native colour phase. The historical unreversed trace has the opposite phase and must be interpreted with its own vertex convention. It is not used to add a sign patch to the native root.

The numerical phase, coupling and Fourier constants multiply as

```text
i * 2*(1-2eps) * (i*CF/4) * 4^3/(4*64) * (-1/2)
    = CF*(1-2eps)/16.
```

The full ordered colour basis uses CA^2=1+2CA CF, hence CF=CA^2 CF-2CA CF^2 before selecting a1. For this tree family, sym=4 and the established crossing extraction yield

```text
G(eps) = exp(-3 EulerGamma eps) Gamma(1-eps) Gamma(1-2eps)
I_g = (1-2eps) G(eps)/8
      * [J_g(z,eps)+z^(3eps) J_g(1/z,eps)]/[1+z^(3eps)].
```

The Euler factor here rewrites the one already in C. The physical I includes these symmetry/crossing operations; ordered J, an IBP block, or a product of two already-normalized I values cannot silently replace it. All44 entries bind this bridge to their own root.

## The eight actual routes

| Evaluation family | Sources | Evaluated quantities used by the reassembly |
| --- | ---: | --- |
| T1 empty |7| Gamma/Beta, F through5, Z through4 |
| T2 empty |7| full-source reflection of T1, fixed z |
| T3 empty |6| endpoint/IBP densities E9,E10,E26 through4 |
| T1 opposite |3| ET and EK through4 |
| T2 opposite |3| full-source reflection of T1 opposite |
| T3 opposite |5| E22/E23 through2 and the separate M product for44 |
| T4 empty |8| raw Q from absolute orders-2 through2, fresh inversion and multiplier through4 |
| T4 opposite |5| integrated global-IBP hard blocks27/28; source reflection29; own factorized41/42 |

[The31 T1/T2/T3 derivation](t123.md) and [the13 T4 derivation](t4.md) give the actual subintegral definitions, domains, boundary terms and order budgets. Original Maple scripts, raw output and logs are retained under snapshots. Python does not launch them. In particular the27/28 raw I text is already an integrated hard result; its inclusion is historical reuse, even though the final selected JSON is excluded from generation. The old T4 relative-order writer is kept as historical code and is not executed; the six adopted r002 corrections are reproduced from original Q instead.

## What a completed function means

Hlog(z,[a1,...,an]) uses kernels dt/(t-a). The empty word is1. The ordinary alphabet is0,-1; the opposite T4 hard finite coefficients additionally retain-1/2,-2. The complete common period order is zeta2,zeta3,zeta5,log2,A4,A6. Standard three-period records embed by adding three zero exponents. Even powers of pi are converted exactly using pi^2=6 zeta2.

M_kj are the original no-factorial moments and N_kj=M_k,j+1-rho M_kj, rho=log z. E and tau labels in a displayed coefficient denote their complete evaluated dictionary, never a new unresolved integral. Every named coefficient is [back-substituted independently](representations.md). Existing evaluated dictionaries remain admitted inputs. Omega6 is separately defined and retained in the T4 note and generated27/28/29 entries. Its cancellation in 2c27,0+c28,0 is checked after preserving the individual coefficients.

All seven positions epsilon^-6 through epsilon^0 are stored, including zero. Fixed z>0, real principal logs and local meromorphic epsilon continuation are the active scope. Archived boundary/complex/numerical evidence keeps its original scope; this batch does not establish a new z=-1 or complex i0 bridge.

## Reproduction and evidence boundaries

Generation parses the admitted raw primitives and performs exact Gamma/Beta, inversion, reflection and shuffle assembly. It is audited to deny final result tables, aggregate coefficient lookups, external processes and network. The result is fixed before final reference export. An actual copied closure then generates byte-identical results before the final-reference directory is added for checking. Archived files may contain historical formulas; this is a generator access separation, not a blinded human calculation or an OS sandbox.

[Independent controls](independent-controls.md) distinguish source identities, finite-I controls and representation checks. They do not confer a new independent original-domain finite integration on every source. Passing regression/hash/management checks is distinct from a proof of the original integral. Formal=not_run; RG is not rerun for this unchanged source selection.
