# The current calculation

## Continuous derivations and individual audit

The [current integral reading layer](../../integrals/README.md) now provides continuous source-to-Laurent derivations for2 one-loop,12 two-loop and39 three-loop exchange entries. The [fixed standard and reusable skill](../../workflow/integral-documentation-standard.md) separate exposition, implementation, replay, mathematical validation and publication. The [39-entry audit](../../workflow/checks/integral-documentation-6e-20261103/review/README.md) accepts10 direct derivations and29 precisely adopted evaluated lemmas. This is not39 new independent CAS integrations. Statements saying "not yet published" inside frozen reports describe their construction-time status; this release publishes the reviewed reading layer.

## P5/P6 publication: 3 November 2026

[The three-loop index](../../integrals/three-loop/README.md) now covers98 source memberships:39 pure exchanges,15 two-point/exchange insertions and44 tree-three-gluon sources. The new immutable packets provide59 own articles,413 scalar Laurent positions,238 external-colour positions and651 independent function back-substitutions. See [publication validation](../../validation/three-loop-p5-p6-2026-11-03/README.md) for actual copied-public replays and the limits of reused primitive integration. The remaining16 insertions and35 four-leg source memberships are subsequent curation batches; the complete four-colour scientific result is unchanged.

Use one entrypoint for the current physical calculation. You do not select a
historical revision or apply a family-dependent sign when reading an integral
or assembling a color channel.

From the repository root, with Python 3.12 and SymPy 1.14:

```text
python -m pip install -r docs/current/requirements.txt
python -m workflow.current_calculation check
python -m workflow.current_calculation integral 2e-Pi1-01
python -m workflow.current_calculation integral 4e-Pi1-01 --definition
python -m workflow.current_calculation result 1
python -m workflow.current_calculation result 2
python -m workflow.current_calculation result 3
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

The replay directory must be new. Inside the checkout, use `replay-output` or
`.wilson`; a new directory outside the checkout is also accepted. The reader
checks its exact data/code/source bindings before use. Missing files, mixed
contracts, extra family signs, duplicate Euler factors, and incorrect orbit
references fail explicitly. No private publication checkout, CAS installation
or Wilson proof database is required for these commands.

## Follow the definitions to the result

| Object | Entry |
| --- | --- |
| Physical sign and normalization contract | [Single convention](../../workflow/current-conventions.md) |
| 163 source orbits and selected integral identities | [Current catalogue](../../benchmarks/current/catalog.json) |
| 13 native full-color prepared C/J packages | [Prepared inputs](../../benchmarks/current/prepared/) |
| Canonical low-loop I tables through positive epsilon orders | [Low-loop tables](../../benchmarks/current/low-loops/integrals/) |
| One-loop and two-loop full assembled series | [One loop](../../benchmarks/current/results/1loop.json), [two loops](../../benchmarks/current/results/2loop.json) |
| Complete four-color three-loop series | [Current result](../../benchmarks/current/results/3loop.json), [reading guide](../../results/2026-11-02-current/README.md) |
| Exact runtime bindings | [Current manifest](../../workflow/current-calculation.json) |
| Actual public-copy reproduction | [Validation](../../validation/current-2026-11-02/README.md) |

There are 2 one-loop, 12 two-loop and 149 three-loop source orbits. The current
catalogue identifies 151 result objects and 12 algebraic zeros; this is not a
claim that 151 separate finite integrals were newly proved. In particular, the
35 four-leg source orbits are consumed through nine joint objects, with their
reflection maps and zero terms retained. They must not be added a second time.

Native prepared C/J preserves the source's declared cube or ordered domain,
measures, maps and i0. It adds Euler compensation once and performs no legacy
sign conversion. Full-color preparation, fixed-basis a1 projection and the
connected physical I are separate stages. The original 47 ordered/a1 solver
roots retain their identity; this publication does not replace their proofs
with a new full-color preparation claim.

## What the reader returns

`integral ID --definition` reads native prepared C/J. `integral ID` returns
explicit whole-I coefficient fields, with actual available powers and channel
identity. A J series, IBP contribution or known partial sum is never silently
presented as the complete I. Missing epsilon orders are not filled with zero.
Gamma channels remain distinct; Pi2 shared geometry still requires its external
operator/color factors. Use `--provenance` only when inspecting the original
source/version bindings.

The low-loop canonical tables contain one-loop coefficients through epsilon3
and two-loop coefficients through epsilon2. Their 96 scalar positions were
checked exactly during normalization. The extra orders are actual evaluated
input data. Three-loop finite counterterms consume one loop through epsilon2
and two loops through epsilon1.

## What the replay computes

The replay has four stages: low-loop assembly, three matter colors, AA, and
the combined scale/counterterm/RG calculation. Af is now assembled directly
from the 31 current native I objects in a common z-HLOG representation. It does
not use an old Af total to construct the result. AA directly reads 83 tree
orbits, 31 insertion orbits and the four-leg joint objects. Nonstandard letters
and periods in the three relevant 5e diagrams are retained until their exact
cancellation is checked. AA's seven MN/zeta coefficients are independently
obtained by exact rational basis conversion and back substitution.

Existing totals are opened after construction as holdouts or exact compact
function-representation targets. The original integral evaluations, dictionaries
and admitted I2344/I1234 research results are reused at their recorded scope.
This is a direct coefficient reconstruction and validation workflow, not a new
integration of every original high-dimensional kernel. `formal=not_run`.

## Dated deliveries and next work

The [upstream](../upstream/README.md), [complete low-loop example](../../results/2026-11-02-low-loops/README.md)
and older matter-channel archives retain their exact dated data and validation
scope. The current reader and assembly commands above are the normal entrypoints.
The [39 three-loop pure-exchange evaluation articles](../../integrals/three-loop/README.md) are now public, with a standalone six-route evaluator and original source/result bindings. Other three-loop families follow. The current four-color numerical/symbolic result and its runtime manifest remain unchanged.
