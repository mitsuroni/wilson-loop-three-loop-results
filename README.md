# Three-loop Wilson-loop research data

## P5/P6 publication: 3 November 2026

[The three-loop index](integrals/three-loop/README.md) now covers98 source memberships:39 pure exchanges,15 two-point/exchange insertions and44 tree-three-gluon sources. The new immutable packets provide59 own articles,413 scalar Laurent positions,238 external-colour positions and651 independent function back-substitutions. See [publication validation](validation/three-loop-p5-p6-2026-11-03/README.md) for actual copied-public replays and the limits of reused primitive integration. The remaining16 insertions and35 four-leg source memberships are subsequent curation batches; the complete four-colour scientific result is unchanged.

Research materials for the fundamental lightlike quadrilateral Wilson-loop exponent in massless QCD, prepared by Dengyao for sharing and review.

## Current calculation: one convention and one entrypoint

Start with the [current guide](docs/current/README.md) and [single convention](workflow/current-conventions.md).
The current catalogue covers all163 source orbits. The directly assembled
[four-color three-loop result](results/2026-11-02-current/README.md) is available
through the finite term, with scale/counterterm/RG checks and explicit research
reuse boundaries. You do not choose an old revision or apply a family sign patch.

```text
python -m pip install -r docs/current/requirements.txt
python -m workflow.current_calculation check
python -m workflow.current_calculation result 3
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

The exact public runtime closure has passed an [isolated complete replay](validation/current-2026-11-02/README.md).
Older dated deliveries remain intact. The [39 three-loop pure-exchange evaluation articles](integrals/three-loop/README.md) now join the insertion/tree-three-gluon packets, covering98 source memberships in total.

## Follow the calculation

| Stage | Public entry | This publication |
|---|---|---|
| 1. Generate and classify diagram orbits | [Homemade generator](packages/qgraf-pipeline/README.md) · [guide](docs/upstream/generation.md) | Source, tests, frozen definitions and provenance |
| 2. Convert diagram definitions to parameter integrals C and J | [Montopar](packages/montopar/README.md) · [guide](docs/upstream/parameterization.md) | Source, tests, integral inputs and versioned contracts |
| 3. Evaluate individual or joint integrals | [One-loop articles](integrals/one-loop/README.md) · [two-loop articles](integrals/two-loop/README.md) · [three-loop 6e](integrals/three-loop/README.md) | English evaluation packets for 2+12 low-loop orbits and98 three-loop source memberships |
| 4. Assemble color channels | [Current four-color calculation](results/2026-11-02-current/README.md) · [low-loop example](results/2026-11-02-low-loops/README.md) | Direct native-current low/matter/AA assembly and exact current results |
| 5. Renormalize, extract finite terms and verify RG | [Low-loop derivation](results/2026-11-02-low-loops/derivation.md) · [checks](results/2026-11-02-low-loops/validation/README.md) | Counterterms, logZ, finite results, RG and complete external finite comparison |

Start with the [complete one/two-loop example](results/2026-11-02-low-loops/README.md) to follow the full calculation, or the [upstream guide](docs/upstream/README.md) to inspect its source packages. A source package, a parameter-integral definition, an evaluated Laurent coefficient, and an independent validation are distinct deliverables.

## Three-loop exchange batch: 3 November 2026

All39 pure-exchange 6e orbits now have [individual English derivations and exact data](integrals/three-loop/README.md), sharing six mathematical methods and one [portable evaluator](evaluations/three-loop-6e/README.md). The packet retains the original ordered roots and full-colour source bindings, all273 Laurent positions and function-family conversions. Ten Gamma calculations are fresh; complex routes explicitly reuse checked original raw CAS output before new assembly. [Actual validation](validation/three-loop-6e-2026-11-03/README.md) separates exact comparison, original-domain controls and copied-closure reproduction. No new all39 independent CAS or formal proof is claimed.

## Complete low-loop pilot: 2 November 2026

All fourteen one/two-loop orbits now have English derivations, exact data and portable calculation entrypoints. The tables extend to one-loop epsilon3 and two-loop epsilon2. Their full color/scales/counterterm assembly reproduces the finite functions and independent external finite result, with explicit RG checks and actual isolated reproduction. See the [data-flow table](results/2026-11-02-low-loops/README.md), [reproduction commands](docs/reproduction.md) and [lessons for three-loop publication](docs/low-loops/curation-lessons.md). Original integrations and reused HLOG identities have separate evidence scopes.

## First upstream snapshot: 2 November 2026

- [Source packages](packages/README.md): the self-authored combinatorial diagram generator and Montopar. Source provenance, publication adaptations, dependencies and actual tests are recorded with each package.
- [163-integral English index](docs/upstream/catalogue.md) and [versioned inputs](data/README.md): frozen diagram definitions and raw C/J packages for one, two and three loops, plus available ordered/prepared inputs. Each integral has an English entry page and exact source pointers.
- [Reproduction levels](docs/reproduction.md): inspect the frozen inputs, run the source packages, then move to evaluated results. Generating C/J does not evaluate J.
- [Actual upstream checks](validation/upstream-2026-11-02/README.md) and the [historical full-tree manifest at its original commit](https://github.com/mitsuroni/wilson-loop-three-loop-results/blob/ed854f02b0e79840f96c7a73c2642c21fef07307/UPSTREAM-PUBLICATION.json). It binds that snapshot, not this evolving navigation.
- [Roadmap and publication status](docs/publication-status.md): what is already public, what is being curated next, and how historical snapshots should be read.

The calculation workspace reached a complete four-color three-loop Laurent-through-finite and RG research milestone on 2 November. The upstream foundation and complete low-loop pilot are now public. The current four-color data and direct reconstruction are now public. The39 pure-exchange three-loop articles are now public; the fifteen two-point/exchange insertions and44tree5e packets are also public; the remaining16 insertions and35 four-leg source memberships are the next curation stages. Older three-color files retain their dated scope.

## Existing result archives

- [23 October matter-channel summary](docs/results/index.md).
- [Exact three-channel coefficients](results/2026-10-23/coefficients.json), [31-orbit ledger](results/2026-10-23/diagram-ledger.json), and [M/N–Hlog dictionary](results/2026-10-23/conversion-dictionary.json).
- [Illustrated end-to-end data flow](docs/data-flow/README.md) · [中文说明](docs/data-flow/reading-guide-zh.md).
- [Original downloadable archives](https://github.com/mitsuroni/wilson-loop-three-loop-results/releases).

The 23 October standalone result replay remains:

```text
python results/2026-10-23/reproduce/reproduce.py --output-dir replay-output
```

Its original dependencies and validation boundaries are preserved in that dated package. The 18 October archive preserves the wider historical catalogue.

## Provenance and reuse

Use the source manifests and stable orbit IDs when comparing versions. Frozen mathematical inputs retain their original bytes, fields, domains and continuation conventions. Current package source snapshots may differ from historical producer versions; each relationship is reported explicitly. The generator's native path is a homemade combinatorial implementation; optional QGRAF interfaces do not make the catalogue a newly rerun native-QGRAF calculation.

This is a research delivery repository. It contains no Wilson proof database, CAS installation, downloaded third-party source tree, credentials, or personal execution environment. Original third-party notices in historical assets remain intact; no new license is assigned to those assets.
