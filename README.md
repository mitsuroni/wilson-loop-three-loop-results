# Three-loop Wilson-loop research data

Research materials for the fundamental lightlike quadrilateral Wilson-loop exponent in massless QCD, prepared by Dengyao for sharing and review.

## Follow the calculation

| Stage | Public entry | This publication |
|---|---|---|
| 1. Generate and classify diagram orbits | [Homemade generator](packages/qgraf-pipeline/README.md) · [guide](docs/upstream/generation.md) | Source, tests, frozen definitions and provenance |
| 2. Convert diagram definitions to parameter integrals C and J | [Montopar](packages/montopar/README.md) · [guide](docs/upstream/parameterization.md) | Source, tests, integral inputs and versioned contracts |
| 3. Evaluate individual or joint integrals | [Integral input catalogue](data/README.md) · [derivation structure](integrals/README.md) | Input records and English entry pages; full evaluation articles will be curated next |
| 4. Assemble color channels | [Result archive](results/README.md) | The dated three-matter-color delivery remains available |
| 5. Renormalize, extract finite terms and verify RG | [Reporting conventions](docs/results/conventions.md) · [validation scope](docs/results/validation-scope.md) | Existing dated evidence remains in its original scope |

Start with the [upstream publication guide](docs/upstream/README.md), then open one integral's input record. The release is being organized in data-flow order. A source package, a parameter-integral definition, an evaluated Laurent coefficient, and an independent validation are distinct deliverables.

## First upstream snapshot: 2 November 2026

- [Source packages](packages/README.md): the self-authored combinatorial diagram generator and Montopar. Source provenance, publication adaptations, dependencies and actual tests are recorded with each package.
- [163-integral English index](docs/upstream/catalogue.md) and [versioned inputs](data/README.md): frozen diagram definitions and raw C/J packages for one, two and three loops, plus available ordered/prepared inputs. Each integral has an English entry page and exact source pointers.
- [Reproduction levels](docs/reproduction.md): inspect the frozen inputs, run the source packages, then move to evaluated results. Generating C/J does not evaluate J.
- [Actual upstream checks](validation/upstream-2026-11-02/README.md) and the [publication manifest](UPSTREAM-PUBLICATION.json).
- [Roadmap and publication status](docs/publication-status.md): what is already public, what is being curated next, and how historical snapshots should be read.

The calculation workspace reached a complete four-color Laurent-through-finite and RG research milestone on 2 November. This repository update deliberately publishes the upstream packages and input data first. The older result files below remain a dated three-color snapshot until the new downstream release is curated; they are not relabeled as a four-color data release.

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
