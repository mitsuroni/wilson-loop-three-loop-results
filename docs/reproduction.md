# Reproduction levels

## Current calculation

The [current guide](current/README.md) is the normal reader and assembly route.
From the repository root:

```text
python -m workflow.current_calculation check
python -m workflow.current_calculation integral 2e-Pi1-01
python -m workflow.current_calculation result 3
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

The complete isolated public-copy replay is [recorded here](../validation/current-2026-11-02/README.md).
The levels below retain the dated source and low-loop research workflows and
their precise integration/reuse scopes.


## Reproduce the39 three-loop pure-exchange integrals

The [English index](../integrals/three-loop/README.md) connects each original source to its seven Laurent coefficients and shared mathematical route. From the repository root, using the versions declared by the packet:

```text
python evaluations/three-loop-6e/scripts/verify_manifest.py
python evaluations/three-loop-6e/scripts/evaluate_batch.py --output replay-output/three-loop-6e-new
python evaluations/three-loop-6e/scripts/verify_batch.py --results replay-output/three-loop-6e-new/generated.json --out replay-output/three-loop-6e-new/verification.json
python evaluations/three-loop-6e/scripts/independent_controls.py --input evaluations/three-loop-6e/engine/benchmarks/2ne/prepared/6e.json --results replay-output/three-loop-6e-new/generated.json --out replay-output/three-loop-6e-new/independent-controls.json
```

No installed Maple/HyperInt or private Wilson workspace is needed. The default command computes ten Gamma routes afresh and reuses exact-profile archived original CAS outputs for the others. It then performs new physical assembly; no final coefficient reference is read during generation. The checker opens final references afterwards. See [actual public-copy validation](../validation/three-loop-6e-2026-11-03/README.md) and the [independent-control coverage](../evaluations/three-loop-6e/methods/independent-controls.md). This does not freshly reintegrate all39 complex sources or provide formal proofs.

## Inspect and verify input data

Use [data/README.md](../data/README.md) to find an orbit's definition, C/J package, available prepared representation, stable ID and source hash. Run the [upstream package audit](../tools/check_upstream.py) from the repository root:

```text
python tools/check_upstream.py --output-dir upstream-audit
```

This checks the published manifests, catalogue pointers, preserved bytes and navigation. It does not evaluate divergent parameter integrals or prove their original physics.

## Run the source packages

The [generator README](../packages/qgraf-pipeline/README.md) and [Montopar README](../packages/montopar/README.md) give exact entrypoints, dependencies and new-output-directory examples. Use the documented non-PDF path if TeX is unavailable. Montopar color-backend requirements depend on the family.

Current-source replays and frozen historical inputs have separate version identities. Compare only the fields and mathematical conventions stated by the replay record. Do not overwrite the checked-in data to make a comparison pass.

## Reassemble evaluated results

The [complete low-loop packet](../results/2026-11-02-low-loops/README.md) rebuilds one/two-loop color sums, scales, counterterms, logZ, finite functions and RG/external finite checks. Run from the public root, with the recorded dependencies:

```text
python tools/check_low_loops.py --output-dir replay-output/audit-new
python results/2026-11-02-low-loops/reproduce.py --output-dir replay-output/totals-new
python results/2026-11-02-low-loops/numeric.py --output-dir replay-output/finite-numeric-new
```

The [23 October result replay](../results/2026-10-23/reproduce.md) still rebuilds its dated three-matter-color expressions from evaluated coefficients. It does not rerun all original integrations. The complete three-loop four-color delivery is already available through [the current guide](current/README.md); its direct assembly and the original-integral evaluation routes above have separate scopes.

## Recompute an integral

The [one-loop](../integrals/one-loop/README.md) and [two-loop](../integrals/two-loop/README.md) articles distinguish actual original/reduced-integral computation from representation verification and reused identities. To execute their declared routes before assembling:

```text
python results/2026-11-02-low-loops/reproduce.py --output-dir replay-output/full-new --regenerate-integrals
python integrals/two-loop/check_integrals.py --output-dir replay-output/kernel-controls-new
```

One-loop primitives and Gamma/M-moment series are freshly computed; the two-loop shared engines expand the documented mother functions and compare complete HLOG coordinates. The original TGV04 HyperInt integration and fixed dictionaries are explicitly reused. No new full CAS integration of fourteen graphs is claimed. Historical CAS source scripts are archival; they retain original installation paths and must not be run as portable output writers.
