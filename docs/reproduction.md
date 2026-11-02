# Reproduction levels

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

The existing [23 October result replay](../results/2026-10-23/reproduce.md) rebuilds its three-matter-color expressions, scales, counterterms and checks from the already evaluated coefficients. It does not rerun all original integrations. The new complete four-color delivery will be curated in the later downstream publication step.

## Recompute an integral

Each future [evaluation article](../integrals/README.md) will distinguish actual original-integral computation from result verification and reused proved identities. It will list the scripts, source data, CAS versions, branch prescriptions, checkpoints and resource requirements for its actual route. No all-integral one-command recomputation is claimed by this first upstream snapshot.
