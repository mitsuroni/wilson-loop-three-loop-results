# Reproduce the results

Use Python3.12 and the declared dependencies. Run commands from the repository root, choosing a new output directory for each calculation.

## Complete current result

```text
python -m pip install -r docs/current/requirements.txt
python -m workflow.current_calculation check
python -m workflow.current_calculation result 3
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

The reader checks its exact runtime bindings. Replay reconstructs low-loop color sums, all four three-loop channels, scale restoration, coupling counterterms, poles, finite functions and RG relations. The output matches the [current analytic data](../benchmarks/current/results/3loop.json). Read the [formula article](../results/current/four-color-results.md) for the mathematical quantities returned.

## Regenerate the result article

```text
python -B tools/render_four_color_results.py --check
```

The renderer checks its96 displayed exact expressions and213 algebraic identities against the current result data. These include the finite slice at L=0, scale restoration, lower-loop counterterm convolution, pole separation and the finite RG relation. To regenerate the article and its check report, run the same command without `--check`.

## Individual integral batches

The per-entry code adapters use shared same-loop or same-family evaluators. Their derivation and data indices identify the actual source integral and evaluated subintegrals.

```text
python -B integrals/one-loop/2e-01/code/reproduce.py --output-dir replay-output/one-loop-new
python -B integrals/two-loop/4e-01/code/reproduce.py --output-dir replay-output/two-loop-new
python -B integrals/three-loop/6e-01/code/reproduce.py --output replay-output/6e-new
```

For the6e comparison and controls after generation:

```text
python -B integrals/families/three-loop-6e/runtime/scripts/verify_batch.py --results replay-output/6e-new/generated.json --out replay-output/6e-new/verification.json
python -B integrals/families/three-loop-6e/runtime/scripts/independent_controls.py --input integrals/families/three-loop-6e/runtime/engine/benchmarks/2ne/prepared/6e.json --results replay-output/6e-new/generated.json --out replay-output/6e-new/controls.json
```

The original CAS integrations, their scripts and raw outputs accompany each family. The supported replay checks and reuses those evaluated objects before physical assembly. The [validation index](../validation/README.md) identifies exact checks, numerical controls and the associated mathematical audits.

Additional family-specific commands are in the [15 insertion packet](../evaluations/three-loop-insertions/README.md), [44 tree-three-gluon packet](../evaluations/three-loop-tree5e/README.md), and [upstream generation guide](upstream/README.md). Original dependency versions and execution conditions are retained with each package.

## Historical reproduction

[Dated releases](archive/README.md) keep the commands for their original inputs and result snapshots. Use the current commands above to reproduce the present result.
