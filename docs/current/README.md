# Current calculation and data

**[Analytic four-color results](../../results/current/four-color-results.md)** · **[Exact three-loop JSON](../../benchmarks/current/results/3loop.json)** · **[Source catalogue](../../benchmarks/current/catalog.json)** · **[Reproduction](../reproduction.md)**

The current interface connects163 source orbits to prepared C/J, selected integral objects, complete color sums, coupling counterterms and finite/RG results. It uses the [single physical convention](../../workflow/current-conventions.md).

| You need | Entry |
|---|---|
| Four-color Laurent coefficients and finite functions | [Result article](../../results/current/four-color-results.md) |
| Exact one/two/three-loop result files | [Current result data](../../benchmarks/current/results/) |
| Graph identity and selected integral | [Current catalogue](../../benchmarks/current/catalog.json) |
| Full-color native parameter inputs | [Prepared packages](../../benchmarks/current/prepared/) |
| Low-loop integral coefficients through positive orders | [Low-loop tables](../../benchmarks/current/low-loops/) |
| Individual mathematical derivation | [Integral reading index](../../integrals/README.md) |
| Reproduction and validation | [Commands](../reproduction.md), [check index](../../validation/README.md) |

## Read an integral or result

```text
python -m pip install -r docs/current/requirements.txt
python -m workflow.current_calculation check
python -m workflow.current_calculation integral 2e-Pi1-01 --definition
python -m workflow.current_calculation integral 2e-Pi1-01
python -m workflow.current_calculation result 1
python -m workflow.current_calculation result 2
python -m workflow.current_calculation result 3
```

`--definition` returns the source's native prepared C/J with its actual domain and normalization. The integral command returns its evaluated coefficient fields and available epsilon powers. Result commands return the complete connected exponent data.

There are2 one-loop,12 two-loop and149 three-loop source orbits. The35 four-leg sources enter through nine joint objects and their explicit source-consumption ledger. One-loop tables extend through epsilon3; two-loop tables extend through epsilon2. The three-loop file contains the four complete Laurent series through epsilon0.

## Reconstruct the calculation

```text
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

The output directory must be new. Replay assembles the low-loop results, three matter channels, AA, and then the combined scale/counterterm/RG result. It reads the exact included source snapshot and evaluated integral data. [Recorded public-copy validation](../../validation/current-2026-11-02/README.md) documents the execution; [integral audits](../../validation/README.md#integral-derivations-and-public-replay) document the individual mathematics and function evaluations.

The [runtime manifest](../../workflow/current-calculation.json) binds the calculation inputs and code. [Historical releases](../archive/README.md) preserve the preceding versions and their own commands.
