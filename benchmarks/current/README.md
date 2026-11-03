# Current source and exact-result data

**[Four-color analytic formulas](../../results/current/four-color-results.md)** · **[Three-loop exact JSON](results/3loop.json)** · **[Source catalogue](catalog.json)** · **[Current commands](../../docs/current/README.md)**

| Object | Data |
|---|---|
| All163 source orbits and selected integral identities | [catalog.json](catalog.json) |
| Native full-color C/J definitions and domains | [prepared/](prepared/) |
| Low-loop I tables through positive epsilon orders | [low-loops/](low-loops/) |
| Complete assembled one/two/three-loop results | [1loop.json](results/1loop.json), [2loop.json](results/2loop.json), [3loop.json](results/3loop.json) |
| Exact runtime data/source/code bindings | [Calculation manifest](../../workflow/current-calculation.json) |

The result fields and all displayed analytic formulas are explained in [the result article](../../results/current/four-color-results.md). Use the [single convention](../../workflow/current-conventions.md) for symbols and normalization, [individual derivations](../../integrals/README.md) for source evaluations, and [validation](../../validation/README.md) for assembly and RG checks.

```text
python -m workflow.current_calculation integral 2e-Pi1-01 --definition
python -m workflow.current_calculation integral 2e-Pi1-01
python -m workflow.current_calculation result 3
```

[Migration provenance](migration.json) records the exact connection of the canonical low-loop tables to their original data. Historical snapshots and release notes remain under [the archive index](../../docs/archive/README.md).
