# Current input and result data

Use the [current calculation guide](../../docs/current/README.md) for the only
normal reader/replay entrypoint and the [single convention](../../workflow/current-conventions.md).

The [catalogue](catalog.json) connects all163 source orbits to native prepared
C/J and adopted whole-integral fields. [Prepared packages](prepared/) preserve
the declared domain representation and add Euler once, without sign conversion.
The [low-loop tables](low-loops/) include the actual positive-order input budget.

Read the [one-loop](results/1loop.json), [two-loop](results/2loop.json), and
[four-color three-loop](results/3loop.json) results with:

```text
python -m workflow.current_calculation result 3
```

The [runtime manifest](../../workflow/current-calculation.json) binds these files
and their exact sources. [Migration provenance](migration.json) identifies the
unchanged scalar data behind the canonical low-loop tables. It is a one-time
data/notation migration record, not a runtime sign patch. Current reading does
not infer missing coefficients, confuse J/IBP fragments with I, or require a
historical revision choice. The151 result identities and12 algebraic zeros are
not a count of newly integrated/proved finite objects.

The [publication validation](../../validation/current-2026-11-02/README.md) records
actual copied-source runs. Original integration and function-dictionary reuse
remain explicit; formal=not_run.
