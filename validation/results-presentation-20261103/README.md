# Four-color formula and presentation validation

The [result article](../../results/current/four-color-results.md) is generated from the exact one/two/three-loop files under `benchmarks/current/results/`.

The [formula report](formula-check.json) records96 displayed exact expressions and213 exact algebraic checks: conversion of pi powers to zeta values, finite slices at L=0, full scale restoration, counterterm convolution from the actual lower-loop series, pole/logZ separation and the four-color finite RG identity. The renderer also compares the complete generated article with the checked-in Markdown.

```text
python -B tools/render_four_color_results.py --check
```

The underlying integral evaluations, color assembly and public-copy execution are recorded in the [validation index](../README.md). This presentation update preserves all exact calculation data and original mathematical records.
