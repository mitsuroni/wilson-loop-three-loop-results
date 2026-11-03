# Calculation validation and mathematical audits

The checks below connect the [current analytic formulas](../results/current/four-color-results.md), exact data, source assembly and integral derivations.

## Current four-color result

| Check | Record |
|---|---|
| Displayed formulas, finite slices, scale factors, lower-loop counterterms and finite RG | [Result-presentation validation](results-presentation-20261103/README.md) |
| Complete current calculation in an isolated public-data closure | [Current calculation validation](current-2026-11-02/README.md) |
| Exact result/source/code bindings | [Current runtime manifest](../workflow/current-calculation.json) |

The current calculation validation covers low-loop assembly, all three matter colors, AA, the complete counterterm and scale relations, and the combined RG result. The new presentation check verifies the formulas actually displayed in the result article against those exact data.

## Integral derivations and public replay

| Material | Record |
|---|---|
| All53 current reading entries,369 coefficient positions, actual public checkout | [Continuous-derivations publication replay](continuous-derivations-2026-11-03/README.md) |
| One/two-loop14-entry scientific reading audit | [Per-entry review](../workflow/checks/integral-documentation-low-loops-20261103/review/README.md) |
| Three-loop39-entry scientific reading audit | [Per-entry review](../workflow/checks/integral-documentation-6e-20261103/review/README.md) |
| Insertion and tree-three-gluon family calculation packages | [P5/P6 public-copy checks](three-loop-p5-p6-2026-11-03/README.md) |
| Original graph/C/J publication | [Upstream-data checks](upstream-2026-11-02/README.md) |

Each record identifies the mathematical objects, actual operations and shared inputs it checks. Original integral evaluation, exact expression comparison, numerical control, source reconstruction and RG provide complementary evidence. The per-entry audits state precisely which evaluated subintegrals are adopted and how they enter the complete coefficients.

[Historical releases](../docs/archive/README.md) retain their original validation records.
