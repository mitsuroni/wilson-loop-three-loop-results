# Integral input data

Use the [current163-source catalogue](../benchmarks/current/catalog.json) and [native prepared C/J packages](../benchmarks/current/prepared/) for the current calculation. The [exact evaluated results](../benchmarks/current/results/) and [four-color formula article](../results/current/four-color-results.md) are their downstream result entrypoints.

| Data layer | Entry |
|---|---|
| Current source identities and selected integral objects | [Current catalogue](../benchmarks/current/catalog.json) |
| Current native full-color prepared C/J | [Prepared inputs](../benchmarks/current/prepared/) |
| Original graph definitions and unmodified producer outputs | [Frozen input catalogue](2026-11-02/README.md) |
| Per-source original graph/raw C/J bindings | [Original machine catalogue](2026-11-02/catalog.json) |
| Source versions and exact field/file hashes | [Original source manifest](2026-11-02/source-manifest.json) |
| Diagram generation and parameterization implementation | [Generator and Montopar](../packages/README.md) |

There are2 one-loop,12 two-loop and149 three-loop source orbits. Original snapshot entries live under `2026-11-02/integrals/<N>loop/<orbit_id>/`; each input index points to the original graph and full C/J record. The current catalogue adds the selected whole-integral identity and its consumption in the assembled result.

The [integral reading layer](../integrals/README.md) explains evaluation from C/J to Laurent coefficients. The [current guide](../docs/current/README.md) explains how to read definitions and evaluated results with the same interface. Dated input snapshots and their original checks are listed under [historical records](../docs/archive/README.md).
