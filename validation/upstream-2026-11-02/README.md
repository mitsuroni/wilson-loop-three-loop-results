# Upstream publication checks

2 November 2026, WILSON-5.20. This evidence concerns source packages and input data, before integral evaluation.

- [Summary and exact coverage](summary.json).
- [Generator regression/portability](../../packages/qgraf-pipeline/VALIDATION.json) and [Montopar regression/portability](../../packages/montopar/VALIDATION.json).
- [Non-PDF adapter check](generator-adapter.json): default-preserving behavior and mathematical exports unchanged.
- [Offline source build and installation](offline-install.json): both packages built and installed into an isolated target; no global installation or network dependency download.
- [Representative CLI regeneration](representative-replay.json): two-loop12 and three-loop149 definitions and the native two-orbit Pi1 parameter package equal frozen products byte for byte; no PDF/CAS runtime used.
- [Independent source/data/navigation audit](independent-audit.json).
- [Frozen data export checks](data-export-checks.json).

Actual run records can contain original local paths as provenance; portable entrypoints are given in the package READMEs. Test counts overlap and are not a count of independent physical proofs. Input hashes fix versions; graph generation and C/J export do not evaluate the original parameter integrals. Optional external color/PDF backends and full Laurent evaluation are outside this package check.

For a new offline integrity audit from a fresh checkout:

```text
python tools/check_upstream.py --output-dir upstream-audit
python tools/build_catalogue_view.py --check
```

These commands do not modify frozen inputs. Original source snapshots and older dated result deliveries keep their recorded scientific scope.
