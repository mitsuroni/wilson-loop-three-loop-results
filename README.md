# Three-loop Wilson-loop research results

Research data for the fundamental lightlike quadrilateral Wilson-loop exponent in massless QCD. Prepared by Dengyao for sharing and review.

## Start here

- [Illustrated end-to-end data flow](docs/data-flow/README.md) · [中文说明](docs/data-flow/reading-guide-zh.md)

- [Research summary and reading order](docs/results/index.md)
- [Color-channel expressions](results/2026-10-23/coefficients.json): bare sums, counterterms, poles of log Z4, finite functions and anomalous dimensions.
- [Individual diagram data](results/2026-10-23/diagram-ledger.json): 31 unique orbits and 35 matter-channel contributions.
- [Function dictionary](docs/results/conventions.md) and [exact M/N–Hlog conversions](results/2026-10-23/conversion-dictionary.json).
- [Validation and limitations](docs/results/validation-scope.md).
- [Downloadable archives](https://github.com/mitsuroni/wilson-loop-three-loop-results/releases): the new three-channel delivery and the historical all-diagram catalogue.

The linked guides and result reading pages render formulas directly on GitHub. For offline reading, download the latest complete research snapshot from Releases and open `docs/data-flow/index.html`. The dated three-channel ZIP remains available with its original data and `index.html` report.

## Results and scope

The CF CA TF nf, CF^2 TF nf and CF (TF nf)^2 channels are available through epsilon^0, including coupling and contour renormalization. Their totals use the original M/N moments, logarithms and zeta constants. All three pass the cusp anomalous-dimension, pole-structure and finite RG checks in the stated conventions.

The mixed CF CA TF nf finite result is calculated, but its independent original-input reconstruction remains open. Independent finite reconstructions for the other two matter colors are preserved historical evidence. A same-definition external complete finite comparison remains open. Computational checks are distinguished from formal proofs.

The direct pure-gluon CF CA^2 calculation still lacks the finite coefficients of double-TGV diagrams 15 and 24 (reflection-related), and the single-pole and finite coefficients of diagram 16. These gaps do not affect the three delivered matter channels. Additional total-level research does not determine those missing individual coefficients.

## Reproduce

With Python 3.12, SymPy 1.14.0 and mpmath 1.3.0, run from the repository root:

```text
python results/2026-10-23/reproduce/reproduce.py --output-dir replay-output
```

The script verifies the bundled source hashes and writes into a new directory. It requires no Maple or Mathematica session. The historical HyperInt mixed-argument summation witness is reused, not newly integrated. The packaged standalone replay passed; see its [record](results/2026-10-23/validation/standalone-reproduction.json).

## Provenance

The versioned delivery contains 185 source snapshots, exact expressions, check records and 41 offline report pages. Sources retain their original bytes and scientific revision labels. [Source manifest](results/2026-10-23/source-manifest.json) and [package manifest](results/2026-10-23/package-manifest.json) identify the snapshot. This repository is a delivery repository, not the Wilson calculation framework or its live selector database.

The October 18 archive in Releases preserves the broader 149-orbit catalogue as a historical snapshot. Its progress labels should be read alongside the newer three-channel report. Original third-party notices are retained under the delivery assets; this repository does not assign a new license to those materials.
