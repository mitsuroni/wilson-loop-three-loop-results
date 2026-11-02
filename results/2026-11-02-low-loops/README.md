# One and two loops: an end-to-end research example

These fourteen orbit packets connect the published diagram and momentum definitions
to parameter integrals, individually evaluated Laurent series, the connected color
sum, renormalization and the finite Wilson loop. They are the publication pilot
before curating the complete three-loop calculation.

Read the [continuous calculation](derivation.md), then follow an orbit through the
table below. The [machine flow map](data-flow.json) retains original graph/raw-source
pointers and checksums. Shared programs belong to a family; each orbit has its own
derivation, parameters and exact result.

| Loops | Orbit | Family | Graph / momentum / C,J | English derivation | Exact normalized integral |
| --- | --- | --- | --- | --- | --- |
| 1 | 2e-01 | exchange | [input](../../data/2026-11-02/integrals/1loop/2e-01/input.json) | [derivation](../../integrals/one-loop/2e-01/derivation.md) | [I](../../integrals/one-loop/2e-01/result.json) |
| 1 | 2e-02 | exchange | [input](../../data/2026-11-02/integrals/1loop/2e-02/input.json) | [derivation](../../integrals/one-loop/2e-02/derivation.md) | [I](../../integrals/one-loop/2e-02/result.json) |
| 2 | 4e-01 | pure | [input](../../data/2026-11-02/integrals/2loop/4e-01/input.json) | [derivation](../../integrals/two-loop/4e-01/derivation.md) | [I](../../integrals/two-loop/4e-01/result.json) |
| 2 | 4e-02 | pure | [input](../../data/2026-11-02/integrals/2loop/4e-02/input.json) | [derivation](../../integrals/two-loop/4e-02/derivation.md) | [I](../../integrals/two-loop/4e-02/result.json) |
| 2 | 4e-03 | pure | [input](../../data/2026-11-02/integrals/2loop/4e-03/input.json) | [derivation](../../integrals/two-loop/4e-03/derivation.md) | [I](../../integrals/two-loop/4e-03/result.json) |
| 2 | 4e-04 | pure | [input](../../data/2026-11-02/integrals/2loop/4e-04/input.json) | [derivation](../../integrals/two-loop/4e-04/derivation.md) | [I](../../integrals/two-loop/4e-04/result.json) |
| 2 | 4e-05 | pure | [input](../../data/2026-11-02/integrals/2loop/4e-05/input.json) | [derivation](../../integrals/two-loop/4e-05/derivation.md) | [I](../../integrals/two-loop/4e-05/result.json) |
| 2 | 4e-06 | pure | [input](../../data/2026-11-02/integrals/2loop/4e-06/input.json) | [derivation](../../integrals/two-loop/4e-06/derivation.md) | [I](../../integrals/two-loop/4e-06/result.json) |
| 2 | 3e-3g-01 | TGV | [input](../../data/2026-11-02/integrals/2loop/3e-3g-01/input.json) | [derivation](../../integrals/two-loop/3e-3g-01/derivation.md) | [I](../../integrals/two-loop/3e-3g-01/result.json) |
| 2 | 3e-3g-02 | TGV | [input](../../data/2026-11-02/integrals/2loop/3e-3g-02/input.json) | [derivation](../../integrals/two-loop/3e-3g-02/derivation.md) | [I](../../integrals/two-loop/3e-3g-02/result.json) |
| 2 | 3e-3g-03 | TGV | [input](../../data/2026-11-02/integrals/2loop/3e-3g-03/input.json) | [derivation](../../integrals/two-loop/3e-3g-03/derivation.md) | [I](../../integrals/two-loop/3e-3g-03/result.json) |
| 2 | 3e-3g-04 | TGV | [input](../../data/2026-11-02/integrals/2loop/3e-3g-04/input.json) | [derivation](../../integrals/two-loop/3e-3g-04/derivation.md) | [I](../../integrals/two-loop/3e-3g-04/result.json) |
| 2 | 2e-Pi1-01 | Pi1 | [input](../../data/2026-11-02/integrals/2loop/2e-Pi1-01/input.json) | [derivation](../../integrals/two-loop/2e-Pi1-01/derivation.md) | [I](../../integrals/two-loop/2e-Pi1-01/result.json) |
| 2 | 2e-Pi1-02 | Pi1 | [input](../../data/2026-11-02/integrals/2loop/2e-Pi1-02/input.json) | [derivation](../../integrals/two-loop/2e-Pi1-02/derivation.md) | [I](../../integrals/two-loop/2e-Pi1-02/result.json) |

## What has actually been delivered

One-loop I coefficients are explicit through epsilon3; all twelve two-loop I
coefficients are explicit through epsilon2. [coefficients.json](coefficients.json)
provides the scale-free color sum U, restored bare exponent, coupling counterterm,
post-counterterm series, pole logZ, finite F, and cusp/noncusp coefficients. The
[ledger](assembly-ledger.json) records every color weight. One-loop epsilon2 and
two-loop epsilon1 have also been [reassembled and matched](three-loop-budget.json)
to the accepted lower-order input for the three-loop finite counterterm.

The finite one/two-loop results have an exact comparison to independently
translated author ancillary expressions, in addition to RG and cusp checks.
The full two-loop epsilon2 series may contain the already defined tau_D4 function;
the physical epsilon0 finite result needs only logarithms, M/N and zeta values.

## Reproduce from a checkout

Use Python 3.12 and the [recorded dependencies](requirements.txt). From the
repository root, always choose new output directories:

```text
python -m pip install -r results/2026-11-02-low-loops/requirements.txt
python tools/check_low_loops.py --output-dir replay-output/audit-new
python results/2026-11-02-low-loops/reproduce.py --output-dir replay-output/totals-new
python results/2026-11-02-low-loops/reproduce.py --output-dir replay-output/full-new --regenerate-integrals
python integrals/two-loop/check_integrals.py --output-dir replay-output/integral-controls-new
python results/2026-11-02-low-loops/numeric.py --output-dir replay-output/finite-numeric-new
```

The first mathematical command consumes the released individual tables and rebuilds
all downstream stages. Adding `--regenerate-integrals` first executes the actual
one-loop primitive/Gamma/moment calculation and the two-loop mother-function/HLOG
engines. Historical HyperInt integrations and fixed HLOG identities are explicitly
reused where stated; no new Maple or Mathematica execution is claimed. The separate
integral controls evaluate the original/reduced kernels listed in the family guide.
The last command evaluates the original M11 integral, the finite functions and their
external comparison at two numerical precisions. To read freshly generated totals,
add `--coefficients replay-output/totals-new/coefficients.json`.

Original third-party ancillary downloads are not bundled. Default comparison uses
the attributed, source-bound [derived fixture](reference/fixture.json). If the
original files are available, repeat their exact translation without executing them:

```text
python results/2026-11-02-low-loops/reproduce.py --output-dir replay-output/reference-new --reference-dir /path/to/ancillary
```

Expected files and hashes are in [reference provenance](reference/provenance.json).
Archival scripts under `sources/` keep their original bytes and paths; use the
portable entrypoints above, rather than running historical output writers in place.

## Validation and its limits

Actual downstream reproduction passed [68 checks including 5 negative controls](checks.json).
An actual full child regeneration is recorded [separately](validation/fresh-regeneration.json).
The [new finite numerical samples](numeric-samples.json) compare 40 and 60 working
digits. The independent [version bridge](validation/version-bridge-checks.json)
checks the Pi1 source-sign migration at the C/J and individual Laurent levels.
See [isolated reproduction and review](validation/README.md) for the final publication
audit. Numerical differences are empirical controls, not certified intervals.

The Fourier/Feynman rules, source geometry and the mathematical identities identified
in each article are explicit research inputs. `formal=not_run`. Neither checksums nor
RG determine an arbitrary L-independent finite shape. The complete author finite
comparison supplies an additional independent constraint in this low-loop example.

[Publication lessons for three loops](../../docs/low-loops/curation-lessons.md) explain
what this pilot found difficult and what to carry forward. Frozen upstream packages,
input data and prior dated deliveries remain unchanged.
