# Wilson-loop diagram generator

This is the published source snapshot of `qgraf-pipeline` 0.1.0. Its default
production path is the project's **combinatorial diagram generator**, followed
by cyclic-orbit classification and momentum-space definition export. The name
does not mean that a QGRAF binary generated the released diagram set.

The supported catalogue at loop orders 1, 2 and 3 contains respectively 2, 12
and 149 orbit representatives. These are counts for the specified finite null
quadrilateral calculation, not a proof of general QCD graph completeness.

## Installation and use

From the public repository root, with Python 3.12 or newer:

```text
python -m pip install -e packages/qgraf-pipeline
qgraf-pipeline 2 --output-dir generated/2loop --no-pdf
```

The original declared minimum for this package alone is Python 3.10. Python
3.12 is recommended when using it together with Montopar. Core Python modules
use the standard library. `--no-pdf` is the single publication adapter: it skips
only the optional `pdflatex` rendering step and preserves all mathematical
generation and exports. The default without this flag retains the original
PDF-rendering behavior and needs a suitable LaTeX installation.

This command produces `integrals.md`, `integrals.json`, `diagrams.tex` and the
downstream `orbit_definitions.json`. Use an output path ending in `<n>loop`, as
shown: the historical definition exporter constructs that loop subdirectory
under the output directory's parent. Always provide a fresh output directory.

The equivalent module entry point after installation is:

```text
python -m qgraf_pipeline.cli 2 --output-dir generated/2loop --no-pdf
```

`--out-file existing.out` parses an already generated QGRAF `.out` file instead
of using the combinatorial generator. Native QGRAF model/control/runner helpers
are retained in the source, but no QGRAF executable or local binary path is
bundled or needed by the default generator. This command does not invoke QGRAF.

`--audit-recursion` writes a distinct topology audit at loop orders 1--3. It is
not the normal production catalogue and must not silently replace it.

## Physical conventions and output boundary

The native source declares `wilson-qcd-action-path-fourier/1.0`: ordered Wilson
insertions, late-to-early color trace, Fourier `exp(-ipX)`, momentum Wick `+i`,
coordinate Wick `-i`, complete self-energy/vertex internal-loop signs and
external Dyson coefficient `+1`. Edge labels, ordered attachment bounds,
numerators, momentum routing and phase conventions belong to the definition.

The output is a graph/orbit and momentum-space specification. It does not
evaluate parameter integrals, compute Laurent coefficients, or prove the
downstream physics. Montopar supplies the next C/J representation.

## Source snapshot and tests

`PUBLICATION-SOURCE.json` records every original and published hash, omitted
directories and the exact `--no-pdf` adapter. Source/tests otherwise preserve
their original bytes. [Original README](README.original.md) and theory docs
retain their dates and historical terminology. Current dirty source bytes are
not automatically identical to the programs that produced an older frozen
research input. Compare source and product manifests before claiming a replay.

Install test dependencies and run the original tests from this directory:

```text
python -m pip install -r requirements-test.txt
python -m pytest -q tests -p no:cacheprovider
```

No personal agent configuration, output directories, compiled binaries or
third-party library downloads are part of this package.
