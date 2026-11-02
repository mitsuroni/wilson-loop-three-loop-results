# Stage 1: homemade diagram generation

The package called `qgraf-pipeline` contains the project's native combinatorial generator and orbit classifier. Its name also reflects optional tools for consuming QGRAF output; the native catalogue route does not require a QGRAF executable.

The generated records describe ordered attachments to the four Wilson-loop segments, internal connectivity, family/loop labels, representatives, cyclic or reflection relations, and associated algebraic metadata. Frozen definitions carry stable orbit IDs used by the input catalogue and later result ledgers.

## Read and run

- [Package README and actual run commands](../../packages/qgraf-pipeline/README.md).
- [Published source](../../packages/qgraf-pipeline/src/qgraf_pipeline/).
- [Original tests](../../packages/qgraf-pipeline/tests/).
- [Frozen definitions and per-integral entries](../../data/README.md).

Use an explicit new output directory. For systems without TeX, use the documented non-PDF route. PDF rendering is a presentation step separate from graph generation and JSON export. The package README records the supported entrypoint and any publication portability adaptation.

## Version boundary

The installed source at publication time and the preserved producer versions of earlier JSON exports are distinct. Per-file manifests identify the actual source bytes, and the data catalogue records the producer/version information that is available for each snapshot. A successful structural replay validates the compared output fields; it is not a new proof of QCD diagram completeness.

Proceed to [parameterization](parameterization.md).
