# Frozen upstream inputs: 2026-11-02

This version publishes the adopted one-, two- and three-loop graph catalogue,
the complete original Montopar C/J packages, and the prepared representations
that were found in the calculation archive. Original payloads are copied
**byte for byte**. Newly generated per-integral JSON files are indexes, not
new definitions of the integrals.

| Coverage | Count |
| --- | ---: |
| One-loop source orbits | 2 |
| Two-loop source orbits | 12 |
| Three-loop source orbits | 149 |
| Families / native raw C/J packages | 13 |
| Orbits with frozen native raw C/J | 163 |
| Default original ordered-s inputs | 47 |
| Orbits with a separate existing prepared payload | 126 |
| Prepared representation occurrences, including distinct versions | 144 |
| Orbits without a separate prepared payload in this export | 37 |

See [catalog.json](catalog.json), [coverage.json](coverage.json),
[source-manifest.json](source-manifest.json), and the
[nine four-leg joint objects](four-leg-joint-map.json).
The catalogue links each orbit's own English page and `input.json`.

Examples:

- [2e-01](integrals/1loop/2e-01/README.md): one-loop exchange and its
  original ordered-s / prepared0.3 alternatives.
- [2e-Pi1-01](integrals/2loop/2e-Pi1-01/README.md): a self-energy operator,
  ordered kernel, full native-sign C and existing prepared versions.
- [4e-3g+3g-15](integrals/3loop/4e-3g+3g-15/README.md): original raw and IBP
  source archives beneath joint I2344; no new finite calculation here.
- [3e-Gamma1-04](integrals/3loop/3e-Gamma1-04/README.md): the full channel/block
  sum, with per-block measures, orders and combined C factors.

## Adopted, historical and current snapshots

The adopted graph and raw C/J files are the project's **frozen native-sign
production snapshots**, bound by `source-generation.json`, `all-native.json`
and `native-check.json`. All three recorded graph hashes and all13 recorded
C/J package hashes were checked before export. Their native convention is
`wilson-qcd-action-path-fourier/1.0`.

The recorded graph producer is the combinatorial qgraf-pipeline
generator/classifier/exporter. That production record explicitly says the
native QGRAF executable was not used. This data exporter does not rerun that
producer or any CAS. Selected new package-level replays are separate checks;
they are not substituted for these original frozen records.

Earlier graph definitions, original ordered/cube files and versioned
mathematical raw/IBP archives remain explicitly labelled alternatives.
These matter because a mathematical derivation may bind an earlier source
file rather than a later regeneration. A shared orbit ID or similar formula
does not silently grant equality between versions.

[current-observations.json](current-observations.json) separately records
read-only comparisons with the external current workspace. Fifteen of the
sixteen observed family/graph files match the adopted frozen bytes. The
current double-TGV family file differs; it is **not** silently promoted to
the historical producer or substituted for the frozen source. Original local
paths inside payloads are provenance, not required runtime locations.

## Original domains and prepared inputs

The three default ordered original packages were found under the original
manifest's `benchmarks/2ne` path base, and all recorded hashes match. Their
47 original s-parameter kernels preserve dependent bounds and the exact
outer-to-inner order. The prepared0.3 versions retain the original C,
complete color coordinates, the a1 projection and one Euler compensation.

The remaining prepared snapshots include the native-contract self-energy
inputs, the versioned mathematical self-energy inputs, the two TGV families,
six nonzero FGV orbits and eight double-TGV representatives. FGV prepared
inputs explicitly use a unit cube; the TGV and selected DT inputs retain
their declared s domains. No domain is relabelled and no integral order is
permuted by this export. Multiple prepared versions are alternatives, not
terms to be added together.

The 37 entries without a separate prepared payload are the12 leg-insertion
orbits, four Gamma1 orbits,18 DT source orbits and three zero-contact FGV
orbits. Their native raw ordered or channel/block inputs remain present.
This is an input-format coverage gap, not a claim that37 integral values
are missing. Gamma1 has a full block sum rather than one scalar integrand;
its per-term `combined_C` already contains the common C.

## Factors that must remain distinct

- Raw native C/J has no extra Euler compensation. Prepared files state
  their own color and Euler operation. A projected a1 coefficient is not
  the full physical color vector.
- Native leg/Gamma signs already follow the shared action/path/Fourier
  convention. Applying the old family minus sign again is an error.
- C4 orbit multiplicity, stabilizer, crossing and prepared color projection
  are distinct fields. Apply each only at its documented normalization layer.
- Raw and IBP forms, or ordered and cube forms, describe alternative
  representations. In IBP/block data, a `combined_C` must not be multiplied
  by the global C a second time.
- The35 FGV/DT source orbits map to nine ordered four-leg objects. Keep all
  zero sources. Normalized nonzero I values add with +1, while raw J kernels
  carry their original unequal prefactors. DT25 already includes its rotated
  second channel. Consume either a joint object or its source components,
  without counting both.

## Limits of this publication

This batch supplies definitions, original kernels, versioned source links
and existing preparation data. It supplies no new graph-completeness proof,
no re-evaluation of all163 amplitudes, no new Laurent derivation, and no
general complex-branch or Fubini theorem. In particular, a historical
`evaluation_status` or `whole_integral_complete` field is retained unchanged;
it is not the current result selector.

The public [source packages](../../packages/README.md) explain executable
generation and its external CAS requirements. The standalone file audit in
[data/README.md](../README.md) checks provenance and pointers only.
`formal=not_run`.
