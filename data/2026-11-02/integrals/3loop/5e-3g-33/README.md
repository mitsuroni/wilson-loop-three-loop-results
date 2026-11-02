# 5e-3g-33: graph and original input

Stable ID `L3/5e-3g-33`; family `5e-3g`; 3 loop(s).
This is a source orbit/canonical representative, not a claim of one microscopic
Feynman graph. Effective self-energy/vertex operators retain their stated
subgraph loop order. [Machine input index](input.json) contains exact file and
JSON-field hashes; paths inside that JSON are relative to `data/2026-11-02`.

## Input definition and source

- [Adopted frozen native graph definition](../../../snapshots/workflow/checks/shared-sign-production/generated/qgraf/3loop/orbit_definitions.json) at `/orbits/98`.
- [Adopted frozen native C/J package](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/5e-3g/montopar-initial.json) at `/orbits/32`.
- Raw C is at `/orbits/32/C`. Its complete combined expression, color
  channels, phase ledger and component semantics are authoritative; components
  must not be blindly multiplied twice.

- `J`: `source-native representation`; [complete source](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/5e-3g/montopar-initial.json) at `/orbits/32/J`. Declared order: `block-specific / see complete payload`.
- `J_raw_s`: `source-native representation`; [complete source](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/5e-3g/montopar-initial.json) at `/orbits/32/J_raw_s`. Declared order: `['s12', 's11', 's22', 's21', 's31']`.
- `J_ibp`: `source-native representation`; [complete source](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/5e-3g/montopar-initial.json) at `/orbits/32/J_ibp`. Declared order: `block-specific / see complete payload`.

## Available prepared representations

- `historical_TGV_calculation_prepared_snapshot`: [payload](../../../snapshots/benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json) at `/orbits/32`. Existing single-canonical-representative calculation input; its original C, projected C, raw/IBP definitions and I map are preserved. Native-sign equivalence uses the separate historical production bridge; no equality is inferred from names.

## Versioned calculation source archives

- [Versioned calculation source](../../../snapshots/benchmarks/tgv/5e-3g/sources/canonical-r001.json) at `/orbits/32` (`canonical-r001`). Its historical sign, raw/IBP stage and normalization are preserved; no version equivalence is inferred.

## Normalization and grouping

The adopted native physical convention is
`wilson-qcd-action-path-fourier/1.0`: late Wilson generators act on the left,
Fourier phases and insertion signs are already part of raw C/J. Raw Euler
count is zero; use the selected prepared payload's explicit Euler/color rules.
Do not add the historical leg/Gamma minus sign again. C4 multiplicity,
crossing and fixed-basis a1 preparation are different operations; their
source fields are retained and must each be applied at the documented layer.
The exported default47 a1 value is not the full physical color vector.

Four-leg joint object: `not applicable`.
If applicable, consult the [joint-source map](../../../four-leg-joint-map.json).
Normalized physical I sources add with +1; their raw J prefactors differ.
DT25 already contains the rotated second channel. Do not add individual
sources again after consuming a joint object. All zero source records stay
in the catalogue.

## Available now and later work

This batch publishes graph definitions and original C/J plus the existing
prepared representations listed above. Laurent coefficients and their
derivations are **not exported by this batch**, and no new integral evaluation
was run. This omission is not a statement that an existing mathematical
result is unknown. Frozen payload `evaluation_status` fields retain their
historical meaning and do not override later research results. Branches,
regulator groups, dependent bounds and block-specific orders are unchanged;
there is no newly asserted Fubini or general complex-continuation theorem.
`formal=not_run`.
