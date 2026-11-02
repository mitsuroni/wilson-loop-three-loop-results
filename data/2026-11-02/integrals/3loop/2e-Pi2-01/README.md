# 2e-Pi2-01: graph and original input

Stable ID `L3/2e-Pi2-01`; family `2e-Pi2`; 3 loop(s).
This is a source orbit/canonical representative, not a claim of one microscopic
Feynman graph. Effective self-energy/vertex operators retain their stated
subgraph loop order. [Machine input index](input.json) contains exact file and
JSON-field hashes; paths inside that JSON are relative to `data/2026-11-02`.

## Input definition and source

- [Adopted frozen native graph definition](../../../snapshots/workflow/checks/shared-sign-production/generated/qgraf/3loop/orbit_definitions.json) at `/orbits/2`.
- [Adopted frozen native C/J package](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/2e-Pi2/montopar-initial.json) at `/orbits/0`.
- Raw C is at `/orbits/0/C`. Its complete combined expression, color
  channels, phase ledger and component semantics are authoritative; components
  must not be blindly multiplied twice.

- `J`: `ordered_s`; [complete source](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/2e-Pi2/montopar-initial.json) at `/orbits/0/J`. Declared order: `['s31', 's41']`.
- `J_raw_s`: `ordered_s`; [complete source](../../../snapshots/workflow/checks/shared-sign-production/generated/montopar/3loop/2e-Pi2/montopar-initial.json) at `/orbits/0/J_raw_s`. Declared order: `['s31', 's41']`.

## Available prepared representations

- `native_contract_prepared_snapshot`: [payload](../../../snapshots/workflow/checks/shared-sign-production/generated/prepared/2e-Pi2.json) at `/orbits/0`. Full native-sign source contract, one Euler compensation, no legacy sign conversion; the existing original C/J and all domain information are retained.
- `existing_mathematical_prepared_r002_snapshot`: [payload](../../../snapshots/benchmarks/self-energy/pi2/prepared/2e-Pi2-r002.json) at `/orbits/0`. Existing versioned mathematical prepared input. Its sign/color projection and I_target must be read from this payload, not inferred from the newer native-contract wrapper.

## Versioned calculation source archives

- [Versioned calculation source](../../../snapshots/benchmarks/self-energy/pi2/sources/2e-Pi2-r002.json) at `/orbits/0` (`r002`). Its historical sign, raw/IBP stage and normalization are preserved; no version equivalence is inferred.
- [Versioned calculation source](../../../snapshots/benchmarks/self-energy/pi2/sources/montopar-2e-Pi2-r001.json) at `/orbits/0` (`montopar-2e-Pi2-r001`). Its historical sign, raw/IBP stage and normalization are preserved; no version equivalence is inferred.

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
