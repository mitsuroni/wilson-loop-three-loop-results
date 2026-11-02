# Double-TGV source contract (2026-09-19)

The schema 1.3 mathematical output is unchanged. Source validation now uses a closed topology contract before constructing C/J.

Supported encodings:

1. Four Wilson vertices, four external propagators, four external phases, and one effective TGV with four ordered external loop names. The first two and last two loops define the oriented vertex pairs.
2. The same external objects and two TGV loop triples sharing one internal loop name. Cyclic rotation puts this connector last at each vertex. The connector is implicit; it must not occur in the external propagator or phase lists.

Explicit connector propagators are unsupported and rejected. Explicit internal-vertex Lorentz/color annotations are also unsupported and rejected; indices at each vertex are inherited from its ordered external loops. Do not supply a detailed tensor graph and expect unused fields to be ignored.

Each Wilson vertex has exactly `line`, `lorentz`, `color`; each external propagator exactly `loop`, `lorentz`, `color`, with two Lorentz and two color endpoints. Each phase has exactly `loop`, `expression`. Each internal vertex has `type`, `loops`, and optionally the descriptive `id`. Other top-level descriptive metadata is preserved. Wilson lines are integer 1..4.

The four propagators must cover the four distinct Wilson endpoints bijectively. External colors must match their endpoints. Internal external-leg indices must be distinct within each index namespace and disjoint from the Wilson indices. Loop and phase sets must match exactly. Additional vertices, propagators, phases, index collisions and unsupported per-propagator attributes are rejected. Existing ordered-domain and phase-geometry validation remains in force.

Array order does not define vertex pairing. Reversing both index arrays of an external propagator is permitted. Cyclic triples and exchanging entire vertex groups preserve the previously documented orientation conventions.

`validate_double_tgv_result` still reconstructs outputs for consistency, but the reconstruction must first pass this source contract. It is not an independent mathematical proof. Current 26 production sources and their C/J are unchanged by this guard repair.

For orbit 04, keep both numerator terms when using a finite regulator. The special simplified expression in the theory document is an analytic/boundary-value simplification, not an identity at finite eta. This repair changes neither the numerator formulas nor the Feynman prescription.
