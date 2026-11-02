# Stage 2: Montopar C/J definitions

Montopar reads the supported orbit-definition schema and applies the family's parameterization formulas, geometry, color rules and optional contour-IBP identities. Consult the [package README](../../packages/montopar/README.md) for precise CLI syntax, Python version and color-backend requirements. Some families use a Wolfram/FeynCalc color calculation; native branches are documented separately. Generating parameter integrals is not a CAS-free operation for every family.

## Preserve the full definition

An integral input consists of its coefficient C, J integrand, measure, variables, ordered bounds, contour/continuation prescription and representative metadata. The outer-to-inner integration order is part of the definition. An ordered-s expression and its unit-cube transform describe two representations of the same object under the stated change of variables; they are not separate terms to add.

For an optional IBP output:

```math
\mathcal A_\rho=\sum_t C_tJ_t,\qquad
C_t=C_{\rho,\mathrm{raw}}r_t.
```

`term.combined_C` already contains the global C. Keep every term's own integrand, variables, bounds, boundary origin and common regulator group. Reapplying the global C would change the definition.

## Normalization is stage-specific

Raw inputs retain their recorded color and Gamma/Fourier conventions. Downstream prepared inputs may perform a fixed-basis color projection and one Euler compensation. The normalized physical I additionally follows its documented orbit/crossing convention. A raw C/J entry, a prepared a1 entry and a full physical color-channel contribution must therefore be identified by their own stage.

Historical versus native sign conventions are preserved as versioned inputs. Do not apply a historical extra sign to a native source that already contains it. The [input catalogue](../../data/README.md) records actual source locations and hashes rather than selecting data solely from an old status label.

## Read the implementation

- [Source modules](../../packages/montopar/src/montopar/).
- [Original tests](../../packages/montopar/tests/).
- [Package documentation](../../packages/montopar/README.md).
- [Frozen input data](../../data/README.md).

The source snapshot and the frozen calculation inputs can have different producer versions. Their documented relationship should be checked before treating a new generation as a reproduction of a particular historical input.
