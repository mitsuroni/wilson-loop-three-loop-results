# Double-TGV parent operator export

`--double-tgv-operator` adds `J_operator` to each double-TGV result and initial
package. Its schema is `double-tgv-operator-1.0`; the existing root schema stays
1.3. The ordinary scalar fields are still generated for compatibility. This
mode preserves a choice for the consumer; it does not yet avoid the producer's
existing scalar computation or promise a runtime improvement.

```powershell
python -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type "4e-3g+3g" --double-tgv-operator --package-initial
```

It can be combined with `--double-tgv-forms` and `--tgv-ibp`. Existing results
can be packaged with `--package-only --package-initial`; an export-generation
flag cannot be used with `--package-only`. Only double-TGV targets accept the
new flag. Other families and the default output remain unchanged.

## Mathematical object

The same raw C multiplies

$$
J_{\rm op}=\frac1{(1-3\epsilon)A_0}
\int_{\mathcal D}d^4s\int d\mu_5\,
\prod_r\alpha_r^{-\epsilon}U^{-1-2\epsilon}
\mathscr R\left[\mathcal O(B_0/A_0+i\eta_{\rm reg}U)^{-1+3\epsilon}\right].
$$

`kernel.definitions` contains an explicit polynomial in six independent
distances d_ij=-(Xi-Xj)^2, together with U. `operator.tensor_terms` gives the
four vertex coefficient tensors: identity and tangent-vector outer products.
The overall minus sign in O=-VL dot VR is included. Partial derivatives act in
independent D-dimensional positions, holding the tangents, A0, alpha, epsilon
and the coordinate auxiliary regulator fixed. The D trace is 4-2 epsilon;
the contour-vector Gram matrix is not the spacetime metric trace.

`restriction.positions` subsequently maps each Xi to its original Wilson
path. `operation_order` requires parent differentiation before restriction.
All six parent distances remain even when a distance is zero on that path.
The block copies the original ordered domain, dependent bounds, outer-to-inner
order and five-alpha simplex. No integral reordering, endpoint subtraction,
Euler compensation, color projection, C4 multiplication or crossing is done.

eta_reg is a dimensionless coordinate auxiliary regulator, kappa=A0*eta_reg.
The chain-rule conversion is valid at a shared finite regulator. Physical
limits and derivative/integral exchanges still require the original common
boundary-value or analytic-continuation argument; this export adds no proof
of those exchanges.

## Scalar consumer

```python
from montopar.double_tgv_operator import scalar_fields, to_scalar

fields = scalar_fields(result["J_operator"])
ordered_scalar = to_scalar(result["J_operator"])
assert ordered_scalar == result["J_raw_s"]
```

The consumer differentiates the exported parent B0 with respect to its
distance variables, computes parent gradients and cross Hessians, contracts
the exported tensors, and only then uses the path map. It retains both P1 and
P2 and the original C/J convention. It does not call `compute_numerators` or
read a saved scalar answer from the block. Conversion returns the existing
ordered scalar representation with the same C; do not multiply a second C.

All geometric fields are computed before color screening. On a color-zero
orbit the effective scalar integrand is zero while its geometric kernel and
unapplied operator remain available for checks.

This is a finite descriptor for the current two-vertex module, not a general
tensor/executable language. The standalone consumer checks derivative scope,
the supported tensor constructors, symbol definitions, simplex, normalization
and regulator conventions. Full result validation additionally rebuilds the
descriptor from the original oriented source and checks its scalar recovery.
Hashes identify the source; they are not mathematical or Lean evidence.
