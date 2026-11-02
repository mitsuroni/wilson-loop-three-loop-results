# Double-TGV ordered two-layer IBP output

2026-09-20. Applies to the four-external-leg, two-TGV source contract already
accepted by the raw exporter. The complete derivation is
[theory r003](4e-3g+3g-ibp-theory-r003.md). The existing raw derivation remains
the reference for `C`, `J_raw_s` and `J_raw`.

```powershell
python -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type "4e-3g+3g" --tgv-ibp --package-initial
```

The Python entry point is `build_double_tgv_result(orbit, ibp_enabled=True)`.
Without this option the old schema 1.3 JSON is byte-identical. With this option
one additional `J_ibp` block is emitted, with its own schema
`double-tgv-ibp-1.0`. There is one source-oriented representative, with no A/B
duplication. Package validation recomputes and compares the entire IBP block.

## How to consume C and J

The total is `Sum[term.combined_C * Integral[term.integrand_s]]`.
Each term also exports a cube `integrand`, its map and its Jacobian. Choose
one of the ordered/cube representations for each term. Do not multiply the
common `C` again after using `combined_C`.

Each term contains:

- `C_ratio_to_raw` and `combined_C`;
- a scalar numerator, denominator polynomial and exponent;
- all four independent alpha integrations, their simplex constraint and weight;
- ordered variables, bounds and explicit outer-to-inner order;
- a composed map from the face back to all four parent path parameters;
- the boundary/vertex sources, including the contact-normal correction;
- one shared regulator group for the complete orbit.

`G_L_edge_basis` and `G_R_edge_basis` are the dimensionful vectors
`R_L Bhat` and `R_R Bhat` in the four-edge basis. In that basis the metric
`dot_product_factor` gives `u_i.u_j/v14`. Consequently `k_over_A0=-k_over_v14/2`,
since `A0=-2*v14`. Final `pnum_s` and `C_ratio_to_raw` are dimensionless.
The D-dimensional Hessian keeps `D=4-2*eps` before any expansion.

The four scalar row types are:

| Kind | Ratio C_term/C_raw | J numerator | Denominator power |
|---|---|---|---|
| double_boundary | -k/(lambda*A0) | 1 | lambda |
| single_boundary | 1 | N_face/A0 | lambda+1 |
| double_recoil_gradient | -(lambda+1) | G_L.G_R/A0 | lambda+2 |
| double_recoil_hessian | 1 | H_RR/A0 | lambda+1 |

Here `lambda=1-3*eps`. The last two rows form one body-remainder group.
All rows include `(product alpha)^(-eps)*U^(-1-2*eps)` in J.
No Euler factor, color projection, crossing or C4 multiplicity is inserted.
For color-zero graphs, effective integrands vanish; geometric expressions
are retained separately so that the algebra can still be audited.

## Regulated equivalence

The construction starts from `V_L=T_L+R_L` and `V_R=T_R+R_R`, applies
contour IBP to TT/TR/RT, and retains RR. On a contact face the remaining
parent derivative has both a tangential component and a normal component.
Derivatives of the parent B are evaluated before restriction to that face.

The exact equality applies at fixed interior alpha and positive dimensionless
eta, using `Bhat+I*eta*U`. `eta=kappa/A0`; kappa is the dimensionful coordinate
regulator (called eta in the historical raw branch metadata).
The optional `regulated_integrand_s_before_color_screening` makes the finite
regulator explicit. Common alpha cutoffs can then be imposed on the full sum.
Take regulator limits and analytic continuation for the grouped expression.
Individual unregulated terms are not claimed to share a convergence strip.
In particular, a boundary with geometric B=0 must not be set to zero or treated
as a separately convergent integral by dropping the regulator.

This is integral equality under IBP; the lower-dimensional face integrands
are not pointwise equal to the original four-path-variable integrand.
The auxiliary-regulator numerical checks do not prove the physical complex
boundary limit or the original eight-dimensional singular integral.

## Validation and scope

All 26 current geometries, including the nine color-zero kernels, are checked
against frozen raw output. Checks include symbolic epsilon-dependent gradient
and Hessian identities, complete facet enumeration, composed maps, ordered
polynomial moments, and numerical contour integration at fixed alpha.
The mutation tests reject missing normal terms, signs, C factors, domains,
regulators and D-trace changes. The 04/06/25 counts are 9/0/0, 7/6/1 and
16/8/1 (double boundaries / single boundaries / body groups).

This exports a representation, not evaluated integrals or Lean proofs.
Stage 3's existing single-J bootstrap still needs an explicit multi-term
consumer before it can process the full `J_ibp.terms` automatically.
The production generator uses Python exact polynomial arithmetic, with no
new CAS runtime dependency. Independent checking uses SymPy and NumPy.
