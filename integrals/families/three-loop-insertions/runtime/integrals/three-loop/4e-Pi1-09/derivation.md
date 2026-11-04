# 4e-Pi1-09: Pi1-P4

This labelled source belongs to the fifteen-integral insertion batch. We derive the connected scalar I with its own source normalization, then retain the complete insertion-operator colour convolution. The [shared method](../../../methods/pi1.md) gives the continuous endpoint/IBP/Gamma argument; [the method index](../../../methods/README.md) fixes the physical and function conventions.

## Original momentum insertion and parameter source

The momentum source contains D0 Pi1 D0 with one bare exchange for the four-endpoint branch, D0 Pi1 D0 Pi1 D0 for a serial chain, or D0 Pi2 D0 for a single two-loop1PI insertion. Here D0 is the Feynman-gauge gluon propagator and Pi is the complete transverse effective self-energy operator, with its original phases. The adopted internal operator is an input; the present packet does not claim a new internal bubble/two-loop-self-energy diagram integration.

[input.json](input.json) binds the actual graph, raw C/J, full-colour native preparation and selected reference separately. The exact original ordered kernel is

```text
(2*s31 - 2*s21*s31)*((s31 - s21*s31)^(-2+2*eps)/(2*(1+eps)))*((-(1+z1))*(1 - s22 - s41 + s22*s41 + s22*s41*z1)^(-1+eps))
```

Its outer-to-inner order is s22 → s21 → s31 → s41.

| Parameter | Lower bound | Upper bound |
| --- | --- | --- |
| s22 | 0 | 1 |
| s21 | 0 | s22 |
| s31 | 0 | 1 |
| s41 | 0 | 1 |

The full transverse numerator is retained in this kernel, including its longitudinal/derivative part. Each polynomial and tensor identity is checked against the source, not inferred from an integer-regulator sample. Polynomial controls then compare this actual order with a separately constructed cube/Jacobian integration; their scope is described below.

The original multiplicative C components are:

```json
{
  "wilson_color": "CF*(CF-CA/2)",
  "self_energy_color": "((10-6*eps)*CA-8*(1-eps)*TF*nf)/(3-2*eps)",
  "coupling": "(alpha_s/Pi)^3",
  "gamma_fourier": "(Pi^eps*Gamma[1-eps])^3",
  "rational": "1/(32*eps*(1-2*eps))",
  "scale": "(-2*v14+I0)^(3*eps)",
  "euler": "exp(-3*EulerGamma*eps)"
}
```

Their phase ledger and full trace remain in the bound native source. Its normalization has assembly_sign=1, one Euler compensation and no native sign conversion. The scalar I uses its documented orbit/crossing extraction and, for the four-endpoint Pi1 branch, the connected fixed-basis a1. The unprojected raw trace is not identified with that projected object.

## Reduction, functions and required orders

The route is `Pi1-P4`, shared with `4e-Pi1-04` while retaining this source's own kernel. Its reduction and actual epsilon budgets are explicitly recorded in [subintegrals.json](subintegrals.json):

```json
{
  "raw_ordered_J": {
    "representation": "ordered_s",
    "variables": [
      "s22",
      "s21",
      "s31",
      "s41"
    ],
    "bounds": {
      "s21": [
        "0",
        "s22"
      ],
      "s22": [
        "0",
        "1"
      ],
      "s31": [
        "0",
        "1"
      ],
      "s41": [
        "0",
        "1"
      ]
    },
    "iterated_order_outer_to_inner": [
      "s22",
      "s21",
      "s31",
      "s41"
    ],
    "source_bounds_array": [
      {
        "var": "s22",
        "low": "0",
        "high": "1"
      },
      {
        "var": "s41",
        "low": "0",
        "high": "1"
      },
      {
        "var": "s21",
        "low": "0",
        "high": "s22"
      },
      {
        "var": "s31",
        "low": "0",
        "high": "1"
      }
    ],
    "order_basis": "ascending contour edges; descending slots; checked against source integrals.md at publication",
    "jacobian": "1",
    "integrand": "(2*s31 - 2*s21*s31)*((s31 - s21*s31)^(-2+2*eps)/(2*(1+eps)))*((-(1+z1))*(1 - s22 - s41 + s22*s41 + s22*s41*z1)^(-1+eps))",
    "tensor_decomposition": {
      "metric": "(2*s31 + 2*eps*s31 - 2*s21*s31 - 2*eps*s21*s31)*((s31 - s21*s31)^(-2+2*eps)/(2*(1+eps)))*((-(1+z1))*(1 - s22 - s41 + s22*s41 + s22*s41*z1)^(-1+eps))",
      "derivative": "(-2*eps*s31 + 2*eps*s21*s31)*((s31 - s21*s31)^(-2+2*eps)/(2*(1+eps)))*((-(1+z1))*(1 - s22 - s41 + s22*s41 + s22*s41*z1)^(-1+eps))",
      "identity": "integrand = metric + derivative; common C"
    },
    "se_numerator": "2*s31 - 2*s21*s31",
    "metric_numerator": "2*s31 + 2*eps*s31 - 2*s21*s31 - 2*eps*s21*s31",
    "derivative_numerator": "-2*eps*s31 + 2*eps*s21*s31",
    "geometry_polynomials": [
      {
        "loop": "l2",
        "kind": "se",
        "edges": [
          2,
          3
        ],
        "s_vars": [
          "s21",
          "s31"
        ],
        "X_coefficients": [
          "0",
          "-1 + s21",
          "-s31",
          "0"
        ],
        "f": "s31 - s21*s31",
        "d": "1",
        "u": "-s31",
        "w": "-1 + s21"
      },
      {
        "loop": "l1",
        "kind": "bare",
        "edges": [
          2,
          4
        ],
        "s_vars": [
          "s22",
          "s41"
        ],
        "X_coefficients": [
          "0",
          "-1 + s22",
          "-1",
          "-s41"
        ],
        "f": "1 - s22 - s41 + s22*s41 + s22*s41*z1",
        "d": "-(1+z1)",
        "u": "-1 + s41 + s41*z1",
        "w": "1 - s22 - s22*z1"
      }
    ],
    "power_prescription": "each f power is the boundary value (f+I0)^power in A_ref>0; continue jointly with C"
  },
  "source_paired_J_jet_expression": "((-M_0_1*eps**2*rho**2/2 + M_0_2*eps**2*rho - M_0_3*eps**2/2 - 3*M_1_1*eps**3*rho**2/4 - M_1_1*eps + 3*M_1_2*eps**3*rho/2 - 19*M_1_3*eps**3/24 - M_2_1*eps**2 - 13*M_3_1*eps**3/24)*cosh(eps*rho) - (-M_0_1*eps**3*rho**3/3 - 2*M_0_1*eps*rho + M_0_2*eps**3*rho**2 + 2*M_0_2*eps - 13*M_0_3*eps**3*rho/12 + 5*M_0_4*eps**3/12 - 3*M_1_1*eps**2*rho + 3*M_1_2*eps**2 - 9*M_2_1*eps**3*rho/4 + 9*M_2_2*eps**3/4)*sinh(eps*rho)/2)/(4*eps**2*(eps + 1)*cosh(3*eps*rho/2))",
  "required_jet_orders": {
    "Gamma": 2,
    "weighted_H": 3
  },
  "actual_conservative_Gamma_jet": 5,
  "external_operator_Taylor_through": 2,
  "primitives": "Existing original M moments and rho_id_0_0_0/1 HLOG evaluations; no final I table read."
}
```

The method derives these finite jets from the original kernel, preserving the full boundary term in the opposite-edge mixed derivative. The evaluator then uses an exact Gamma/rational jet and the admitted original M/HLOG dictionary. Existing evaluated moment data is reused as such; no final integral table is used for generation.

rho means log(z), not the physical average scale L. M_kj denotes the original no-factorial moment; N_kj=M_k,j+1-rho M_kj. The chain display abbreviations M01,M11,M21,M03 mean M_0,1,M_1,1,M_2,1,M_0,3. The separate source labels rho_id_0_0_0 and rho_id_0_0_1 denote the already evaluated T coefficient functions, rather than rho itself: the first is M_0,1 squared/6; both complete HLOG definitions are in the [bound structural dictionary](../../../snapshots/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json). These are explicit evaluated functions, not names for remaining unsolved integrals.

## Complete Laurent coefficients and external colours

At fixed positive z, the scalar series is sum(c_p(z)epsilon^p,p=-6,...,0)+O(epsilon). In this entry the scalar is the connected scalar I. All seven positions are explicit, including established zeros.

| Power | Exact scalar coefficient | HLOG coordinates |
| ---: | --- | ---: |
| -6 | $0$ | 0 |
| -5 | $0$ | 0 |
| -4 | $0$ | 0 |
| -3 | $0$ | 0 |
| -2 | $\frac{M_{1 1}}{128}$ | 5 |
| -1 | $- \frac{M_{0 1} \rho^{2} - M_{0 3} - 2 M_{1 1} - 2 M_{2 1}}{256}$ | 16 |
| 0 | $- \frac{12 M_{0 1} \rho^{2} - 12 M_{0 3} + 33 M_{1 1} \rho^{2} - 72 M_{1 1} - 6 \pi^{2} M_{1 1} - 19 M_{1 3} - 24 M_{2 1} - 13 M_{3 1}}{3072}$ | 40 |

The complete per-order coloured tables are in [result.json](result.json). Each row below gives a monomial outside its scalar multiplier; the last column is its complete finite coefficient after epsilon convolution.

| External monomial | Exact multiplier of the scalar I | Coloured finite coefficient |
| --- | --- | --- |
| `CA^2*CF` | `(10 - 6*eps)/(3 - 2*eps)` | $- \frac{M_{0 1} \rho^{2}}{72} + \frac{M_{0 3}}{72} - \frac{55 M_{1 1} \rho^{2}}{1536} + \frac{5 \pi^{2} M_{1 1}}{768} + \frac{35 M_{1 1}}{432} + \frac{95 M_{1 3}}{4608} + \frac{M_{2 1}}{36} + \frac{65 M_{3 1}}{4608}$ |
| `CA*CF*TF*nf` | `(8*eps - 8)/(3 - 2*eps)` | $\frac{M_{0 1} \rho^{2}}{144} - \frac{M_{0 3}}{144} + \frac{11 M_{1 1} \rho^{2}}{384} - \frac{\pi^{2} M_{1 1}}{192} - \frac{11 M_{1 1}}{216} - \frac{19 M_{1 3}}{1152} - \frac{M_{2 1}}{72} - \frac{13 M_{3 1}}{1152}$ |

The series of each external multiplier is retained through the actual scalar pole order before multiplication. For Pi2 its Gamma ratio is expanded to the higher intermediate order needed before cancellation; the full native A+B/U operator and the direct Gamma route agree. A serial-Pi1 operator is not substituted for it. Other connected colour channels absent from the complete polynomial are zero in the stated projected/connected convention; this does not erase the unprojected raw-source trace.

## Evidence and reproduction

The candidate was frozen before the final reference tables were opened. [The per-orbit reference check](checks.json) compares every scalar coefficient with the adopted current object and, where stored, each coloured coefficient. The whole batch also independently back-substitutes all343 scalar/colour function expressions into complete HLOG/period coordinates with rational shuffle arithmetic.

[Independent controls](../../../methods/independent-controls.md) cover all15 tensors/operators/colour convolutions and30 paired original-domain polynomial identities. Six adjacent Gamma structures additionally have new two-radius Cauchy coefficient controls. Integer-epsilon checks do not determine finite terms of divergent sources; the other nine integrals do not gain a new independent original-multidimensional finite check through those controls. Numerical residuals are empirical, not strict enclosures.

Use [reproduce.md](reproduce.md) for the shared batch. Original sources, Gamma rules and evaluated functions retain their own evidence scopes. Positive real z/local meromorphic epsilon continuation is the scope; no new complex i0 or z=-1 connection, internal Pi2 literature proof or formal integral proof is claimed. Formal=not_run.
