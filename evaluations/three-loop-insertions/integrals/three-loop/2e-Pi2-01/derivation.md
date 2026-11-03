# 2e-Pi2-01: Pi2-adjacent

This labelled source belongs to the fifteen-integral insertion batch. We derive the shared Pi2 geometry I_geom, before its separate external operator with its own source normalization, then retain the complete insertion-operator colour convolution. The [shared method](../../../methods/chains-pi2.md) gives the continuous endpoint/IBP/Gamma argument; [the method index](../../../methods/README.md) fixes the physical and function conventions.

## Original momentum insertion and parameter source

The momentum source contains D0 Pi1 D0 with one bare exchange for the four-endpoint branch, D0 Pi1 D0 Pi1 D0 for a serial chain, or D0 Pi2 D0 for a single two-loop1PI insertion. Here D0 is the Feynman-gauge gluon propagator and Pi is the complete transverse effective self-energy operator, with its original phases. The adopted internal operator is an input; the present packet does not claim a new internal bubble/two-loop-self-energy diagram integration.

[input.json](input.json) binds the actual graph, raw C/J, full-colour native preparation and selected reference separately. The exact original ordered kernel is

```text
(2*s41*z1^2 + eps*s41*z1^2 - 2*s31*s41*z1^2 - eps*s31*s41*z1^2)*(s41*z1 - s31*s41*z1)^(-2+3*eps)/(2*(1+2*eps))
```

Its outer-to-inner order is s31 → s41.

| Parameter | Lower bound | Upper bound |
| --- | --- | --- |
| s31 | 0 | 1 |
| s41 | 0 | 1 |

The full transverse numerator is retained in this kernel, including its longitudinal/derivative part. Each polynomial and tensor identity is checked against the source, not inferred from an integer-regulator sample. Polynomial controls then compare this actual order with a separately constructed cube/Jacobian integration; their scope is described below.

The original multiplicative C components are:

```json
{
  "wilson_color": "CF",
  "self_energy_scalar": "(CA^2)*(((2*eps^2-4*eps-1)/(2*eps^3*(2*eps-3)))+((2*eps-1)*(6*eps^5-eps^4-44*eps^3+79*eps^2-30*eps-2)/(2*eps^3*(eps+1)*(2*eps-3)*(3*eps-2)*(3*eps-1)))*(Gamma[1+2*eps]*Gamma[1-2*eps]^2/(Gamma[1+eps]^2*Gamma[1-eps]*Gamma[1-3*eps])))+(CA*TF*nf)*((2*(eps-1)*(2*eps^2-3*eps+2)/(eps^3*(2*eps-3)))+(-2*(eps-1)*(2*eps-1)*(2*eps^4-11*eps^3-10*eps^2+17*eps-4)/(eps^3*(eps+1)*(2*eps-3)*(3*eps-2)*(3*eps-1)))*(Gamma[1+2*eps]*Gamma[1-2*eps]^2/(Gamma[1+eps]^2*Gamma[1-eps]*Gamma[1-3*eps])))+(CF*TF*nf)*((-4*(eps-1)*(2*eps^2-eps+2)/(eps^3*(2*eps-3)))+(8*(eps-1)*(2*eps-1)^2*(eps^2-2*eps+2)/(eps^3*(2*eps-3)*(3*eps-2)*(3*eps-1)))*(Gamma[1+2*eps]*Gamma[1-2*eps]^2/(Gamma[1+eps]^2*Gamma[1-eps]*Gamma[1-3*eps])))",
  "coupling": "(alpha_s/Pi)^3",
  "gamma_fourier": "(Pi^eps*Gamma[1-eps])^3",
  "gamma_residual": "Gamma[1+eps]^2*Gamma[1-eps]*Gamma[1-3*eps]/(Gamma[1-2*eps]^2*Gamma[1+2*eps])",
  "rational": "1/(32*(1-2*eps)^2)",
  "scale": "(-2*v14+I0)^(3*eps)",
  "euler": "exp(-3*EulerGamma*eps)"
}
```

Their phase ledger and full trace remain in the bound native source. Its normalization has assembly_sign=1, one Euler compensation and no native sign conversion. The scalar I uses its documented orbit/crossing extraction and, for the four-endpoint Pi1 branch, the connected fixed-basis a1. The unprojected raw trace is not identified with that projected object.

## Reduction, functions and required orders

The route is `Pi2-adjacent`, shared with `2e-Pi1+Pi1-01` while retaining this source's own kernel. Its reduction and actual epsilon budgets are explicitly recorded in [subintegrals.json](subintegrals.json):

```json
{
  "adjacent": "Two original endpoint integrals 1/(3eps), multiplied in prescribed order.",
  "opposite": "H0(3eps,z) plus complete mixed-derivative boundary (1+z^(3eps))/(6eps(1+2eps)).",
  "H0_series": "M01+(3eps/2)M11+9eps^2(M21/8+M03/24)+O(eps^3)",
  "regular_geometry_budget": 4,
  "operator_multiplier_budget": 4,
  "Pi2_Gamma_ratio_budget_before_cancellation": 5,
  "H0_jet_order": 2,
  "evaluated_M_reuse": [
    {
      "path": "workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_1.txt",
      "sha256": "a5f017e08433739b41a21c014c158eeeff300af551d7143096e8f519a2802457"
    },
    {
      "path": "workflow/checks/single-opposite/runs/20260916T090538777022Z/M_1_1.txt",
      "sha256": "baacaae02ee6de07ea9dd30e3f7a47efbab82aeb57bcfdb5aed548de4a1e1b6b"
    },
    {
      "path": "workflow/checks/single-opposite/runs/20260916T090538777022Z/M_2_1.txt",
      "sha256": "e455df6fe255b79c0fa8ac6b9b01e4a20d813ed6ad41bd951e0c46473c757ca6"
    },
    {
      "path": "workflow/checks/single-opposite/runs/20260916T090538777022Z/M_0_3.txt",
      "sha256": "d3b8c7cd93aeeeb25fe1a8e6f7b4b0e0141142b9e133d3fe9070492371827269"
    }
  ]
}
```

The method derives these finite jets from the original kernel, preserving the full boundary term in the opposite-edge mixed derivative. The evaluator then uses an exact Gamma/rational jet and the admitted original M/HLOG dictionary. Existing evaluated moment data is reused as such; no final integral table is used for generation.

rho means log(z), not the physical average scale L. M_kj denotes the original no-factorial moment; N_kj=M_k,j+1-rho M_kj. The chain display abbreviations M01,M11,M21,M03 mean M_0,1,M_1,1,M_2,1,M_0,3. The separate source labels rho_id_0_0_0 and rho_id_0_0_1 denote the already evaluated T coefficient functions, rather than rho itself: the first is M_0,1 squared/6; both complete HLOG definitions are in the [bound structural dictionary](../../../snapshots/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json). These are explicit evaluated functions, not names for remaining unsolved integrals.

## Complete Laurent coefficients and external colours

At fixed positive z, the scalar series is sum(c_p(z)epsilon^p,p=-6,...,0)+O(epsilon). In this entry the scalar is the shared Pi2 geometry I_geom, before its separate external operator. All seven positions are explicit, including established zeros.

| Power | Exact scalar coefficient | HLOG coordinates |
| ---: | --- | ---: |
| -6 | $0$ | 0 |
| -5 | $0$ | 0 |
| -4 | $\frac{1}{576}$ | 1 |
| -3 | $\frac{5}{1152}$ | 1 |
| -2 | $\frac{\pi^{2}}{2304} + \frac{1}{64}$ | 2 |
| -1 | $\frac{5 \pi^{2}}{4608} + \frac{7 \zeta\left(3\right)}{576} + \frac{5}{144}$ | 3 |
| 0 | $\frac{67 \pi^{4}}{276480} + \frac{35 \zeta\left(3\right)}{1152} + \frac{\pi^{2}}{256} + \frac{7}{72}$ | 4 |

The complete per-order coloured tables are in [result.json](result.json). Each row below gives a monomial outside its scalar multiplier; the last column is its complete finite coefficient after epsilon convolution.

| External monomial | Exact multiplier of the scalar I | Coloured finite coefficient |
| --- | --- | --- |
| `CA^2*CF` | `4*eps**2*((2*eps**2 - 4*eps - 1)/(2*eps**3*(2*eps - 3)) + (2*eps - 1)*(6*eps**5 - eps**4 - 44*eps**3 + 79*eps**2 - 30*eps - 2)*gamma(1 - 2*eps)**2*gamma(2*eps + 1)/(2*eps**3*(eps + 1)*(2*eps - 3)*(3*eps - 2)*(3*eps - 1)*gamma(1 - 3*eps)*gamma(1 - eps)*gamma(eps + 1)**2))` | $- \frac{1187401}{559872} - \frac{2941 \pi^{2}}{62208} + \frac{\pi^{2} \zeta\left(3\right)}{576} + \frac{71 \pi^{4}}{276480} + \frac{7 \zeta\left(5\right)}{144} + \frac{605 \zeta\left(3\right)}{5184}$ |
| `CA*CF*TF*nf` | `4*eps**2*((2 - 2*eps)*(2*eps - 1)*(2*eps**4 - 11*eps**3 - 10*eps**2 + 17*eps - 4)*gamma(1 - 2*eps)**2*gamma(2*eps + 1)/(eps**3*(eps + 1)*(2*eps - 3)*(3*eps - 2)*(3*eps - 1)*gamma(1 - 3*eps)*gamma(1 - eps)*gamma(eps + 1)**2) + (2*eps - 2)*(2*eps**2 - 3*eps + 2)/(eps**3*(2*eps - 3)))` | $\frac{223 \pi^{4}}{207360} + \frac{\pi^{2} \zeta\left(3\right)}{72} + \frac{199 \pi^{2}}{7776} + \frac{7 \zeta\left(5\right)}{18} + \frac{73 \zeta\left(3\right)}{216} + \frac{97655}{139968}$ |
| `CF^2*TF*nf` | `4*eps**2*((4 - 4*eps)*(2*eps**2 - eps + 2)/(eps**3*(2*eps - 3)) + (2*eps - 1)**2*(8*eps - 8)*(eps**2 - 2*eps + 2)*gamma(1 - 2*eps)**2*gamma(2*eps + 1)/(eps**3*(2*eps - 3)*(3*eps - 2)*(3*eps - 1)*gamma(1 - 3*eps)*gamma(1 - eps)*gamma(eps + 1)**2))` | $- \frac{557 \zeta\left(3\right)}{648} - \frac{7 \zeta\left(5\right)}{9} - \frac{\pi^{2} \zeta\left(3\right)}{36} - \frac{\pi^{4}}{324} + \frac{23 \pi^{2}}{864} + \frac{7993}{3888}$ |

The series of each external multiplier is retained through the actual scalar pole order before multiplication. For Pi2 its Gamma ratio is expanded to the higher intermediate order needed before cancellation; the full native A+B/U operator and the direct Gamma route agree. A serial-Pi1 operator is not substituted for it. Other connected colour channels absent from the complete polynomial are zero in the stated projected/connected convention; this does not erase the unprojected raw-source trace.

## Evidence and reproduction

The candidate was frozen before the final reference tables were opened. [The per-orbit reference check](checks.json) compares every scalar coefficient with the adopted current object and, where stored, each coloured coefficient. The whole batch also independently back-substitutes all343 scalar/colour function expressions into complete HLOG/period coordinates with rational shuffle arithmetic.

[Independent controls](../../../methods/independent-controls.md) cover all15 tensors/operators/colour convolutions and30 paired original-domain polynomial identities. Six adjacent Gamma structures additionally have new two-radius Cauchy coefficient controls. Integer-epsilon checks do not determine finite terms of divergent sources; the other nine integrals do not gain a new independent original-multidimensional finite check through those controls. Numerical residuals are empirical, not strict enclosures.

Use [reproduce.md](reproduce.md) for the shared batch. Original sources, Gamma rules and evaluated functions retain their own evidence scopes. Positive real z/local meromorphic epsilon continuation is the scope; no new complex i0 or z=-1 connection, internal Pi2 literature proof or formal integral proof is claimed. Formal=not_run.
