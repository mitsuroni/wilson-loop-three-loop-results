# 5e-3g-19: t1-empty

This source is one of the44 tree-three-gluon entries. Its evaluation family is `t1-empty` and its representative is `5e-3g-19`. The complete scalar physical I is reconstructed through epsilon^0 from the admitted evaluated primitives below. Read the continuous [family derivation](../../../methods/t123.md) together with the [common momentum/C/J/I conventions](../../../methods/README.md).

## Momentum source and the actual domain

The five Wilson insertions, four gluon propagators and one three-gluon vertex are bound by [input.json](input.json). The original Lorentz/colour slots, incidence and Fourier phases are reproduced in [source.json](source.json); they determine the momentum integrand, rather than a diagram identifier standing in for it. The field/Fourier/path-order rules are the byte-bound convention in the common method. Each edge's ordered Wilson integrations follow the source bounds.

The present native C/J record uses the stored cube map. Its complete parameter integrand is

```text
(1) * (x3^2*x4) * f1^(-1+eps) * ((alpha1*alpha2*alpha3)^(-eps)) * (1/4*alpha1*alpha2*z1^2 - 1/4*alpha1*alpha3*z1^2 - 1/4*x2*alpha1*alpha2*z1^2 + 1/4*x2*alpha1*alpha3*z1^2) * (x3*alpha1*alpha3*z1 - x2*x3*alpha1*alpha3*z1 + x3*x4*alpha1*alpha2*z1 - x2*x3*x4*alpha1*alpha2*z1)^(-(2 - 2*eps))
```

The x variables range over their stored unit intervals. The alpha measure is the ordinary simplex delta(1-alpha1-alpha2-alpha3), or alpha1 in[0,1], alpha2 in[0,1-alpha1], alpha3=1-alpha1-alpha2. The factor (alpha1 alpha2 alpha3)^(-epsilon) is retained. The stored map and Jacobian connect this cube representation to the original `J_raw_s` ordered domain; both full kernels and bounds remain in [source.json](source.json). No unrecorded permutation of an iterated divergent integral is used.

The native multiplicative C is

```json
{
  "color": "I/4*CF",
  "coupling": "(4*Pi*alphaS)^3",
  "ft_prefactor": "Gamma[1-eps]^1*Gamma[1-2*eps]/((4*Pi^(2-eps))^1*64*Pi^(4-2*eps))",
  "derivative": "2*(1-2*eps)",
  "v14_scale": "(-2*v14+I0)^(3*eps)",
  "i_factor": "I",
  "twoscale": "(-2)^(-1)",
  "euler": "exp(-3*EulerGamma*eps)"
}
```

It carries one Euler compensation. Its phase product and complete colour remain in the source. The current assembly sign is1. Under CA^2=1+2CA CF the complete colour is reduced in the ordered basis before a1 is selected. The family bridge is

```text
(1-2eps)*G/8*(J(z)+z^(3eps)*J(1/z))/(1+z^(3eps))
G = exp(-3*EulerGamma*eps)*Gamma(1-eps)*Gamma(1-2eps)
```

Thus this I includes this orbit's C4/crossing extraction; it is neither ordered J nor an individual boundary/IBP block. The source checks compare full native C/J with raw C/J and the accepted legacy bridge, without altering source signs.

## Reduction and what is actually recalculated

The route is `exact Gamma/Beta jet and shuffle arithmetic from saved evaluated primitives`. The family note supplies the original reduction, boundary terms and primitive definitions. This member consumes:

- [Z1.txt](../../../snapshots/workflow/checks/tgv-t1-empty/runs/20260917T115009990802Z/Z1.txt)
- [Z2.txt](../../../snapshots/workflow/checks/tgv-t1-empty/runs/20260917T115009990802Z/Z2.txt)
- [Z3.txt](../../../snapshots/workflow/checks/tgv-t1-empty/runs/20260917T115009990802Z/Z3.txt)
- [Z4.txt](../../../snapshots/workflow/checks/tgv-t1-empty/runs/20260917T115009990802Z/Z4.txt)

The saved files are actual historical CAS outputs, not fresh integrations. The present execution parses them and recomputes Gamma/Beta factors, reflection, products and order convolutions. For27/28/29 the saved global-IBP output is already the integrated hard coefficient; reparsing it is explicitly historical integration reuse. The final selected I JSON is held until the new candidate is fixed.

The precise order budget is:

```json
{
  "gamma_and_beta_through": 6,
  "F_through": null,
  "Z_through": 4,
  "explicit_denominator_powers": [
    6,
    5,
    4
  ],
  "remainder": "O(eps) at fixed z>0"
}
```

The complete reduction and intermediate expansions are machine-readable in [subintegrals.json](subintegrals.json). Family/reflection relationships for this member are:

```json
{
  "scope": "This representative is assembled directly from its evaluated primitives."
}
```

For the six corrected T4 members the generator uses absolute Q orders-2 through2, and prefactor orders through4. It does not execute the old relative-budget writer. Historical corrected I text is a later witness, not the raw-Q generator input.

## Complete Laurent functions

At fixed z>0, I=sum(c_p(z)epsilon^p,p=-6,...,0)+O(epsilon). Each coefficient, including every zero, is stored as complete finite HLOG coordinates in [result.json](result.json). The following named forms are the fixed reference display, admitted only after freezing and exact comparison; their dictionaries are retained and separately back-substituted.

| Power | Verified function expression or representation | HLOG coordinates |
| ---: | --- | ---: |
| -6 | `0` | 0 |
| -5 | `0` | 0 |
| -4 | `0` | 0 |
| -3 | `-(-9*zeta(3) + pi**2)/192` | 2 |
| -2 | `-(-2160*zeta(3) + 7*pi**4 + 240*pi**2)/23040` | 3 |
| -1 | `(-540*ell**2*zeta(3) + 60*pi**2*ell**2 - 37*pi**4 - 240*pi**2 + 15*pi**2*zeta(3) + 2160*zeta(3) + 3780*zeta(5))/11520` | 7 |
| 0 | `(-181440*ell**2*zeta(3) + 588*pi**4*ell**2 + 20160*pi**2*ell**2 - 12432*pi**4 - 80640*pi**2 - 737*pi**6 - 25200*pi**2*zeta(3) + 362880*zeta(3)**2 + 725760*zeta(3) + 1270080*zeta(5))/1935360` | 10 |

The HLOG endpoint is z and its letter a means dt/(t-a). Period coordinates are rational with the declared generator ordering. M uses the original no-factorial convention; rho=log z. The individual27/28/29 coefficients keep letters-1/2,-2 and periods log2,A4,A6 through Omega6; they are not forced into the old alphabet. In the exact finite combination 2c27+c28 the extra coordinates cancel. That identity alone is not a separate proof of the individual integrals.

## Checks, limits and reproduction

[checks.json](checks.json) records all seven exact fixed-current comparisons for this source. [Independent controls](../../../methods/independent-controls.md) audit original source tensors, path-ordered colour, domains and reflections, and four Gamma-only full-I finite controls. Other family primitives and numerical records retain their historical scopes. Polynomial source moments do not prove the finite coefficient of a regulated divergent integral.

The whole batch was also executed from a real copied closure with the final-reference directory absent during generation. This uses the same installed Python/SymPy and admitted historical primitives; it is not an independent CAS integration or an OS sandbox. Positive real z and local meromorphic epsilon continuation are the present scope; no new negative-one/i0 connection or formal proof is claimed. Formal=not_run.

Use [reproduce.md](reproduce.md) for the batch commands. Original local machine paths inside archived logs are provenance only; active scripts resolve from this packet.
