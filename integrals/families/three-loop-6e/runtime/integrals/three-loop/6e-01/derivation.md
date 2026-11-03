# 6e-01: the empty three-loop exchange integral

This article follows this labelled graph from its ordered C/J root to all seven Laurent positions of the physical I. Fresh original-chain evaluation and Gamma expansion; the completed P1 calculation is retained as prior evidence. The shared [empty derivation](../../../methods/empty.md) supplies the integration and endpoint argument; [common normalization](../../../methods/common-normalization.md) explains colour projection, Fourier convention and crossing.

## Original graph and kernel

The source has six ordered Wilson insertions and three freely integrated gluon momenta. Each momentum-space factor is the source eikonal contraction of the Feynman-gauge propagator with exp(-i l.X); its actual labelled separation is retained below. The Fourier/source conventions are reused, and no new FeynCalc calculation is asserted.

| Momentum | Endpoint Lorentz labels | Source separation/phase argument |
| --- | --- | --- |
| l1 | mu1, mu4 | `(1-s31)*v3 + s41*v4` |
| l2 | mu2, mu5 | `(1-s32)*v3 + s42*v4` |
| l3 | mu3, mu6 | `(1-s33)*v3 + s43*v4` |

The exact full-colour source, graph, native cube presentation and original ordered root are separately bound in [input.json](input.json). The original J is the integral of

```text
(z1^3) * ((1-s31)*s41*z1)^(-1+eps) * ((1-s32)*s42*z1)^(-1+eps) * ((1-s33)*s43*z1)^(-1+eps)
```

The prescribed outer-to-inner order is s33 → s43 → s32 → s42 → s31 → s41. This is part of the root interpretation.

| Variable | Lower bound | Upper bound |
| --- | --- | --- |
| s33 | 0 | 1 |
| s43 | 0 | 1 |
| s32 | 0 | s33 |
| s42 | 0 | s43 |
| s31 | 0 | s32 |
| s41 | 0 | s42 |

For positive z all propagator bases are positive in the interior. The family derivation states the convergence/continuation argument allowing its grouped integrations; the source-to-cube Jacobian is independently reconstructed in [checks.json](checks.json). It is not silently absorbed into the ordered root.

## Reduction and physical assembly

This orbit has a1=1/2, orbit multiplicity=4, and kappa=-1/8. Its original full colour is `((1 + CA^2)*CF)/(4*CA^2)`; all three fixed-basis coordinates are `['1/2', '-3/2', '1']`. This is the fixed-basis reconstruction, with Euler compensation once. It is not a large-N projection. The common physical convention is wilson-qcd-action-path-fourier/1.0.

The original orbit-specific reduction is recorded explicitly below. The symbols and mother definitions are given by the shared method; active primitive requests, original J and complete consumed order inventory are also available in [subintegrals.json](subintegrals.json).

```json
{
  "r": 3,
  "chains": [
    {
      "variables_low_to_high": [
        "s31",
        "s32",
        "s33"
      ],
      "word": "BBB"
    },
    {
      "variables_low_to_high": [
        "s41",
        "s42",
        "s43"
      ],
      "word": "AAA"
    }
  ],
  "K": {
    "language": "WolframLanguage/InputForm",
    "expression": "(1/36)/eps^6",
    "latex": "\\frac{1}{36 \\varepsilon^{6}}"
  },
  "I_compact": {
    "language": "WolframLanguage/InputForm",
    "expression": "-1/288*G^3/eps^6",
    "latex": "- \\frac{G^{3}}{288 \\varepsilon^{6}}"
  }
}
```

The active implementation recognizes the original endpoint chains, integrates or replays the required regular moments and convolves the actual pole, Euler/Gamma and crossing factors. It never uses the final reference coefficient table to generate the candidate. A fixed epsilon jet is truncated only after the physical order budget is satisfied.

## Complete Laurent coefficients

Write I=sum(c_p(z) epsilon^p,p=-6,...,0)+O(epsilon) at fixed positive z. M and N are the original unnormalised moments; N_kj=M_k,j+1-rho M_kj with rho=log(z). All generator HLOG definitions, including tau_D,4 and E when present, are in [the common dictionary](../../../references/basis-and-conversions.json). No factorial-normalised M is introduced.

| epsilon power | Function-family coefficient | HLOG coordinate count |
| ---: | --- | ---: |
| -6 | $-\frac{1}{288}$ | 1 |
| -5 | $0$ | 0 |
| -4 | $-\frac{1}{192} \zeta_{2}$ | 1 |
| -3 | $-\frac{1}{288} \zeta_{3}$ | 1 |
| -2 | $-\frac{19}{3840} \zeta_{2} \zeta_{2}$ | 1 |
| -1 | $-\frac{1}{192} \zeta_{2} \zeta_{3}-\frac{1}{480} \zeta_{5}$ | 2 |
| 0 | $-\frac{1}{576} \zeta_{3} \zeta_{3}-\frac{631}{161280} \zeta_{2} \zeta_{2} \zeta_{2}$ | 2 |

[result.json](result.json) gives every rational/zeta/HLOG coordinate and this function-family representation. HLOG uses dt/(t-a), base point0, endpointz, alphabet0,-1. Products of generators are expanded with a separately implemented rational shuffle, and every representation back-substitution has zero residual. The finite coefficient contains no unsolved parameter integral.

## Sharing, checks and reproduction

The recorded sharing statement is `own_orbit_shared_family_method` with representative `6e-01` and factor `1`. Its exact scope is: complete physical I under the documented positive-real-z normalization and convergence/continuation domain No equality/proportionality of the original J is inferred from a physical-I relation. Full-epsilon sharing is reused from the family proof, rather than inferred from seven matching coefficients. The six batch alias cases are described in [the method index](../../../methods/README.md).

All seven positions, full-colour reconstruction, the native/ordered polynomial bridge and C/I normalization pass the actual batch checks. [The independent-control matrix](../../../methods/independent-controls.md) separates representative polynomial/Jacobian checks from direct finite-source numerical checks; its scope is not enlarged to untested finite directions. Formal state remains not_run.

Run the common command in [reproduce.md](reproduce.md). Historical source methods are preserved byte-for-byte as evidence; their absolute paths and old status labels do not become current commands or proofs.
