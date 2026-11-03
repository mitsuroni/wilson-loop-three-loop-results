# 6e-35: the empty three-loop exchange integral

This article follows this labelled graph from its ordered C/J root to all seven Laurent positions of the physical I. Fresh original-chain evaluation and Gamma expansion; the completed P1 calculation is retained as prior evidence. The shared [empty derivation](../../../methods/empty.md) supplies the integration and endpoint argument; [common normalization](../../../methods/common-normalization.md) explains colour projection, Fourier convention and crossing.

## Original graph and kernel

The source has six ordered Wilson insertions and three freely integrated gluon momenta. Each momentum-space factor is the source eikonal contraction of the Feynman-gauge propagator with exp(-i l.X); its actual labelled separation is retained below. The Fourier/source conventions are reused, and no new FeynCalc calculation is asserted.

| Momentum | Endpoint Lorentz labels | Source separation/phase argument |
| --- | --- | --- |
| l1 | mu1, mu5 | `-((1-s41)*v4 + s11*v1)` |
| l2 | mu2, mu4 | `(1-s21)*v2 + s32*v3` |
| l3 | mu3, mu6 | `(1-s31)*v3 + s42*v4` |

The exact full-colour source, graph, native cube presentation and original ordered root are separately bound in [input.json](input.json). The original J is the integral of

```text
(1^2*z1) * ((1-s41)*s11*1)^(-1+eps) * ((1-s21)*s32*1)^(-1+eps) * ((1-s31)*s42*z1)^(-1+eps)
```

The prescribed outer-to-inner order is s11 → s21 → s32 → s42 → s31 → s41. This is part of the root interpretation.

| Variable | Lower bound | Upper bound |
| --- | --- | --- |
| s11 | 0 | 1 |
| s21 | 0 | 1 |
| s32 | 0 | 1 |
| s42 | 0 | 1 |
| s31 | 0 | s32 |
| s41 | 0 | s42 |

For positive z all propagator bases are positive in the interior. The family derivation states the convergence/continuation argument allowing its grouped integrations; the source-to-cube Jacobian is independently reconstructed in [checks.json](checks.json). It is not silently absorbed into the ordered root.

## Reduction and physical assembly

This orbit has a1=1/4, orbit multiplicity=4, and kappa=-1/16. Its original full colour is `CF/(4*CA^2)`; all three fixed-basis coordinates are `['1/4', '-1', '1']`. This is the fixed-basis reconstruction, with Euler compensation once. It is not a large-N projection. The common physical convention is wilson-qcd-action-path-fourier/1.0.

The original orbit-specific reduction is recorded explicitly below. The symbols and mother definitions are given by the shared method; active primitive requests, original J and complete consumed order inventory are also available in [subintegrals.json](subintegrals.json).

```json
{
  "r": 1,
  "chains": [
    {
      "variables_low_to_high": [
        "s11"
      ],
      "word": "A"
    },
    {
      "variables_low_to_high": [
        "s21"
      ],
      "word": "B"
    },
    {
      "variables_low_to_high": [
        "s31",
        "s32"
      ],
      "word": "BA"
    },
    {
      "variables_low_to_high": [
        "s41",
        "s42"
      ],
      "word": "BA"
    }
  ],
  "K": {
    "language": "WolframLanguage/InputForm",
    "expression": "(R2 - 1)^2/eps^6",
    "latex": "\\frac{\\left(R_{2} - 1\\right)^{2}}{\\varepsilon^{6}}"
  },
  "I_compact": {
    "language": "WolframLanguage/InputForm",
    "expression": "-1/16*G^3*(R2 - 1)^2*Cosh[(1/2)*ell*eps]/(eps^6*Cosh[(3/2)*ell*eps])",
    "latex": "- \\frac{G^{3} \\left(R_{2} - 1\\right)^{2} \\cosh{\\left(\\frac{\\ell \\varepsilon}{2} \\right)}}{16 \\varepsilon^{6} \\cosh{\\left(\\frac{3 \\ell \\varepsilon}{2} \\right)}}"
  }
}
```

The active implementation recognizes the original endpoint chains, integrates or replays the required regular moments and convolves the actual pole, Euler/Gamma and crossing factors. It never uses the final reference coefficient table to generate the candidate. A fixed epsilon jet is truncated only after the physical order budget is satisfied.

## Complete Laurent coefficients

Write I=sum(c_p(z) epsilon^p,p=-6,...,0)+O(epsilon) at fixed positive z. M and N are the original unnormalised moments; N_kj=M_k,j+1-rho M_kj with rho=log(z). All generator HLOG definitions, including tau_D,4 and E when present, are in [the common dictionary](../../../references/basis-and-conversions.json). No factorial-normalised M is introduced.

| epsilon power | Function-family coefficient | HLOG coordinate count |
| ---: | --- | ---: |
| -6 | $0$ | 0 |
| -5 | $0$ | 0 |
| -4 | $0$ | 0 |
| -3 | $0$ | 0 |
| -2 | $-\frac{1}{16} \zeta_{2} \zeta_{2}$ | 1 |
| -1 | $\frac{1}{4} \zeta_{2} \zeta_{3}$ | 1 |
| 0 | $-\frac{1}{4} \zeta_{3} \zeta_{3}-\frac{1}{8} \zeta_{2} \zeta_{2} M_{0,1}(z)-\frac{93}{160} \zeta_{2} \zeta_{2} \zeta_{2}$ | 3 |

[result.json](result.json) gives every rational/zeta/HLOG coordinate and this function-family representation. HLOG uses dt/(t-a), base point0, endpointz, alphabet0,-1. Products of generators are expanded with a separately implemented rational shuffle, and every representation back-substitution has zero residual. The finite coefficient contains no unsolved parameter integral.

## Sharing, checks and reproduction

The recorded sharing statement is `own_orbit_shared_family_method` with representative `6e-35` and factor `1`. Its exact scope is: complete physical I under the documented positive-real-z normalization and convergence/continuation domain No equality/proportionality of the original J is inferred from a physical-I relation. Full-epsilon sharing is reused from the family proof, rather than inferred from seven matching coefficients. The six batch alias cases are described in [the method index](../../../methods/README.md).

All seven positions, full-colour reconstruction, the native/ordered polynomial bridge and C/I normalization pass the actual batch checks. [The independent-control matrix](../../../methods/independent-controls.md) separates representative polynomial/Jacobian checks from direct finite-source numerical checks; its scope is not enlarged to untested finite directions. Formal state remains not_run.

Run the common command in [reproduce.md](reproduce.md). Historical source methods are preserved byte-for-byte as evidence; their absolute paths and old status labels do not become current commands or proofs.
