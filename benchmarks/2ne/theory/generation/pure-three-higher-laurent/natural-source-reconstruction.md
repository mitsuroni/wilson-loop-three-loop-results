# Natural-source coordinates through physical epsilon cubed

The [machine-readable relations](coefficients-through-eps3-natural-sources.json) express all 16 independent physical/bare coefficients of I20, I21=I22, W and Lambda at orders zero through three in independently defined source moments and their products. Candidate construction uses no stored three-kernel coefficient and no old E/tau generator. The complete physical normalization is retained.

The new source is

\[
\mathcal B_{s;tu}(a,b,c;z)=\int_0^1\frac{dx}{x}
(1-X_sx)^aJ_b(X_tx)J_c(X_ux),\qquad X_0=1,\ X_1=-z.
\]

It comes from the [complete-sector Lambda identity](lambda-one-dimensional-insertions.md), with independent elementary and J exponents. Its coefficient rule is a finite shuffle and its full parameter dependence has a single-index recurrence with fixed boundary. All 158 source moments through weight nine were generated before the target comparison.

The other sources are M/N, registered periods and normalization logarithms, independent Gauss/Euler and Mellin moments, and ordinary independent double-Beta/3F2 moments. Definitions and factorial conventions are explicit in the JSON. A plus/minus suffix always denotes a half-sum/half-difference under z to 1/z.

To read a coefficient, select for example

```text
relations["I20_3"]["terms"]
```

Each row means `rational coefficient * product(factors)`. The `generators` dictionary specifies weight, parity and source moment names. Primitive elimination-basis aliases have already been expanded to actual source names. Relations use raw parameter coefficients for F3 and BJ, while E moments include their stated derivative factorials.

| Physical coefficient | Weight | Product terms after source-alias expansion |
|---|---:|---:|
| I20/I21 order zero | 6 | 36 each |
| I20/I21 order one | 7 | 85 each |
| I20/I21 order two | 8 | 165 each |
| I20/I21 order three | 9 | 366 each |

These counts and the 303 source-coordinate entries used across the relation file describe this finite representation. They do not count independent generating functions or master integrals. The coordinate representation is not particularly compact and has large rational coefficients.

All relations were independently replayed by multiplying the actual source Hlog vectors, after expanding every basis alias. The candidate even ranks at weights six through nine are 49,99,199,399; adding either target pair at its weight leaves the rank unchanged. Without the elementary endpoint power, the tested J-product subspace still has two missing weight-nine directions. [Actual unseeded comparison](../../../../../workflow/checks/pure-three-higher-laurent-20261016/bilinear-j-runs/20261016T040428164992Z/unseeded/comparison.json), [actual source-product replay](../../../../../workflow/checks/pure-three-higher-laurent-20261016/published-natural-sources.json).

These are exact finite relations in registered rational Hlog/formal-period coordinates, without a transcendence or all-order minimality assertion. They do not establish a whole-epsilon formula for I20/I21 using BJ alone. The [separate order-four result](coefficients-eps4-gpl.md) retains the fixed complete generator and has independent original-integral validation. Formal proof status is `not_run`.
