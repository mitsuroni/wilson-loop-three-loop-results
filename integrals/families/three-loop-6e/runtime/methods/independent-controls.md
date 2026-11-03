# What the independent controls establish

The batch replays existing analytic primitives and reconstructs each physical
integral. An agreement with the historical table checks that reconstruction,
but is not a second derivation of the primitives. The independent controls in
[`independent_controls.py`](../scripts/independent_controls.py) therefore start
from the original ordered propagator product and its dependent bounds. They do
not import the Wilson evaluator, its parser, its HLOG arithmetic, or any raw CAS
output. Only the frozen candidate supplied on the command line is compared.

For every one of the 39 sources, a separate restricted arithmetic reader checks
the propagator product against the expanded integrand and verifies that each
bound refers only to available outer variables. It reconstructs the supplied
full color vector under `CA^2=1+2 CA CF`, checks the single Euler compensation,
and obtains

\[
\kappa=\frac{\mathrm{multiplicity}\,a_1}{2(-2)^3}.
\]

This source-derived value is compared with the candidate normalization. These
checks connect metadata and arithmetic to the original C/J; they are not a
formal proof of the initial momentum-space/Fourier rules.

## Two exact integration methods on representative sources

At positive integer epsilon, the original product is a polynomial. The script
integrates that polynomial in the declared inner-to-outer order, retaining every
dependent upper bound. Independently, it constructs a triangular unit-cube map
directly from those bounds,

\[
s_i=\rho_i(s_{<i})+
  \bigl(h_i(s_{<i})-\rho_i(s_{<i})\bigr)x_i,
\]

checks the full determinant against the product of interval lengths, expands the
transformed polynomial, and uses the exact monomial moments
`integral_0^1 x^n dx=1/(n+1)`. This second method does not trust the archived
`input_export.s_to_x` map or Jacobian. The resulting symbolic polynomials in z
are compared exactly.

| Mathematical route | Representative | Exact controls |
| --- | --- | --- |
| Empty / adjacent-only | 6e-01 | epsilon=1 and 2 |
| Single-opposite | 6e-04 | epsilon=1 and 2 |
| Same-channel-two | 6e-10 | epsilon=1 and 2 |
| Same-channel-three | 6e-20 | epsilon=1 and 2 |
| Double-opposite | 6e-32 | epsilon=1 and 2 |
| Mixed-opposite | 6e-37 and 6e-38 | epsilon=1 and 2 for each |

There are fourteen paired exact checks. They test domain and polynomial
manipulations at the two specified regulator values; they do not determine the
Laurent coefficients near epsilon=0. Coverage of one representative is not an
independent calculation of the other members' finite coefficients.

## Finite terms that can be sampled directly

For the holomorphic all-opposite sources 6e-20/21/22/37/38, the accepted endpoint
analysis permits epsilon=0 in the original integral. The script confirms that
each denominator has one of the two opposite-channel forms and samples the
original six-dimensional domain by randomized, shifted Halton points. It uses
both reciprocal channels and the source-derived normalization:

\[
c_0(z)=\frac{\kappa}{2}
\int_{\text{original ordered domain}}
       \bigl(j(s,0,z)+j(s,0,1/z)\bigr)\,d^6s.
\]

Each of these five comparisons uses z=1/2, eight shifts, and, by default, 16384
points per shift. The candidate HLOG expression is evaluated with a separate
log-polynomial Taylor-series implementation of the kernels `dt/t` and
`dt/(1+t)`. Its elementary weight-one and dilogarithm values are checked against
mpmath. This comparison starts from the actual source, rather than from the
truncated Laurent expansion at a positive integer epsilon.

The QMC tolerance is `6 empirical standard errors + 0.005 abs(candidate)`.
Both the standard error and observed series convergence are numerical controls,
not strict error enclosures or guaranteed accuracy. No original-domain finite
comparison is claimed here for the divergent adjacent/single/two-opposite
routes. Their finite-term evidence continues to include the explicit analytic
primitive replay, normalization checks and exact historical-table comparison.

Six negative controls reject an unavailable dependent bound, duplicated Euler
compensation, an incorrect projected color, direct epsilon-zero integration of
a pole source, deletion of an ordered upper bound, and promotion of the 37/38
finite equality to an all-epsilon identity. The last test cancels the common
Gamma factor before evaluating the exact rational source integrals at epsilon=2;
it never numerically evaluates Gamma at a pole.

From the package directory, run against a **newly generated and frozen** batch:

```powershell
D:/py/python.exe -X utf8 -B scripts/independent_controls.py --input engine/benchmarks/2ne/prepared/6e.json --results PATH_TO_NEW_GENERATED_JSON --out PATH_TO_NEW_CONTROL_JSON
```

The output records the input/candidate hashes, every actual result, the coverage
matrix and limitations. It refuses to overwrite an existing report.
`formal_state=not_run` applies throughout.
