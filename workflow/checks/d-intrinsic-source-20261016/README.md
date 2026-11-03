# Intrinsic rectangle source for C and D

2026-10-16. Mathematical result: [continuous Chinese derivation](../../../benchmarks/2ne/theory/generation/d-intrinsic-source-20261016/derivation.md) and [independent endpoint proof](../../../benchmarks/2ne/theory/generation/d-intrinsic-source-20261016/independent-kernel.md).

The result is a structural reorganization of the existing ordered integral,
using an independently evaluated single-propagator rectangle. Its complete
boundary equals the original M generator. C and D use the same rectangle;
D has indexed, crossing-even two-dimensional moments with rational constant
assembly coefficients and no external log(rho) convolution. The remaining
outer integration dimension has not decreased. Minimality, reduction to a
small extension of M/N/P/Q/R, and a new canonical Hlog D4 are not claimed.

## Checks actually completed

- [check_algebra.py](check_algebra.py) and [actual result](algebra-verification.json): 40 exact SymPy checks and five negative controls. These cover the coordinate/Jacobian connections, endpoint products, crossing, removable internal lines, twelve integer-parameter rectangle integrals, and the complete D0–D4 moment coefficient rule. The initial interchanged line-sign assertions and their correction are documented in [attempt-notes.md](attempt-notes.md).
- [Native Wolfbook source](wolfbook/algebra.wl) and [session export](wolfbook/session-report.md): actual one-off MCP evaluation, 19 positive checks and four negatives, ASSERT PASS and no kernel messages. This is a separate algebra check, not a formal proof or a full integration. The report was exported by Wolfbook and copied byte-identically; no success response was reconstructed by hand.
- [Numerical route](numeric/README.md) and [actual result](numeric/verification.json): seven unexpanded D comparisons, including epsilon zero and a negative epsilon; ten D0–D4 coefficient comparisons; twelve generic rectangle comparisons and six integer rectangle controls; fourteen omitted-scale negatives. The full J-based two-dimensional route is compared with the original four-dimensional ordered route. The coefficient checks use separately assembled intrinsic logs but share the earlier three-dimensional interval grid, so they are normalization checks rather than wholly independent integration domains.
- [Interior-line controls](numeric/line-verification.json): thirty exact/near-line numerical comparisons, with finite extensions of both divided differences.

The original noninteger target values were frozen first in
[original-reference.json](numeric/original-reference.json). Existing original
D0–D4 values are explicitly reused with their source hashes. Numerical
refinements give empirical agreement, not certified intervals. The new
two-dimensional representation does not alter any old coefficient file.

[Individual-moment controls](numeric/moments-verification.json) additionally
evaluate all 22 index positions needed through D4 at z=.7 and its reciprocal:
22 individual crossing checks, ten coefficient reconstructions, six M-boundary
checks and two shared-C comparisons pass. These are 22 evaluated moments,
not a claim of 22 linearly independent masters. Two wrong-density C controls
are rejected. The C0 original-coordinate check uses its finite lower-node
grid and is weaker than the comparison to the exact old identity H0^3/18;
see the numerical README for that distinction and preserved endpoint failures.

The final [audit](verification.json) checks the reports and recorded bindings.

## Reproduction and scope

Use `D:/py/python.exe -X utf8`. Run `check_algebra.py` for the exact local
identities, and see the numerical README for its separate scripts. Reports
are immutable outputs: use a fresh revision directory to rerun producers.
`verify.py` is a read-only evidence/links audit except for writing its new
audit JSON. It checks the recorded hashes rather than repeating integrations.

No Wilson/Lean implementation, accepted input root, color projection,
physical I normalization, existing result or historical receipt was changed.
No formal chain was admitted. `formal=not_run`. No Git commit was requested.
