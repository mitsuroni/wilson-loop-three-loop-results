# Intrinsic D source: independent numerical checks

The target is the bare geometric `D = R231`, with no Gamma, color, physical
cosh factor, or external rho power. The new source is implemented directly in
[intrinsic.py](intrinsic.py). All previous research files are read-only.

The two factors

`w_s = L+V = log(h/(1-s))`,
`w_t = L-V = log(z*t/h)`

are divided separately by `s-alpha` and `t-alpha`, using their continuous
limits. This evaluates the proposed cosh difference without subtracting two
nearly equal coshes across the internal lines.

## Completed checks

[original-reference.json](original-reference.json) freezes four fresh values
of the original ordered integral at epsilon `.13/.31` and z `.4/1.4` before
the candidate was numerically evaluated. The original calculation retains
four ordered coordinates after two analytic interval primitives. It also
binds and explicitly reuses the earlier independent D0 through D4 values at
z `.4/.7`; those coefficients were not integrated again in this task.

[verification.json](verification.json) records the actual run of [check.py](check.py):

- Seven full-epsilon comparisons of the new two-dimensional J expression
  with the original four-coordinate integral, including epsilon `-.05` and
  `0`, have maximum absolute difference `1.93e-9`.
- Ten coefficient comparisons, D0 through D4 at z `.4/.7`, have maximum
  absolute difference `1.37e-9` against the frozen original 512-node values.
  The intrinsic source was evaluated at 96, 160, and 256 nodes.
- Twelve generic rectangle J formulas agree with a separate one-dimensional
  exact interval primitive. These include negative, zero, and positive
  epsilon and check the eta reflection symmetry.
- Six integer-epsilon rectangle integrals reproduce the exact epsilon-one
  and epsilon-two polynomials.
- Fourteen negative controls omit half of `log((1-t)*s)` from the intrinsic
  exponent; all are rejected.

[line-verification.json](line-verification.json) records thirty further
pointwise comparisons at and near both removable alpha lines. The maximum
absolute difference is `9.90e-13`. That check intentionally shares the
rectangle value and tests only the two power differences and their scale.

The separate bounded extension [moments.py](moments.py) and
[moments-verification.json](moments-verification.json) retain every one of the
22 individual `X_(m,n,j)` moments entering D0 through D4 at z `0.7` and its
reciprocal, at 96/160/256 nodes. Here the n-index is the raw n-th epsilon
derivative of B, so reconstruction uses `2/(m! n! (2j)!)`.

- All 22 individual crossing comparisons pass; the largest absolute
  difference at 256 nodes is `5.62e-9`.
- Both reconstructions through D4 agree with the frozen original values;
  the largest absolute difference is `1.63e-9`.
- Six noninteger one-loop boundary values of B agree with the independent
  one-dimensional M generating kernel to `1.69e-12`.
- The same B block reconstructs C at epsilon `.13`, z `.7`, agreeing with
  the untouched original ordered route to `7.88e-11`. C0 agrees with the
  known `H0^3/18` to `1.07e-14`. Its coarse original 20/24-node comparison
  is less accurate, differing by `1.46e-7` within the observed refinement.
- Two omitted-C-density negative controls are rejected.

For example, at z `0.7`, the 256-node values are

| Raw moment | Value |
| --- | ---: |
| `X_(0,0,1)` | -4.214101729444072 |
| `X_(0,1,1)` | 3.1113970146073395 |
| `X_(1,0,1)` | 11.205875712167886 |
| `X_(0,0,2)` | -2.5706661580062926 |

The first coefficient is `D1=X_(0,1,1)+X_(1,0,1)`. Full tables, empirical
refinements, and reconstruction factors are in the JSON. Two numerical
endpoint representation failures were retained in
[moments-attempt-notes.md](moments-attempt-notes.md); no old code or successful
report was changed to repair them.

## Coefficient extraction and independence

The universal rectangle is integrated in one coordinate first. Its exact
interval factor is

`I0 * exp(epsilon*center) * sinhc(epsilon*width/2)`.

Together with the two outer sinhc factors, every coefficient through order
four follows from the sum of the three squared widths and fourth powers.
The code expands this expression before quadrature; it does not fit epsilon
samples, import a GPL answer, or convolve a power of rho.

After the analytic rectangle coordinate is integrated, this coefficient
route uses the same three outer variables and tensor grid as the earlier
original interval route. It independently assembles the intrinsic moments
and checks the scale and coefficient transformation; it is not an independent
partition of the integration domain. The full-epsilon J expression versus
the original four-coordinate integral is the separate numerical route.

All comparisons use floating-point quadrature and observed refinements,
which are not rigorous error bounds. These finite checks do not prove
all-parameter identities or all-weight closure. `formal = not_run`.

Reproduction scripts refuse to overwrite their output reports. Run them
from a new copied check directory, preserving the relative workspace layout,
or select a new output filename rather than altering the frozen reports:

```powershell
D:/py/python.exe -X utf8 workflow/checks/d-intrinsic-source-20261016/numeric/freeze_original.py
D:/py/python.exe -X utf8 workflow/checks/d-intrinsic-source-20261016/numeric/check.py
D:/py/python.exe -X utf8 workflow/checks/d-intrinsic-source-20261016/numeric/check_lines.py
D:/py/python.exe -X utf8 workflow/checks/d-intrinsic-source-20261016/numeric/moments.py
```
