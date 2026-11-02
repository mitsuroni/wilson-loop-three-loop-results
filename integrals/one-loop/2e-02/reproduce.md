# Reproduce 2e-02

Run the same small family calculation from the repository root:

```text
python integrals/one-loop/reproduce.py --output-dir one-loop-replay --order 3 --dps 40
```

The new directory contains `2e-02/result.json` and the shared source,
symbolic and numerical records. It must be outside the frozen scientific
input/result trees. The Python dependencies are pinned in
[requirements.txt](../requirements.txt). No archived solver, external CAS,
network operation or original Wilson workspace is invoked.

The replay checks the actual source phase, eikonal numerator, full opposite
denominator, original order and native/prepared C connection. It
differentiates the actual affine-denominator primitive, verifies the source
reflection and moment variable change, then generates the exact Gamma/M
series. M01 and M03 compactification comes from an explicit Mellin generating
formula, not from fitting coefficients. The new result is persisted before
opening the original r002 coefficient table.

Numerical controls independently evaluate:

- the original ordered two-dimensional J at a noninteger epsilon;
- the source-derived one-dimensional removable difference quotient;
- the standard all-epsilon Meijer G answer;
- the original M01, M03 and M11 integrals, including their removable point;
- epsilon-Cauchy coefficients of the complete normalized I through epsilon3.

The two shape points2/3 and3/2 form one reflection orbit. We do not count
them as two independent shapes. The actual source-reflection check is also
symbolic. These controls separate integral evaluation, function evaluation
and exact coefficient algebra; their empirical discrepancies are not strict
interval bounds.

The [published check record](checks.json) and
[complete run record](../validation/checks.json) retain the actual successful
execution. The original r002 result has coefficients only through finite
order. Its comparison is therefore a post-computation holdout at those
orders, not the source of the newly generated positive-order coefficients.

The normalized I in [result.json](result.json) includes its m=2 orbit,
crossing and one Euler compensation, and removes CF, lambda and S1. Its
moment powers have no factorial normalization. The source-domain theorem
is for positive z and local epsilon continuation; the full Mellin-Barnes
answer is initially identified in its common convergence strip. No general
physical complex-kinematics extension or new formal integral proof is
asserted. `formal=not_run`.
