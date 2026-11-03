# Bounded moment check: boundary evaluation correction

The first execution of `moments.py` on 2026-10-16 completed both sets of
96/160/256-node moments and the reconstruction/crossing assertions, then
failed the one-loop boundary comparison at epsilon `-0.05`, z `1.4`.
The computed rectangle value was `-3.5146779169876092`, versus
`-5.284030833734253` from the one-dimensional M kernel.

This was a numerical endpoint representation error. On this boundary the
two J arguments are exactly `-z` and `-1/z`, so their product is exactly one.
Reconstructing them from rounded exponentials of eta and kappa pushed the
product just above one; the reused real J evaluator does not support that
out-of-domain value. The original strict-interior rectangle checks were
unaffected and their files remain unchanged.

The new boundary check evaluates the analytically exact expression

`-kappa^(-epsilon) / epsilon * [J_e(-z)+J_e(-1/z)-2 J_e(1)]`.

No successful report existed at the failed run, and no report was modified.
The full bounded script is rerun after this correction; its actual successful
output, if obtained, is `moments-verification.json`.

The second execution passed all six M-boundary checks, then encountered an
infinite value in the unchanged original C0 comparator. Its nested endpoint
`v+(1-v)*u` rounds to exactly one at grids of 28 nodes or more, and the original
epsilon-zero primitive then evaluates `log(0)`. Direct bounded diagnostics
gave finite C0 values at 16, 20, and 24 nodes, and infinite values at 28 and
32 nodes. The final script therefore compares C0 using the original 20/24-node
grids with their observed refinement, and independently checks the candidate
against the exact known `H0^3/18`. The nonzero epsilon original C primitive
supports that endpoint and continues to use 28/48 nodes. No old code changed.
