# Reproduce 2e-01

From the repository root, with the versions in
[requirements.txt](../requirements.txt), run:

```text
python integrals/one-loop/reproduce.py --output-dir one-loop-replay --order 3 --dps 40
```

Use a fresh output directory outside data/, packages/, results/ and
integrals/. The command writes `one-loop-replay/2e-01/result.json`, together
with the opposite-exchange result, the shared checks and numerical controls.
It creates no files inside this packet and does not execute source-generation
code or any archived historical script. It needs neither Mathematica nor
Maple and never downloads a dependency or a reference.

The actual computation reads the byte-bound graph, native C/J, original
ordered input and prepared0.3 input referenced by [input.json](input.json).
It derives the geometry from the source phase and null Gram matrix, checks
native C against original C, checks the single Euler preparation, verifies
the endpoint primitive, evaluates the two endpoint factors and regenerates
the Gamma series through the required regular order5. Its regular
epsilon-Cauchy readout independently checks all coefficients of normalized I
from epsilon^-2 through epsilon3.

The file `generated-before-holdout.json` is written before the old r002
coefficient table is opened. The latter supplies only the existing pole and
finite comparison targets; it supplies no positive-order answer to the new
computation. Do not edit a reference table to make the replay pass.

The accepted run is retained in [checks.json](checks.json) and the shared
[original run record](../validation/checks.json). The full family has80
checks, including four concrete negative controls, and27 numerical controls.
Checks overlap; this count is not a count of independent physics proofs.

Compare the new `2e-01/result.json` semantically with the published
[result.json](result.json). The exact coefficients and normalization are
deterministic; run time and working-precision readouts in the check record
may vary. The public result contains the C4/crossing/Euler-normalized I,
not raw J or a scale-restored Wilson coefficient. Restore CF lambda S1 only
at the physical assembly layer.

For another finite truncation the reproducer accepts order0..3. The exact
all-epsilon Gamma expression is available independently of that truncation.
No general complex-z continuation or new integral formal proof is claimed;
numeric agreement is empirical, not a certified interval. `formal=not_run`.
