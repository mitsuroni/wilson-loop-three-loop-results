# Reproduce 2e-01

Use the repository root as working directory. Python3 with the adopted
[requirements](../../families/low-loop-runtime/integrals/one-loop/requirements.txt)
is sufficient; the local configured interpreter is D:/py/python.exe.
The entry [code/reproduce.py](code/reproduce.py) delegates to the byte-preserved
same-loop [runtime](../../families/low-loop-runtime/integrals/one-loop/reproduce.py).
It executes the whole two-entry one-loop batch; select the member 2e-01 in its output.
There are not two separate cold CAS integrations hidden behind these entry points.

```powershell
D:/py/python.exe integrals/one-loop/2e-01/code/reproduce.py --output-dir D:/wilson-agent/replay-output/2e-01-trial-new
```

Choose a genuinely new output directory outside the immutable runtime and
scientific input/result trees. To rerun, change the final directory name.
The wrapper does not overwrite source objects, accepted result files or proofs.
Do not run an old writer whose absolute output path points into a frozen packet.

## What each stage actually does

1. Input and binding: read the original ordered parameter source, prepared C/J
   and the source manifest. Preserve the order and color/normalization stage
   described in the derivation; the manifest is a binding record, not a proof.
2. Adopted mathematics: Two convergent cusp endpoint primitives give J=z^epsilon/epsilon². The C/J and crossing bridge gives I=-Z/epsilon²; the Gamma jet is required through5.
3. Generation: reconstruct the actual source geometry, endpoint primitive, Gamma/M coefficients and numerical controls.
   The final reference table is not used as a substitute for this coefficient
   generation. generated-before-holdout.json is written before the historical comparison.
4. Validation: compare the generated function to recorded results and run the
   stated controls. A successful comparison does not independently reprove the
   Fourier rule, accepted prepared input or adopted HLOG dictionary.

The outputs are the batch checks.json, numeric.json, generated-before-holdout.json and per-entry result.json.
The data and check indices identify stable adopted source/primitives separately
from the new trial report. The accepted result.json in this entry is kept
unchanged while the calculation writes its fresh result to the new directory.

## Language and reuse boundary

The entry script is an executable Python adapter, not an integration algorithm.
The actual source-derived one-loop calculation is in the shared Python runtime.
Wolfram/native data supplies the original C/J where recorded; the one-loop
Fourier rule and standard Euler/Beta identities are declared accepted inputs.
No Wolfram process or Maple process is invoked by this command.
No unresolved parameter subintegral. The Fourier propagator is accepted upstream physics; the current derivation is not a formal QFT proof.



For a long cold integration the native checkpoints must be preserved, but this
finite algebra replay needs no fabricated checkpoint/resume protocol. A failed
run leaves its output for inspection; use a new directory on the next attempt.
