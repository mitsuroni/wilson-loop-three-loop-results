# Reproduce 4e-03

Use the repository root as working directory. Python3 with the adopted
[requirements](../../families/low-loop-runtime/integrals/one-loop/requirements.txt)
is sufficient; the local configured interpreter is D:/py/python.exe.
The entry [code/reproduce.py](code/reproduce.py) delegates to the byte-preserved
same-loop [runtime](../../families/low-loop-runtime/integrals/two-loop/reproduce.py).
It executes the whole twelve-entry two-loop batch; select the member 4e-03 in its output.
There are not twelve separate cold CAS integrations hidden behind these entry points.

```powershell
D:/py/python.exe integrals/two-loop/4e-03/code/reproduce.py --output-dir D:/wilson-agent/replay-output/4e-03-trial-new
```

Choose a genuinely new output directory outside the immutable runtime and
scientific input/result trees. To rerun, change the final directory name.
The wrapper does not overwrite source objects, accepted result files or proofs.
Do not run an old writer whose absolute output path points into a frozen packet.

## What each stage actually does

1. Input and binding: read the original ordered parameter source, prepared C/J
   and the source manifest. Preserve the order and color/normalization stage
   described in the derivation; the manifest is a binding record, not a proof.
2. Adopted mathematics: Two independent endpoints and one mixed ordered Beta chain give J=z^epsilon(1-R2)/epsilon⁴. The Gamma-ratio zero removes two apparent poles; R2 through6 is needed.
3. Generation: use family_engines.py to generate the declared exact mother-function coefficients and complete Gamma/cosh convolution, then convert adopted moments to exact whole-HLOG coordinates.
   The final reference table is not used as a substitute for this coefficient
   generation. Historical complete-I tables are opened for post-generation comparison and selection of an equivalent compact representative.
4. Validation: compare the generated function to recorded results and run the
   stated controls. A successful comparison does not independently reprove the
   Fourier rule, accepted prepared input or adopted HLOG dictionary.

The outputs are the batch catalog.json, verification.json and per-entry result.json.
The data and check indices identify stable adopted source/primitives separately
from the new trial report. The accepted result.json in this entry is kept
unchanged while the calculation writes its fresh result to the new directory.

## Language and reuse boundary

The entry script is an executable Python adapter, not an integration algorithm.
The actual finite exact reconstruction is in the shared Python family engine.
Wolfram/native data supplies the original C/J where recorded; the one-loop
Fourier rule and standard Euler/Beta identities are declared accepted inputs.
No Wolfram process or Maple process is invoked by this command.
No unresolved parameter integral; only the elementary identity M01=-(ell²+pi²)/2 is used to compactify logs.



For a long cold integration the native checkpoints must be preserved, but this
finite algebra replay needs no fabricated checkpoint/resume protocol. A failed
run leaves its output for inspection; use a new directory on the next attempt.

The separate original-kernel numerical controls are not called by the entry
wrapper. Run their actual implementation explicitly if needed:

```powershell
D:/py/python.exe integrals/families/low-loop-runtime/integrals/two-loop/check_integrals.py --output-dir D:/wilson-agent/replay-output/two-loop-integral-controls-new
```

These controls include individual graph roots, weighted real moments and4e04
subtracted two-dimensional coefficients. They share accepted source inputs and
do not independently generate a new complete complex-domain solution.
