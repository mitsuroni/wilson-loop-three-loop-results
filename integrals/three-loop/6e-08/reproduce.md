# Reproduce 6e-08

From the repository root, use Python3 and the shared
[requirements](../../families/three-loop-6e/runtime/requirements.txt).
The local configured interpreter is D:/py/python.exe.

```powershell
D:/py/python.exe integrals/three-loop/6e-08/code/reproduce.py --output D:/wilson-agent/replay-output/6e-08-6e-new
```

The [adapter](code/reproduce.py) runs the complete39-root same-family batch;
select this orbit's own generated member. The output must be a NEW directory
separate from the immutable inputs, primitives and scientific results. To
rerun, choose a different final name. Do not execute historical scripts against
their old absolute output paths.

The actual [batch evaluator](../../families/three-loop-6e/runtime/scripts/evaluate_batch.py)
loads the copied engine and original ordered prepared roots, verifies their
bindings, reexecutes structural recognition and physical coefficient assembly,
and writes generated.json, primitive-replay.json and run.json. It opens no old
final-I reference during generation and launches no external CAS or network
process. Frozen final tables enter only the subsequent verifier.

For this graph: Actual ordered adjacent endpoints -> explicit evaluated Beta/Gamma word product -> physical crossing -> order6 finite convolution.
The downstream calculation performs a direct exact Gamma/Beta expansion; no archived whole-I value supplies its answer.

The [archive](../../families/three-loop-6e/runtime/primitives/objects/cas-runs/)
contains real historical Maple inputs, logs and evaluated coefficients. The
[shared derivation](../../families/three-loop-6e/adjacent-chains.md) describes their actual
mathematical role. These are not the present command, and Python does not
pretend to rerun Maple. No Wolfram process is required merely to perform this
curation replay. Diagram generation, Fourier transformation and color projection
are accepted upstream stages.

The [data index](data/index.json) identifies the concrete primitive, assembly
and comparison roles. The [current report](../../../workflow/checks/integral-documentation-6e-20261103/replay.json) gives the actual
execution and output paths. No checkpoint/resume mechanism is fabricated for
this finite replay: retain failed output for diagnosis and rerun into a new
output directory. The original CAS checkpoints remain part of their native
historical records.
