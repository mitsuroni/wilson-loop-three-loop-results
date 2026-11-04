# Reproduce 2e-Pi1+Pi1-02

From the repository root, use the installed Python3.12 environment with the shared [requirements](../../families/three-loop-insertions/runtime/requirements.txt). The configured interpreter is D:/py/python.exe. Choose a new output directory:

```powershell
D:/py/python.exe integrals/three-loop/2e-Pi1+Pi1-02/code/reproduce.py --output D:/wilson-agent/replay-output/2e-Pi1+Pi1-02-new
```

The own [adapter](code/reproduce.py) runs the same supported fifteen-entry batch; select this member in generated.json. It does not execute fifteen independent cold integrations. A failed or completed output directory is retained; use a fresh name for another run.

The actual [chains_pi2.py](../../families/three-loop-insertions/runtime/scripts/chains_pi2.py) reads the exact native prepared sources, reconstructs the full numerator/distance, checks the original ordered domain, performs the endpoint or mixed-derivative reduction, constructs Gamma jets and repeats the full colour convolution. The current [shared mathematical derivation](../../families/three-loop-insertions/chains-pi2.md) gives each stage's explicit equations. In particular the opposite branch retains all four corners, and Pi2 regular multipliers are generated from their native A/B rational functions with U and inverse U through5, not read from the old operator-expansion table.

Four original M text outputs supply already evaluated HLOG coordinates; the internal bubble rule, Pi2 full-D 1PI source where applicable, and the Fourier rule are accepted inputs. Python re-executes the downstream source identities, Gamma/rational expansion and complete colour assembly. It does not rerun a historical Maple/Wolfram integration. The byte-preserved old writer scripts may contain absolute installation/output paths and are historical evidence, not current batch commands.

Generation writes generated.json, per-entry outputs and run metadata before final-reference comparison. Exact source, scalar/colour comparisons and rational shuffle checks run afterwards. The [data index](data/index.json) identifies adopted functions, actual code and comparison-only whole-I tables by object and version.

The [current execution report](../../../workflow/checks/integral-documentation-insertions-20261104/replay.json) records the actually executed task batch and copied-closure check. Its scientific scope is separately described in [validation](validation.md). No automatic Git action, scientific selector change, or formal proof is performed by the command.
