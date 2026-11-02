# Two-loop evaluated integrals

Twelve source-connected orbit contributions through epsilon². Start with [conventions and evaluated mother functions](METHODS.md). The [machine catalogue](catalog.json) fixes native physical signs and color weights.

| Integral | Family | Route |
|---|---|---|
| [4e-01](4e-01/derivation.md) | pure | ordered elementary powers |
| [4e-02](4e-02/derivation.md) | pure | one weighted opposite kernel |
| [4e-03](4e-03/derivation.md) | pure | ordered Beta difference |
| [4e-04](4e-04/derivation.md) | pure | actual Q/K evaluation and M/N sources |
| [4e-05](4e-05/derivation.md) | pure | own ordered weight and crossing bridge |
| [4e-06](4e-06/derivation.md) | pure | two independent opposite kernels |
| [3e-3g-01](3e-3g-01/derivation.md) | TGV | one-dimensional Beta difference |
| [3e-3g-02](3e-3g-02/derivation.md) | TGV | explicit ordered reflection then Beta |
| [3e-3g-03](3e-3g-03/derivation.md) | TGV | contact-preserving IBP and single kernel |
| [3e-3g-04](3e-3g-04/derivation.md) | TGV | six faces + subtracted body, HLOG and tau |
| [2e-Pi1-01](2e-Pi1-01/derivation.md) | Pi1 | full transverse operator, elementary endpoints |
| [2e-Pi1-02](2e-Pi1-02/derivation.md) | Pi1 | full transverse operator, four corners and V |

## Reproduce

```text
python integrals/two-loop/reproduce.py --output-dir replay-output/two-loop-new
python integrals/two-loop/check_integrals.py --output-dir replay-output/two-loop-controls-new
```

The coefficient replay performs149 exact comparisons and two negative controls; separate original-kernel/moment controls perform68 checks. These counts overlap and are not a count of independent integral proofs. Both scripts require a new output directory. Historical HyperInt integrations and the M/N HLOG dictionary are reused explicitly; this publication does not silently rerun original Maple/Mathematica writers.

The nativePi1 r002 convention is used: assembly_sign=+1, weights+B1(epsilon). Original legacy sign history remains preserved in sources; never reapply its-1. Four-color assembly, coupling renormalization and finite RG are subsequent observable-level stages.
