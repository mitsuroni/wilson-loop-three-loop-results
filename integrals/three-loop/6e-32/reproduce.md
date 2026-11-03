# Reproduce 6e-32

From the repository root, with the installed Python 3.12 environment and the declared SymPy/mpmath/jsonschema libraries, choose a fresh output directory:

```text
python -B integrals/three-loop/6e-32/code/reproduce.py --output replay-output/6e-32-fresh
```

This adapter runs the full 39-entry same-batch generator. It does not execute only this member. The output `generated.json` contains this source's coefficients; the accompanying run and primitive-replay receipts state the actual consumed orders and subintegral files. The shared [locked requirements](../../families/three-loop-6e/runtime/requirements-lock.txt) record the tested versions; no online installation is needed for the existing local environment.

The mathematical stages are the explicit original ordered root, the source-specific endpoint/reduction steps in [derivation.md](derivation.md), the adopted already integrated HLOG primitives printed in [the appendix](appendices/evaluated-jets.md), and the new physical Gamma/crossing convolution. In this route the historical primitive quantity is the uncentered regular gap kernel through orders 0,1,2, before its external double pole and scalar phases. It is never an already normalized whole-I reference table. The runtime validates the exact primitive request, regenerated script, original execution/raw output and parser before consuming it.

For a separate exact comparison after generation, run the byte-preserved batch verifier:

```text
python -B integrals/families/three-loop-6e/runtime/scripts/verify_batch.py --results replay-output/6e-32-fresh/generated.json --out replay-output/6e-32-fresh/verification.json
```

The output directory must be new for generation. Historical raw scripts retain their installation paths and old writers as immutable evidence; this current entry does not execute those absolute paths. No Maple/FeynCalc process, network fetch or fresh cold CAS calculation is part of this replay. The generator does not read old final-I reference tables; verification occurs after generation.

The actual current batch evidence is linked in [validation](validation.md): all 39 roots,273 coefficient positions and273 function substitutions passed. Root input/version checks, primitive reuse, physical assembly and independent controls have distinct scopes. `formal=not_run` is unchanged.

To repeat the separate check of the new mathematical convolution for all 14 complex entries, after the primary evidence is available:

```text
python -B workflow/checks/integral-documentation-6e-20261103/author-checks/complex_assembly.py --out replay-output/complex-prose-convolution.json
```

It compares23 coefficient formulas using shared raw primitives and the primary generated candidate, with an independent finite-convolution implementation. It performs no new CAS integration.
