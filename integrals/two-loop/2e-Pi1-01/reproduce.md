# Reproduce 2e-Pi1-01

Run from the repository root with Python, SymPy, mpmath and NumPy:

```text
python integrals/two-loop/2e-Pi1-01/code/reproduce.py --output-dir replay-output/2e-Pi1-01-fresh
```

Use a new empty output location. This own-entry adapter runs the same-loop batch runtime, reconstructs the generic Gamma/moment-family coefficient jets, adopts the explicitly identified already assembled TGV04 scalar coefficients, and performs the post-computation reference comparisons. Other two-loop entries therefore occur in its batch outputs. Its result for this entry must include all seven powers -4 through 2.

The stages are parameter root in [input.json](input.json), the explicit reduction in [derivation.md](derivation.md), adopted integrated source data indexed under [data](data/), finite-series assembly by the shared runtime, then exact comparison and the separately described numerical controls indexed under [checks](checks/). Reference answer tables belong to the comparison stage; they are not inputs for the generic Gamma/moment generation. TGV04 specifically assigns the accepted closed formulas for I at orders -2,-1,0 and reads `compact.json` for the already assembled normalized I coefficients at orders 1 and 2. It then compares their complete historical HLOG coordinates. It does not live-convolve the D/H/O subintegral jets. The D/H/O evaluation and its convolution are the explicitly exposed historical mathematics, not operations newly executed by this Python replay.

Historical Wolfram/FeynCalc work established the tensor/parameter root; historical Maple/HyperInt scripts evaluated the stated rational-logarithmic kernels. Python performs exact coefficient convolution, function-coordinate comparison and numerical checks. No requirement to execute all three languages is implied. The original CAS scripts retain their old installation paths and output writers and are not safe portable entry points without an explicit adapter.

The current new batch execution passed 149 assembly/source checks and two negative controls, and the separate original-kernel/function control run passed 68 checks. [Validation](validation.md) links the exact receipts and distinguishes those batch totals from per-entry coverage. The own-entry adapter can be rerun in the fresh directory specified above. A repeat comparison shares the adopted moment dictionary and is not independent verification of that dictionary. No Maple, Wolfram or Lean run is implied by this command; `formal=not_run`.
