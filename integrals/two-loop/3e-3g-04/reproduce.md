# Reproduce 3e-3g-04

Run from the repository root with Python, SymPy, mpmath and NumPy:

```text
python integrals/two-loop/3e-3g-04/code/reproduce.py --output-dir replay-output/3e-3g-04-fresh
```

Use a new empty output location. This own-entry adapter runs the same-loop batch runtime, reconstructs the generic Gamma/moment-family coefficient jets, adopts the explicitly identified already assembled TGV04 scalar coefficients, and performs the post-computation reference comparisons. Other two-loop entries therefore occur in its batch outputs. Its result for this entry must include all seven powers -4 through 2.

The historical mathematical stages are the original parameter root, six-face IBP, cusp/Beta evaluation, D/H/O coefficient integration, and the finite convolution exposed explicitly in [derivation.md](derivation.md). Its exact integrated subresults are printed in the mathematical appendices and indexed under [data](data/). The current executable has a narrower role: it assigns the accepted closed formulas for normalized I at orders -2,-1,0 and reads `compact.json` for the already assembled normalized I coefficients at orders 1 and 2. It compares the complete historical HLOG coordinates after substituting the adopted function dictionary. It does not live-convolve the D/H/O jets and does not newly solve the five-dimensional graph. The [result table](result.json) is bound for the later reference-comparison stage; the already assembled adopted scalar source is separately identified rather than presented as new integration input.

Historical Wolfram/FeynCalc work established the tensor/parameter root; historical Maple/HyperInt scripts evaluated the stated rational-logarithmic kernels. Python performs exact coefficient convolution, function-coordinate comparison and numerical checks. No requirement to execute all three languages is implied. The original CAS scripts retain their old installation paths and output writers and are not safe portable entry points without an explicit adapter.

The current new batch execution passed 149 assembly/source checks and two negative controls, and the separate original-kernel/function control run passed 68 checks. [Validation](validation.md) links the exact receipts and distinguishes those batch totals from per-entry coverage. The own-entry adapter can be rerun in the fresh directory specified above. A repeat comparison shares the adopted moment dictionary and is not independent verification of that dictionary. No Maple, Wolfram or Lean run is implied by this command; `formal=not_run`.
