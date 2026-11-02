# 2e-01: the adjacent one-loop exchange

The representative joins edges3 and4. Its four cyclic images are the four
adjacent edge pairs of the null quadrilateral. The original ordered kernel
factorizes into two endpoint integrals and gives a double epsilon pole.

- [Continuous derivation](derivation.md): momentum/path contribution,
  raw/prepared C/J, actual parameter integration and the physical I map.
- [Input and source index](input.json): complete frozen-source file and
  field hashes, preserved order and normalization.
- [All-epsilon expression and Laurent coefficients through epsilon3](result.json).
- [Fresh source, symbolic and numerical checks](checks.json).
- [Reproduction commands and scope](reproduce.md).

The result is I=-exp(-EulerGamma epsilon)Gamma(1-epsilon)/epsilon².
It already includes m=4, crossing and one Euler compensation. CF, lambda and
the common S1 scale factor are removed. Nothing is reconstructed by fitting
the finite term to an RG equation.
