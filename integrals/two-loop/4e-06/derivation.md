# 4e-06: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/4e-06/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

Each of the four lightlike edges carries one attachment. The original
order s11,s31,s21,s41 has four independent unit intervals; the two propagators
are Dbar(s11,s31)=z(1-s11)(1-s31)+s11*s31 and D(s21,s41), with numerator (1+z)².

The two pairs contain disjoint variables, so the integral really factorizes.
The change (s11,s31)->(1-s11,1-s31) maps Dbar to D with unit Jacobian.
Centering each opposite-edge integral gives

$$V(\epsilon)=-(1+z)z^{-\epsilon/2}
\int_0^1du\int_0^1dv\,D(u,v)^{\epsilon-1},\qquad
J=z^\epsilon V(\epsilon)^2.$$

Here there is no ordered-chain factor1/2, because the pairs are independent.
The orbit size is only one, giving kappa=-1/16 rather than the coefficient
used for 4e-02:

$$I_{4e-06}=-\frac{Z^2V(\epsilon)^2}{16c_2}.$$

The centered source is regular at epsilon0. Thus every negative Laurent
coefficient vanishes and a0=-M01²/16. Its positive orders follow by actually
multiplying the V series, external Z² and reciprocal c2, not by treating the
finite square as an all-epsilon result. Original four-variable integer-power
checks independently confirm the factorization and multiplicity.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=0.$$

$$a_{0}(z)=- \frac{M_{0,1}^{2}}{16}.$$

$$a_{1}(z)=- \frac{M_{0,1} M_{1,1}}{16}.$$

$$a_{2}(z)=- \frac{11 M_{0,1}^{3}}{192} - \frac{\pi^{2} M_{0,1}^{2}}{24} - \frac{M_{0,1} M_{2,1}}{64} - \frac{M_{1,1}^{2}}{64}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/2ne/results/4e-06/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
