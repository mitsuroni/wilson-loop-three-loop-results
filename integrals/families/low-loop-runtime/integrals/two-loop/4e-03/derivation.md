# 4e-03: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/4e-03/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

This representative has attachments (2,3,3,4) and original order
s21,s32,s41,s31 with s31<s32. Its integrand is

$$z^\epsilon s_{32}^{\epsilon-1}s_{41}^{\epsilon-1}
(1-s_{21})^{\epsilon-1}(1-s_{31})^{\epsilon-1}.$$

The endpoints s21 and s41 each give 1/epsilon. The remaining mixed ordered
chain is evaluated, not renamed:

$$\int_0^1dv\,v^{\epsilon-1}\int_0^vdu(1-u)^{\epsilon-1}
=\frac{1-R_2}{\epsilon^2}.$$

The incomplete inner primitive supplies [1-(1-v)^epsilon]/epsilon;
integrating its two terms gives the stated Beta ratio. Hence

$$J=-\frac{z^\epsilon(R_2-1)}{\epsilon^4},\qquad
I_{4e-03}=\frac{Z^2(R_2-1)}{4\epsilon^4c_2}.$$

The Gamma ratio begins at epsilon². This explains why the apparent fourth
pole is actually a double pole. Crossing removes the centered z^epsilon
factor but leaves 1/c2, which contributes to the finite and positive orders.
The fresh control compares both the original ordered polynomial J and this
closed expression at integer epsilon; the gamma factor itself is removed
there to avoid evaluating Gamma(1-epsilon) at its poles.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=- \frac{\pi^{2}}{24}.$$

$$a_{-1}(z)=\frac{\zeta\left(3\right)}{2}.$$

$$a_{0}(z)=- \frac{\pi^{2} M_{0,1}}{24} - \frac{49 \pi^{4}}{1440}.$$

$$a_{1}(z)=\frac{M_{0,1} \zeta\left(3\right)}{2} + \frac{3 \zeta\left(5\right)}{2} + \frac{2 \pi^{2} \zeta\left(3\right)}{9}.$$

$$a_{2}(z)=- \frac{5 \pi^{2} M_{0,1}^{2}}{144} - \frac{23 \pi^{4} M_{0,1}}{480} - \frac{223 \pi^{6}}{12096} + \frac{5 \zeta^{2}\left(3\right)}{6}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/2ne/results/4e-03/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
