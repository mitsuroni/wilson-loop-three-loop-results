# 4e-05: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/4e-05/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

This distinct representative has attachments (2,2,3,4), m=4, original
order s31,s22,s41,s21 and s21<s22. The propagators are
(1-s21)s31 and D(s22,s41), with numerator -(1+z).

Integrate s31 independently, obtaining 1/epsilon. At fixed t=s22 the ordered
s21 primitive is [1-(1-t)^epsilon]/epsilon, so

$$J=-\frac{1+z}{\epsilon^2}\int_0^1dt\int_0^1ds\,
D(t,s)^{\epsilon-1}[1-(1-t)^\epsilon].$$

This is a different bare weight from 4e-02: no equality of the unrotated J
has been assumed. Reflect the endpoints, exchange the kernel arguments,
and track the z^epsilon factors in the cyclic partner. With the same
kappa=-1/4 the full normalized result becomes

$$I_{4e-05}=\frac{Z^2}{4\epsilon^2c_2}
[V(2\epsilon)-c_1V(\epsilon)]=I_{4e-02}.$$

The equality is an evaluated full-epsilon relation with identical physical
normalization. It preserves two physical orbit contributions in assembly.
The new polynomial controls integrate this representative's own ordered
root, including its different dependent bound, before comparing the common
mother expression. The same coefficient engine is reused only at this final
stage.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=\frac{M_{1,1}}{8}.$$

$$a_{0}(z)=\frac{M_{0,1}^{2}}{32} + \frac{\pi^{2} M_{0,1}}{32} + \frac{3 M_{2,1}}{32}.$$

$$a_{1}(z)=\frac{5 M_{0,1} M_{1,1}}{32} + \frac{19 \pi^{2} M_{1,1}}{192} + \frac{7 M_{1,3}}{192} + \frac{7 M_{3,1}}{192}.$$

$$a_{2}(z)=\frac{5 M_{0,1}^{3}}{192} + \frac{37 \pi^{2} M_{0,1}^{2}}{768} + \frac{13 M_{0,1} M_{2,1}}{128} + \frac{31 \pi^{4} M_{0,1}}{1536} + \frac{M_{0,5}}{512} + \frac{M_{1,1} \zeta\left(3\right)}{12} + \frac{17 \pi^{2} M_{2,1}}{256} + \frac{5 M_{2,3}}{256} + \frac{5 M_{4,1}}{512}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/2ne/results/4e-05/r001.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
