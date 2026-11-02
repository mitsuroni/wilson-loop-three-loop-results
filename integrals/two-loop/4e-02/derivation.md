# 4e-02: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/4e-02/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The attachments are (2,3,4,4), with the original outer-to-inner order
s21,s31,s42,s41 and s41<s42. Write D(t,s)=(1-t)(1-s)+zts. The numerator is
-z(1+z), and the propagators are D(s21,s41), z(1-s31)s42.

The endpoint s31 is independent and integrates to 1/epsilon. At fixed s41,
integrating the ordered upper endpoint s42 from s41 to1 gives
(1-s41^epsilon)/epsilon. Thus the same J becomes

$$J=-\frac{(1+z)z^\epsilon}{\epsilon^2}
\int_0^1dt\int_0^1ds\,D(t,s)^{\epsilon-1}(1-s^\epsilon).$$

This is a weighted opposite-edge kernel, not a replacement of the original
domain by four independent endpoint integrals. Pair the reflected weight and
retain the z^(2epsilon) crossing factor. The weighted Euler source then gives

$$I_{4e-02}=\frac{Z^2}{4\epsilon^2c_2}
[V(2\epsilon)-c_1V(\epsilon)].$$

The source's constant terms cancel before division by epsilon². Its first
nonzero term is proportional to epsilon, leaving only a simple pole,
a_-1=M11/8. The coefficients follow from the explicit V source, including
the c2 reciprocal and Gamma mixed orders. This calculation is shared with
4e-05 only after their two original weights and cyclic images are related;
the two orbit IDs remain separate in the color sum.

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
- [Byte-preserved current individual result](../sources/benchmarks/2ne/results/4e-02/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
