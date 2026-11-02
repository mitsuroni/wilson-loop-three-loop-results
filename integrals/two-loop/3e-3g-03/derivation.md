# 3e-3g-03: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/3e-3g-03/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The legs lie on (2,4,4). With (t,x,y)=(s21,s41,s42), 0<x<y<1,

$$S=a[bD(t,x)+cD(t,y)],\qquad
P=(1+z)a(b-c)[1-(1+z)t]/4.$$

The numerator equals (1+z)(partial_y-partial_x)S/4. Integrate the a scale,
then integrate by parts on the original triangle. Its diagonal x=y is a
contact face and occurs with the correct two orientations. The resulting
weight is evaluated by matching its derivative and zero average:

$$F_\epsilon(u)=\int_0^u db\,b^{-\epsilon}(1-b)^{-1-\epsilon}
+\int_u^1db\,b^{-1-\epsilon}(1-b)^{-\epsilon}
-2B_\epsilon=\frac{[u(1-u)]^{-\epsilon}-B_\epsilon}{\epsilon}.$$

Here B_epsilon=B(1-epsilon,1-epsilon); the subtraction contains the contact
face. The remaining single-kernel integral gives

$$I_{03}=\frac{G_2}{8\epsilon^2c_2}
[\mathscr C(\epsilon,z)-B_\epsilon V(2\epsilon,z)],$$
$$\mathscr C=\frac{\pi}{\epsilon\sin(\pi\epsilon)}
[\cos(\pi\epsilon)-\cosh(\epsilon\ell)].$$

Both terms begin with M01; the constant cancels and only a simple pole
remains. Expanding the two different moment sources, the Beta ratio and the
full reciprocal c2 produces the coefficients. The original contact face is
why replacing Fepsilon by its first unsubtracted integral would fail. The
source proof and HLOG dictionary are reused, while the coefficient/Gamma
convolution and current individual comparisons are rerun in this packet.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=- \frac{M_{0,1}}{4} - \frac{M_{1,1}}{8}.$$

$$a_{0}(z)=- \frac{M_{0,1}}{2} + \frac{\pi^{2} M_{0,1}}{48} - \frac{M_{1,1}}{4} - \frac{M_{2,1}}{16}.$$

$$a_{1}(z)=- \frac{5 M_{0,1}^{2}}{24} - \frac{M_{0,1} M_{1,1}}{8} - \frac{\pi^{2} M_{0,1}}{6} - M_{0,1} + \frac{M_{0,1} \zeta\left(3\right)}{4} - \frac{\pi^{2} M_{1,1}}{12} - \frac{M_{1,1}}{2} - \frac{M_{1,3}}{48} - \frac{M_{2,1}}{8} - \frac{M_{3,1}}{48}.$$

$$a_{2}(z)=- \frac{5 M_{0,1}^{2}}{12} + \frac{5 \pi^{2} M_{0,1}^{2}}{288} - \frac{M_{0,1} M_{1,1}}{4} - \frac{M_{0,1} M_{2,1}}{16} - \frac{\pi^{2} M_{0,1}}{3} - 2 M_{0,1} - \frac{M_{0,1} \zeta\left(3\right)}{6} + \frac{59 \pi^{4} M_{0,1}}{2880} - \frac{\pi^{2} M_{1,1}}{6} - M_{1,1} - \frac{M_{1,1} \zeta\left(3\right)}{12} - \frac{M_{1,3}}{24} - \frac{\pi^{2} M_{2,1}}{24} - \frac{M_{2,1}}{4} - \frac{M_{2,3}}{96} - \frac{M_{3,1}}{24} - \frac{M_{4,1}}{192}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/tgv/3e-3g/results/3e-3g-03/r001.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
