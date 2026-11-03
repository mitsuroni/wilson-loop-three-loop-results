# 3e-3g-02: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/3e-3g-02/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

This source has legs (3,3,4), not the same attachment list as01.
Let (x,y,t)=(s31,s32,s41), retain 0<x<y<1 and outer order y,t,x.
Its numerator and star polynomial are

$$S=czt[a(1-x)+b(1-y)],\qquad P=cz^2t(b-a)/4.$$

Reflect the oriented ordered domain, reversing the endpoint order together
with the edge reflection, and exchange the corresponding simplex labels
a and c. The simplex measure is invariant and the full numerator is carried
to the01 kernel, including its sign. The prepared C and orbit size4 agree.
This is why the full evaluated integrals coincide, not a deletion of the
second diagram.

The surviving ordered calculation is therefore the same Beta difference
as01, giving

$$I_{02}=\frac{G_2}{32\epsilon^4}
[2\Gamma(1+\epsilon)\Gamma(1-\epsilon)-B(1-\epsilon,1-\epsilon)].$$

The tabulated coefficients are independently bound to the02 result record.
The portable algebra checks that record at every known negative and finite
order; it does not merely duplicate the01 JSON without checking the second
source. The original reflection proof is reused in this release; no fresh
five-dimensional02 numerical integration is claimed.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=\frac{1}{32}.$$

$$a_{-3}(z)=- \frac{1}{16}.$$

$$a_{-2}(z)=- \frac{1}{8} + \frac{5 \pi^{2}}{192}.$$

$$a_{-1}(z)=- \frac{1}{4} - \frac{\pi^{2}}{96} + \frac{7 \zeta\left(3\right)}{48}.$$

$$a_{0}(z)=- \frac{1}{2} - \frac{\pi^{2}}{48} - \frac{\zeta\left(3\right)}{24} + \frac{119 \pi^{4}}{11520}.$$

$$a_{1}(z)=-1 - \frac{\pi^{2}}{24} - \frac{7 \pi^{4}}{5760} - \frac{\zeta\left(3\right)}{12} + \frac{31 \zeta\left(5\right)}{80} + \frac{23 \pi^{2} \zeta\left(3\right)}{288}.$$

$$a_{2}(z)=-2 - \frac{\pi^{2}}{12} - \frac{7 \pi^{4}}{2880} - \frac{\zeta\left(3\right)}{6} - \frac{\pi^{2} \zeta\left(3\right)}{144} - \frac{\zeta\left(5\right)}{40} + \frac{31 \zeta^{2}\left(3\right)}{144} + \frac{569 \pi^{6}}{161280}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/tgv/3e-3g/results/3e-3g-02/r001.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
