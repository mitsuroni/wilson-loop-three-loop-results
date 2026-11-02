# 3e-3g-01: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/3e-3g-01/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The three vertex legs attach to edges (3,4,4). Put
(t,x,y)=(s31,s41,s42), with outer order t,y,x and 0<x<y<1.
For a+b+c=1 the complete kernel is

$$S=az(1-t)(bx+cy),\qquad P=az^2(1-t)(b-c)/4.$$

The alpha-simplex weight (abc)^(-epsilon) is part of J. Integrate the
isolated simplex scale and two elementary path scales while retaining the
ordered same-edge pair. The genuine surviving one-dimensional integral is

$$\frac{J}{z^{2\epsilon}}=\frac1{16\epsilon^3}
\int_0^1db\,b^{-\epsilon-1}(1-b)^{-\epsilon}(2b-1)
\frac{1-(1-b)^{-1+2\epsilon}}{-1+2\epsilon}.$$

The endpoint difference must remain together. Beta evaluation and
meromorphic continuation yield

$$J=\frac{z^{2\epsilon}}{16(1-2\epsilon)\epsilon^4}
[2\Gamma(1+\epsilon)\Gamma(1-\epsilon)-B(1-\epsilon,1-\epsilon)],$$
$$I_{01}=\frac{G_2}{32\epsilon^4}
[2\Gamma(1+\epsilon)\Gamma(1-\epsilon)-B(1-\epsilon,1-\epsilon)].$$

The cyclic combination removes z, but it must be applied to J proportional
to z^(2epsilon); the self-crossing formula used for03/04 is not a starting
assumption here. The new replay expands this exact Gamma combination to
epsilon², and the new real checks integrate the surviving original Beta
kernel at two noninteger epsilon values. Its apparent epsilon=1/2 poles in
intermediate formulas are handled by the complete combination's limit.

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
- [Byte-preserved current individual result](../sources/benchmarks/tgv/3e-3g/results/3e-3g-01/r001.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
