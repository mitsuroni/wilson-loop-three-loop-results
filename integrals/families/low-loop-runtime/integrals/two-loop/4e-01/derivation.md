# 4e-01: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/4e-01/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The attachment order is (3,3,4,4). Both exchanges join adjacent edges;
their small parameters obey s31<s32 and s41<s42. Define A(u,v)=(1-u)v.
The entire ordered J is

$$J=z^2\int_0^1ds_{32}\int_0^1ds_{42}\int_0^{s_{32}}ds_{31}
\int_0^{s_{42}}ds_{41}\,[zA(s_{31},s_{41})]^{\epsilon-1}
[zA(s_{32},s_{42})]^{\epsilon-1}.$$

Group the two B-type endpoint factors (1-s)^(epsilon-1) and the two A-type
endpoint factors s^(epsilon-1). Each two-point ordered chain equals half the
square of its one-point integral, because its two equal factors are symmetric
on the square. In the positive-epsilon convergence domain this is an ordinary
Fubini calculation, not deletion of the ordering. Each chain gives
1/(2epsilon²); therefore J=z^(2epsilon)/(4epsilon⁴).

Applying the four rotations and extracting the symmetric external scale gives

$$\mathcal B_2[J]=\frac1{4\epsilon^4},\qquad
I_{4e-01}=-\frac{Z^2}{16\epsilon^4}.$$

The cancellation of z is a result of the cyclic normalization. Expanding the
Gamma logarithm through order six gives every displayed coefficient through
epsilon². In particular there is no epsilon^-3 term: the single Euler
compensation has cancelled the linear Euler constant. The fresh original-kernel
controls integrate the four-variable polynomial at positive integer epsilon
and compare the same Gamma-stripped physical expression.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=- \frac{1}{16}.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=- \frac{\pi^{2}}{96}.$$

$$a_{-1}(z)=- \frac{\zeta\left(3\right)}{24}.$$

$$a_{0}(z)=- \frac{7 \pi^{4}}{5760}.$$

$$a_{1}(z)=- \frac{\pi^{2} \zeta\left(3\right)}{144} - \frac{\zeta\left(5\right)}{40}.$$

$$a_{2}(z)=- \frac{31 \pi^{6}}{241920} - \frac{\zeta^{2}\left(3\right)}{72}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/2ne/results/4e-01/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
