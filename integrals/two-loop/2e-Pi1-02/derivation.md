# 2e-Pi1-02: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/2e-Pi1-02/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The opposite-edge insertion connects edges2 and4 with multiplicity2.
Set s=s21,t=s41 and f=(1-s)(1-t)+zst. Retaining both terms of the transverse
operator gives

$$J_O=\int_0^1ds\int_0^1dt\,
\frac{-2(1+z)f+(1-2\epsilon)z}{2(1+\epsilon)}f^{2\epsilon-2}.$$

The same mixed-derivative numerator identity as the adjacent case is checked
with d=-(1+z). Its derivative part integrates to four actual corners in the
positive-epsilon domain. Define the already evaluated centered opposite-edge
source H(nu,z)=-(1+z)z^(-nu/2)integral f^(nu-1)dsdt=V(nu,z). Then

$$J_O=z^\epsilon V(2\epsilon,z)
+\frac{1+z^{2\epsilon}}{4\epsilon(1+\epsilon)},$$
$$I_O=-\frac{Z^2}{16\epsilon(1-2\epsilon)}
\left[\frac{V(2\epsilon)}{c_2}+\frac1{2\epsilon(1+\epsilon)}\right].$$

The second term is a genuine corner contribution. Dropping it changes the
double pole and every lower result. The source V(2epsilon), not V(epsilon),
has expansion M01+epsilon*M11+epsilon²(M21/2+M03/6)+epsilon³(M31+M13)/6.
This fixes exactly which moments are needed through epsilon² after the
external simple pole; the separate corner rational term has a double pole.

As for the adjacent insertion the released scalar is native physical r002,
already with its negative sign. Assembly therefore uses +B1. Its two color
channels share this scalar integral, but their epsilon-dependent rational
weights are convolved separately. The new checks compare the own individual
source at every stored order and independently integrate the polynomial
kernel/corners, not just the sum of the two Pi1 diagrams.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=- \frac{1}{32}.$$

$$a_{-1}(z)=- \frac{M_{0,1}}{16} - \frac{1}{32}.$$

$$a_{0}(z)=- \frac{M_{0,1}}{8} - \frac{M_{1,1}}{16} - \frac{3}{32} - \frac{\pi^{2}}{192}.$$

$$a_{1}(z)=- \frac{5 M_{0,1}^{2}}{96} - \frac{\pi^{2} M_{0,1}}{24} - \frac{M_{0,1}}{4} - \frac{M_{1,1}}{8} - \frac{M_{2,1}}{32} - \frac{5}{32} - \frac{\pi^{2}}{192} - \frac{\zeta\left(3\right)}{48}.$$

$$a_{2}(z)=- \frac{5 M_{0,1}^{2}}{48} - \frac{M_{0,1} M_{1,1}}{16} - \frac{\pi^{2} M_{0,1}}{12} - \frac{M_{0,1}}{2} - \frac{M_{0,1} \zeta\left(3\right)}{24} - \frac{\pi^{2} M_{1,1}}{24} - \frac{M_{1,1}}{4} - \frac{M_{1,3}}{96} - \frac{M_{2,1}}{16} - \frac{M_{3,1}}{96} - \frac{11}{32} - \frac{\pi^{2}}{64} - \frac{7 \pi^{4}}{11520} - \frac{\zeta\left(3\right)}{48}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/self-energy/pi1/results/2e-Pi1-02/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
