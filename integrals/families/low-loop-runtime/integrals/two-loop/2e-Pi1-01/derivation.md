# 2e-Pi1-01: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/2e-Pi1-01/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

This adjacent-edge insertion connects edges3 and4 with multiplicity4.
Integrating the internal massless gluon/ghost/quark bubble gives the transverse
two-point tensor and the full B1(epsilon) color weight. The outer Fourier
transform produces both metric and longitudinal derivative terms. One cannot
discard the latter merely because the operator is transverse.

With s=s31,t=s41, f=z(1-s)t, the complete prepared kernel reduces to

$$J_A=\int_0^1ds\int_0^1dt\,
\frac{z}{1+\epsilon}f^{2\epsilon-1}
=\frac{z^{2\epsilon}}{4\epsilon^2(1+\epsilon)}.$$

The derivative reduction is explicit: for d=z, the numerator
2(1+epsilon)*d*f+f*f_st+(2epsilon-1)*f_s*f_t equals2*z*f.
The source polynomial identity and integer-power integrations are rerun.
Let R1=-Z²/[16epsilon(1-2epsilon)]. The prepared C is
lambda_b²*CF*B1*S14*R1, where S14=[pi*(-2v14+i0)]^(2epsilon), before the later
MSbar coupling/scale conversion. Applying m/2 times the crossing average gives

$$I_A=-\frac{Z^2}{32\epsilon^3(1-2\epsilon)(1+\epsilon)}.$$

This is the current physical negative scalar I. The original r001 source used
the opposite I convention; that legacy-to-native multiplier is -1 once.
Since this packet already contains native r002, the color contribution is
+B1*I_A, not -B1*I_A. Each coefficient is checked against the actual individual
r002 record. The highest pole is triple, so Z² and the rational prefactor
must be kept through epsilon5 to reach epsilon². The gluonic/quark rational
weights must also be convolved with these poles, rather than frozen at D=4.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=- \frac{1}{32}.$$

$$a_{-2}(z)=- \frac{1}{32}.$$

$$a_{-1}(z)=- \frac{3}{32} - \frac{\pi^{2}}{192}.$$

$$a_{0}(z)=- \frac{5}{32} - \frac{\pi^{2}}{192} - \frac{\zeta\left(3\right)}{48}.$$

$$a_{1}(z)=- \frac{11}{32} - \frac{\pi^{2}}{64} - \frac{7 \pi^{4}}{11520} - \frac{\zeta\left(3\right)}{48}.$$

$$a_{2}(z)=- \frac{21}{32} - \frac{5 \pi^{2}}{192} - \frac{\zeta\left(3\right)}{16} - \frac{7 \pi^{4}}{11520} - \frac{\pi^{2} \zeta\left(3\right)}{288} - \frac{\zeta\left(5\right)}{80}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/self-energy/pi1/results/2e-Pi1-01/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
