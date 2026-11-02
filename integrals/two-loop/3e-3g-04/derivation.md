# 3e-3g-04: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/3e-3g-04/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The three legs are on different edges (2,3,4), so x=s21,y=s31,w=s41
are independent. With a+b+c=1,

$$S=ab\,y(1-x)+acD(x,w)+zbc\,w(1-y),$$
$$4P=ab[-(1+z)(1-x)+zy]
+ac[1+3z-(1+z)x+z(1+z)w]+zbc[1-y-(1+z)w].$$

This is the difficult two-loop vertex calculation, not a Gamma-only case.
Write the original numerator as an exact path derivative plus its body:

$$P=\frac{-z\partial_xS-(1+z)\partial_yS+\partial_wS}{4}+P_{rec},$$
$$P_{rec}=\frac{1+z}{2}c[a(1-x+zw)-zbw].$$

All six endpoint faces and the recoil integral survive. If F_i,b is the
face integral of S^(-1+2epsilon) with the original simplex measure, define
D=epsilon(F_2,0-F_2,1), H=F_1,0 and O=the P_rec body. The exact organized J is

$$J=\frac1{4(1-2\epsilon)}\left[(1+z^{2\epsilon})C_{cusp}
-\frac{1+z}{\epsilon}D-zH(\epsilon,z)
-z^{2\epsilon-1}H(\epsilon,1/z)\right]+O,$$
$$C_{cusp}=\frac{\Gamma(1+\epsilon)\Gamma(1-\epsilon)
-\Gamma(1-\epsilon)^2/\Gamma(1-2\epsilon)}{4\epsilon^4}.$$

The original computation first normalized the endpoint-singular body, then
expanded integrable logarithmic kernels. To obtain the physical epsilon¹
coefficient one must integrate D2,H1,O1; epsilon² needs D3,H2,O2. These six
subintegrals were genuinely evaluated in separate HyperInt sessions with
divergence checks enabled. Their complete scripts and HLOG outputs are
preserved in this packet; the current publication does not rerun those CAS
sessions. Naming D,H,O is an intermediate decomposition, not a solved answer.

For example, after subtraction the recoil density is (A/epsilon+B)exp(epsilon Lb),
so O1=integral(B*Lb+A*Lb²/2), O2=integral(B*Lb²/2+A*Lb³/6). The H1/H2 kernels
are the original boundary density times Lh and Lh²/2. The archived easy/hard/
recoil scripts specify every rational factor and domain. The original raw
kernel has the proved crossing relation J(epsilon,1/z)=z^(-2epsilon)J(epsilon,z).
Consequently the two terms in the common crossing numerator both equal J:
the factor2 cancels the common prefactor1/2. Only at this point does the
physical formula reduce to

$$I_{04}=\frac{e^{-2\gamma_E\epsilon}\Gamma(2-2\epsilon)}
{1+z^{2\epsilon}}J_{04}(\epsilon,z).$$

Its epsilon² coefficient needs this prefactor through order4 because J
starts at epsilon^-2. No inversion term or factor1/2 was dropped.

The complete coefficients start with zeta2/(8epsilon²),
(-M11+zeta3)/(8epsilon), and M01²/48-M21/8+9zeta2²/10. At epsilon¹ the complete
answer is M/N/zeta. At epsilon² the existing tau_D4 direction enters with
coefficient1/8; its original two-variable definition is in METHODS, not left
undefined. The corresponding complete HLOGs have23 and48 nonzero coordinates.

This replay freshly verifies the complete integrated HLOG and M/N/tau
equality, all lower coefficients against the individual source, and new real
M/N/tau quadrature against frozen HLOG values. It does not represent those
tests as a new raw five-dimensional integration. The old original-chart QMC
and independent pySecDec checks are historical evidence with their original
scope, not rerun claims. No new all-epsilon finite special-function closure
or minimal-function-space theorem is asserted.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=\frac{\pi^{2}}{48}.$$

$$a_{-1}(z)=- \frac{M_{1,1}}{8} + \frac{\zeta\left(3\right)}{8}.$$

$$a_{0}(z)=\frac{M_{0,1}^{2}}{48} - \frac{M_{2,1}}{8} + \frac{\pi^{4}}{40}.$$

$$a_{1}(z)=- \frac{M_{0,1} M_{1,1}}{16} + \frac{3 M_{0,1} \zeta\left(3\right)}{4} - \frac{\pi^{2} M_{1,1}}{32} - \frac{M_{1,3}}{96} - \frac{5 M_{3,1}}{96} - \frac{N_{1,2}}{16} - \frac{37 \pi^{2} \zeta\left(3\right)}{144} + \frac{5 \zeta\left(5\right)}{8}.$$

$$a_{2}(z)=\frac{223 M_{0,1}^{3}}{11520} - \frac{13 \pi^{2} M_{0,1}^{2}}{11520} - \frac{5 M_{0,1} M_{2,1}}{64} + \frac{169 \pi^{4} M_{0,1}}{11520} + \frac{M_{1,1}^{2}}{32} + \frac{5 M_{1,1} \zeta\left(3\right)}{3} - \frac{7 \pi^{2} M_{2,1}}{128} + \frac{25 M_{2,3}}{1536} - \frac{39 M_{4,1}}{1024} - \frac{N_{2,2}}{32} + \frac{\tau^+_{D,4}}{8} - \frac{19 \zeta^{2}\left(3\right)}{24} + \frac{7381 \pi^{6}}{362880}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/tgv/3e-3g/results/3e-3g-04/r001.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
