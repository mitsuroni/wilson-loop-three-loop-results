# One opposite propagator in each channel: 6e-32, 6e-33 and 6e-34

These three graphs have D_z, the reversed-channel kernel Dtilde_z, and one adjacent propagator. Their ordered endpoints permit the adjacent pair to be integrated first and the two channel factors to separate. This separation is checked from the actual bounds; it would fail if an additional constraint coupled endpoints from the two factors.

The four possible elementary adjacent placements relative to an opposite endpoint s are

| Ordered word | Adjacent endpoint integral |
| --- | --- |
| AO | s^epsilon/epsilon |
| OA | (1-s^epsilon)/epsilon |
| BO | [1-(1-s)^epsilon]/epsilon |
| OB | (1-s)^epsilon/epsilon |

Re epsilon>0 supplies integrability. With distinct propagator endpoints the unrestricted cube bound factors into integrable adjacent densities and two integrable opposite kernels, justifying this regrouping by absolute convergence. The explicit two factors of 1/epsilon must remain until their numerator cancellations have been computed.

Let A=H0 and B=H1, using the centered definition (6) in [the one-opposite method](single-opposite.md). Reversing the channel uses

\[
-(1+z)\widetilde D_z^{\epsilon-1}
=z^\epsilon[-(1+1/z)D_{1/z}^{\epsilon-1}].
\]

Reflection of its two endpoints then converts its weight to the appropriate H channel. This step transports both the kernel scale and the endpoint weight; retaining only one would change the answer. The actual results are

\[
J_{32}=J_{34}=\frac{z^{2\epsilon}}{\epsilon^2}B(A-B),\qquad
J_{33}=\frac{z^{2\epsilon}}{\epsilon^2}(A-B)^2.
\tag{19}
\]

For 32 and 34 the independently integrated factors exchange positions. Both have kappa=-1/16; 33 has kappa=-1/8. Therefore 32=34 holds for original J and for complete physical I. This is a specific proved factor-exchange relation; it is not a rule that physical-I aliases generally identify their J roots.

## Finite physical assembly

Write F=B(A-B), K=(A-B)^2 and ell=log z. Substitution in [the normalization](common-normalization.md) gives

\[
I_{32}=I_{34}=-\frac{\mathcal G^3}{32\epsilon^2\cosh(3\epsilon\rho/2)}
\left[e^{\epsilon\rho/2}F(z)+e^{-\epsilon\rho/2}F(1/z)\right],
\]

\[
I_{33}=-\frac{\mathcal G^3}{16\epsilon^2\cosh(3\epsilon\rho/2)}
\left[e^{\epsilon\rho/2}K(z)+e^{-\epsilon\rho/2}K(1/z)\right].
\tag{20}
\]

The H0/H1 jets q=0,1,2 and the scalar factors through q=2 suffice for c0. Since A(0)=B(0), (19) contains cancellations. They remove the double pole in 32/34 and all poles in 33. Let t_q=tau_(0,q); the exact result is

\[
c_{-1}^{32}=c_{-1}^{34}=\frac{t_0t_1}{16},\qquad
c_0^{32}=c_0^{34}=\frac{3t_0t_2}{16}+\frac{t_1^2}{8}
+\frac{7t_0^3}{576}+\frac{7\zeta_2t_0^2}{192},
\]

\[
c_0^{33}=-\frac{t_1^2}{8}+\frac{t_0^3}{144}+\frac{\zeta_2t_0^2}{48}.
\tag{21}
\]

All other negative-power coefficients are zero. The known conversions

\[
t_0=M_{0,1},\qquad t_1=M_{1,1}/2,\qquad
t_2=M_{2,1}/8-M_{0,1}^2/24
\]

give the compact M/zeta result. Complete HLOG values remain available to check every conversion.

The original family calculation explicitly reused previously evaluated H0/H1 jets, rather than running another HyperInt calculation for the same moments. The current engine's route is a direct ordered-gap reduction through q=2 that keeps the original J hash and bounds. P2 replays that exact archived gap output and assembles I again. These are two organizations of the same root mathematics, not two independent new integrations performed by this batch. The per-orbit records distinguish their execution histories.

See the preserved [historical method](../evidence/historical/benchmarks/2ne/families/double-opposite/method.md) for the full product reduction. Positive-real-z continuation and formal=not_run have the [common scope](common-normalization.md).
