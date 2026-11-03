# 4e-04: from its original parameter kernel to the complete I

[Original graph, momentum and C/J definition](../../../data/2026-11-02/integrals/2loop/4e-04/README.md) · [Exact coefficients](result.json) · [Common conventions](../METHODS.md)

The difficult exchange representative has attachments (2,2,4,4), m=2
and independent ordered chains 0<t1<t2<1, 0<s1<s2<1:

$$J=(1+z)^2\int_\Omega
D(t_1,s_1)^{\epsilon-1}D(t_2,s_2)^{\epsilon-1}d^4s,\qquad
Q=z^{-\epsilon}J,\qquad I=-\frac{Z^2Q}{8c_2}.$$

The two propagators couple the low endpoints to each other and the high
endpoints to each other. Squaring one independent kernel would erase this
ordering. For positive epsilon an integrable product-kernel bound permits
the particular reorderings used below.

First integrate t1 from0 to u=t2, and s2 from v=s1 to1. The exact result is

$$Q=\frac{(1+z)^2z^{-\epsilon}}{\epsilon^2}
\int_0^1du\int_0^1dv\,
\frac{[D^\epsilon-(1-v)^\epsilon][(zu)^\epsilon-D^\epsilon]}
{[(1+z)v-1][(1+z)u-1]}.$$

Both difference quotients have removable internal points. Their terms cannot
be separated into unrelated principal-value integrals. With p=(1+z)u-1 and
q=(1+z)v-1, set A=(1+p)^epsilon, B=(1-p/z)^epsilon and
C=(1+pq/z)^epsilon. An antisymmetric integrable addition has zero integral
on the common square [-1,z]², turning the numerator into
(A+B)(C-1)-(C²-1)-(A-1)(B(q)-1).

Use the actual primitive J_a(x)=integral_0^x((1-t)^a-1)dt/t to integrate q.
One integration by parts then leaves endpoint products and the unique
antisymmetric ordered integral

$$K(a,b;z)=\int_{0<x<y<1}\frac{[(1+zx)^a-1][(1-y)^b-1]
-[(1-x)^b-1][(1+zy)^a-1]}{xy}\,dxdy.$$

K was explicitly evaluated in the source derivation by a z derivative,
integration by parts and Beta summation, with the boundary constant fixed.
It equals 2F(a,b;-z)-2F(a,0;-z)-J_b(1)J_a(-z), where

$$F(a,b;x)=-\frac{ax}{1+b}\,{}_4F_3(1,1,1,1-a;2,2,2+b;x).$$

Putting it back into the same Q, then converting the two reflected Euler
sources to the fixed moments, gives the completely evaluated mother formula

$$Q=\frac12V(\epsilon)^2-2V_N(\epsilon)
+\frac{2\chi(\epsilon)}\epsilon V(2\epsilon).$$

There is no unknown K in this answer. The endpoint constant chi/epsilon is
regular and essential. For a separate coefficient-integral route set
DeltaL=log(D)-log(1-v), DeltaU=log(zu)-log(D),
Phi=log(D)+[log(1-v)+log(u)-ell]/2 and
K0=(1+z)² DeltaL DeltaU/([(1+z)v-1][(1+z)u-1]). Then

$$Q_0=\int K_0,\quad Q_1=\int K_0\Phi,\quad
Q_2=\int K_0\left[\Phi^2/2+(\Delta_L^2+\Delta_U^2)/24\right].$$

The original HyperInt scripts integrated these genuine two-variable kernels.
This publication reuses their HLOG outputs and freshly checks
a0=-Q0/8, a1=-Q1/8 and a2=-(Q2+(zeta2-ell²/2)Q0)/8. Independent new real
quadrature also compares all three original kernels to the M/N formulas.
The first term is a0=-M01²/48 and all poles vanish. All displayed positive
orders require M/N but no tau_D4. The source proof and full original scripts
are preserved in the packet; Maple was not rerun here.

## Every Laurent coefficient through epsilon²

The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=0.$$

$$a_{0}(z)=- \frac{M_{0,1}^{2}}{48}.$$

$$a_{1}(z)=- \frac{M_{0,1} M_{1,1}}{16} + M_{0,1} \zeta\left(3\right) - \frac{\pi^{2} M_{1,1}}{12} + \frac{N_{1,2}}{8}.$$

$$a_{2}(z)=- \frac{41 M_{0,1}^{3}}{576} - \frac{M_{0,1} M_{2,1}}{64} + \frac{M_{0,5}}{30} - \frac{M_{1,1}^{2}}{64} + M_{1,1} \zeta\left(3\right) - \frac{\pi^{2} M_{2,1}}{24} - \frac{N_{0,4}}{32} + \frac{N_{2,2}}{16}.$$

## Calculation and verification files

- [Family coefficient engine](../family_engines.py) and [portable replay](../reproduce.py).
- [Original-kernel and independent real-moment checks](../check_integrals.py).
- [Byte-preserved current individual result](../sources/benchmarks/2ne/results/4e-04/r002.json).
- [Source manifest and precise historical reuse](../source-manifest.json).

The current publication runs the stated algebra/integral controls, preserves the original integration proof at its recorded scope, and does not claim a new formal proof or complete complex-domain continuation. See [verification](../verification.json) and [independent integral controls](../integral-checks.json).
