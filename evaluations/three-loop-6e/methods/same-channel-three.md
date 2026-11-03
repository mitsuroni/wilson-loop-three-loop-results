# Three opposite propagators in one channel: 6e-20, 6e-21 and 6e-22

There is no adjacent endpoint to generate an epsilon pole. There are still two ordered edges, so the three propagators do not factor into three independent one-opposite integrals. The remaining information is their pairing permutation.

Let

\[
\Omega_3=\{0<t_1<t_2<t_3<1,\quad0<s_1<s_2<s_3<1\}
\]

and define the centered ordered mother

\[
R_\pi(\epsilon,z)=-(1+z)^3z^{-3\epsilon/2}
\int_{\Omega_3}\prod_{i=1}^3D_z(t_i,s_{\pi(i)})^{\epsilon-1}\,d^3t\,d^3s.
\tag{13}
\]

The actual roots give J_pi=z^(3 epsilon/2) R_pi. Their permutations and kappa values are (123,-1/16) for 20, (132,-1/32) for 21, and (213,-1/32) for 22. The permutation is reconstructed from propagators and edge bounds.

Each endpoint belongs to one D kernel. At either zero corner a kernel is comparable to x+y, and integral_(0,delta)^2 (x+y)^(-1-eta) dxdy is finite for eta<1. Enlarging the ordered region to the cube yields a product of three such bounds. This gives absolute convergence and a locally uniform bound for epsilon near zero, including its finite logarithmic derivatives. Therefore J and I are holomorphic at zero, and c_-6,...,c_-1 vanish. Only the ordinary original integral at epsilon=0 is needed for c0; no discarded finite-part prescription is introduced.

Exchanging the supporting edges sends pi to pi^(-1). Simultaneously reflecting and reversing both edges sends pi to r pi r, r=321, and gives R_pi(1/z)=R_(r pi r)(z). In particular R123 is crossing even, while R132 and R213 exchange. With ell=log z,

\[
I_{20}=-\frac{\mathcal G^3R_{123}}{16\cosh(3\epsilon\rho/2)},\qquad
I_{21}=I_{22}=-\frac{\mathcal G^3(R_{132}+R_{213})}{64\cosh(3\epsilon\rho/2)}.
\tag{14}
\]

This is a complete physical-I equality for 21/22. The original J21 and J22 generally differ.

## Reducing the finite integral before HyperInt

Set a=1/(1+z), d_a(u,v)=a(1-u)(1-v)+(1-a)uv, and

\[
L_-(u,v)=\log d_a(u,v)-\log a-\log(1-v),
\]
\[
L_+(u,v)=\log(1-a)+\log v-\log d_a(u,v).
\]

Integrating the smallest and largest t endpoints gives

\[
R_{123}(0,z)=-\int_0^1du\int_0^1dw\int_0^w dv\int_w^1dr
\frac{L_-(u,v)L_+(u,r)}{(v-a)(r-a)d_a(u,w)}.
\tag{15}
\]

For the 132 ordering another endpoint is elementary:

\[
R_{132}(0,z)=-\int_0^1du\int_0^1dw\int_0^w dv
\frac{L_-(u,v)L_+(u,w)L_+(w,u)}{(v-a)(w-a)(u-a)}.
\tag{16}
\]

All apparent internal poles at a are removable in their complete logarithmic differences. The historical family computation integrated these finite combinations and transformed the answers to z-HLOG. The current original-kernel engine has a separate archived q=0 calculation that emits the root's own six ordered bounds directly. P2 reuses and checks the latter raw output, then reconstructs (14); it does not relabel either old execution as a new CAS run.

## The complete finite answer and the existing E function

Put T3=R123(0,z). The complete HLOG calculation gives

\[
R_{132}^+(0,z)=\frac{M_{0,1}^3}{36}-T_3,
\]

and hence

\[
c_0^{20}=-T_3/16,\qquad
c_0^{21}=c_0^{22}=T_3/32-M_{0,1}^3/1152.
\tag{17}
\]

The compact existing-family value is

\[
T_3=2E-\tfrac12 M_{0,1}M_{2,1}+\tfrac12 M_{1,1}^2
+\frac{\pi^2}{9}M_{0,3}-\frac1{12}M_{0,5}.
\tag{18}
\]

E is the already defined double-kernel insertion direction, not a name invented for an uncomputed part of T3. For a self-contained definition, let

\[
Q_{\mathrm{id};0,1}(\epsilon,z)=(1+z)^2z^{-\epsilon}
\int_{\Omega_2}D_z(t_1,s_1)^{\epsilon-1}D_z(t_2,s_2)^{\epsilon-1}s_2^\epsilon\,d^4s,
\]

\[
E=\log z\,[\epsilon^1]\frac{Q_{\mathrm{id};0,1}(\epsilon,z)-Q_{\mathrm{id};0,1}(\epsilon,1/z)}2.
\]

Equivalently, E is log(z)(1+z)^2/2 times the Omega2 integral of log[s2/(1-s1)]/(D1 D2). Its complete saved HLOG is available with the representation data. In particular E(1)=0. Equations (17)--(18) are identities in full word/period coordinates; a fit at sample points would not establish them.

The three source entries retain their own seven Laurent slots. No claim is made here that the complete epsilon-dependent R functions reduce to the finite basis in (18). Historical endpoint and numerical controls remain separately labelled, and formal=not_run.

The original [method](../evidence/historical/benchmarks/2ne/families/same-channel-three/method.md) and [basis definitions](../evidence/historical/benchmarks/2ne/families/same-channel-three/basis-functions.md) are preserved for provenance.
