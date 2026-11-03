# Mixed opposite channels: retain the distinct37 and 38 kernels

The two roots have one reverse-channel propagator with independent endpoints $a=s_{11},b=s_{31}$ and two ordered forward-channel propagators on
$\Omega_2=\{0<t_1<t_2<1,0<s_1<s_2<1\}$. Their complete roots are
\[
J_\pi=-(1+z)^3\int_{[0,1]^2\times\Omega_2}
\widetilde D_z(a,b)^{e-1}D_z(t_1,s_{\pi(1)})^{e-1}D_z(t_2,s_{\pi(2)})^{e-1}d^6s,
\]
with id pairing for 37 and swap for 38. There is no adjacent density and no explicit regulator pole. Each own article prints its original nested order; neither paired kernel is replaced by the other's input.

The three propagators have six distinct endpoints. Their zero corners have the model $(x+y)^{e-1}$. Enlarging the actual ordered domain to the cube gives a product of three finite two-dimensional corner bounds, locally uniform for $|e|<\delta<1$, with slightly stronger powers for any fixed number of logarithmic derivatives. Therefore Fubini is justified and both original $J$ and normalized $I$ are holomorphic at zero. All $c_{-6},\ldots,c_{-1}$ vanish, and $c_0$ is the ordinary original $e=0$ integral, without any isolated finite-part replacement.

Define the already evaluated centered kernels
\[
H_0=-(1+z)z^{-e/2}\int_{[0,1]^2}D_z(u,v)^{e-1}du\,dv,
\quad
Q_\pi=(1+z)^2z^{-e}\int_{\Omega_2}D_1^{e-1}D_2^{e-1}d^4s.
\]
The independent square permits actual factorization. Reflection $(a,b)\mapsto(1-a,1-b)$ has absolute Jacobian one and sends $\widetilde D_z$ to $D_z$; the original numerator sign combines with the negative sign in $H_0$ to give
\[
J_\pi=z^{3e/2}H_0Q_\pi.
\]
This does not factor $Q_\pi$ into independent kernels: the two supporting edges are still ordered. Simultaneous reflection and reversal of these two edges preserves id and swap separately, proving $Q_\pi(e,1/z)=Q_\pi(e,z)$, as well as $H_0(e,1/z)=H_0(e,z)$.

The checked $\kappa$ values differ, $-1/16$ for 37 and $-1/32$ for 38. The exact physical functions are consequently
\[
I_{37}=-\frac{Z^3H_0Q_{\rm id}}{16\cosh(3e\ell/2)},\qquad
I_{38}=-\frac{Z^3H_0Q_{\rm swap}}{32\cosh(3e\ell/2)}.
\]
Only the zero-order kernels are required for this task. The adopted evaluations, on the stated real branch and ordinary domains, are
\[
H_0(0,z)=M_{0,1}=-(\ell^2+\pi^2)/2,\quad
Q_{\rm id}(0,z)=M_{0,1}^2/6,\quad Q_{\rm swap}(0,z)=M_{0,1}^2/3.
\]
For $H_0$ the elementary first path primitive and Beta reflection are explicit in [the one-opposite method](single-opposite.md). For $Q$ the first endpoint primitive is $\log[D_z(u,s)/(1-s)]/[(1+z)s-1]$ on $0<t_1<u$; the complete three-variable divided-difference kernel, ordinary endpoints and the prior evaluated two-kernel source are documented in [the two-kernel method](same-channel-two.md) and the [two-loop04 derivation](../../two-loop/4e-04/derivation.md). The ratio1:2 is a property of their evaluated zero-order HLOG coefficients, not an all-$e$ pairing rule.

There is an especially direct complete evaluation route in the current engine: it adopts each distinct original six-variable $J_0$, with its actual original inner-to-outer order. Both full raw and reflected finite HLOG expressions are printed in [37's appendix](../../three-loop/6e-37/appendices/evaluated-jets.md) and [38's appendix](../../three-loop/6e-38/appendices/evaluated-jets.md), not just linked as hidden data. They equal $M_{0,1}^3/6$ and $M_{0,1}^3/3$, respectively. Those full original-$J_0$ tables and the exact scalar normalization alone determine the answer, even without the optional $H_0Q$ factorization. The adopted computational lemma specifies the actual original kernel, its ordinary convergent domain, integration order and complete endpoint periods. Historical HyperInt checked divergences, integrated successive rational/HLOG kernels, and fibred both channels to $z$. The new Python runtime verifies and parses those archived subintegrals and performs the physical assembly; it does not launch another six-dimensional CAS calculation.

At zero, $Z(0)=1$ and $\cosh(0)=1$. Both scalar finite terms are therefore
\[
c_0^{37}=c_0^{38}=-M_{0,1}^3/96=(\ell^2+\pi^2)^3/768,
\quad I_{37,38}=c_0+O(e).
\]
The equality is **only through the requested finite order**. It neither identifies the two original $J$ roots nor the complete physical $I$ functions. For a direct countercheck, $Q_{\rm id}(1,z)=Q_{\rm swap}(1,z)=(1+z)^2/(4z)$ from the ordered-domain volume. Full-$I$ equality would require $Q_{\rm swap}=2Q_{\rm id}$ and is therefore false. This argument evaluates regular kernels at an integer regulator, not Gamma functions at their poles. The current independent controls also give the exact reduced physical difference $(I_{37}/Z^3-I_{38}/Z^3)|_{e=2,z=1/2}=1161/266240$.

The valid full-kernel identity is instead $Q_{\rm id}+Q_{\rm swap}=H_0^2/2$. The product of two unordered squares partitions into four order sectors; exchange of the two kernels pairs these as two id and two swap sectors. Hence
\[
I_{37}+2I_{38}=-\frac{Z^3H_0^3}{32\cosh(3e\ell/2)}.
\]
It is an assembly control and cannot replace either own record. The [original research method](../../../benchmarks/2ne/families/mixed-opposite/method.md) preserves the full source factorization and checks. No new function family, complex continuation or formal proof is claimed; `formal=not_run`.
