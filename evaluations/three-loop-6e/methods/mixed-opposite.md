# One and two opposite propagators: 6e-37 and 6e-38

Both roots contain one reversed-channel propagator with independent endpoints a,b, and two ordered same-channel propagators. Let Omega2={0<t1<t2<1,0<s1<s2<1}. Their actual source is

\[
J_\pi=-(1+z)^3\int_{[0,1]^2\times\Omega_2}
\widetilde D_z(a,b)^{\epsilon-1}
D_z(t_1,s_{\pi(1)})^{\epsilon-1}
D_z(t_2,s_{\pi(2)})^{\epsilon-1}\,d^6s.
\tag{22}
\]

The pairing is id for 37 and swap for 38. There are no adjacent endpoints and no explicit epsilon^-2 factor. The independent a,b domain is what permits a product reduction; the two ordered same-channel propagators remain coupled.

At fixed z>0, each opposite kernel has only the integrable corner model (x+y)^(epsilon-1). A cube bound factors into three two-dimensional corner bounds, uniformly for epsilon in a sufficiently small disk around zero. Thus J and I are holomorphic there, all six negative-power coefficients vanish, and c0 is the ordinary original integral evaluated at epsilon=0.

Define H0 as in [the single-opposite method](single-opposite.md), and

\[
Q_\pi=(1+z)^2z^{-\epsilon}
\int_{\Omega_2}D_z(t_1,s_{\pi(1)})^{\epsilon-1}
D_z(t_2,s_{\pi(2)})^{\epsilon-1}\,d^4s.
\tag{23}
\]

Reflection (a,b)->(1-a,1-b) converts Dtilde to D, with unit absolute Jacobian. Factorization of the independent domain now gives

\[
J_\pi=z^{3\epsilon/2}H_0Q_\pi.
\]

H0 and each zero-insertion Q_pi are individually crossing even. The latter follows by simultaneous reflection and reversal on its two ordered edges; id and swap each remain themselves. The checked kappa factors differ: -1/16 for 37 and -1/32 for 38. Therefore

\[
I_{37}=-\frac{\mathcal G^3H_0Q_{\mathrm{id}}}{16\cosh(3\epsilon\rho/2)},\qquad
I_{38}=-\frac{\mathcal G^3H_0Q_{\mathrm{swap}}}{32\cosh(3\epsilon\rho/2)}.
\tag{24}
\]

## The finite coefficient and the scope of its equality

The already evaluated zero-order mothers are

\[
H_0(0,z)=M_{0,1},\qquad Q_{\mathrm{id}}(0,z)=M_{0,1}^2/6,
\qquad Q_{\mathrm{swap}}(0,z)=M_{0,1}^2/3.
\]

Only these q=0 values are needed, since the prefactors in (24) are regular. The full finite answer is

\[
c_0^{37}=c_0^{38}=-\frac{M_{0,1}^3}{96}
=\frac{(\log^2z+\pi^2)^3}{768}.
\tag{25}
\]

This equality is restricted to c0. The complete functions in (24) generally differ. To see the distinction without evaluating singular physical-I factors at an integer regulator, compare the regular Q kernels: at epsilon=1 both equal (1+z)^2/(4z), whereas equality of the complete I functions would require Q_swap=2 Q_id. The contradiction concerns the regular kernel identity, and does not evaluate Gamma(1-epsilon) at its pole.

Partitioning the product domain supplies the valid complete-kernel relation

\[
Q_{\mathrm{id}}+Q_{\mathrm{swap}}=H_0^2/2,
\]

or I37+2 I38=-Gcal^3 H0^3/[32 cosh(3 epsilon ell/2)]. This is a useful assembly control but cannot replace the separate source and coefficient records.

The original family computation reused existing H0 and Q zero-order HyperInt outputs. P2 instead replays the current engine's archived direct epsilon=0 integration of each ordered root and then repeats its physical normalization. Both are documented evaluation routes. Neither is called a new cold CAS computation in this batch. Formula (25), tau_(0,0)^3 and complete HLOG coordinates are equivalent representations of the same solved coefficient, with no new function required.

The preserved [historical method](../evidence/historical/benchmarks/2ne/families/mixed-opposite/method.md) gives the full source factorization and the existing all-epsilon kernel representation. The [common scope](common-normalization.md) keeps physical complex continuation and formal proof separate.
