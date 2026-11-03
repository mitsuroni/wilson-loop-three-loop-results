# Three ordered opposite propagators in one channel

The three roots 20/21/22 have no adjacent propagator. Write $t_i=s_{2i}$, $s_i=s_{4i}$, and
\[
\Omega_3=\{0<t_1<t_2<t_3<1,\ 0<s_1<s_2<s_3<1\},\quad
J_\pi=-(1+z)^3\int_{\Omega_3}\prod_{i=1}^3D_z(t_i,s_{\pi(i)})^{e-1}d^3t\,d^3s.
\]
The permutations are respectively123,132,213; $\kappa=-1/16,-1/32,-1/32$. Each own article retains its initial nested variable order and all dependent bounds.

Every endpoint occurs in exactly one kernel. For fixed $z>0$, each zero corner is comparable to $x+y$, and $\int_0^\delta\int_0^\delta(x+y)^{-1-\eta}dxdy<\infty$ for $\eta<1$. The cube enlargement factors into three such bounds, with slightly stronger exponents to dominate logarithmic derivatives. This proves a locally uniform majorant for $e$ near zero and therefore holomorphy, as well as the Fubini bridge used in the reductions. No finite-part prescription replaces the original integral. All six negative Laurent coefficients are zero; the only needed primitive is the ordinary original $J_0$.

To separate scaling from ordering, define $R_\pi=z^{-3e/2}J_\pi$. Exchanging the two edges gives $R_\pi=R_{\pi^{-1}}$. Reflect and reverse both ordered edges,
\[
t_i\mapsto1-t_{4-i},\qquad s_i\mapsto1-s_{4-i},\quad r(i)=4-i\quad\text{(one-line permutation }(3,2,1)\text{)}.
\]
Their Jacobians have absolute value one and $D_z(1-t,1-s)=zD_{1/z}(t,s)$, so
$R_\pi(1/z)=R_{r\pi r}(z)$. Thus123 is self-crossing and132/213 exchange. The exact physical functions are
\[
I_{20}=-\frac{Z^3R_{123}}{16\cosh(3e\ell/2)},\qquad
I_{21}=I_{22}=-\frac{Z^3(R_{132}+R_{213})}{64\cosh(3e\ell/2)}.
\]
This is a full-$e$ equality of the two physical $I$ functions and does not assert $J_{21}=J_{22}$.

For a readable actual reduction of the finite integral, put $a=1/(1+z)$ and $d_a(u,v)=a(1-u)(1-v)+(1-a)uv$. At $e=0$, $(1+z)^3/(D_1D_2D_3)=1/(d_1d_2d_3)$. The actual smallest/largest endpoint primitives are
\[
\int_0^u\frac{dt}{d_a(t,v)}=\frac{L_-(u,v)}{v-a},\qquad
\int_u^1\frac{dt}{d_a(t,v)}=\frac{L_+(u,v)}{v-a},
\]
\[
L_-=\log d_a(u,v)-\log a-\log(1-v),\qquad
L_+=\log(1-a)+\log v-\log d_a(u,v).
\]
At $v=a$ both formulas mean the ordinary continuous divided difference. Writing $u=t_2$ and $(v,w,r)=(s_1,s_2,s_3)$ gives the actual four-dimensional finite kernel
\[
R_{123}(0,z)=-\int_0^1du\int_0^1dw\int_0^w dv\int_w^1dr\,
\frac{L_-(u,v)L_+(u,r)}{(v-a)(r-a)d_a(u,w)}.
\]
For132 the additional endpoint integration yields
\[
R_{132}(0,z)=-\int_0^1du\int_0^1dw\int_0^w dv\,
\frac{L_-(u,v)L_+(u,w)L_+(w,u)}{(v-a)(w-a)(u-a)}.
\]
Every apparent interior pole is removable in its complete logarithmic difference. The historical family sessions integrated these complete kernels, not separated divergent terms. The current engine has a distinct archived six-variable original-$J_0$ integration with the original inner-to-outer order. The new batch consumes the latter exact request/script/raw output, which is printed separately for each entry in its [evaluated-jet appendix](../../three-loop/6e-20/appendices/evaluated-jets.md) and the21/22 counterparts. These tables give **original $J_0$**, not an already normalized $I$ and not an assumed all-$e$ $R$. At zero $R_0=J_0$, which is the exact bridge between the two presentations.

The adopted computational lemma consists of each stated original kernel, its ordinary domain and order, and its complete finite raw and reflected HLOG tables. The historical implementation checks divergence, partial-fractions and integrates hyperlogarithms successively, combines ordinary endpoint limits, then fibers to $z$ and reflects on the same positive-real branch. Zeta period rows in the printed tables fix the endpoint constants. The new replay validates/reads those primitives and applies $c_0=\kappa[J_0(z)+J_0(1/z)]/2$; it does not cold-run HyperInt.

In compact existing coordinates set $T_3=R_{123}(0,z)$. The fully evaluated finite tables imply, by exact shuffle-coordinate subtraction,
\[
R_{132}^+(0,z)=M_{0,1}^3/36-T_3,\quad
T_3=2E-\tfrac12M_{0,1}M_{2,1}+\tfrac12M_{1,1}^2+\tfrac{\pi^2}{9}M_{0,3}-\tfrac1{12}M_{0,5}.
\]
This names an already printed finite HLOG expression, not an unknown three-kernel function. Consequently
\[
c_0^{20}=-T_3/16,\qquad c_0^{21}=c_0^{22}=T_3/32-M_{0,1}^3/1152.
\]
Here $E$ is the previously evaluated two-kernel endpoint insertion. Its independent structural definition is
\[
E=\frac{\ell(1+z)^2}{2}\int_{\Omega_2}
\frac{\log[s_2/(1-s_1)]}{D_z(t_1,s_1)D_z(t_2,s_2)}d^4s.
\]
Equivalently $E=\ell[e](Q_{\mathrm{id};0,1}(e,z)-Q_{\mathrm{id};0,1}(e,1/z))/2$, where $Q_{\mathrm{id};0,1}=(1+z)^2z^{-e}\int_{\Omega_2}D_1^{e-1}D_2^{e-1}s_2^e$. Its adopted HLOG evaluation and the exact conversion are part of the common function dictionary. The whole requested answer is already fixed by the explicit $J_0$ tables even before this optional compact notation is used. No claim about an all-$e$ reduction to this finite $M/E$ space follows.

For the original detailed calculation and its source-specific convergence/representation records see the [research method](../../../benchmarks/2ne/families/same-channel-three/method.md). Complex kinematics, $z=-1$ boundary continuation and integral-level formal proof are separate; `formal=not_run`.
