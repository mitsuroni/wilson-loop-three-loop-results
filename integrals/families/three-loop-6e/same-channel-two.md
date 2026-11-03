# Two ordered opposite propagators in one channel

This route treats 6e-10/11/12/23/24/25. The initial object is each entry's actual ordered $J$, not a product-measure interpretation imposed on its source. Put $e=\epsilon$, $D_z(t,s)=(1-t)(1-s)+zts$, $z>0$, and $\ell=\log z$. The [common normalization](common-normalization.md) fixes the unique Euler/color convention.

A neighboring propagator with source scale $q=1$ or $z$ contributes $q^e A^{e-1}B^{e-1}$ (the canonical 10/11/12 scale is $q=z$), with $A$ its ordinary endpoint and $B$ the distance of its other endpoint from one. On $0<\Re e<1/8$, its primitives are
\[
\int_l^u x^{e-1}dx=(u^e-l^e)/e,\quad
\int_l^u(1-x)^{e-1}dx=[(1-l)^e-(1-u)^e]/e.
\]
The values $0^e=0$ are taken in this initial convergence domain before continuation. Applying the actual bounds, one independent endpoint gives $1/e$, while an OOA chain gives $(1-w^e)/e$ and an OAO chain gives $(w^e-v^e)/e$. Their endpoints and substitutions are spelled out individually in the six derivations.

Let $\Omega_2=\{0<t_1<t_2<1,0<v<w<1\}$. For id pair $(t_1,v),(t_2,w)$, and for swap pair $(t_1,w),(t_2,v)$. The centered integral is
\[
V_{\pi,X}(e,z)=(1+z)^2z^{-e}\int_{\Omega_2}D_1^{e-1}D_2^{e-1}W_X(e)\,d^4s,
\quad W_A=1-w^e,\quad W_B=w^e-v^e.
\]
The canonical reductions are $J_{10}=z^{2e}V_{\mathrm{id},A}/e^2$, $J_{11}=z^{2e}V_{\mathrm{id},B}/e^2$, and $J_{12}=z^{2e}V_{\mathrm{swap},A}/e^2$.

The **actual current runtime primitive has a different normalization**. Its $M(e,z)$ is the uncentered four-variable integral of the two original propagators and the complete normalized gap weight, with no $(1+z)^2$, no $z^{-e}$ and no $1/e^2$. Thus, for these three canonical cases,
\[
V=(1+z)^2z^{-e}M,\quad J=(1+z)^2z^eM/e^2.
\]
The reflected 23/24/25 roots use their own original variables and gap weights, with $J=(1+z)^2M/e^2$. Their raw $M$ tables cannot be silently identified with the canonical $V$ tables. Each own [evaluated-jet appendix](../../three-loop/6e-10/appendices/evaluated-jets.md) and its counterparts prints all $M_0,M_1,M_2$ and reflected values, including endpoint constants.

To justify the expansion, $D$ vanishes only at two corners, with local model $x+y$. For $|e|\le\delta<1/8$, use $D^{-1-\delta}$ and
\[
|(1-w^e)/e|\le|\log w|w^{-\delta},\qquad
|(w^e-v^e)/e|\le|\log v|v^{-\delta}.
\]
After enlarging the ordered domain to a cube, each independent corner kernel integrates first to at worst an endpoint power $s^{-\delta}$; the extra weight contributes at most another $s^{-\delta}$ and logarithms. This is integrable for the stated small disk. Slightly stronger powers also bound the finitely many regulator derivatives. Fubini is justified first on the initial strip and the reduced complete kernels are holomorphic near zero. Since the complete gap weight starts at $e$, $M_0=V_0=0$ and the explicit double pole reduces to a simple pole. The zero is checked as a complete coefficient, not assumed from the graph's appearance.

The historical dimension reduction illustrates the actual integration. Put $u=t_2$, take $(s,r)=(v,w)$ for id or $(w,v)$ for swap, and set $D_1=D_z(u,s)$, $D_2=D_z(u,r)$, $B=1-s$, $d=(1+z)s-1$. The first path primitive is
\[
\int_0^uD_z(t,s)^{e-1}dt=\frac{D_1^e-B^e}{ed};
\]
its $d=0$ limit is $uB^{e-1}$. Keeping that complete divided difference, write $V=eU_1+e^2U_2+O(e^3)$. The actual coefficient kernels are
\[
K_0=(1+z)^2\frac{\log(D_1/B)}{dD_2},\quad K_{\log}=\log D_2+\tfrac12(\log D_1+\log B)-\ell,
\]
\[
U_1=\int_0^1du\int_0^1dw\int_0^w dv\,K_0W_1,\quad
U_2=\int_0^1du\int_0^1dw\int_0^w dv\,K_0(W_2+W_1K_{\log}).
\]
For A, $(W_1,W_2)=(-\log w,-\log^2w/2)$; for B they are $(\log(w/v),(\log^2w-\log^2v)/2)$. Historical HyperInt evaluated these three-dimensional kernels in order $v,u,w$ after $a=1/(1+z)$. The current archived engine route instead integrates the original four-variable logarithmic jets in the exact inner-to-outer order printed for each entry. It is the latter raw output that the new Python runtime consumes. These are two organizations of the same reduction, not two new independent CAS calculations.

There is a complete explicit bridge to the solved values: let $M_j$ be the finite coordinate tables in an entry's appendix. For a canonical case,
\[
U_1=(1+z)^2M_1,\qquad U_2=(1+z)^2(M_2-\ell M_1).
\]
Thus the required $U$ jets are fully specified finite HLOG expressions by those tables; no $U$ is left as an unknown named integral. HyperInt's historical mechanism is partial fractions of the current rational kernel, shuffle expansion of logarithmic products, the primitive $G(a,\mathbf w;q)$ for $G(\mathbf w;q)/(q-a)$, IBP for repeated poles, and the combined ordinary endpoint limits. Fibration to $z$ fixes the alphabet $\{0,-1\}$; the printed empty-word rows retain the endpoint zeta periods. The exact computational lemma is the stated kernel/domain/order equals the full table, rather than a claim that CAS success itself proves the graph.

With $U_j^\pm=[U_j(z)\pm U_j(1/z)]/2$, the physical coefficients of a canonical case are
\[
c_{-6}=\cdots=c_{-2}=0,\qquad
c_{-1}=\kappa U_1^+,\quad c_0=\kappa(U_2^++\ell U_1^-/2).
\]
$Z^3$ and $\operatorname{sech}(3e\ell/2)$ have no linear term, but the crossing-odd jet does contribute to $c_0$.

For an exact formula matching the actual uncentered primitive, let $a$ be the adjacent $z^e$ power ($a=1$ for 10/11/12 and $a=0$ for 23/24/25), $A_j=(1+z)^2M_j(z)$ and $\bar A_j=(1+1/z)^2M_j(1/z)$. Then
\[
c_{-1}=\tfrac\kappa2(A_1+\bar A_1),\qquad
c_0=\tfrac\kappa2[A_2+\bar A_2+(a-3/2)\ell(A_1-\bar A_1)].\tag{S}
\]
This formula and the finite tables determine the whole requested Laurent expansion without a final reference table.

The complete source map, not seven matched coefficients, proves the aliases. For 23 use the canonical 10 coordinates
\[
(t_1,t_2;v,w;b,a)_{10}=(1-s_{42},1-s_{41};1-s_{23},1-s_{22};1-s_{31},1-s_{21})_{23}
\]
at $1/z$; for 24 use the same map with canonical swap case 12; for 25 use canonical 11 with $(v,w;a)=(1-s_{23},1-s_{21};1-s_{22})$. Bounds and absolute Jacobians are preserved. Every propagator becomes its 23/24/25 counterpart divided by $z$, while the canonical numerator at $1/z$ is the target numerator divided by $z^3$. Hence
\[
J_{23}=z^{3e}J_{10}(1/z),\quad J_{24}=z^{3e}J_{12}(1/z),\quad J_{25}=z^{3e}J_{11}(1/z).
\]
The paired $\kappa$ values agree, so the normalized physical $I$ functions are equal pairwise. This does not equate the unreflected $J$ roots.

The original [research method](../../../benchmarks/2ne/families/same-channel-two/method.md) supplies the recorded three-dimensional route, convergence and representation comparisons. The new own-entry bodies and explicit finite appendices supply the present reading chain. The function-family formulas use previously solved $M,N,\tau,E$ as defined in the [common convention](common-normalization.md); their HLOG equality is exact coordinate substitution. The restricted historical basis ranks are not universal nonexistence theorems. Positive-real $z$ and local regulator continuation are the scope; formal remains `not_run`.
