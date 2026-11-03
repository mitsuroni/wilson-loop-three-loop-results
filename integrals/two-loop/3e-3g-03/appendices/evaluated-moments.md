# Evaluated real moments and their HLOG convention

This appendix is an exact adopted function lemma used by the present low-loop entries, rather than a new unknown subintegral. It fixes the definition, branch, source order and normalization needed for its application.

For fixed $z>0$, put
\[
\ell=\log z,\quad r=\frac z{1+z},\quad W=\log[\beta(1-\beta)],\quad R=\log\frac{1-\beta}{\beta},\quad L=\ell+R,
\]
\[
M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-r}d\beta,\qquad
N_{k,j}=M_{k,j+1}-\ell M_{k,j}\quad(k\ge0,j\ge1).
\]
Powers are raw, without factorials. Since $L(r)=0$, the interior singularity is removable; endpoint logarithms are integrable. No new principal-value prescription is attached to these functions. The exact generator and its finite coefficient rule are
\[
\mathcal V(A,B;z)=\int_0^1\frac{e^{AW}\sinh(BL)}{B(\beta-r)}d\beta,
\qquad [x^q]\mathcal V(x,x)=\sum_{j=0}^{\lfloor q/2\rfloor}\frac{M_{q-2j,2j+1}}{(q-2j)!(2j+1)!}.
\]
The expression at $B=0$ is the continuous limit. The notation in the shared runtime is $V(\delta)=\mathcal V(\delta/2,\delta/2)$, so $V(2e)=\mathcal V(e,e)$.

The adopted HLOG convention is
\[
G(a,\mathbf w;z)=\int_0^z\frac{G(\mathbf w;t)}{t-a}dt,\quad G(;z)=1,\quad G(0^n;z)=\frac{\ell^n}{n!}.
\]
All the moments needed here have a finite rational combination of words in $\{0,-1\}$ and zeta periods. Multiplication uses shuffle including repeated-word multiplicity. For example, with $G_{\mathbf w}=G(\mathbf w;z)$,
\[
M_{0,1}=-G_{0,0}-3\zeta_2,
\quad M_{1,1}=2G_{-1,0,0}+6\zeta_2G_{-1}+2\zeta_3-3\zeta_2G_0-G_{0,0,0},
\]
\[
\begin{aligned}
M_{2,1}={}&-\frac{33}{5}\zeta_2^2-8\zeta_3G_{-1}-24\zeta_2G_{-1,-1}-8G_{-1,-1,0,0}
+12\zeta_2G_{-1,0}+4G_{-1,0,0,0}\\
&+4\zeta_3G_0+12\zeta_2G_{0,-1}+4G_{0,-1,0,0}-4\zeta_2G_{0,0}-2G_{0,0,0,0}.
\end{aligned}
\]
The higher moments appearing in the complete coefficient tables are **the same evaluated family**, not additional unevaluated graph kernels. Their original exact integration and weight-five/six dictionary are in [single-opposite high-weight](../../../../benchmarks/2ne/theory/generation/single-opposite-high-weight-20261016/dictionary.md), [weight-six coefficients](../../../../benchmarks/2ne/theory/generation/single-opposite-weight-six-20261016/hlog-coefficients.md), and the [original function-source derivation](../../../../benchmarks/2ne/theory/generation/mothers.md). The data index binds the corresponding exact native HLOG maps used by replay. Adopting these evaluations is part of the declared input to the present calculation, not a claim that the replay independently integrates them.

Their integration rule is concrete: expand the logarithm products by shuffle, partial-fraction each rational kernel, integrate a simple pole by
$\int G(\mathbf w;\beta)d\beta/(\beta-a)=G(a,\mathbf w;\beta)$, reduce higher poles by IBP and combine the full endpoint expansions before the regularized endpoint value is taken. The removable point $\beta=r$ is combined before endpoint limiting. Fibration to $z$ then yields the displayed alphabet. The source proof fixes the endpoint zeta periods; a symbol without these constants would not determine this lemma. This entry uses the exact finite dictionary, not an ansatz that all such moments must lie in it.

For 03 through $e^2$, only $M_{0,1}$, $M_{1,1}$, $M_{2,1}$, $M_{3,1}$, $M_{4,1}$, $M_{1,3}$, $M_{2,3}$ and the elementary $M_{0,3},M_{0,5}$ can occur before reduction. For Pi1 opposite through $e^2$, the source stops one degree earlier. Low-index reductions include $M_{0,3}=-M_{0,1}^2$ and $\ell^2=-2M_{0,1}-6\zeta_2$. These identities account for the compact moment representatives in the final tables.
