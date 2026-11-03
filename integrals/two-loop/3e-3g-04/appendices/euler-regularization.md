# Why the recoil Euler subtraction preserves the original integral

This proof concerns the chart variables $r,u$ of the body recoil, with $s,t,v$ held as the other integration variables. It is unrelated to nonzero Wilson-path endpoint fluxes in the six-face IBP, which have already been kept explicitly.

Use the exact $R,b_e,F,T,X$ defined in the main text and $E=r\partial_r-u\partial_u$. Relative to the flat positive-orthant measure $dr\,du\,ds\,dt\,dv$, this vector field has divergence $1-1=0$. The claim is
\[
\int E(Rb_e)=0,\qquad
\mathcal O(e,z)=\int Rb_e=\int b_e\left[R-\frac{ER+R E\log b_e}{e}\right].\tag{A}
\]
It is first proved on $1/2<\Re e<1$, with $z>0$, and then meromorphically continued. No pointwise cutoff at an overlapping corner is silently dropped.

Here is an absolute-integrability control on that strip. Set $p=2-2\Re e$, so $0<p<1$. In the original cube/simplex variables, replace the signed recoil numerator by the sum of its absolute $ac$ and $bc$ monomials. For an $ac$ term use $S\ge ac D_z(x,w)$ and, where needed, $D_z\ge(1-x)(1-w)$; for a $bc$ term use $S\ge zbcw(1-y)$. Each paired simplex power is $1-p-\Re e=\Re e-1>-1$, the remaining power is $-\Re e>-1$, and each path endpoint power is at worst $-p>-1$. Polynomial path factors only improve this estimate. The absolute original recoil is therefore integrable. This bound also controls the positive monomial sum in the chart.

The rational quantities $r/T$, $u/(1+u)$ and $EF/F$ are bounded on the positive orthant: $F$ is a positive polynomial and $E$ multiplies each monomial by its finite integer scaling degree. Thus $E\log b_e$ is bounded on compact regulator subsets. Also $ER$ is bounded by a fixed multiple of the positive monomial sum just controlled. It follows that both $Rb_e$ and $E(Rb_e)$ are absolutely integrable on this same strip.

Choose smooth cutoffs $\chi_L(\log r,\log u)$ supported on compact rectangles in logarithmic coordinates, tending pointwise to one, with uniformly bounded $E\chi_L$ tending pointwise to zero. Ordinary compact-support integration by parts gives
\[
\int\chi_L E(Rb_e)=-\int Rb_e E\chi_L.
\]
Dominated convergence applies on both sides by the preceding bounds. Its limit is zero on the right, proving (A). This handles $r=0,\infty$ and $u=0,\infty$, including their intersections, together. It is stronger than checking fixed-variable endpoint limits alone.

Expanding the explicit logarithmic derivative gives
\[
E\log b_e=-e+(-1-e)r/T+(1+2e)u/(1+u)+(-2+2e)EF/F=d_0+ed_1.
\]
The main text's $A=-ER-Rd_0$, $B=R(1-d_1)$ follow algebraically, and the regularized density has precisely a simple explicit $1/e$ pole. The dangerous mixed scaling $r\mapsto\rho r$, $u\mapsto u/\rho$ gives $R\sim\rho^{-1}$, $b_e\sim\rho^{1+e}$; the Euler subtraction cancels that leading density because $E(Rb_e)\sim eRb_e$. To expand near $e=0$, convergence of the resulting coefficient kernels is an additional issue, not a consequence of strip convergence alone. The historical `findDivergences` checks reported no remaining coefficient divergence with regulator expansion at zero; HyperInt ran with divergence checking and abort enabled on the exact kernels used for $O_{-1},O_0,O_1,O_2$. These are computational convergence/integration controls, not a Lean theorem or a rigorous complex-domain estimate.

The preserved original implementation is [source-step10_optimized_hlog.mpl](../../../../workflow/checks/tgv04-high-order/runs/positive-epsilon/source-step10_optimized_hlog.mpl). It constructs `Rreg`, verifies the remaining divergence set, builds the explicit pole/finite kernels, and integrates them. The later `recoil1.mpl` and `recoil2.mpl` use the same $A,B,L_b$ to form precisely the higher jets in equation (7). The implementation and the cutoff proof have complementary roles: the proof justifies the equal integral, while the computational checks cover the adopted coefficient kernels and their evaluated jets.
