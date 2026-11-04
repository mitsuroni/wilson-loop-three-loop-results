# 4e-Pi1-07: its transverse root to the complete scalar and coloured series

Work at D=4-2e, e=epsilon, fixed positive real z and rho=log z, with the inherited source i0. The accepted starting object is this graph's original prepared C/J, including its connected projection and exactly one Euler compensation. The target is every scalar and connected-colour coefficient at powers -6 through0, with O(e).

## Its actual ordered integral

The insertion is the line with endpoints `s21,s41`; the bare exchange has endpoints `s22,s42`. In the adopted momentum expression these are the Fourier factors D0 Pi1 D0 and D0, respectively. Their actual source separations in the null tangent basis are

```text
Phi_A = (-1 + s21)v_2 + (-1)v_3 + (-s41)v_4
```

```text
Phi_B = (-1 + s22)v_2 + (-1)v_3 + (-s42)v_4
```

The [shared momentum/Fourier root](../../families/three-loop-insertions/pi1.md) fixes the full tensor and its physical sign. The exact parameter integral is
\[
J=\int_{\mathcal D}^{\mathrm{ordered}}d^4s\;
\frac{N}{2(1+e)}F^{-2+2e}\,d_0F_0^{-1+e},
\]

$$ F=s_{21} s_{41} z + s_{21} s_{41} - s_{21} - s_{41} + 1. $$

$$ d=- z - 1. $$

$$ F_0=s_{22} s_{42} z + s_{22} s_{42} - s_{22} - s_{42} + 1. $$

$$ d_0=- z - 1. $$

$$N=-2(1+z)F+(1-2e)z.$$

This compact numerator is the complete transverse source, not a metric-only replacement. It satisfies identity(2) of the shared derivation. The original outer-to-inner order is `s22,s21,s42,s41`, with these actual dependent bounds:

| Variable | Lower | Upper |
|---|---|---|
| $s_{22}$ | $0$ | $1$ |
| $s_{21}$ | $0$ | $s_{22}$ |
| $s_{42}$ | $0$ | $1$ |
| $s_{41}$ | $0$ | $s_{42}$ |

No old cube Jacobian is inserted. On the starting positive-regulator strip the full transverse corner is bounded by D^(-2+2Re e), whose radial integral converges for Re e>0; the specific regroupings below follow from the two-line product bound. Local meromorphic continuation is taken after evaluation, rather than asserting holomorphy of the original divergent source.

## Native normalization and its concrete reduction

This graph has m=2, a1=-1/2. With Z=exp(-gamma_E e)Gamma(1-e), the actual native rational factor is 1/[32e(1-2e)] and
\[
I=\frac{ma_1}{64}\frac{Z^3}{e(1-2e)}K,
\quad K=\frac{J(z)+z^{3e}J(1/z)}{1+z^{3e}},\quad \frac{ma_1}{64}=- \frac{1}{64}.
\]
The unprojected trace is C_F(C_F-C_A/2)B1. The adopted connected projection and the complete B1 multiplier are distinct stages; neither the orbit nor the sign is reinserted after this equation.

Put $(s_1,s_2,t_1,t_2)=(s_{21},s_{22},s_{41},s_{42})$. The insertion has the lower endpoints and the bare line the upper ones. The metric body is $z^{3e/2}T(2e,e;z)$. Integrate the insertion derivative over its *actual rectangle* $[0,s_2]\times[0,t_2]$: the signed boundary mass is
\[
D_z(s_2,t_2)^{2e}-(1-t_2)^{2e}-(1-s_2)^{2e}+1.
\]
Consequently its own raw root is exactly
\[
J_{07}=z^{3e/2}T(2e,e)-\frac{1+z}{4e(1+e)}\int_{[0,1]^2}D_z(s_2,t_2)^{e-1}
[D_z(s_2,t_2)^{2e}-(1-t_2)^{2e}-(1-s_2)^{2e}+1]ds_2dt_2.
\]
Every one of these four square integrals is an evaluated weighted H. To form the physical reciprocal combination, the simultaneous ordered reflection $(s_1,t_1,s_2,t_2)\mapsto(1-s_2,1-t_2,1-s_1,1-t_1)$ sends07 to08 with full scaling $J_{08}=z^{3e}J_{07}(1/z)$. The shared derivation prints both rectangles and evaluates all eight corners into A before assembly.

## Evaluate all remaining kernels before expansion

Define $H_0^{(a)}=H(ae,0;z)$, $H_2=H(e,2e;z)$, $\bar H_2=H(e,2e;1/z)$ and
\[
A=H_0^{(3)}+\cosh(e\rho)(H_0^{(1)}-H_2-\bar H_2)-\sinh(e\rho)(H_2-\bar H_2).
\]
After the explicit eight-corner evaluation, the *whole* paired root is
\[
K=\frac{T(2e,e)+T(e,2e)+A/[2e(1+e)]}{2\cosh(3e\rho/2)}.
\]
The needed mixed-body coefficients are $2\rho_0+3e\rho_1+O(e^2)$, with $\rho_0=M_{0,1}^2/6$. The thirteen-term full HLOG expression of rho1, its exact four-variable definition and its reduced logarithmic integration kernel are printed in [the moment appendix](../../families/three-loop-insertions/appendices/pi1-moments.md). It is an adopted solved function, not an unresolved remainder. Substitution gives the explicit finite assembly
\[
c_0=-\rho^2M_{0,1}/512-M_{0,1}^2/192+M_{0,3}/384+(M_{1,1}+M_{2,1})/256-3\rho_1/128.
\]
This uses only the indicated first derivative of the mixed body; it does not replace its general-e value by a same-exponent body.

The actual finite budget is: **H through2 and the mixed T sum through1; the Gamma constant suffices**. The shared equation(5) expands H into a finite polynomial in the original no-factorial M moments; every needed HLOG value, including endpoint constants and rho1 where needed, is explicit. No code-success assertion or new symbol supplies a missing mathematical integration.

## The highest coefficient and the regular colour operator

Let k_r=[e^r]K for the explicit evaluated K above, and let a_n=[e^n]Z^3/(1-2e). The shared Gamma recurrence fixes a_n exactly. The full physical rule is c_p=(ma1/64)sum(a_n k_(p+1-n)), restricted to existing source orders. In particular, for this graph the finite term is the explicit finite sum
\[
c_0=\frac{ma_1}{64}\sum_{n=0}^{1}a_n k_{1-n}.
\]
This is coefficient extraction from the solved functions just given; it contains no unknown integral. It retains the pole-times-positive-order terms. The c0 expression printed below is the resulting exact moment form, with the complete HLOG form available by the explicit shared dictionary.

The two connected colour multipliers are phi_A=(10-6e)/(3-2e) and phi_f=-8(1-e)/(3-2e), multiplying CA^2 CF and CA CF TF nf. Their constants are10/3 and-8/3; for n>=1 their Taylor coefficients are (2/9)(2/3)^(n-1) and (8/9)(2/3)^(n-1). Therefore the coloured finite term is sum(phi_X,n c_-n,n=0,...,1), not phi_X,0 c0 alone. The other two connected channels vanish in this projection; this does not delete the unprojected source trace.

## Complete coefficients

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=0.$$

$$c_{-2}(z)=0.$$

$$c_{-1}(z)=\frac{M_{1,1}}{256} - \frac{\rho_0}{64}.$$

$$c_{0}(z)=- \frac{M_{0,1} \rho^{2}}{512} + \frac{M_{0,3}}{384} + \frac{M_{1,1}}{256} + \frac{M_{2,1}}{256} - \frac{\rho_0}{32} - \frac{3 \rho_1}{128}.$$

The remainder is O(e). Both complete coloured finite coefficients are:

| External monomial | Complete coefficient at e^0 |
|---|---|
| `CA^2*CF` | $- \frac{5 M_{0,1} \rho^{2}}{768} + \frac{5 M_{0,3}}{576} + \frac{M_{1,1}}{72} + \frac{5 M_{2,1}}{384} - \frac{31 \rho_0}{288} - \frac{5 \rho_1}{64}$ |
| `CA*CF*TF*nf` | $\frac{M_{0,1} \rho^{2}}{192} - \frac{M_{0,3}}{144} - \frac{M_{1,1}}{144} - \frac{M_{2,1}}{96} + \frac{5 \rho_0}{72} + \frac{\rho_1}{16}$ |

[result.json](result.json) gives every coloured order, including structural zeros; the finite-convolution rule above reconstructs them. M/N retain the original definitions, N_kj=M_k,j+1-rho M_kj. The indexed rho0/rho1 are the explicitly evaluated ordered-body Taylor coefficients and differ from the shape logarithm rho. [Reproduce](reproduce.md) and [validation](validation.md) state precisely which computations are rerun and which source/function evaluations are adopted. The claim is restricted to the stated positive-real branch and local regulator continuation; formal=not_run.
