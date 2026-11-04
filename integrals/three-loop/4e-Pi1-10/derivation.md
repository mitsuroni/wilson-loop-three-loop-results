# 4e-Pi1-10: its transverse root to the complete scalar and coloured series

Work at D=4-2e, e=epsilon, fixed positive real z and rho=log z, with the inherited source i0. The accepted starting object is this graph's original prepared C/J, including its connected projection and exactly one Euler compensation. The target is every scalar and connected-colour coefficient at powers -6 through0, with O(e).

## Its actual ordered integral

The insertion is the line with endpoints `s22,s41`; the bare exchange has endpoints `s21,s31`. In the adopted momentum expression these are the Fourier factors D0 Pi1 D0 and D0, respectively. Their actual source separations in the null tangent basis are

```text
Phi_A = (-1 + s22)v_2 + (-1)v_3 + (-s41)v_4
```

```text
Phi_B = (-1 + s21)v_2 + (-s31)v_3
```

The [shared momentum/Fourier root](../../families/three-loop-insertions/pi1.md) fixes the full tensor and its physical sign. The exact parameter integral is
\[
J=\int_{\mathcal D}^{\mathrm{ordered}}d^4s\;
\frac{N}{2(1+e)}F^{-2+2e}\,d_0F_0^{-1+e},
\]

$$ F=s_{22} s_{41} z + s_{22} s_{41} - s_{22} - s_{41} + 1. $$

$$ d=- z - 1. $$

$$ F_0=- s_{21} s_{31} + s_{31}. $$

$$ d_0=1. $$

$$N=-2(1+z)F+(1-2e)z.$$

This compact numerator is the complete transverse source, not a metric-only replacement. It satisfies identity(2) of the shared derivation. The original outer-to-inner order is `s22,s21,s31,s41`, with these actual dependent bounds:

| Variable | Lower | Upper |
|---|---|---|
| $s_{22}$ | $0$ | $1$ |
| $s_{21}$ | $0$ | $s_{22}$ |
| $s_{31}$ | $0$ | $1$ |
| $s_{41}$ | $0$ | $1$ |

No old cube Jacobian is inserted. On the starting positive-regulator strip the full transverse corner is bounded by D^(-2+2Re e), whose radial integral converges for Re e>0; the specific regroupings below follow from the two-line product bound. Local meromorphic continuation is taken after evaluation, rather than asserting holomorphy of the original divergent source.

## Native normalization and its concrete reduction

This graph has m=4, a1=-1/2. With Z=exp(-gamma_E e)Gamma(1-e), the actual native rational factor is 1/[32e(1-2e)] and
\[
I=\frac{ma_1}{64}\frac{Z^3}{e(1-2e)}K,
\quad K=\frac{J(z)+z^{3e}J(1/z)}{1+z^{3e}},\quad \frac{ma_1}{64}=- \frac{1}{32}.
\]
The unprojected trace is C_F(C_F-C_A/2)B1. The adopted connected projection and the complete B1 multiplier are distinct stages; neither the orbit nor the sign is reinserted after this equation.

Keep the insertion endpoints $x=s_{22}$, $t=s_{41}$. The bare variables $s_{31}\in[0,1]$ and $s_{21}\in[0,x]$ supply $[1-(1-x)^e]/e^2$; the full insertion density is the same explicit transverse D kernel. This gives its own reduced root
\[
J=\frac1{e^2}\int_{[0,1]^2}[1-(1-x)^e]
\left[-(1+z)D_z^{2e-1}+\frac{\partial_x\partial_tD_z^{2e}}{4e(1+e)}\right]dxdt.
\]
For the concrete reflection, set
\[
(s_{21},s_{41},s_{42},s_{31})_{03}=(1-s_{41},1-s_{22},1-s_{21},1-s_{31})_{10}
\]
at the reciprocal shape. Its absolute Jacobian is1, its bounds preserve $s_{21}<s_{22}$, and insertion/bare labels stay attached to their lines. The insertion and bare Fourier kernels scale by $z^{-2e}$ and $z^{-e}$ respectively, proving $J_{10}=z^{3e}J_{03}(1/z)$. Therefore the physical crossing K, and hence I, equals03; no equality of the unreflected bare J is asserted. The nonzero03 boundary mass remains in that mapped object.

## Evaluate all remaining kernels before expansion

With $C=\cosh(3e\rho/2)$, $H_0^{(2)}=H(2e,0;z)$, $H_h=H(2e,e;z)$ and $\bar H_h=H(2e,e;1/z)$, the exact physical paired root is
\[
K=\frac{\cosh(e\rho/2)[H_0^{(2)}-(H_h+\bar H_h)/2]-\sinh(e\rho/2)(H_h-\bar H_h)/2}{e^2 C}
+\frac{1+3(1-Q)\cosh(e\rho/2)/C}{12e^3(1+e)}.
\]
The second summand is the complete derivative-boundary contribution. The H functions are evaluated by(5) of the shared derivation; their required ten M values are fully printed in its appendix.

The actual finite budget is: **weighted H through3, Q through4, Gamma/rational through4**. The shared equation(5) expands H into a finite polynomial in the original no-factorial M moments; every needed HLOG value, including endpoint constants and rho1 where needed, is explicit. No code-success assertion or new symbol supplies a missing mathematical integration.

## The highest coefficient and the regular colour operator

Let k_r=[e^r]K for the explicit evaluated K above, and let a_n=[e^n]Z^3/(1-2e). The shared Gamma recurrence fixes a_n exactly. The full physical rule is c_p=(ma1/64)sum(a_n k_(p+1-n)), restricted to existing source orders. In particular, for this graph the finite term is the explicit finite sum
\[
c_0=\frac{ma_1}{64}\sum_{n=0}^{4}a_n k_{1-n}.
\]
This is coefficient extraction from the solved functions just given; it contains no unknown integral. It retains the pole-times-positive-order terms. The c0 expression printed below is the resulting exact moment form, with the complete HLOG form available by the explicit shared dictionary.

The two connected colour multipliers are phi_A=(10-6e)/(3-2e) and phi_f=-8(1-e)/(3-2e), multiplying CA^2 CF and CA CF TF nf. Their constants are10/3 and-8/3; for n>=1 their Taylor coefficients are (2/9)(2/3)^(n-1) and (8/9)(2/3)^(n-1). Therefore the coloured finite term is sum(phi_X,n c_-n,n=0,...,4), not phi_X,0 c0 alone. The other two connected channels vanish in this projection; this does not delete the unprojected source trace.

## Complete coefficients

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=- \frac{1}{384}.$$

$$c_{-3}(z)=- \frac{1}{384}.$$

$$c_{-2}(z)=\frac{M_{1,1}}{64} - \frac{5 \pi^{2}}{1536} - \frac{1}{128}.$$

$$c_{-1}(z)=- \frac{M_{0,1} \rho^{2}}{256} + \frac{M_{0,3}}{256} + \frac{M_{1,1}}{32} + \frac{5 M_{2,1}}{256} - \frac{5 \pi^{2}}{1536} - \frac{5}{384} + \frac{17 \zeta\left(3\right)}{384}.$$

$$c_{0}(z)=- \frac{M_{0,1} \rho^{2}}{128} + \frac{M_{0,3}}{128} - \frac{11 M_{1,1} \rho^{2}}{512} + \frac{\pi^{2} M_{1,1}}{256} + \frac{M_{1,1}}{16} + \frac{13 M_{1,3}}{1536} + \frac{5 M_{2,1}}{128} + \frac{19 M_{3,1}}{1536} + \frac{\pi^{2} \rho^{2}}{384} - \frac{7 \pi^{4}}{4096} - \frac{5 \pi^{2}}{512} - \frac{11}{384} + \frac{17 \zeta\left(3\right)}{384}.$$

The remainder is O(e). Both complete coloured finite coefficients are:

| External monomial | Complete coefficient at e^0 |
|---|---|
| `CA^2*CF` | $- \frac{31 M_{0,1} \rho^{2}}{1152} + \frac{31 M_{0,3}}{1152} - \frac{55 M_{1,1} \rho^{2}}{768} + \frac{5 \pi^{2} M_{1,1}}{384} + \frac{47 M_{1,1}}{216} + \frac{65 M_{1,3}}{2304} + \frac{155 M_{2,1}}{1152} + \frac{95 M_{3,1}}{2304} + \frac{5 \pi^{2} \rho^{2}}{576} - \frac{35 \pi^{4}}{6144} - \frac{175 \pi^{2}}{5184} - \frac{583}{5832} + \frac{17 \zeta\left(3\right)}{108}$ |
| `CA*CF*TF*nf` | $\frac{5 M_{0,1} \rho^{2}}{288} - \frac{5 M_{0,3}}{288} + \frac{11 M_{1,1} \rho^{2}}{192} - \frac{7 M_{1,1}}{54} - \frac{\pi^{2} M_{1,1}}{96} - \frac{13 M_{1,3}}{576} - \frac{25 M_{2,1}}{288} - \frac{19 M_{3,1}}{576} - \frac{\pi^{2} \rho^{2}}{144} - \frac{17 \zeta\left(3\right)}{216} + \frac{341}{5832} + \frac{55 \pi^{2}}{2592} + \frac{7 \pi^{4}}{1536}$ |

[result.json](result.json) gives every coloured order, including structural zeros; the finite-convolution rule above reconstructs them. M/N retain the original definitions, N_kj=M_k,j+1-rho M_kj. The indexed rho0/rho1 are the explicitly evaluated ordered-body Taylor coefficients and differ from the shape logarithm rho. [Reproduce](reproduce.md) and [validation](validation.md) state precisely which computations are rerun and which source/function evaluations are adopted. The claim is restricted to the stated positive-real branch and local regulator continuation; formal=not_run.
