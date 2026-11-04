# 4e-Pi1-09: its transverse root to the complete scalar and coloured series

Work at D=4-2e, e=epsilon, fixed positive real z and rho=log z, with the inherited source i0. The accepted starting object is this graph's original prepared C/J, including its connected projection and exactly one Euler compensation. The target is every scalar and connected-colour coefficient at powers -6 through0, with O(e).

## Its actual ordered integral

The insertion is the line with endpoints `s21,s31`; the bare exchange has endpoints `s22,s41`. In the adopted momentum expression these are the Fourier factors D0 Pi1 D0 and D0, respectively. Their actual source separations in the null tangent basis are

```text
Phi_A = (-1 + s21)v_2 + (-s31)v_3
```

```text
Phi_B = (-1 + s22)v_2 + (-1)v_3 + (-s41)v_4
```

The [shared momentum/Fourier root](../../families/three-loop-insertions/pi1.md) fixes the full tensor and its physical sign. The exact parameter integral is
\[
J=\int_{\mathcal D}^{\mathrm{ordered}}d^4s\;
\frac{N}{2(1+e)}F^{-2+2e}\,d_0F_0^{-1+e},
\]

$$ F=- s_{21} s_{31} + s_{31}. $$

$$ d=1. $$

$$ F_0=s_{22} s_{41} z + s_{22} s_{41} - s_{22} - s_{41} + 1. $$

$$ d_0=- z - 1. $$

$$N=2dF.$$

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

The insertion endpoints are $s_{21},s_{31}$, with the former below the bare endpoint $x=s_{22}$; the other bare endpoint is $t=s_{41}$. The full adjacent insertion gives the free $s_{31}$ primitive $1/(2e)$ and the ordered primitive $[1-(1-x)^{2e}]/(2e)$. Thus its own reduced root is
\[
J=-\frac{1+z}{4e^2(1+e)}\int_{[0,1]^2}[1-(1-x)^{2e}]D_z(x,t)^{e-1}dxdt.
\]
The concrete map
\[
(s_{21},s_{41},s_{42},s_{31})_{04}=(1-s_{41},1-s_{22},1-s_{21},1-s_{31})_{09}
\]
at reciprocal z preserves the actual ordered bound, insertion labels and unit absolute Jacobian. Its line scaling gives $J_{09}=z^{3e}J_{04}(1/z)$, hence equal complete physical I after the checked equal m,a1. This proves an I identity, not equality of the unreflected J.

## Evaluate all remaining kernels before expansion

Put $C=\cosh(3e\rho/2)$, $H_0^{(1)}=H(e,0;z)$, $H_2=H(e,2e;z)$ and $\bar H_2=H(e,2e;1/z)$. The actual paired root is
\[
K=\frac{\cosh(e\rho)[H_0^{(1)}-(H_2+\bar H_2)/2]-\sinh(e\rho)(H_2-\bar H_2)/2}{4e^2(1+e)C}.
\]
Both the even and odd H channels remain before coefficient extraction. The evaluated finite-M rule and complete HLOG values supply all required terms, including the order3 terms multiplied by the explicit poles.

The actual finite budget is: **weighted H through3, remaining regular factors through2**. The shared equation(5) expands H into a finite polynomial in the original no-factorial M moments; every needed HLOG value, including endpoint constants and rho1 where needed, is explicit. No code-success assertion or new symbol supplies a missing mathematical integration.

## The highest coefficient and the regular colour operator

Let k_r=[e^r]K for the explicit evaluated K above, and let a_n=[e^n]Z^3/(1-2e). The shared Gamma recurrence fixes a_n exactly. The full physical rule is c_p=(ma1/64)sum(a_n k_(p+1-n)), restricted to existing source orders. In particular, for this graph the finite term is the explicit finite sum
\[
c_0=\frac{ma_1}{64}\sum_{n=0}^{2}a_n k_{1-n}.
\]
This is coefficient extraction from the solved functions just given; it contains no unknown integral. It retains the pole-times-positive-order terms. The c0 expression printed below is the resulting exact moment form, with the complete HLOG form available by the explicit shared dictionary.

The two connected colour multipliers are phi_A=(10-6e)/(3-2e) and phi_f=-8(1-e)/(3-2e), multiplying CA^2 CF and CA CF TF nf. Their constants are10/3 and-8/3; for n>=1 their Taylor coefficients are (2/9)(2/3)^(n-1) and (8/9)(2/3)^(n-1). Therefore the coloured finite term is sum(phi_X,n c_-n,n=0,...,2), not phi_X,0 c0 alone. The other two connected channels vanish in this projection; this does not delete the unprojected source trace.

## Complete coefficients

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=0.$$

$$c_{-2}(z)=\frac{M_{1,1}}{128}.$$

$$c_{-1}(z)=- \frac{M_{0,1} \rho^{2}}{256} + \frac{M_{0,3}}{256} + \frac{M_{1,1}}{128} + \frac{M_{2,1}}{128}.$$

$$c_{0}(z)=- \frac{M_{0,1} \rho^{2}}{256} + \frac{M_{0,3}}{256} - \frac{11 M_{1,1} \rho^{2}}{1024} + \frac{\pi^{2} M_{1,1}}{512} + \frac{3 M_{1,1}}{128} + \frac{19 M_{1,3}}{3072} + \frac{M_{2,1}}{128} + \frac{13 M_{3,1}}{3072}.$$

The remainder is O(e). Both complete coloured finite coefficients are:

| External monomial | Complete coefficient at e^0 |
|---|---|
| `CA^2*CF` | $- \frac{M_{0,1} \rho^{2}}{72} + \frac{M_{0,3}}{72} - \frac{55 M_{1,1} \rho^{2}}{1536} + \frac{5 \pi^{2} M_{1,1}}{768} + \frac{35 M_{1,1}}{432} + \frac{95 M_{1,3}}{4608} + \frac{M_{2,1}}{36} + \frac{65 M_{3,1}}{4608}$ |
| `CA*CF*TF*nf` | $\frac{M_{0,1} \rho^{2}}{144} - \frac{M_{0,3}}{144} + \frac{11 M_{1,1} \rho^{2}}{384} - \frac{\pi^{2} M_{1,1}}{192} - \frac{11 M_{1,1}}{216} - \frac{19 M_{1,3}}{1152} - \frac{M_{2,1}}{72} - \frac{13 M_{3,1}}{1152}$ |

[result.json](result.json) gives every coloured order, including structural zeros; the finite-convolution rule above reconstructs them. M/N retain the original definitions, N_kj=M_k,j+1-rho M_kj. The indexed rho0/rho1 are the explicitly evaluated ordered-body Taylor coefficients and differ from the shape logarithm rho. [Reproduce](reproduce.md) and [validation](validation.md) state precisely which computations are rerun and which source/function evaluations are adopted. The claim is restricted to the stated positive-real branch and local regulator continuation; formal=not_run.
