# 4e-Pi1-05: its transverse root to the complete scalar and coloured series

Work at D=4-2e, e=epsilon, fixed positive real z and rho=log z, with the inherited source i0. The accepted starting object is this graph's original prepared C/J, including its connected projection and exactly one Euler compensation. The target is every scalar and connected-colour coefficient at powers -6 through0, with O(e).

## Its actual ordered integral

The insertion is the line with endpoints `s21,s32`; the bare exchange has endpoints `s31,s41`. In the adopted momentum expression these are the Fourier factors D0 Pi1 D0 and D0, respectively. Their actual source separations in the null tangent basis are

```text
Phi_A = (-1 + s21)v_2 + (-s32)v_3
```

```text
Phi_B = (-1 + s31)v_3 + (-s41)v_4
```

The [shared momentum/Fourier root](../../families/three-loop-insertions/pi1.md) fixes the full tensor and its physical sign. The exact parameter integral is
\[
J=\int_{\mathcal D}^{\mathrm{ordered}}d^4s\;
\frac{N}{2(1+e)}F^{-2+2e}\,d_0F_0^{-1+e},
\]

$$ F=- s_{21} s_{32} + s_{32}. $$

$$ d=1. $$

$$ F_0=- s_{31} s_{41} z + s_{41} z. $$

$$ d_0=z. $$

$$N=2dF.$$

This compact numerator is the complete transverse source, not a metric-only replacement. It satisfies identity(2) of the shared derivation. The original outer-to-inner order is `s21,s32,s31,s41`, with these actual dependent bounds:

| Variable | Lower | Upper |
|---|---|---|
| $s_{21}$ | $0$ | $1$ |
| $s_{32}$ | $0$ | $1$ |
| $s_{31}$ | $0$ | $s_{32}$ |
| $s_{41}$ | $0$ | $1$ |

No old cube Jacobian is inserted. On the starting positive-regulator strip the actual adjacent factors have powers 2Re(e)-1 and Re(e)-1 at their respective endpoints, both integrable for Re(e)>0. The complete ordered-domain kernel is bounded by the corresponding full-cube product, justifying the specific regroupings below. Local meromorphic continuation is taken after the Beta integrations; the original source has poles at zero.

## Native normalization and its concrete reduction

This graph has m=4, a1=-1/2. With Z=exp(-gamma_E e)Gamma(1-e), the actual native rational factor is 1/[32e(1-2e)] and
\[
I=\frac{ma_1}{64}\frac{Z^3}{e(1-2e)}K,
\quad K=\frac{J(z)+z^{3e}J(1/z)}{1+z^{3e}},\quad \frac{ma_1}{64}=- \frac{1}{32}.
\]
The unprojected trace is C_F(C_F-C_A/2)B1. The adopted connected projection and the complete B1 multiplier are distinct stages; neither the orbit nor the sign is reinserted after this equation.

Use $(a,u,b,t)=(s_{21},s_{32},s_{31},s_{41})$, with $0<b<u<1$, $a,t\in[0,1]$. The complete adjacent insertion and bare factors are $(1-a)^{2e-1}u^{2e-1}/(1+e)$ and $z^e(1-b)^{e-1}t^{e-1}$. The first two independent primitives give $1/(2e)$ and $1/e$; the ordered $b$ integration gives $[1-(1-u)^e]/e$. Hence
\[
J=\frac{z^e}{2e^3(1+e)}\int_0^1u^{2e-1}[1-(1-u)^e]du
=\frac{z^e(1-Q)}{4e^4(1+e)},\quad Q=\frac{\Gamma(1+e)\Gamma(1+2e)}{\Gamma(1+3e)}.
\]
The last identity is $1/(2e)-B(2e,1+e)=(1-Q)/(2e)$. Its complete paired root is
\[
K=\frac{(1-Q)\cosh(e\rho/2)}{4e^4(1+e)\cosh(3e\rho/2)}.
\]
The cancellation $1-Q=O(e^2)$ is a source endpoint cancellation, lowering the actual I pole to3.

The actual finite budget is: **Q through5 and Gamma/rational through5 conservatively; the actual pole is3**. The displayed Beta/Gamma formulas directly evaluate every required source coefficient; this adjacent route needs no adopted M or ordered-double-kernel dictionary.

## The highest coefficient and the regular colour operator

Let k_r=[e^r]K for the explicit evaluated K above, and let a_n=[e^n]Z^3/(1-2e). The shared Gamma recurrence fixes a_n exactly. The full physical rule is c_p=(ma1/64)sum(a_n k_(p+1-n)), restricted to existing source orders. In particular, for this graph the finite term is the explicit finite sum
\[
c_0=\frac{ma_1}{64}\sum_{n=0}^{3}a_n k_{1-n}.
\]
This is coefficient extraction from the solved functions just given; it contains no unknown integral. It retains the pole-times-positive-order terms. The c0 expression printed below is the resulting exact moment form, with the complete HLOG form available by the explicit shared dictionary.

The two connected colour multipliers are phi_A=(10-6e)/(3-2e) and phi_f=-8(1-e)/(3-2e), multiplying CA^2 CF and CA CF TF nf. Their constants are10/3 and-8/3; for n>=1 their Taylor coefficients are (2/9)(2/3)^(n-1) and (8/9)(2/3)^(n-1). Therefore the coloured finite term is sum(phi_X,n c_-n,n=0,...,3), not phi_X,0 c0 alone. The other two connected channels vanish in this projection; this does not delete the unprojected source trace.

## Complete coefficients

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=0.$$

$$c_{-4}(z)=0.$$

$$c_{-3}(z)=- \frac{\pi^{2}}{384}.$$

$$c_{-2}(z)=- \frac{\pi^{2}}{384} + \frac{3 \zeta\left(3\right)}{64}.$$

$$c_{-1}(z)=\frac{\pi^{2} \rho^{2}}{384} - \frac{37 \pi^{4}}{23040} - \frac{\pi^{2}}{128} + \frac{3 \zeta\left(3\right)}{64}.$$

$$c_{0}(z)=- \frac{3 \rho^{2} \zeta\left(3\right)}{64} + \frac{\pi^{2} \rho^{2}}{384} - \frac{37 \pi^{4}}{23040} - \frac{5 \pi^{2}}{384} - \frac{5 \pi^{2} \zeta\left(3\right)}{768} + \frac{9 \zeta\left(3\right)}{64} + \frac{21 \zeta\left(5\right)}{64}.$$

The remainder is O(e). Both complete coloured finite coefficients are:

| External monomial | Complete coefficient at e^0 |
|---|---|
| `CA^2*CF` | $- \frac{5 \rho^{2} \zeta\left(3\right)}{32} + \frac{\pi^{2} \rho^{2}}{108} - \frac{37 \pi^{4}}{6480} - \frac{89 \pi^{2}}{1944} - \frac{25 \pi^{2} \zeta\left(3\right)}{1152} + \frac{35 \zeta\left(3\right)}{72} + \frac{35 \zeta\left(5\right)}{32}$ |
| `CA*CF*TF*nf` | $- \frac{\pi^{2} \rho^{2}}{216} + \frac{\rho^{2} \zeta\left(3\right)}{8} - \frac{7 \zeta\left(5\right)}{8} - \frac{11 \zeta\left(3\right)}{36} + \frac{5 \pi^{2} \zeta\left(3\right)}{288} + \frac{49 \pi^{2}}{1944} + \frac{37 \pi^{4}}{12960}$ |

[result.json](result.json) gives every coloured order, including structural zeros; the finite-convolution rule above reconstructs them. [Reproduce](reproduce.md) and [validation](validation.md) state precisely which computations are rerun and which source/function evaluations are adopted. The claim is restricted to the stated positive-real branch and local regulator continuation; formal=not_run.
