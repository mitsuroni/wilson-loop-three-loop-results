# 4e-Pi1-02: its transverse root to the complete scalar and coloured series

Work at D=4-2e, e=epsilon, fixed positive real z and rho=log z, with the inherited source i0. The accepted starting object is this graph's original prepared C/J, including its connected projection and exactly one Euler compensation. The target is every scalar and connected-colour coefficient at powers -6 through0, with O(e).

## Its actual ordered integral

The insertion is the line with endpoints `s32,s42`; the bare exchange has endpoints `s31,s41`. In the adopted momentum expression these are the Fourier factors D0 Pi1 D0 and D0, respectively. Their actual source separations in the null tangent basis are

```text
Phi_A = (-1 + s32)v_3 + (-s42)v_4
```

```text
Phi_B = (-1 + s31)v_3 + (-s41)v_4
```

The [shared momentum/Fourier root](../../families/three-loop-insertions/pi1.md) fixes the full tensor and its physical sign. The exact parameter integral is
\[
J=\int_{\mathcal D}^{\mathrm{ordered}}d^4s\;
\frac{N}{2(1+e)}F^{-2+2e}\,d_0F_0^{-1+e},
\]

$$ F=- s_{32} s_{42} z + s_{42} z. $$

$$ d=z. $$

$$ F_0=- s_{31} s_{41} z + s_{41} z. $$

$$ d_0=z. $$

$$N=2dF.$$

This compact numerator is the complete transverse source, not a metric-only replacement. It satisfies identity(2) of the shared derivation. The original outer-to-inner order is `s32,s31,s42,s41`, with these actual dependent bounds:

| Variable | Lower | Upper |
|---|---|---|
| $s_{32}$ | $0$ | $1$ |
| $s_{31}$ | $0$ | $s_{32}$ |
| $s_{42}$ | $0$ | $1$ |
| $s_{41}$ | $0$ | $s_{42}$ |

No old cube Jacobian is inserted. On the starting positive-regulator strip the actual adjacent factors have powers 2Re(e)-1 and Re(e)-1 at their respective endpoints, both integrable for Re(e)>0. The complete ordered-domain kernel is bounded by the corresponding full-cube product, justifying the specific regroupings below. Local meromorphic continuation is taken after the Beta integrations; the original source has poles at zero.

## Native normalization and its concrete reduction

This graph has m=4, a1=-1/2. With Z=exp(-gamma_E e)Gamma(1-e), the actual native rational factor is 1/[32e(1-2e)] and
\[
I=\frac{ma_1}{64}\frac{Z^3}{e(1-2e)}K,
\quad K=\frac{J(z)+z^{3e}J(1/z)}{1+z^{3e}},\quad \frac{ma_1}{64}=- \frac{1}{32}.
\]
The unprojected trace is C_F(C_F-C_A/2)B1. The adopted connected projection and the complete B1 multiplier are distinct stages; neither the orbit nor the sign is reinserted after this equation.

The same domain map $x=1-s_{31}$, $y=1-s_{32}$, $a=s_{41}$, $b=s_{42}$ has $y<x$, $a<b$ and Jacobian1. Here the *upper* labelled line carries the insertion: the full adjacent transverse density is still $dF^{-1+2e}/(1+e)$ but its powers are on $y,b$. Thus
\[
J=\frac{z^{3e}}{1+e}\left[\int_{0<y<x<1}x^{e-1}y^{2e-1}dxdy\right]
\left[\int_{0<a<b<1}a^{e-1}b^{2e-1}dadb\right]
=\frac{z^{3e}}{18e^4(1+e)}.
\]
The two factors are now $1/(6e^2)$ and $1/(3e^2)$. This explicitly proves the present source evaluation, rather than asserting equality from two matching final tables. Its physical crossing gives
\[
K=\frac1{18e^4(1+e)},\qquad I=-\frac{Z^3}{576e^5(1-2e)(1+e)}.
\]

The actual finite budget is: **Gamma/rational through5; no M or ordered-body input is mathematically needed**. The displayed Beta/Gamma formulas directly evaluate every required source coefficient; this adjacent route needs no adopted M or ordered-double-kernel dictionary.

## The highest coefficient and the regular colour operator

Let k_r=[e^r]K for the explicit evaluated K above, and let a_n=[e^n]Z^3/(1-2e). The shared Gamma recurrence fixes a_n exactly. The full physical rule is c_p=(ma1/64)sum(a_n k_(p+1-n)), restricted to existing source orders. In particular, for this graph the finite term is the explicit finite sum
\[
c_0=\frac{ma_1}{64}\sum_{n=0}^{5}a_n k_{1-n}.
\]
This is coefficient extraction from the solved functions just given; it contains no unknown integral. It retains the pole-times-positive-order terms. The c0 expression printed below is the resulting exact moment form, with the complete HLOG form available by the explicit shared dictionary.

The two connected colour multipliers are phi_A=(10-6e)/(3-2e) and phi_f=-8(1-e)/(3-2e), multiplying CA^2 CF and CA CF TF nf. Their constants are10/3 and-8/3; for n>=1 their Taylor coefficients are (2/9)(2/3)^(n-1) and (8/9)(2/3)^(n-1). Therefore the coloured finite term is sum(phi_X,n c_-n,n=0,...,5), not phi_X,0 c0 alone. The other two connected channels vanish in this projection; this does not delete the unprojected source trace.

## Complete coefficients

$$c_{-6}(z)=0.$$

$$c_{-5}(z)=- \frac{1}{576}.$$

$$c_{-4}(z)=- \frac{1}{576}.$$

$$c_{-3}(z)=- \frac{1}{192} - \frac{\pi^{2}}{2304}.$$

$$c_{-2}(z)=- \frac{5}{576} - \frac{\pi^{2}}{2304} - \frac{\zeta\left(3\right)}{576}.$$

$$c_{-1}(z)=- \frac{11}{576} - \frac{\pi^{2}}{768} - \frac{19 \pi^{4}}{276480} - \frac{\zeta\left(3\right)}{576}.$$

$$c_{0}(z)=- \frac{7}{192} - \frac{5 \pi^{2}}{2304} - \frac{19 \pi^{4}}{276480} - \frac{\zeta\left(3\right)}{192} - \frac{\pi^{2} \zeta\left(3\right)}{2304} - \frac{\zeta\left(5\right)}{960}.$$

The remainder is O(e). Both complete coloured finite coefficients are:

| External monomial | Complete coefficient at e^0 |
|---|---|
| `CA^2*CF` | $- \frac{3353}{26244} - \frac{89 \pi^{2}}{11664} - \frac{19 \pi^{4}}{77760} - \frac{35 \zeta\left(3\right)}{1944} - \frac{5 \pi^{2} \zeta\left(3\right)}{3456} - \frac{\zeta\left(5\right)}{288}$ |
| `CA*CF*TF*nf` | $\frac{\zeta\left(5\right)}{360} + \frac{19 \pi^{4}}{155520} + \frac{11 \zeta\left(3\right)}{972} + \frac{\pi^{2} \zeta\left(3\right)}{864} + \frac{49 \pi^{2}}{11664} + \frac{1897}{26244}$ |

[result.json](result.json) gives every coloured order, including structural zeros; the finite-convolution rule above reconstructs them. [Reproduce](reproduce.md) and [validation](validation.md) state precisely which computations are rerun and which source/function evaluations are adopted. The claim is restricted to the stated positive-real branch and local regulator continuation; formal=not_run.
