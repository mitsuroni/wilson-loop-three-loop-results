# How an adjacent exchange produces the double pole

The first question is where the familiar one-loop double pole enters the
calculation. We can answer it before any expansion: one end of an adjacent
exchange approaches the cusp from each of its two incident edges. Once the
Fourier propagator has been integrated, those two independent distances
become two elementary endpoint integrals. This example fixes the factors
that must not be counted again in higher-loop assembly.

We work in D=4-2 epsilon, with null edge vectors satisfying sum_i v_i=0.
Choose x1=0, x2=v1, x3=v1+v2 and x4=-v4; a point on edge i is
X_i(u)=x_i+u v_i. Put s=2 v1.v2, t=2 v1.v4 and z=s/t>0. The starting
region s,t<0 is a region of analytic invariants; it need not be represented
by real Minkowski vectors. The complete input is linked by
[input.json](input.json) and the [frozen upstream entry](../../../data/2026-11-02/integrals/1loop/2e-01/README.md).

## Start with the momentum and path contribution

For the ordered representative joining edges3 and4, use the source Fourier
phase Phi=(1-s31)v3+s41 v4 and the common exp(-i l.X) convention. Before
the momentum integral the contribution is

$$
\mathcal A_{34}=\frac1{N_c}\operatorname{tr}
\int_0^1ds_{31}\int_0^1ds_{41}
(ig_0v_4^\nu T^b)(ig_0v_3^\mu T^a)
\int\frac{d^Dl}{(2\pi)^D}
\frac{-ig_{\mu\nu}\delta^{ab}}{l^2+i0}
e^{il\cdot\Phi}.
$$

Late Wilson generators are on the left, as fixed by the native action/path
contract. The color contraction is CF. We use the accepted coordinate-space
propagator

$$
D_{\mu\nu}^{ab}(X)=-\delta^{ab}g_{\mu\nu}
\frac{\Gamma(1-\epsilon)}{4\pi^{2-\epsilon}}
(-X^2+i0)^{-1+\epsilon}.
$$

This Fourier/QFT rule and its phase convention are reused inputs; the new
Python calculation does not claim to derive the QCD path integral. It does
check their actual geometry and scalar prefactor against both source
versions. The Wilson i² and coordinate-propagator minus sign cancel.
Since Phi²=2(v3.v4)(1-s31)s41, define b=v1.v4 and obtain

$$
f=\frac{\Phi^2}{2b}=z(1-s_{31})s_{41},\qquad
d=\frac{v_3\cdot v_4}{b}=z.
$$

The full null-quadrilateral Gram matrix, the actual source phase and both
polynomials are replayed symbolically. This tests the path-to-kernel
connection rather than simply recognizing a previously stored coefficient.

## Isolate C and keep the original J

Let lambda_tag=alphaS/pi denote the source coupling label and set A=-2b>0
in the initial real region. The native raw scalar factor is

$$
C_{\rm raw}=-\frac12 C_F\lambda_{\rm tag}
\Gamma(1-\epsilon)(\pi A)^{\epsilon},\qquad
J_A=z\int_0^1ds_{31}\int_0^1ds_{41}
[z(1-s_{31})s_{41}]^{\epsilon-1}.
$$

The minus one half comes from extracting -2b from the distance polynomial;
it is not a sign selected by RG. The native C.combined equals the older
original-ordered C.combined exactly. The prepared0.3 C adds one, and only
one, exp(-EulerGamma epsilon). Its one-loop color coordinate is a1=1, so
no physical color information is lost in this particular projection.
The cube export is only x1=s31, x2=s41 with Jacobian1; the actual replay
checks that substitution. The original outer-to-inner order stays
s31 then s41.

Now perform the two integrals in that order. For Re epsilon>0 the positive
bases may be split, both endpoint integrals converge absolutely, and

$$
J_A=z^\epsilon
\int_0^1(1-s_{31})^{\epsilon-1}ds_{31}
\int_0^1s_{41}^{\epsilon-1}ds_{41}
=\frac{z^\epsilon}{\epsilon^2}.
$$

The primitive t^epsilon/epsilon is differentiated back in the new script;
the noninteger endpoint integrals are also evaluated numerically. We do
not insert epsilon=0 into the original divergent integrand.

## Why the normalized I has no z dependence

The C4 multiplicity is m=4. It counts the four adjacent edge pairs, whereas
the source metadata n_diagrams=1 counts the canonical source record. Let

$$G(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon).$$

The project's single physical-I convention follows by matching the
orbit/crossing sum to its two extracted scales:

$$
\frac m2 C_{\rm prep}
\left[J(z,\epsilon)+z^\epsilon J(1/z,\epsilon)\right]
=C_F\lambda_{\rm tag}(\pi A)^\epsilon
(1+z^\epsilon)I(z,\epsilon).
$$

Thus the normalized object is

$$
I=-\frac m4G(\epsilon)
\frac{J(z,\epsilon)+z^\epsilon J(1/z,\epsilon)}{1+z^\epsilon}.
$$

For the adjacent graph the numerator is (1+z^epsilon)/epsilon². Consequently

$$\boxed{I_A=-\frac{G(\epsilon)}{\epsilon^2}.}$$

Orbit multiplicity, crossing and Euler compensation are already included.
CF, the coupling and the common one-loop scale are removed. There is no
additional m to apply when this I is added to the opposite exchange.
Using the bare-to-MS-bar dictionary on the raw amplitude gives the common
physical scale mu_x²=e^(2 gamma_E)mu_MS²/4. The combined coefficient is
CF lambda S1(I_A+I_O), with
S1=exp(epsilon L_s)+exp(epsilon L_t). Restoring that physical scale does
not insert another Euler factor into I.

## Expand only the exact meromorphic answer

The Euler term in log Gamma cancels the explicit compensation, leaving

$$
\log G(\epsilon)=\sum_{n=2}^{\infty}
\frac{\zeta_n}{n}\epsilon^n,\qquad |\epsilon|<1.
$$

The new reproducer exponentiates this series and independently compares it
with the direct SymPy series of exp(-gamma_E epsilon)Gamma(1-epsilon).
Since I has a double pole, coefficients through epsilon3 require the
regular G series through epsilon5. The actual result is

$$
\boxed{\begin{aligned}
I_A={}&-\epsilon^{-2}-\frac{\pi^2}{12}
-\frac{\zeta_3}{3}\epsilon-\frac{\pi^4}{160}\epsilon^2\\
&-\left(\frac{\pi^2\zeta_3}{36}+\frac{\zeta_5}{5}\right)
\epsilon^3+O(\epsilon^4).
\end{aligned}}
$$

The simple-pole coefficient is exactly zero. All higher coefficients are
fixed by the same all-epsilon Gamma expression; this release explicitly
checks and exports only through epsilon3.

## What was executed and what was reused

The new run evaluates the source/Gram and C bridges, the actual endpoint
primitive, the Gamma coefficient recurrence and noninteger original
endpoint integrals. An independent epsilon-Cauchy readout of the exact
regular Gamma function checks every exported coefficient. Wrong orbit
doubling and repeated Euler compensation are detected by concrete negative
controls. The old r002 table is read only after the new predictions have
been saved, and is compared only at its existing orders.

The [check record](checks.json), [reproduction notes](reproduce.md) and
[archived original derivation](../sources/benchmarks/2ne/workspaces/2e-01/derivation.md)
separate those actions from the reused source contract and Fourier rule.
The real-domain computation and local epsilon continuation are research
results; no new general physical complex-kinematics continuation or integral
formal proof is supplied. `formal=not_run`.
