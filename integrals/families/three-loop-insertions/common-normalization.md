# Physical normalization and finite colour convolution

Use the project's [single physical convention](../../../workflow/current-conventions.md), with
\(D=4-2e\), \(e=\epsilon\), \(z>0\), \(\rho=\log z\), and the inherited
Feynman boundary value. The accepted initial objects are the explicitly bound
native \(C/J\) sources. Graph generation and the internal effective self-energy
operators are upstream inputs. Each derivation exposes the subsequent tensor,
parameter and coefficient calculation.

Let \(m\) be the actual cyclic orbit size and
\(S_{14}=(-2\pi v_{14}+i0)^{3e}\), \(S_{12}=z^{3e}S_{14}\).
After the external colour/coupling operator \(\mathcal N\) has been extracted,
write the prepared coefficient as \(C=\mathcal N S_{14}R(e)\). The orbit relation
and the scalar definition are
\[
\frac m2 C[J(z,e)+z^{3e}J(1/z,e)]
=\mathcal N(S_{14}+S_{12})I(z,e),\qquad
I=\frac m2R(e)\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.\tag{N1}
\]
For the single-Pi1 four-endpoint branch its explicitly stated connected
fixed-basis coordinate \(a_1\) multiplies this expression; that projection does
not erase the unprojected source trace. The chain and Pi2 geometries have a
Wilson trace \(C_F\), with the complete insertion operator kept external.
In particular the Pi2 metadata's geometric `a1=1` is not a projection that
discards its real \(C_F^2T_Fn_f\) channel.

Each family gives its actual \(R\) from its own source. It is not fixed by a
family-wide guessed sign. Put
\[
\mathcal G(e)=e^{-\gamma_Ee}\Gamma(1-e),\qquad
\log\mathcal G=\sum_{k\ge2}\frac{\zeta_k}{k}e^k.
\]
The prepared input already supplies one \(e^{-3\gamma_Ee}\) factor. Writing
\(\mathcal G^3\) consumes it exactly once, irrespective of the number of
explicit exchanged momenta remaining after internal loops were integrated.

## A pole fixes the positive-order budget

For \(I=\sum_{p=-P}^{0}c_pe^p+O(e)\) and an external regular scalar
\(\phi=\sum_{j\ge0}\phi_je^j\), the actual coloured coefficient is
\[
[e^k](\phi I)=\sum_{j=0}^{k+P}\phi_jc_{k-j},\qquad
[e^0](\phi I)=\sum_{j=0}^{P}\phi_jc_{-j}.\tag{N2}
\]
Thus replacing an insertion operator by its four-dimensional value loses
finite terms. If intermediate pieces of \(\phi\) have poles, their cancellation
must be carried out at a sufficient absolute order before this rule is used.
The Pi2 example in [chains-pi2](chains-pi2.md) requires a fifth-order Gamma
ratio even though the final adjacent multiplier is needed only through four.

For a regular factor \(F(e)=\exp\sum_{k\ge1}\ell_ke^k\), its exact finite jet
obeys
\[
F_0=1,\qquad F_n=\frac1n\sum_{k=1}^nk\ell_kF_{n-k}.\tag{N3}
\]
Rational factors use polynomial division; products use the finite Cauchy
sum. These operations generate coefficients from the evaluated mathematics,
with historical whole-I tables reserved for subsequent comparisons.

## Evaluated functions and branch

All finite HLOG expressions here use
\[
G(a,\mathbf w;z)=\int_0^z\frac{G(\mathbf w;t)}{t-a}\,dt,
\qquad G(\varnothing;z)=1,\quad G(0;z)=\log z,
\]
with letters \(0,-1\), real logarithms and exact zeta periods. Repeated zero
words mean \(G(0^n;z)=\rho^n/n!\). The original unnormalised moments are
\[
r=\frac z{1+z},\quad W=\log[\beta(1-\beta)],\quad
L=\rho+\log\frac{1-\beta}{\beta},\qquad
M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta\quad(j\ge1),
\]
\[
N_{k,j}=M_{k,j+1}-\rho M_{k,j}.
\]
The zero of \(L\) at \(\beta=r\) removes the apparent internal singularity.
There are no factorials hidden in \(M\). The four moments needed by the chain
and Pi2 geometries are fully evaluated in [their shared derivation](chains-pi2.md).
The additional functions used by the six single-Pi1 routes are defined and
evaluated in [pi1](pi1.md).

The statement is at fixed positive real \(z\), initially in the regulator
strip justified for each kernel, followed by local meromorphic continuation.
It does not establish an additional complex-kinematics or \(z=-1\) branch
theorem. A precise accepted source/operator or evaluated-function lemma keeps
its original scope; exposition, actual replay and formal proof are distinct.
