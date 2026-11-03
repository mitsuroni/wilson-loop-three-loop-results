# Three-loop exchange normalization and finite coefficient assembly

This page applies the project's [single physical interface](../../../workflow/current-conventions.md) to the39 individual6e source entries. Its formulas derive the current accepted C/J normalization and do not introduce another convention. The accepted root of each calculation is its original ordered
prepared C/J. Graph generation, Fourier propagation and color projection are
upstream accepted facts. The integral derivations preserve the concrete J,
its dependent bounds and outer-to-inner order. Work at fixed real z>0,
ell=log z and D=4-2epsilon, initially Re epsilon>0. The inherited source
prescription is(-2v14+i0); the positive-real-z branch uses real logs.

## Extract the physical I once

At three loops let a1 be the first coordinate in the full fixed color basis,
and m the C4 orbit multiplicity. With the common color/coupling/scale removed,
the prepared prefactor is

$$C=\frac{a_1}{(-2)^3}\mathcal N S_{14}Z^3,\qquad
Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad
\mathcal N=C_A^2C_F(\alpha_s/\pi)^3,\quad
S_{14}=(-2\pi v_{14}+i0)^{3\epsilon}.$$

The symmetric external scale is S14(1+z^(3epsilon)). Matching the actual
orbit sum gives

$$\frac m2 C[J(z)+z^{3\epsilon}J(1/z)]
=\mathcal N S_{14}(1+z^{3\epsilon})I,$$
$$\boxed{I=\kappa Z^3\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}},\quad
\kappa=\frac{ma_1}{2(-2)^3}=-\frac{ma_1}{16}.}\tag{N1}$$

The Euler compensation, orbit size, color coordinate and crossing are all
already present in this I. A graph count is not another orbit factor; no
additional color sign or Euler term is subsequently inserted. Each individual
article states its own m,a1,kappa rather than inferring them from a family name.

For a centered source J=z^(3epsilon/2)B(epsilon,z), equation(N1) becomes

$$I=\frac{\kappa Z^3}{\cosh(3\epsilon\ell/2)}
\frac{B(\epsilon,z)+B(\epsilon,1/z)}2.\tag{N2}$$

If a further phase remains in B, it stays inside this paired sum. It cannot
be removed by declaring the original unsymmetrized J crossing even.

## The finite recurrence used by the engine

Write Z³=sum g_n epsilon^n. Its logarithm is

$$\log Z^3=\sum_{k\ge2}\frac{3\zeta_k}{k}\epsilon^k,\qquad
g_0=1,\quad g_n=\frac3n\sum_{k=2}^n\zeta_k g_{n-k}.\tag{N3}$$

For s_n=[epsilon^n]sech(3epsilon ell/2),

$$s_0=1,\quad s_n=-\sum_{j=1}^{\lfloor n/2\rfloor}
\frac{(3\ell/2)^{2j}}{(2j)!}s_{n-2j},\qquad
p_n=\sum_{i=0}^ng_i s_{n-i}.\tag{N4}$$

Thus p0=1,p1=0 and p2=3zeta2/2-9ell²/8. For any known regular centered source with an external epsilon^(-P), let
b_q be the coefficient of its complete crossing-paired average
[B(epsilon,z)+B(epsilon,1/z)]/2, including any remaining phase weights.
Then all exported coefficients
are the finite sums

$$c_k=\kappa\sum_{n=0}^{k+P}p_n b_{k+P-n}\quad(k\ge-P),\qquad
c_k=0\ (k<-P).\tag{N5}$$

The target for each article is every c_-6,...,c0 and remainder O(epsilon),
at fixed z. The adjacent-only route takes P=6; the one-opposite route P=4.
Source zeros may lower the actual pole; these zeros are calculated before
coefficient extraction. For c0, the full regular product must be retained to
P, including all source, Gamma, Beta and crossing terms. A finite source
coefficient alone is not the finite physical I.

Use these identities after convolution when compact M/N coefficients are
wanted:

$$M_{k,j}=\int_0^1\frac{W^kL^j}{\beta-r}d\beta,\quad
N_{k,j}=M_{k,j+1}-\ell M_{k,j},\quad
r=\frac z{1+z},\ W=\log[\beta(1-\beta)],\
L=\ell+\log\frac{1-\beta}{\beta},$$
$$M_{0,1}=-\frac{\ell^2+\pi^2}{2},\qquad M_{0,3}=-M_{0,1}^2.$$

The j>=1 numerators vanish at beta=r, so these are ordinary removable
integrals, with raw powers and no factorial normalization. M/N and the tau
coefficients defined in the single-opposite article are adopted established
function representations with exact HLOG dictionaries, not uncomputed
remainders. HLOG means

$$G(a,\mathbf w;z)=\int_0^z\frac{G(\mathbf w;t)}{t-a}dt,\quad
G(\varnothing;z)=1,\quad G(0;z)=\log z.$$

The present solved dictionaries use letters0,-1 and exact rational zeta
periods. Neither basis minimality nor a general complex-z or z=-1 continuation
is established here. The derivations and executed CAS controls are research
mathematics; no new formal integral proof is claimed.


The optional compact symbol E(z) is the established evaluated endpoint
insertion, defined independently on
Omega2={0<t1<t2<1,0<s1<s2<1} by

$$E(z)=\frac{\ell(1+z)^2}{2}\int_{\Omega_2}
\frac{\log[s_2/(1-s_1)]}{D_z(t_1,s_1)D_z(t_2,s_2)}d^4s.$$

Its actual connection to the finite three-kernel result and its adopted
HLOG evaluation are stated in [same-channel-three](same-channel-three.md).
The exact evaluated mapping is in the [function dictionary](runtime/references/basis-and-conversions.json),
with full-HLOG coefficient alternatives in each relevant entry. Thus E in a
compact table does not stand for an unknown remainder. Likewise tau_D4 is
an ordinary epsilon⁴ coefficient of HD's crossing-even part,

$$\tau_{D,4}(z)=\frac{h_{D,4}(z)+h_{D,4}(1/z)}2,$$

with no extra4!; both complete evaluated functions are supplied by the
[single-kernel jet appendix](appendices/single-opposite-jets.md).
