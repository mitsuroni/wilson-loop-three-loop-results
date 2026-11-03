# Two-loop parameter integrals: conventions and evaluated sources

This packet follows the twelve orbit contributions of the connected Wilson-loop
exponent. It keeps the six exchange representatives, four tree three-gluon
vertices and two one-loop self-energy insertions distinct. Equality of two
evaluated integrals is not a reason to delete either orbit.

We use fixed real z>0, ell=log(z), D=4-2epsilon and ordinary coefficients
$I=\sum_q a_q\epsilon^q$. All twelve tables extend from epsilon^-4 through
epsilon², with remainder O(epsilon³) at fixed z. More negative displayed zeros
in historical records are consistent pole-bound zeros, not additional poles.

## What the source and prefactor mean

The public [upstream catalogue](../../data/2026-11-02/catalog.json) stores each
graph, momentum-space contribution, raw C/J and available prepared alternative.
A raw C/J record and its prepared counterpart are different stages. The exact
ordered domain and any alpha-simplex measure are part of J. A field named
`integrand` is not an evaluated integral. Current raw native data contains no
downstream Euler compensation or a1 projection; cyclic factors must be read
from the relevant family contract, not added from metadata a second time.

For the exchange family let m=(4,4,4,2,4,1) and

$$
Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad
\mathcal B_2[J]=\frac{J(z)+z^{2\epsilon}J(1/z)}{1+z^{2\epsilon}},\quad
I_{4e-r}=-\frac{m_r}{16}Z^2\mathcal B_2[J_r].
$$

The factor -m/16 is m*a1/[2*(-2)²] with a1=-1/2. It is already in I. The
prepared common C is

$$
C_{4e}=-\frac{C_AC_F}{2}e^{-2\gamma_E\epsilon}(4\pi\alpha_s)^2
\frac{\Gamma(1-\epsilon)^2}{(4\pi^{2-\epsilon})^2}
\frac{(-2v_{14}+i0)^{2\epsilon}}4.
$$

For TGV orbits all four cyclic multiplicities equal four, but the vertex
prefactor is different:

$$
I_r=\frac{e^{-2\gamma_E\epsilon}\Gamma(2-2\epsilon)}2\mathcal B_2[J_r],
\quad J_r=\int_{\mathcal D_r}d^3s\int[d\alpha]_3
\frac{(\alpha_1\alpha_2\alpha_3)^{-\epsilon}P_r}{S_r^{2-2\epsilon}}.
$$

Here $[d\alpha]_3=d^3\alpha\,\delta(1-\sum\alpha)$, and the prepared
representative C is the raw vertex prefactor
$C_AC_Fg_0^4\Gamma(2-2\epsilon)(-2v_{14}+i0)^{2\epsilon}/(64\pi^{4-2\epsilon})$
times the single Euler compensation. The original loop measure and phase
definitions remain in the raw source. Replacing this Gamma by Z² is wrong.

The Pi1 pages use native physical r002 scalar I. Their external color weight
is $\mathcal B_1=a(\epsilon)C_A+b(\epsilon)T_Fn_f$, where
$a=(10-6\epsilon)/(3-2\epsilon)$ and $b=-8(1-\epsilon)/(3-2\epsilon)$.
The old r001 convention needed an additional -1; r002 has absorbed it. The
released tables therefore use assembly_sign=+1 and are consumed with +B1.

## One real moment definition, fixed across every graph

$$
r=\frac z{1+z},\quad W=\log[\beta(1-\beta)],\quad
R=\log\frac{1-\beta}{\beta},\quad L_\beta=\ell+R,
$$
$$
M_{k,j}=\int_0^1\frac{W^kL_\beta^j}{\beta-r}d\beta,\qquad
N_{k,j}=\int_0^1\frac{W^kRL_\beta^j}{\beta-r}d\beta
=M_{k,j+1}-\ell M_{k,j}.
$$

These are raw powers, without factorials. At beta=r, Lbeta=0 and the kernel
has a removable limit; no principal value is added. The HLOG convention is
$G(a,\mathbf w;z)=\int_0^zG(\mathbf w;t)dt/(t-a)$, so G(0;z)=log(z).
All stored HLOG words in this packet have letters 0,-1 and rational zeta-period
coefficients. M/N are convenient representatives, not a claim of a minimal
algebraically independent basis.

## The evaluated mother functions used by the actual engine

$$
V(\epsilon)=\int_0^1\frac{e^{\epsilon W/2}}{\beta-r}
\frac{\sinh(\epsilon L_\beta/2)}{\epsilon/2}d\beta
=\sum_{q\ge0}\epsilon^q\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,2j+1}}{2^q(q-2j)!(2j+1)!},
$$
$$
V_N(\epsilon)=\int_0^1\frac{e^{\epsilon W}R}{\beta-r}
\frac{\cosh(\epsilon L_\beta)-1}{\epsilon^2}d\beta
=\sum_{q\ge0}\epsilon^q\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,2j+2}}{(q-2j)!(2j+2)!}.
$$

V(2epsilon) doubles the source parameter, not the external Z². Also set

$$
c_1=\cosh(\epsilon\ell/2),\quad c_2=\cosh(\epsilon\ell),\quad
R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\quad
\chi=\gamma_E+\psi(1+2\epsilon),\quad
G_2=e^{-2\gamma_E\epsilon}\Gamma(1-2\epsilon).
$$

The code constructs the Gamma logarithms and their finite exponential
recurrences, differentiates the fixed V/VN sources by their coefficient
rules, then performs the complete prefactor convolution. For example Z²
must be expanded to epsilon⁶ for an epsilon^-4 source to reach epsilon².
Only after rebuilding the coefficients are historical tables read for
comparison. Equivalent compact representatives are published after their
full rational HLOG coordinate difference is zero.

## A genuinely reused extra function at TGV04 epsilon²

The symbol tau_D4 is the previously defined double-endpoint moment, not a name
assigned to an unknown remainder:

$$
D=(1-t)(1-s)+zts,\quad
X_+=\log\frac{Dts}{\sqrt z},\quad
X_-=\log\frac{D(1-t)(1-s)}{\sqrt z},\quad
\tau^+_{D,4}=-\frac{1+z}{48}\int_0^1\!\int_0^1
\frac{X_+^4+X_-^4}{D}dt\,ds.
$$

Its complete stored HLOG dictionary is included. The coefficient is 1/8 in
TGV04 epsilon². No new E function or non-0/-1 letter is needed here. The
restricted M/N non-membership statement in the original study is not an
impossibility theorem for every conceivable function representation.

## Reproduction and its precise scope

From the public repository root, with Python 3.12, SymPy, mpmath and NumPy:

```text
python integrals/two-loop/reproduce.py --output-dir replay-output/two-loop-new
python integrals/two-loop/check_integrals.py --output-dir replay-output/two-loop-controls-new
```

Both commands reject a nonempty output directory and reject immutable public
science/docs/package trees as output locations. Inside the public repository,
only replay-output/<new-name> is allowed. The first hashes the bundled
sources and rebuilds all84 coefficients and complete HLOG comparisons. The
second performs new original ordered polynomial integrations, Pi1 tensor/corner
identities, reduced TGV01 Beta integrals, direct M/N/tau quadrature and the
4e04 two-variable coefficient integrals. Numeric mesh errors are empirical,
not rigorous intervals.

The historical HyperInt integrations and fixed HLOG moment dictionary are
explicitly reused. Archival CAS scripts under sources retain their original
bytes, installation paths and old output writers. They are not the portable
execution entry, and were not rerun in this publication task. No Maple,
Mathematica, QGRAF or TeX executable is needed for the two Python commands.
The packet proves neither a new global complex-kinematics theorem nor a Lean
integral chain; formal=not_run. Source binding and a successful replay are
evidence components, not substitutes for the integral derivations below.
