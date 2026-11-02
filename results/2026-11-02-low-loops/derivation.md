# From an ordered graph integral to a finite Wilson loop

A useful low-loop example must show where the physics enters, where an integral is
actually evaluated, and where the evaluated answer is tested. Following only the
final polynomial hides a sign or an orbit factor that can change several diagrams
while leaving part of a total unchanged. We therefore keep all fourteen orbits until
their normalized contributions have been connected to the exponent.

## Fix the original object before changing its variables

The [generator](../../packages/qgraf-pipeline/README.md) records ordered Wilson-line
attachments, propagators, vertices, color vectors and orbit metadata. The
[Montopar stage](../../packages/montopar/README.md) applies the family Feynman rules
to the momentum-space contribution and returns C times J. Its frozen per-orbit
momentum expressions and raw/prepared parameter kernels can be opened from the
[input table](README.md). This publication does not assert a new theorem that the
enumerated graph set is complete; its precise source-generation tests remain in
the upstream package records.

The dependent parameter bounds, iterated order, simplex delta measures, i0,
numerators and C factor are part of the integral definition. Raw and prepared
objects differ: preparation includes the recorded fixed-basis color projection
and exactly one Euler compensation. The original full color vector is retained.
A change to cube variables is legitimate only with the explicit family map and
Jacobian; a folder called `prepared` does not supply that argument.

We set D=4-2epsilon, z>0, rho=log(z), lambda=alpha_s_MSbar/pi and T=TF*nf. Shape
logarithms in historical kernels are often named ell or L; in this article L is
always the average *physical scale* logarithm. Thus Ls=L+rho/2 and Lt=L-rho/2.
The published I already contains its orbit/crossing combination, a1/web projection
where applicable, and Euler factor. CF, the coupling and the common external scale
sum are removed. These conventions are established per family, not by applying a
universal factor to every raw C.

## Evaluate the scalar contributions while keeping their common conventions

For the adjacent one-loop exchange the original ordered integral is
J=z^epsilon/epsilon². Its normalized answer is

$$I_{2e-01}=-\frac{e^{-\gamma_E\epsilon}\Gamma(1-\epsilon)}{\epsilon^2}.$$

The opposite graph has D(u,v)=(1-u)(1-v)+zuv. Integrating one parameter gives a
removable divided difference. A beta-variable change reduces the remaining
integral to the M generating function. The [opposite-graph derivation](../../integrals/one-loop/2e-02/derivation.md)
records the all-epsilon form and the expansion; its new checked table was frozen
before the old coefficients were read. Expanding the Gamma logarithm and this
moment source produces the explicit one-loop series through epsilon3.

For two loops, the six exchange diagrams, four three-gluon vertices and two
self-energy insertions use different prefactors. Their [shared derivation](../../integrals/two-loop/METHODS.md)
gives the exact Gamma/Beta and moment mother functions. The actual code expands
those functions and convolves their Gamma factors, then compares every coefficient
in complete rational HLOG coordinates. TGV04's historical HyperInt integration is
reused explicitly, including the defined tau_D4 at epsilon2 and its HLOG dictionary.
The per-orbit articles state their actual source integrations and new controls.

All moments in this packet use *raw powers*, with no hidden factorials:

$$
r=\frac z{1+z},\quad W=\log[\beta(1-\beta)],\quad
R=\log\frac{1-\beta}{\beta},\quad L_\beta=\rho+R,
$$
$$
M_{k,j}=\int_0^1\frac{W^kL_\beta^j}{\beta-r}d\beta,\qquad
N_{k,j}=\int_0^1\frac{W^kRL_\beta^j}{\beta-r}d\beta
=M_{k,j+1}-\rho M_{k,j}.
$$

Only the convergent moments required by the articles are used. For j>=1 the
interior apparent pole is removable because Lbeta(r)=0; no principal value is
silently added. HLOG means G(a,w;z)=integral_0^z G(w;t)dt/(t-a). The stored words
use 0,-1. Useful exact reductions are M01=-(rho²+pi²)/2, M03=-M01² and
N02=M01²/3+4*zeta(2)*M01. M/N here are convenient functions, not a claim that
the displayed monomials form a minimal independent basis.

## Reconstruct the connected color exponent

At two loops a raw exchange trace is CF(CF-CA/2). The CF² disconnected part
must be related to exponentiation before assigning the projected I to w2. For
fixed edge coordinates, the product domain partitions into the 24 endpoint
orderings. The three pairings occur eight times each: 8=2!*(2!)² accounts for
permuting the propagators and reversing their two endpoints. The standard Abelian
shuffle identity is applied to the full contour in a common absolute-convergence
domain, then meromorphically continued. That identity and the original identical
propagator factors give the usual W1²/2 contribution, so subtracting it leaves
the crossed web color -CF*CA/2. This is the ordered-domain shuffle for this
two-loop object, not a rule for deleting CF² terms from an arbitrary graph.
The enumerated pairing count in the driver is a combinatorial control; the
Abelian exponentiation/shuffle identity is explicitly reused mathematical input.
The already connected, projected I tables must not undergo this subtraction again.

The scalar one-loop sum is U1=I_2e-01+I_2e-02. For two loops it is

$$
U_2=C_A\left(\sum_{r=1}^6I_{4e-r}+\sum_{r=1}^4I_{TGV-r}\right)
+\frac{(10-6\epsilon)C_A-8(1-\epsilon)T}{3-2\epsilon}
\left(I_{Pi1-01}+I_{Pi1-02}\right).
$$

In both equations the external CF is removed. The Pi1 entries use native r002,
whose scalar sign already contains the physical operator sign. Legacy r001 used
the opposite scalar I and an extra minus in assembly. The separate
[version bridge](validation/version-bridge-checks.json) checks that J, its domain
and tensor operator are unchanged, C's rational sign changes, and each scalar
Laurent/HLOG coefficient flips. One must consume native I with the plus weight
above. Agreement of only the final sum would be insufficient to infer this bridge.

Restore the scale sum once, at its own loop order:

$$
S_n=e^{n\epsilon(L+\rho/2)}+e^{n\epsilon(L-\rho/2)},\qquad
\frac{w_n^{bare}}{C_F}=S_n U_n. 
$$

In particular S_n[0]=2. Neither the orbit multiplicity nor Euler compensation
is applied again. The full machine ledger retains all orbits even when two
normalized scalar functions happen to be equal.

## Introduce coupling renormalization before taking a finite part

Write b0=(11CA-4T)/12 in the lambda convention. To second order the connected
post-coupling series divided by CF is

$$
P_1=S_1U_1,\qquad P_2=S_2U_2-\frac{b_0}{\epsilon}S_1U_1.
$$

The one-loop scale remains S1 inside the counterterm. Replacing it by S2 changes
the finite answer. Define logZ_n/CF to be the negative-epsilon part of P_n and F_n
its epsilon0 coefficient, so log(W_R)/CF=lambda*F1+lambda²*F2+O(lambda³).
Positive coefficients of P2 are supplied for downstream work; they are not extra
terms in the physical epsilon0 finite function.

The computation gives

$$
F_1=-L^2+\frac{\rho^2}{4}+\frac{\pi^2}{3},\qquad
\frac{\log Z_1}{C_F}=-\frac2{\epsilon^2}-\frac{2L}{\epsilon},
$$
$$
\frac{\log Z_2}{C_F}=\frac{3b_0}{2\epsilon^3}
+\frac{b_0L-k_2/2}{\epsilon^2}
+\frac{-k_2L-h_2/2}{\epsilon},
$$

where k2 and h2 below also have CF removed. They are extracted from the assembled
simple pole and then compared to the independent external coefficients:

$$
k_2=C_A\left(\frac{67}{36}-\frac{\pi^2}{12}\right)-\frac59T,
$$
$$
h_2=C_A\left(\frac{101}{27}-\frac72\zeta_3-\frac{11\pi^2}{72}\right)
+T\left(-\frac{28}{27}+\frac{\pi^2}{18}\right).
$$

The full finite two-loop result is F2=CA*fA+T*ft, with

$$
\begin{aligned}
f_A={}&-\frac{11}{36}L^3+\left(\frac{\pi^2}{12}-\frac{67}{36}\right)L^2
+\left(\frac{11\rho^2}{48}-\frac{101}{27}+\frac72\zeta_3+\frac{11\pi^2}{24}\right)L\\
&-\frac{11}{24}M_{1,1}+\left(\frac{67}{144}-\frac{\pi^2}{48}\right)\rho^2
-\frac{607}{162}-\frac{\pi^4}{48}+\frac{11}{36}\zeta_3+\frac{335\pi^2}{432},\\
f_t={}&\frac19L^3+\frac59L^2+\left(-\frac{\rho^2}{12}-\frac{\pi^2}{6}+\frac{28}{27}\right)L
+\frac16M_{1,1}-\frac{5\rho^2}{36}-\frac{25\pi^2}{108}-\frac19\zeta_3+\frac{82}{81}.
\end{aligned}
$$

These are coefficients of log W_R. For the ordinary Wilson-loop expansion the
second-order coefficient divided by CF is F2+CF*F1²/2. Both objects are recorded,
so a comparison to ordinary W2 cannot be confused with one to w2.

## Test the result at several independent levels

The pole template above, the D-dimensional derivative of logZ, crossing and the
finite equations

$$\partial_LF_1=-2L,\qquad \partial_LF_2-b_0F_1=-2k_2L-h_2$$

all pass exactly by color. RG permits F2+Phi(z) for any L-independent Phi. This
blind direction is explicitly demonstrated rather than treated as a correctness
test of the finite shape.

For that shape the [author reference translation](reference_extract.py) starts
from an independent ordinary Wilson-loop result, converts alpha_s/(4pi) to
lambda, the author's fundamental TR convention to T, and their logarithms to
L,rho. It computes W2-W1²/2, with a common branch for the paired trilogarithms.
The complete finite result, not just its derivative, agrees exactly. Original
source hashes and attribution are in [reference provenance](reference/provenance.json).
Fresh M11 quadrature and 40/60-digit finite samples check the function-reading
conventions. They are not a new multidimensional numerical evaluation of each
two-loop graph.

Finally the counterterm fixes the required epsilon budget. One-loop epsilon3 is
needed for the supplied two-loop P2 through epsilon2. For three-loop finite
assembly, one-loop epsilon2 and two-loop epsilon1 suffice. Those lower-loop bare
coefficients are explicitly reassembled here and match the accepted three-loop
inputs. Publishing only finite one/two-loop answers would omit necessary data.

The evidence supports the stated research calculation and reused identities;
no new formal Lean chain or global complex-z continuation is claimed.
