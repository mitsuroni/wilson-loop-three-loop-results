# Turn the opposite exchange into an evaluated generating function

This current reading-layer trial adopts the source calculation below. Its descriptions
of executed checks refer to the preserved publication run. The new per-entry
trial replay is recorded separately in [validation](validation.md); it does not
redo a cold CAS integration. The starting momentum rule is an accepted physical
input, and the parameter-to-coefficient equalities are the mathematical content
being exposed here.


The opposite exchange complements the adjacent cusp integral: the overall
momentum prefactor is the same, yet its parameter integral is finite at
epsilon=0. We want to understand that difference directly from the path
geometry, then obtain the positive epsilon orders needed by higher-loop
counterterms. A single actual integration reveals both the all-epsilon
answer and the M-moment organization of its coefficients.

Use the null quadrilateral and definitions s=2v1.v2, t=2v1.v4, z=s/t>0,
D=4-2epsilon from the [adjacent example](../2e-01/derivation.md). The original
parameters are u=s21 and v=s41, with outer-to-inner order u then v, each in
[0,1]. The [source index](input.json) binds the graph, raw native package,
original ordered input and prepared0.3 input without changing their bytes.

## The path fixes the denominator and numerator together

The representative momentum contribution is

$$
\mathcal A_{24}=\frac1{N_c}\operatorname{tr}
\int_0^1du\int_0^1dv
(ig_0v_4^\nu T^b)(ig_0v_2^\mu T^a)
\int\frac{d^Dl}{(2\pi)^D}
\frac{-ig_{\mu\nu}\delta^{ab}}{l^2+i0}
e^{il\cdot\Phi},
\quad \Phi=-[v_1+uv_2+(1-v)v_4].
$$

The source uses late Wilson generators on the left. Contracting color gives
CF; using the same Fourier propagator as the adjacent exchange and writing
b=v1.v4 gives

$$
D_z(u,v)=\frac{\Phi^2}{2b}=(1-u)(1-v)+zuv,
\qquad \frac{v_2\cdot v_4}{b}=-(1+z).
$$

The adopted runtime script derives these two polynomials from the source phase and the
null Gram matrix. The original integral and its raw scalar factor are

$$
J_O(z,\epsilon)=-(1+z)\int_0^1du\int_0^1dv\,
D_z(u,v)^{\epsilon-1},\qquad
C_{\rm raw}=-\frac12C_F\lambda_{\rm tag}
\Gamma(1-\epsilon)(\pi[-2b])^\epsilon.
$$

Native and older original-ordered C.combined agree exactly; prepared C
multiplies them by exp(-gamma_E epsilon) exactly once. The source cube is
only x1=u, x2=v, with Jacobian1. There is no inherited leg/Gamma sign change
in this pure exchange.

The only zeroes of D in the closed square occur at (u,v)=(0,1),(1,0).
Near either corner D is comparable to the sum of two nonnegative distances.
The radial integral is therefore of the form integral r^(Re epsilon)dr:
J is absolutely convergent for Re epsilon>-1. On a compact positive-z set,
the same bound with any fixed logarithm power justifies local Taylor
expansion near epsilon=0. It also justifies the reflection and changes of
integration order used below. This is a special proof for the present
kernel, not a general permission to permute prepared input orders.

## Integrate one parameter before expanding

For fixed u set a=zu, d=1-u. The denominator is d+(a-d)v, so the primitive
and definite integral are

$$
\int_0^1[d+(a-d)v]^{\epsilon-1}dv
=\frac{a^\epsilon-d^\epsilon}{\epsilon(a-d)}.
$$

Its derivative is replayed exactly. At a=d the expression is interpreted by
continuity; the difference quotient has no internal pole. With beta=1-u
and r=z/(1+z), the original ordered integral becomes

$$
\boxed{J_O(z,\epsilon)=\frac1\epsilon\int_0^1
\frac{z^\epsilon(1-\beta)^\epsilon-\beta^\epsilon}{\beta-r}\,d\beta.}
$$

The numerator vanishes at beta=r. It must remain a complete numerator;
splitting it into two separately singular integrals would change the task.

Define rho=log z, W=log[beta(1-beta)] and
L_beta=rho+log[(1-beta)/beta]. The two endpoint powers now give

$$
z^\epsilon(1-\beta)^\epsilon-\beta^\epsilon
=2z^{\epsilon/2}e^{\epsilon W/2}
\sinh(\epsilon L_\beta/2).
$$

Removing the common scale power defines a crossing-even function

$$
H_\epsilon=z^{-\epsilon/2}J_O
=\int_0^1\frac{e^{\epsilon W/2}}{\beta-r}
\frac{\sinh(\epsilon L_\beta/2)}{\epsilon/2}\,d\beta.
$$

This is an exact generating representation, not yet the final normalized I.
Introduce raw moments without factorial normalization,

$$
M_{k,j}(z)=\int_0^1\frac{W^kL_\beta^j}{\beta-r}\,d\beta,
\qquad k\ge0,\ j\ge1.
$$

Their internal point is removable because L_beta(r)=0. In particular the
j=1 limiting density is -W(r)^k/[r(1-r)], while higher j give zero there.
Expanding the proven local analytic generating representation yields

$$
\boxed{H_\epsilon=\sum_{k,j\ge0}
\frac{\epsilon^{k+2j}M_{k,2j+1}}
{2^{k+2j}k!(2j+1)!}.}
$$

Thus h0=M01, h1=M11/2, h2=M21/8+M03/24 and
h3=(M31+M13)/48. This obtains the positive orders from the source kernel;
there is no coefficient fit or old-table input.

## The all-epsilon integral is a standard Meijer G

To remove the remaining parameter integral entirely, apply the Mellin-Barnes
identity to zuv+(1-u)(1-v), then perform both Euler Beta integrals. Initially
take real 0<epsilon<1 and epsilon-1<Re w<0. The result is

$$
J_O=-\frac{1+z}{\Gamma(1-\epsilon)\Gamma(1+\epsilon)^2}
\frac1{2\pi i}\int dw\,
\Gamma(-w)\Gamma(\epsilon-w)^2
\Gamma(1-\epsilon+w)\Gamma(1+w)^2z^w.
$$

The strip makes the original Beta integrals convergent and separates the
left and right Gamma poles. Matching the defining Mellin-Barnes kernel of
the standard Meijer G gives

$$
\boxed{J_O=-\frac{1+z}{\Gamma(1-\epsilon)\Gamma(1+\epsilon)^2}
G_{3,3}^{3,3}\!\left(z\;\middle|\;
\begin{matrix}0,0,\epsilon\\0,\epsilon,\epsilon\end{matrix}\right).}
$$

This continues meromorphically in epsilon from the common strip. It is the
same closed special-function answer found in the archived research; here
the complete source connection and normalization are written explicitly.
The historical publication run also evaluates the original difference quotient and this
Meijer G at two noninteger parameter points.

## Restore the orbit and crossing only once

The actual source reflection is

$$D_{1/z}(1-u,1-v)=D_z(u,v)/z,$$

so J_O(1/z,epsilon)=z^(-epsilon)J_O(z,epsilon). The C4 multiplicity is m=2.
Using the same C/I identity as the adjacent example, with
G=exp(-gamma_E epsilon)Gamma(1-epsilon), we therefore get

$$
\begin{aligned}
I_O&=-\frac{G}{2}
\frac{J_O(z,\epsilon)+z^\epsilon J_O(1/z,\epsilon)}{1+z^\epsilon}\\
&=-\frac{GJ_O(z,\epsilon)}{1+z^\epsilon}
=-\frac{G}{2\cosh(\epsilon\rho/2)}H_\epsilon.
\end{aligned}
$$

Equivalently the full evaluated answer is

$$
\boxed{I_O=\frac{(1+z)z^{-\epsilon/2}e^{-\gamma_E\epsilon}}
{2\cosh(\epsilon\rho/2)\Gamma(1+\epsilon)^2}
G_{3,3}^{3,3}\!\left(z\;\middle|\;
\begin{matrix}0,0,\epsilon\\0,\epsilon,\epsilon\end{matrix}\right).}
$$

This I includes m=2, crossing and one Euler compensation. CF, lambda and
the common S1 scale are removed; neither an extra m nor another W-to-log
subtraction is applied. The normalization is fixed from C/J, not from RG.

## Two elementary moments make the coefficients compact

The relations for M01 and M03 can themselves be obtained from the original
moment kernel. Set beta=z/(z+y); reversing its endpoints gives

$$
\frac{-d\beta/dy}{\beta-r}=\frac1{z+y}+\frac1{1-y}.
$$

For -1<Re a<1 the complete subtracted integral

$$
T(a)=\int_0^\infty(y^a-1)
\left(\frac1{z+y}+\frac1{1-y}\right)dy
$$

is ordinary: the subtraction removes y=1, and the combined kernel controls
both outer endpoints. Standard Euler/Mellin evaluation, or first evaluating
the two common principal-value terms and then subtracting their value at
a=0, gives

$$
T(a)=\pi\left[\cot(\pi a)-\frac{e^{\rho a}}{\sin(\pi a)}\right]+\rho
=\sum_{j\ge1}\frac{a^j}{j!}M_{0,j}.
$$

The final subtracted integral and every M moment require no principal-value
prescription. The adopted symbolic replay checks the variable-change identity
and both required Taylor coefficients, giving

$$M_{0,1}=-\frac{\rho^2+\pi^2}{2},\qquad M_{0,3}=-M_{0,1}^2.$$

These are derived identities, not chosen numerical fits. The same original
moments are independently integrated in the numerical controls.

## Laurent coefficients through epsilon3

Expand only the analytic generating representation and the exact prefactor
G/cosh(epsilon*rho/2). Before compactification the epsilon² coefficient is
-M21/16-M03/48-pi² M01/24+rho² M01/16. Applying the two elementary moment
identities gives the constant-coefficient M representation

$$
\boxed{\begin{aligned}
I_O={}&-\frac12M_{0,1}-\frac14M_{1,1}\epsilon\\
&-\left(\frac{M_{2,1}}{16}
+\frac5{48}M_{0,1}^2+\frac{5\pi^2}{48}M_{0,1}\right)\epsilon^2\\
&-\left(\frac{M_{3,1}+M_{1,3}}{96}
+\frac{M_{0,1}M_{1,1}}{16}
+\frac{\zeta_3M_{0,1}}6+\frac{5\pi^2M_{1,1}}{96}\right)
\epsilon^3+O(\epsilon^4).
\end{aligned}}
$$

Both pole coefficients are zero. At epsilon=0 this reduces to
(rho²+pi²)/4, whereas the adjacent exchange retained its double pole.
The finite-order sum generated from H determines any fixed higher order;
only through epsilon3 is explicitly exported and checked in this release.


## A finite coefficient rule, including the highest requested order

Let h_q be the finite moment sum above and write
P(epsilon)=G(epsilon)sech(epsilon ell/2)=sum p_n epsilon^n.
Generate g_n=[epsilon^n]G from g_0=1 and

$$g_n=\frac1n\sum_{k=2}^n\zeta_k g_{n-k}.$$

For s_n=[epsilon^n]sech(epsilon ell/2), use s_0=1 and

$$s_n=-\sum_{j=1}^{\lfloor n/2\rfloor}
\frac{(\ell/2)^{2j}}{(2j)!}s_{n-2j},\qquad
p_n=\sum_{i=0}^ng_i s_{n-i}.$$

Then a_q=-sum_{r=0}^q p_r h_{q-r}/2. In particular

$$a_3=-\frac12\left[h_3+
\left(\frac{\zeta_2}{2}-\frac{\ell^2}{8}\right)h_1
+\frac{\zeta_3}{3}h_0\right].$$

Only h_0,...,h_3 and p_0,...,p_3 are needed because this opposite exchange
is regular. Substituting the explicit h_q, M03=-M01² and
ell²=-2M01-pi² produces the displayed epsilon³ coefficient.
This is a finite reconstruction rule for every exported order, with no
uncomputed subintegral beyond the defined, adopted M family.

The adopted M dictionary uses G(a,w;z)=integral_0^z G(w;t)dt/(t-a),
G(empty;z)=1 and G(0;z)=log z, with letters0,-1 and real values for z>0.
The M definition and the displayed finite moment rule are the exact function
representation here; its HLOG dictionary is a reused evaluated subintegral,
not a new numerical fit. Arbitrary complex-z or z=-1 boundary continuation
requires separate branch work.

## Checks and limits

The successful [record](validation.md) contains source geometry checks from the historical publication run, C/sign,
order, primitive, reflection, moment compactification and coefficient checks.
The historical numerical calculations compare the original two-dimensional J with its
one-dimensional primitive, the primitive with full Meijer G, and the newly
generated coefficients with epsilon-Cauchy readouts. M01, M03 and M11 are
also integrated from their original definitions.40 working digits and the
observed errors are precision controls, not rigorous interval certificates.

Existing r002 coefficients are opened only after saving the new prediction.
The old table has no epsilon1..3 entries, so it cannot supply these new
positive-order computations. The source QFT Fourier rule and standard
Euler/Mellin identities are reused mathematics; archived scripts are not
silently rerun or relabelled as new work. See [reproduce.md](reproduce.md)
and the [historical source manifest](../../families/low-loop-runtime/integrals/one-loop/source-manifest.json).
The result is a research computation for positive z with local epsilon
continuation. No arbitrary physical complex-z continuation or integral
formal proof is claimed. `formal=not_run`.
