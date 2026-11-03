# 4e-04: from the ordered source to evaluated K and M/N coefficients

## Domain, accepted root and physical normalization

The starting object is this representative's accepted prepared C/J, rather
than a fresh rederivation of the QCD path integral. Set D=4-2epsilon, z>0 and
ell=log z. The inherited distance prescription is (-2v14+i0); for positive
real z all bases below are nonnegative and their powers use the real log.
The source is bound by [input.json](input.json). Its exact ordered bounds
are retained below. The Fourier propagator and color projection are upstream
accepted facts; they are not established by a Python replay.

The common prepared exchange prefactor is

$$C_{4e}=-\frac{C_AC_F}{2}e^{-2\gamma_E\epsilon}(4\pi\alpha_s)^2
\frac{\Gamma(1-\epsilon)^2}{(4\pi^{2-\epsilon})^2}
\frac{(-2v_{14}+i0)^{2\epsilon}}4.$$

The projection is a1=-1/2 in the project's fixed color basis. Removing the
common coupling, color and symmetric scale, while including the orbit size m,
gives exactly

$$Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad
\mathcal B_2[J]=\frac{J(z)+z^{2\epsilon}J(1/z)}{1+z^{2\epsilon}},\quad
I=-\frac m{16}Z^2\mathcal B_2[J].$$

Thus m, crossing and the single two-loop Euler compensation are already in I.
In particular neither the source's canonical diagram count nor a second Euler
factor is multiplied into the coefficient table. The target is every ordinary
Laurent coefficient a_q in I=sum_{q=-4}^2 a_q epsilon^q+O(epsilon³).

Write D_z(t,s)=(1-t)(1-s)+zts, c1=cosh(epsilon ell/2),
c2=cosh(epsilon ell). For the sources with cusp endpoints start with
Re epsilon>0. Each original integrand is bounded by an integrable product
of the appropriate adjacent endpoint powers and opposite kernels in that
region; this justifies the specific separations below. The resulting
meromorphic expressions define the local epsilon continuation. No general
Fubini bridge for unrelated prepared inputs is claimed.



## The ordered source and its two exact integrations

This representative has attachments (2,2,4,4), m=2 and source order
t2,s2,t1,s1, with 0<t1<t2<1 and 0<s1<s2<1. Its actual J is

$$J=(1+z)^2\int_0^1dt_2\int_0^1ds_2\int_0^{t_2}dt_1
\int_0^{s_2}ds_1\,
D_z(t_1,s_1)^{\epsilon-1}D_z(t_2,s_2)^{\epsilon-1},\qquad Q=z^{-\epsilon}J.$$

This same-pair ordering is essential. Replacing it by the square of the
one-kernel integral would change the answer. Reflection
(t_i,s_i)->(1-t_{3-i},1-s_{3-i}) preserves the two ordered chains with unit
absolute Jacobian. Using D_{1/z}(1-t,1-s)=D_z(t,s)/z gives
J(1/z)=z^(-2epsilon)J(z), Q(1/z)=Q(z). Thus the physical object will be
I=-Z²Q/(8c2).

Initially take0<Re epsilon<1/4. A product of the integrable opposite kernels
on the full four-cube bounds the ordered domain, so the following specific
reordering is legitimate. Put u=t2,v=s1 (a renaming), and integrate t1 first,
then s2. The linear primitives give

$$\int_0^uD_z(t_1,v)^{\epsilon-1}dt_1=
\frac{D_z(u,v)^\epsilon-(1-v)^\epsilon}{\epsilon[(1+z)v-1]},$$
$$\int_v^1D_z(u,s_2)^{\epsilon-1}ds_2=
\frac{(zu)^\epsilon-D_z(u,v)^\epsilon}{\epsilon[(1+z)u-1]}.$$

Writing D=D_z(u,v), these evaluated subintegrals yield

$$Q=\frac{(1+z)^2z^{-\epsilon}}{\epsilon^2}
\int_{[0,1]^2}\frac{[D^\epsilon-(1-v)^\epsilon]
[(zu)^\epsilon-D^\epsilon]}
{[(1+z)v-1][(1+z)u-1]}du dv.\tag{1}$$

Both internal difference quotients are removable. They must stay whole;
the four expanded summands need not be independently integrable.

## A polynomial change exposes the one remaining ordered subintegral

Let p=(1+z)u-1,q=(1+z)v-1; both new bounds are[-1,z] and
du dv=dp dq/(1+z)². Use sigma=z/(1+z)², distinct from ell=log z, and

$$A(p)=(1+p)^\epsilon,\quad B(p)=(1-p/z)^\epsilon,\quad
C(p,q)=(1+pq/z)^\epsilon.$$

The transformed source is

$$Q=\frac{\sigma^\epsilon}{\epsilon^2}
\int_{[-1,z]^2}\frac{[C-B(q)][A(p)-C]}{pq}dp dq.\tag{2}$$

Add the term (B(p)-B(q))(C-1)/(pq), whose integral is zero by p<->q antisymmetry.
The complete addition is integrable across both coordinate axes. The new
numerator is exactly

$$(A(p)+B(p))(C-1)-(C^2-1)-(A(p)-1)(B(q)-1).\tag{3}$$

Define actual endpoint primitives

$$J_a(x)=\int_0^x\frac{(1-t)^a-1}{t}dt,\qquad
T_a(x)=\int_0^x\frac{J_a(t)}t dt.$$

Then integral_{-1}^z(C_a-1)dq/q=J_a(-p)-J_a(p/z). Set
f=J_epsilon(-p),g=J_epsilon(p/z), j-=J_epsilon(-z),
j0=J_epsilon(1),j+=J_epsilon(-1/z), and
C_a=T_a(-z)+T_a(-1/z)-2T_a(1). Applying this primitive to(3), then using
f'=(A-1)/p,g'=(B-1)/p, leaves

$$\begin{aligned}Q=\frac{\sigma^\epsilon}{\epsilon^2}\bigg[
&\frac12[f^2-g^2]_{-1}^z+\int_{-1}^z(fg'-gf')dp
+2C_\epsilon-C_{2\epsilon}\\
&-(j_--j_0)(j_0-j_+)\bigg].\end{aligned}\tag{4}$$

Here (f'+g')(f-g)=d[(f²-g²)/2]/dp+fg'-gf'. All integration-by-parts
boundaries are explicitly present. Split the remaining integral at p=0.
On its positive half put p=zy, and on its negative half p=-y, retaining the
reversed bounds. One obtains

$$\int_{-1}^z(fg'-gf')dp=K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z),$$
$$K(a,b;z)=\int_{0<x<y<1}\frac{
[(1+zx)^a-1][(1-y)^b-1]-[(1-x)^b-1][(1+zy)^a-1]}
{xy}dxdy.$$

The products at p=-1,z in(4) now give

$$Q=\frac{\sigma^\epsilon}{\epsilon^2}\left[
K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z)
+\frac{(j_-+j_+)^2}{2}-j_0(j_-+j_+)+2C_\epsilon-C_{2\epsilon}\right].\tag{5}$$

## Evaluate K and put its value back into this same Q

Write A_a(y)=(1+zy)^a-1,B_b(y)=(1-y)^b-1,
f(y)=J_a(-zy),g(y)=J_b(y). Then K=integral_0^1(fg'-gf')dy.
The Euler derivative theta=z partial_z gives theta f=A_a(y). Integrating
only the second term by parts therefore removes the nested primitive:

$$\theta K=2P(a,b;z)-J_b(1)[(1+z)^a-1],\qquad
P=\int_0^1\frac{[(1+zy)^a-1][(1-y)^b-1]}y dy.$$

The y=0 boundary is zero; the y=1 boundary displayed above is indispensable.
For|z|<1 expand the first power and integrate each term with its Beta integral:

$$P=\sum_{n\ge1}\frac{(-a)_n(-z)^n}{n!}
\left[B(n,1+b)-\frac1n\right]
=\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{(-z)^n}n-J_a(-z).$$

One integration in dz/z, with K(a,b;0)=0, thus gives a complete evaluation

$$K(a,b;z)=2F(a,b;-z)-2F(a,0;-z)-J_b(1)J_a(-z),$$
$$F(a,b;x)=\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{x^n}{n^2}
=-\frac{ax}{1+b}{}_4F_3\!\left(\begin{matrix}1,1,1,1-a\\2,2,2+b\end{matrix};x\right).$$

Analytic continuation along the negative real hypergeometric argument extends
this derivation to all z>0. It is not a fitted formula. Define
Psi(a,b)=F(a,b;-z)+F(a,b;-1/z)-2F(a,b;1).
The one-kernel result has the endpoint form
V(epsilon)=-sigma^(epsilon/2)(j-+j+-2j0)/epsilon. Also
F(epsilon,epsilon;1)=T_epsilon(1)+j0²/2. Indeed, continue the evaluated K to z=-1 within its convergent endpoint domain: its two weights then coincide, so its antisymmetric integrand vanishes, and0=2F(epsilon,epsilon;1)-2F(epsilon,0;1)-j0². Since F(a,0)=T_a, the endpoint identity follows. Substituting K into(5) cancels
all single endpoint products into

$$Q=\frac12V(\epsilon)^2+
\frac{\sigma^\epsilon}{\epsilon^2}
[2\Psi(\epsilon,\epsilon)-\Psi(2\epsilon,0)].\tag{6}$$

This is already an evaluated standard-function answer to the original source.
The conversion below expresses its coefficients in the adopted M/N family.

## Why the M/N generator, including its endpoint term, is the same answer

Use the moment definitions below and introduce only temporary integrals

$$Y=\int_0^1\frac{e^{\epsilon W}L}{\beta-r}d\beta,\quad
V_0=\int_0^1\frac{e^{\epsilon W}-\sigma^\epsilon}{\beta-r}d\beta,\quad
Y_c=\int_0^1\frac{e^{\epsilon W}L\cosh(\epsilon L)}{\beta-r}d\beta,$$
$$E=2\int_0^1\frac{e^{\epsilon W}[\cosh(\epsilon L)-1]}{\beta-r}d\beta.$$

Let d=J_{2epsilon}(1),J-=J_{2epsilon}(-z),J+=J_{2epsilon}(-1/z).
The [conversion appendix](appendices/endpoint-conversion.md) proves the
necessary derivative identities and fixes their analytic constants:

$$\sigma^\epsilon\Psi(\epsilon,\epsilon)
=Y-\ell V_0-\sigma^\epsilon M_{0,1},\tag{7}$$
$$\sigma^\epsilon[\Psi(2\epsilon,0)-\ell(J_--J_+)+2M_{0,1}]
=2Y_c+2\epsilon dV(2\epsilon),\tag{8}$$
$$E=\sigma^\epsilon(J_+-J_-)-2V_0.\tag{9}$$

These are exact small-epsilon identities, with their full definitions and
boundary constants; none is an assumption chosen to match the finite table.
Put(7),(8) in(6): the two M01 boundaries cancel. Then(9) gives

$$Q=\frac12V(\epsilon)^2+
\frac{2Y-2Y_c+\ell E-2\epsilon dV(2\epsilon)}{\epsilon^2}.$$

All three temporary integrals have the same denominator. Combining their
numerators, and using R=L-ell, gives

$$2Y-2Y_c+\ell E=-2\int_0^1
\frac{e^{\epsilon W}R[\cosh(\epsilon L)-1]}{\beta-r}d\beta.$$

Finally the endpoint Beta derivative is evaluated,
d=integral_0^1[(1-t)^(2epsilon)-1]dt/t=-gamma_E-psi(1+2epsilon).
Define chi=gamma_E+psi(1+2epsilon). We have derived the complete mother

$$\boxed{Q=\frac12V(\epsilon)^2-2V_N(\epsilon)
+\frac{2\chi(\epsilon)}\epsilon V(2\epsilon)},\qquad
V_N=\int_0^1\frac{e^{\epsilon W}R}{\beta-r}
\frac{\cosh(\epsilon L)-1}{\epsilon^2}d\beta.\tag{10}$$

No K,Psi,Y or unknown remainder survives. The essential chi/epsilon is
regular: chi/epsilon=2zeta2-4zeta3 epsilon+8zeta4 epsilon²+O(epsilon³).
The source is analytic at epsilon0. Thus every negative Laurent coefficient
is zero, rather than canceled by a postulated counterterm.

## The evaluated one-kernel family and its finite moment rule

The subintegrals adopted here have one common, fixed definition. Put

$$r=\frac z{1+z},\quad W=\log[\beta(1-\beta)],\quad
R=\log\frac{1-\beta}{\beta},\quad L=\ell+R,$$
$$M_{k,j}=\int_0^1\frac{W^kL^j}{\beta-r}d\beta,\qquad
N_{k,j}=\int_0^1\frac{W^kRL^j}{\beta-r}d\beta
=M_{k,j+1}-\ell M_{k,j}.$$

These are raw powers, with no factorials. Since L(r)=0, each j>=1 density
has a removable internal point and the integrals require no principal value.
The complete opposite kernel is explicitly integrated once:

$$\int_0^1D_z(t,s)^{\epsilon-1}ds=
\frac{(zt)^\epsilon-(1-t)^\epsilon}
{\epsilon[(1+z)t-1]}.$$

Use beta=1-t (absolute Jacobian1, endpoints reverse) and retain the full
difference quotient. Its numerator is
2z^(epsilon/2) exp(epsilon W/2)sinh(epsilon L/2). Therefore

$$V(\epsilon)=-(1+z)z^{-\epsilon/2}
\int_{[0,1]^2}D_z(t,s)^{\epsilon-1}dt ds
=\int_0^1\frac{e^{\epsilon W/2}}{\beta-r}
\frac{\sinh(\epsilon L/2)}{\epsilon/2}d\beta
=\sum_{q\ge0}v_q\epsilon^q,$$
$$v_q=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,2j+1}}{2^q(q-2j)!(2j+1)!}.$$

Thus v0=M01, v1=M11/2, v2=M21/8+M03/24,
v3=(M31+M13)/48, v4=M41/384+M23/192+M05/1920.
This defines a fully evaluated member of the project's established function
family, rather than naming an unknown remainder. The adopted M/N HLOG
dictionary and the source moment integrals retain their original validation
scope. The convention is G(a,w;z)=integral_0^z G(w;t)dt/(t-a), G(empty;z)=1;
its released words have letters0,-1 and zeta periods. Real z>0 fixes the
branch. These statements do not supply a z=-1 or general complex-z theorem.

For the reductions used in the final compact table,

$$M_{0,1}=-\frac{\ell^2+\pi^2}{2},\quad
M_{0,3}=-M_{0,1}^2,\quad
N_{0,2}=\frac{M_{0,1}^2}{3}+4\zeta_2M_{0,1}.$$

They are fixed moment identities, not fitted source coefficients. For example
beta=z/(z+y) gives the combined removable kernel
1/(z+y)+1/(1-y), so

$$\int_0^\infty(y^a-1)\left[\frac1{z+y}+\frac1{1-y}\right]dy
=\pi\left[\cot(\pi a)-\frac{e^{a\ell}}{\sin(\pi a)}\right]+\ell
=\sum_{j\ge1}\frac{M_{0,j}}{j!}a^j.$$

The complete integral converges for -1<Re a<1. Its first and third Taylor
coefficients give M01 and M03; N02 follows from N02=M03-ell M02 with the
same second Taylor coefficient. This supplies the low-weight constant
terms without a numerical ansatz.



## How prefactors and positive orders are assembled

All coefficient extraction uses ordinary power series. If
E(epsilon)=exp(sum_{k>=1} l_k epsilon^k)=sum e_n epsilon^n, use

$$e_0=1,\qquad e_n=\frac1n\sum_{k=1}^n k l_k e_{n-k}.$$

For Z², l1=0 and l_k=2 zeta_k/k (k>=2). For
s_n=[epsilon^n]sech(epsilon ell),

$$s_0=1,\qquad s_n=-\sum_{j=1}^{\lfloor n/2\rfloor}
\frac{\ell^{2j}}{(2j)!}s_{n-2j},\qquad
p_n=[\epsilon^n](Z^2/c2)=\sum_{i=0}^ne_i s_{n-i}.$$

In particular p0=1,p1=0,p2=zeta2-ell²/2. These finite recurrences and the
source coefficient formulas are the actual coefficient engine's operations.
They include pole times positive-order products. Even powers ell² may finally
be replaced by -2M01-pi²; the replacement is applied after convolution.
The table below uses that compact M/N form. The adopted dictionary's exact
whole-HLOG comparison checks the same function, not only selected samples.



## The complete finite-order assembly, including epsilon²

Expanding the now evaluated VN produces

$$V_N=\sum_{q\ge0}n_q\epsilon^q,\qquad
n_q=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,2j+2}}{(q-2j)!(2j+2)!}.$$

Here n0=N02/2,n1=N12/2,n2=N22/4+N04/24. With v_q as defined above,
equation(10) yields the explicit finite coefficient rule

$$Q_q=\frac12\sum_{i=0}^qv_i v_{q-i}-2n_q+
2^{q+2}\sum_{k=0}^q(-1)^k\zeta_{k+2}v_{q-k}.$$

The factor2^(q+2) combines the2epsilon in chi with that in V(2epsilon).
It would be incorrect to apply the scale factor only to one of them. To
obtain a0,a1,a2, the only needed subintegrals are V and VN through2:
M01,M11,M21,M03,N02,N12,N22,N04, plus the displayed zeta constants.

The physical prefactor is P=Z²/c2=1+p2 epsilon²+O(epsilon³),
p2=zeta2-ell²/2. Therefore

$$a_0=-Q_0/8,\qquad a_1=-Q_1/8,\qquad
a_2=-[Q_2+(\zeta_2-\ell^2/2)Q_0]/8,$$
$$Q_0=\frac12M01^2-N02+4\zeta_2M01=\frac16M01^2,$$
$$Q_1=\frac12M01M11-N12+4\zeta_2M11-8\zeta_3M01,$$
$$Q_2=M01(M21/8+M03/24)+M11^2/8-N22/2-N04/12
+2\zeta_2M21+\frac23\zeta_2M03-8\zeta_3M11+16\zeta_4M01.$$

The highest-order line explicitly includes both required source moments and
the prefactor correction. Its compact table uses M03=-M01² and the next
Euler moment identity

$$M05=-\frac{\ell^6}{6}-\frac{5\pi^2\ell^4}{6}
-\frac{7\pi^4\ell^2}{6}-\frac{\pi^6}{2},$$

which is the fifth Taylor coefficient of the same T(a) below. Its fourth coefficient is M04=-ell⁵/5-2pi²ell³/3-7pi⁴ell/15. Together with N04=M05-ell M04 they imply the exact conversion

$$N04=\frac45M05-\frac43M01^3+\frac8{15}\pi^4M01.$$

Applying this identity to Q2 and the prefactor term gives the printed compact table (which retains an equivalent, dependent N04). This identity
explains the alternative M05 term in the compact representative; it does not
add a new independent source integral. These finite formulas reconstruct the
entire table before any comparison with a historical result.

An alternative historical HyperInt route expands the actual two-dimensional
kernel(1): with DeltaL=logD-log(1-v), DeltaU=log(zu)-logD,
Phi=logD+[log(1-v)+logu-ell]/2 and
K0=(1+z)²DeltaL DeltaU/([(1+z)v-1][(1+z)u-1]), it integrates
Q0=integral K0,Q1=integral K0 Phi,Q2=integral K0[Phi²/2+(DeltaL²+DeltaU²)/24].
Those integrated HLOG files are reused by the replay as a separate Q-to-I
comparison. Maple itself is not run by this trial. The analytic route(1)--(10)
above supplies the continuous source bridge independently of their filenames.

## Every Laurent coefficient through epsilon²


The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=0.$$

$$a_{0}(z)=- \frac{M_{0,1}^{2}}{48}.$$

$$a_{1}(z)=- \frac{M_{0,1} M_{1,1}}{16} + M_{0,1} \zeta\left(3\right) - \frac{\pi^{2} M_{1,1}}{12} + \frac{N_{1,2}}{8}.$$

$$a_{2}(z)=- \frac{41 M_{0,1}^{3}}{576} - \frac{M_{0,1} M_{2,1}}{64} + \frac{M_{0,5}}{30} - \frac{M_{1,1}^{2}}{64} + M_{1,1} \zeta\left(3\right) - \frac{\pi^{2} M_{2,1}}{24} - \frac{N_{0,4}}{32} + \frac{N_{2,2}}{16}.$$


The remainder is O(epsilon³) at fixed z>0. [Reproduction](reproduce.md) and [validation](validation.md) distinguish adopted mathematical integration from the current replay.
