# 4e-06: a continuous source-to-Laurent derivation

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


## Factorization follows from disjoint variables in this source

The attachments are (1,2,3,4), m=1. The original outer-to-inner order is
s11,s31,s21,s41, all independent unit intervals. The numerator is (1+z)²:

$$J=(1+z)^2\int_{[0,1]^4}
\widetilde D_z(s_{11},s_{31})^{\epsilon-1}
D_z(s_{21},s_{41})^{\epsilon-1}d^4s,\quad
\widetilde D_z(u,v)=z(1-u)(1-v)+uv.$$

The transformation (s11,s31)->(1-s11,1-s31) has absolute Jacobian1,
interchanges0 and1 in both coordinates, and sends Dtilde to D_z. The other
two variables do not change. Because the resulting two pairs are disjoint,
and the opposite kernel is absolutely integrable for Re epsilon>-1,

$$J=(1+z)^2\left[\int_{[0,1]^2}D_z^{\epsilon-1}\right]^2
=z^\epsilon V(\epsilon)^2.$$

There is no ordering factor1/2. Reflection of both independent pairs gives
B2[J]=V²/c2. Inserting this source's actual m=1 yields

$$I=-\frac{Z^2V(\epsilon)^2}{16c2}.$$

The source and prefactor are analytic at epsilon0, so all four negative
powers vanish. If t_q=sum_{i=0}^q v_i v_{q-i}, then

$$a_q=-\frac1{16}\sum_{r=0}^q p_r t_{q-r}\quad(q=0,1,2).$$

Only V and p through epsilon2 are needed. In particular

$$a_2=-\frac1{16}[2v_0v_2+v_1^2+(\zeta_2-\ell^2/2)v_0^2].$$

The supplied low-moment identities reduce this to the coefficient table.
This source-level disjoint-variable proof is what licenses the product of
one-loop kernels; the factorization is not inferred just from its leg label.


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

## Every Laurent coefficient through epsilon²


The following are coefficients of the complete normalized I, not bare J. Unlisted orders below the stated leading pole are exactly zero; the displayed table includes every power from-4 through2.

$$a_{-4}(z)=0.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=0.$$

$$a_{-1}(z)=0.$$

$$a_{0}(z)=- \frac{M_{0,1}^{2}}{16}.$$

$$a_{1}(z)=- \frac{M_{0,1} M_{1,1}}{16}.$$

$$a_{2}(z)=- \frac{11 M_{0,1}^{3}}{192} - \frac{\pi^{2} M_{0,1}^{2}}{24} - \frac{M_{0,1} M_{2,1}}{64} - \frac{M_{1,1}^{2}}{64}.$$


The remainder is O(epsilon³) at fixed real z>0. See [reproduction](reproduce.md) and [validation](validation.md) for the distinction between adopted integration and the current replay.
