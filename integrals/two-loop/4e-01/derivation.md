# 4e-01: a continuous source-to-Laurent derivation

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


## Integrate the actual two ordered adjacent chains

Attachments are (3,3,4,4), m=4. The original outer-to-inner order is
s32,s42,s31,s41, with 0<s31<s32<1 and 0<s41<s42<1. The full numerator
and propagators are visible in

$$J=z^2\int_0^1ds_{32}\int_0^1ds_{42}
\int_0^{s_{32}}ds_{31}\int_0^{s_{42}}ds_{41}
[z(1-s_{31})s_{41}]^{\epsilon-1}
[z(1-s_{32})s_{42}]^{\epsilon-1}.$$

There is no nontrivial coordinate substitution; Jacobian1. Each edge has two
identical factors. For any integrable f,
integral_{0<u<v<1}f(u)f(v)du dv=(integral_0^1f)^2/2, because the two
triangles partition the square up to its null diagonal and the product is
symmetric. For f=(1-u)^(epsilon-1), respectively u^(epsilon-1), the chain
therefore equals1/(2epsilon²). This proves, without erasing the ordering,

$$J=\frac{z^{2\epsilon}}{4\epsilon^4},\qquad
\mathcal B_2[J]=\frac1{4\epsilon^4},\qquad
I=-\frac{Z^2}{16\epsilon^4}.$$

The z cancellation is the explicit crossing sum. There are no remaining
parameter subintegrals. With e_n=[epsilon^n]Z², every required coefficient
is a_q=-e_{q+4}/16. The budget is e0...e6: using only the finite part of Z²
would miss the requested epsilon² term. For example

$$e_6=\frac{\zeta_6}{3}+\frac{\zeta_2\zeta_4}{2}
+\frac{2\zeta_3^2}{9}+\frac{\zeta_2^3}{6},\qquad a_2=-e_6/16.$$

The absence of e1 explains the exact zero at epsilon^-3. This result uses
only zeta constants, so no M/N subintegral is needed for this entry.


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

$$a_{-4}(z)=- \frac{1}{16}.$$

$$a_{-3}(z)=0.$$

$$a_{-2}(z)=- \frac{\pi^{2}}{96}.$$

$$a_{-1}(z)=- \frac{\zeta\left(3\right)}{24}.$$

$$a_{0}(z)=- \frac{7 \pi^{4}}{5760}.$$

$$a_{1}(z)=- \frac{\pi^{2} \zeta\left(3\right)}{144} - \frac{\zeta\left(5\right)}{40}.$$

$$a_{2}(z)=- \frac{31 \pi^{6}}{241920} - \frac{\zeta^{2}\left(3\right)}{72}.$$


The remainder is O(epsilon³) at fixed real z>0. See [reproduction](reproduce.md) and [validation](validation.md) for the distinction between adopted integration and the current replay.
