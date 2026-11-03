# 4e-03: a continuous source-to-Laurent derivation

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


## Evaluate the mixed ordered chain rather than treating it as a square

Attachments are (2,3,3,4), m=4; the original outer-to-inner order is
s21,s32,s41,s31, with s31<s32. The full accepted J is

$$J=z^\epsilon\int_0^1ds_{21}\int_0^1ds_{32}\int_0^1ds_{41}
\int_0^{s_{32}}ds_{31}\,
s_{32}^{\epsilon-1}s_{41}^{\epsilon-1}
(1-s_{21})^{\epsilon-1}(1-s_{31})^{\epsilon-1}.$$

This expression includes the tangent numerator already combined with the
adjacent distance powers. No change of variables is required; Jacobian1.
The independent s21 and s41 integrations each give1/epsilon. Put v=s32,
u=s31 only as a renaming. The remaining integral is

$$T=\int_0^1v^{\epsilon-1}dv\int_0^v(1-u)^{\epsilon-1}du
=\frac1\epsilon\int_0^1v^{\epsilon-1}[1-(1-v)^\epsilon]dv.$$

For Re epsilon>0 both terms are integrable. The first is1/epsilon and the
second is B(epsilon,1+epsilon)=R2/epsilon, where
R2=Gamma(1+epsilon)²/Gamma(1+2epsilon). Consequently

$$T=\frac{1-R2}{\epsilon^2},\quad
J=\frac{z^\epsilon(1-R2)}{\epsilon^4},\quad
\mathcal B_2[J]=\frac{1-R2}{\epsilon^4c2},\quad
I=\frac{Z^2(R2-1)}{4\epsilon^4c2}.$$

The centered z^epsilon power becomes1/c2, not1. Its effect survives in the
finite and positive Laurent terms. To expand, use

$$\log R2=\sum_{k=2}^\infty
\frac{(-1)^k(2-2^k)\zeta_k}{k}\epsilon^k.$$

The same exponential recurrence gives the ordinary coefficients b_n of R2-1, with b0=b1=0.
Thus a_q=sum_{r=0}^{q+4}p_r b_{q+4-r}/4. The budget is R2 through6 and
p through4; the two vanished b coefficients remove the apparent fourth and
third poles. In particular

$$a_2=\frac14[b_6+p_2b_4+p_3b_3+p_4b_2],$$

with every b and p explicitly fixed by the displayed recurrence. No new
multi-variable special function is introduced. The final occurrence of M01
is merely the exact replacement ell²=-2M01-pi².


## The one moment used to compactify powers of log z

Only M01 is needed in this graph. In the project's fixed M convention,

$$r=\frac z{1+z},\qquad
M01=\int_0^1\frac{\ell+\log[(1-\beta)/\beta]}{\beta-r}d\beta
=-\frac{\ell^2+\pi^2}{2}.$$

The internal point is removable. To evaluate it set beta=z/(z+y) and
use the combined kernel1/(z+y)+1/(1-y). Its subtracted Euler integral is
T(a)=pi[cot(pi a)-exp(a ell)/sin(pi a)]+ell; differentiating at a=0 gives
the displayed M01. This is an elementary evaluated function, not an extra
source integral. The fixed real branch z>0 is used throughout.

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

$$a_{-2}(z)=- \frac{\pi^{2}}{24}.$$

$$a_{-1}(z)=\frac{\zeta\left(3\right)}{2}.$$

$$a_{0}(z)=- \frac{\pi^{2} M_{0,1}}{24} - \frac{49 \pi^{4}}{1440}.$$

$$a_{1}(z)=\frac{M_{0,1} \zeta\left(3\right)}{2} + \frac{3 \zeta\left(5\right)}{2} + \frac{2 \pi^{2} \zeta\left(3\right)}{9}.$$

$$a_{2}(z)=- \frac{5 \pi^{2} M_{0,1}^{2}}{144} - \frac{23 \pi^{4} M_{0,1}}{480} - \frac{223 \pi^{6}}{12096} + \frac{5 \zeta^{2}\left(3\right)}{6}.$$


The remainder is O(epsilon³) at fixed real z>0. See [reproduction](reproduce.md) and [validation](validation.md) for the distinction between adopted integration and the current replay.
