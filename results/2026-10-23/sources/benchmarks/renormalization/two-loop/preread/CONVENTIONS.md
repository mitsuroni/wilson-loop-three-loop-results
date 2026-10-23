# Conventions and function definitions

## Observable and expansion

The calculation concerns massless QCD in Feynman gauge on a non-self-intersecting lightlike quadrilateral. Write v_i = x_(i+1) - x_i, with v_i squared equal to zero and the sum of v_i equal to zero. The scalar products v_ij mean v_i dot v_j. Define

\[
D=4-2\varepsilon,\quad s=2v_{12}=2v_{34},\quad t=2v_{14},\quad
z=s/t>0,\quad \ell=\log z,\quad
\lambda=\alpha_s^{\overline{\mathrm{MS}}}/\pi,\quad T=T_F n_f,\quad T_F=1/2.
\]

The analytic starting region is s,t < 0. The scale powers retain the prescribed +i0. The coefficient functions and comparisons use fixed real z > 0 and compatible real polylogarithm branches. The input name `z1` denotes the same ratio as `z`.

For every representative,

\[
I_r(\varepsilon,z)=\sum_{p=-4}^{2}a_{r,p}(z)\varepsilon^p+O(\varepsilon^3).
\]

These are ordinary Laurent coefficients, with no factorial multiplier. All seven positions are explicit in the machine-readable files. Equality of two normalized results, such as I_(4e-02) = I_(4e-05), does not identify their source orbits.

## From prepared C,J to normalized I

The `inputs/` files retain the original prepared prefactors, parameter kernels, dependent bounds and integration order. Their family-specific schemas and fields are preserved. These snapshots document the source of the results; the uniform `results.json` and `results.wl` are the interfaces for using the evaluated coefficients.

The normalized I values include the cyclic orbit factors, the scale-channel combination and exactly one Euler compensation factor appropriate to two loops. The exchange entries use the exponent color projection. The full C factors in the snapshots belong to C times J and are not additional factors to multiply into an already normalized I. The passage between these objects is given in report Section 3 and the appendix input dictionaries.

Let I_E be the sum of the six exchange representatives, I_V the sum of the four vertex representatives and I_Pi the sum of the two scalar self-energy representatives. Define

\[
a(\varepsilon)=\frac{10-6\varepsilon}{3-2\varepsilon},\qquad
b(\varepsilon)=-\frac{8(1-\varepsilon)}{3-2\varepsilon},
\]

\[
\mathcal S_A=I_E+I_V+a I_{\Pi},\qquad
\mathcal S_f=b I_{\Pi},\qquad
\widehat w_2=C_F S_2(\varepsilon)\,[C_A\mathcal S_A+T\mathcal S_f].
\]

The stored `CA_weight` and `T_weight` series are a I_r and b I_r for each self-energy graph. Use either those weighted series or the rational weights times the scalar series. Applying both would count the weights twice. The epsilon dependence of a and b is retained before truncation. No extra orbit factor, Euler factor or subtraction of W_1 squared over two is applied to this exponent sum.

## Scale and renormalization

\[
\mu_x^2=\frac{e^{2\gamma_E}}4\mu_{\overline{\mathrm{MS}}}^2,\quad
X_s=\mu_x^2(-s+i0),\quad X_t=\mu_x^2(-t+i0),\quad
L=\tfrac12(\log X_s+\log X_t),
\]

\[
S_2=X_s^{2\varepsilon}+X_t^{2\varepsilon}
=2e^{2\varepsilon L}\cosh(\varepsilon\ell),\qquad
\beta_0=(11C_A-4T)/3.
\]

The two-loop coefficient after coupling renormalization is B_2 = w_2-hat - beta_0 w_1-hat/(4 epsilon). Multiplicative contour renormalization uses W_B = Z_4 W_R and strict pole subtraction from log W_B. The symbols Z_(4,1) and Z_(4,2) are coefficients of log Z_4. The finite expansion is

\[
\log W_R=\lambda F_1+\lambda^2 F_2+O(\lambda^3),\qquad
W_R=1+\lambda F_1+\lambda^2(F_2+F_1^2/2)+O(\lambda^3).
\]

The supplied physical coefficient is F_2 at epsilon = 0. Its extraction uses the positive first order of the one-loop input. Extending the entire post-counterterm expression through epsilon squared would require the one-loop input through epsilon cubed; that extended expression is not the physical target of this supplement.

## Canonical symbols and special functions

| JSON and CSV symbol | Wolfram expression | Meaning |
| --- | --- | --- |
| `eps`, `z`, `ell`, `L` | `WilsonTwoLoop` context symbols with the same names | Regulator, ratio, shape logarithm and scale logarithm |
| `CA`, `CF`, `T` | Same names in `WilsonTwoLoop` | Color invariants and T_F n_f |
| `M_k_j` | ``WilsonTwoLoop`MomentM[k,j]`` | M moment defined below |
| `N_k_j` | ``WilsonTwoLoop`MomentN[k,j]`` | N moment defined below |
| `tauD4` | ``WilsonTwoLoop`TauD4[WilsonTwoLoop`z]`` | The weight-six function tau^+_(D,4) |
| `zeta(n)`, `pi` | `Zeta[n]`, `Pi` | Riemann zeta and pi |

Set r = z/(1+z), W(beta) = log[beta(1-beta)] and Q(beta,z) = ell + log[(1-beta)/beta]. Then

\[
M_{k,j}(z)=\int_0^1\frac{W(\beta)^k Q(\beta,z)^j}{\beta-r}\,d\beta,
\qquad k\ge0,\quad j\ge1,
\]

\[
N_{k,j}=M_{k,j+1}-\ell M_{k,j}.
\]

The apparent singularity at beta = r is removable because Q(r,z) = 0. Useful identities are M_(0,1) = -(ell squared + pi squared)/2 and M_(0,3) = -M_(0,1) squared. The data also include the classical-polylogarithm expression for M_(1,1).

The additional function in the epsilon-squared coefficient of 3e-3g-04 is fixed by

\[
D_z(t,u)=(1-t)(1-u)+ztu,\quad
X_+=\log\frac{D_z(t,u)tu}{\sqrt z},\quad
X_-=\log\frac{D_z(t,u)(1-t)(1-u)}{\sqrt z},
\]

\[
\tau^+_{D,4}(z)=-\frac{1+z}{48}
\int_0^1\!\int_0^1\frac{X_+^4+X_-^4}{D_z(t,u)}\,dt\,du.
\]

An explicit finite GPL expansion is included under `functions.tauD4_hlog_terms` and, in the Wolfram file, `functions.tauD4_GPL`. Each term consists of a word in 0 and -1, a rational coefficient, and powers of zeta(2), zeta(3), zeta(5), in that order. Terms with the same word are separate when their zeta monomials differ. The convention is

\[
G(a_1,\ldots,a_n;z)=\int_0^z\frac{dt}{t-a_1}G(a_2,\ldots,a_n;t),
\quad G(;z)=1,\quad G(0,\ldots,0;z)=\log^n(z)/n!.
\]

The Wolfram head `G[word,z]` is symbolic and uses this convention. It does not call an installed GPL package. Appendix A contains the same function definitions and grouped GPL coefficients.
