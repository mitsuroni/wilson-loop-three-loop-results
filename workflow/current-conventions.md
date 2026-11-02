# Single current physical convention

The unique convention ID is `wilson-qcd-action-path-fourier/1.0`. The current
reader rejects legacy sign modes and mixed contracts. This publication exposes
one physical convention; historical revisions are evidence identities.

The original action/path/Fourier derivation is preserved
[byte for byte](../sources/current-20261102/benchmarks/theory/montopar/shared-sign-conventions.md).
That archival Chinese document retains its original research-workspace links;
this page and the [current guide](../docs/current/README.md) provide portable
English navigation. Its [provenance](../sources/current-20261102/provenance.json)
records the exact source identity.

## Determine the sign at the physical source

Use Minkowski metric (+---), D=4-2epsilon, Hermitian generators,
[t^a,t^b]=i f^{abc}t^c and Tr(t^a t^b)=TF delta^{ab}, TF=1/2. The field Fourier
factor is exp(-ipX). The Wilson transport equation is U'=i g0 A U, with later
insertions on the left. Momentum Wick measure contributes +i; an internal
coordinate Wick measure contributes -i. Complete Pi/Gamma insertion operators
already include their defining i and internal Grassmann-loop signs; the external
Dyson coefficient is +1. Vertex labels, momentum directions and ordered domains
follow this same action and transport rule.

These signs enter C/J before parameter integration. Assembly does not decide
signs by diagram names, expected answers or RG residuals. Genuine color traces,
epsilon-dependent color factors, vertex phases, reflection multiplicities and
coupling counterterms may be negative. They are retained as physical factors.
For example the native Pi1 factor is

$$B_1(\epsilon)=\frac{(10-6\epsilon)C_A-8(1-\epsilon)T}{3-2\epsilon},
\qquad T=T_Fn_f.$$

No extra historical family minus sign accompanies native I.

## Keep preparation and physical assembly distinct

Raw C/J has zero Euler compensation. Native preparation checks multiplicative
C, Gamma/Euler weight and the contract, then adds exp(-n gamma_E epsilon) once
without changing J or its sign. Its full-color output is not automatically an
a1 projection or a derivation of connected I.

Each adopted physical I includes the orbit/C4, crossing and Euler factors
specified by its family. Do not restore them again. Four-leg contact and
exchange sources are grouped into joint objects with one consumption per
source. Connected exponent inputs do not undergo a second W-to-log W subtraction.

## Fix notation and scale conventions

| Symbol | Meaning |
| --- | --- |
| eps / epsilon | D=4-2epsilon |
| z, rho | Shape parameter, rho=log(z); fixed real z>0 in the general real calculation |
| L | (Ls+Lt)/2, with Ls=L+rho/2 and Lt=L-rho/2 |
| lambda | alpha_s_MSbar/pi |
| T | TF nf |
| M_kj, N_kj | Raw powers without factorials; N_kj=M_k,j+1-rho M_kj |
| S_n | exp[n epsilon(L+rho/2)]+exp[n epsilon(L-rho/2)] |
| U_n, B_n, CT_n, P_n | Scale-stripped sum, B_n=S_n U_n, coupling CT, P_n=B_n+CT_n |
| logZ, F_n | Negative-epsilon and epsilon0 parts of P_n; F_n is the coefficient in log W_R |

The local moment variable r_M=z/(1+z) and an older leg variable
r_leg=1/(1+z) have different definitions. Current total expressions use common
z-HLOG or the M/N family, rho and L. When reading a source representation,
retain that representation's declared auxiliary variable.

One-loop epsilon2 and two-loop epsilon1 enter three-loop finite counterterms.
Two-loop post-counterterm epsilon2 additionally requires one-loop epsilon3.
Positive epsilon coefficients are downstream input, not extra physical finite
terms. A statement that a finite result has no tau does not extend to every
positive order: the two-loop epsilon2 series contains a defined tau_D4 direction.

For ordinary W=1+lambda W1+lambda^2 W2+..., w1=W1 and w2=W2-W1^2/2.
Already connected I inputs assemble w_n directly. Each lower-loop scale sum
remains at its own loop order inside a counterterm.

## Evidence boundary

Current hashes and refusal checks establish exact version/contract identities.
They do not replace the original integral derivations. Historical results and
proofs retain their definitions and are not rebound to new files. Original
integrations, exact function identities and the admitted research solutions
are reused explicitly. General complex-z continuation is claimed only where
the original source proves it; no new whole-integral Lean proof is claimed.
RG is blind to arbitrary L-independent finite shapes. `formal=not_run`.
