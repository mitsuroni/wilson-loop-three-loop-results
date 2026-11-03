# Independent controls and their scope

The portable checker is [`independent_controls.py`](../scripts/independent_controls.py).
It reads the frozen source snapshots and, when supplied, the newly generated
candidate. It does not read the historical final-`I` records, import the P6
producer adapters, call Maple/Mathematica, or require the original workstation.
The result is an executable scientific control report, not a formal proof.

Run from the P6 working directory containing `package/` and `runs/`:

```powershell
python -X utf8 -B package/scripts/independent_controls.py `
  --results runs/generation-r001/generated.json `
  --output runs/generation-r001/independent-controls.json
```

Adjust those paths to the package/run locations in the installed delivery.
Python, SymPy and mpmath are sufficient. The source-only mode omits `--results`;
it explicitly records that the frozen-candidate checks were not run.

## Coverage contract

| Control | Actual coverage | What it establishes |
| --- | --- | --- |
| Gram geometry, cyclic numerator, path-ordered Fierz color | All 44 sources | Source algebra, color order and the declared kernel |
| Native/raw/prepared C/J, Euler count, C4 weight, dependent bounds | All 44 sources | Root and normalization continuity |
| Exact ordered-domain/projective polynomial moments | One representative in each of eight evaluation families | Domain/Jacobian implementation on nontrivial polynomial probes |
| Complete source reflection/inversion | 21 pairs | The claimed full-`J` relation, followed by the accepted crossing definition of `I` |
| Gamma/Cauchy extraction compared to candidate through epsilon zero | 01, 08, 17, 31 | Independent numerical check including finite `I` in these four cases |
| Bare-`F` Gamma check | Six raw moments shared by T1 | An independent check of this primitive |
| T4 recoil, complementary partition, factorization weight, absolute Laurent budget | Exact structural identities | Specific reduction and assembly requirements |
| Raw Omega defining derivative, crossing, boundary and individual wide-letter sectors | Omega6 and 27/28/29 | The explicit named function and its individual representation |

Replaying an original HyperInt primitive is not a second integration. The
44/308 generated Laurent positions and a successful comparison to a historical
table must not be reported as 44 independently integrated finite coefficients.
The checker performs no new original-domain finite quadrature for the other 40
orbits. RG, complex continuation and Lean integral proofs are outside its scope.

## Sources, ordering and normalization

The checker constructs the four-edge Gram matrix with
`v12=z`, `v13=-1-z`, `v14=1` and the corresponding cyclic entries. Positions
along an edge determine every squared difference. It reconstructs

\[
S=\sum_{i<j}\alpha_i\alpha_j (x_i-x_j)^2/2
\]

and the three cyclic ambient derivatives that give the complete numerator `P`.
The pure-propagator squared distance and tangent contraction are checked
separately. Their cube images are compared with the stored cube polynomials.
The absolute determinant is calculated from the actual `s_to_x` map. Three
exact interior-point comparisons, at epsilon 0, 1 and 2, verify the assembled
integrand on both sides. These are pointwise algebraic checks only.

The current convention is fixed by the staged
[`current-conventions.md`](../snapshots/workflow/current-conventions.md) and
[`shared-sign-conventions.md`](../snapshots/benchmarks/theory/montopar/shared-sign-conventions.md).
Slots are early to late; `U'=ig A U` puts the later generator on the left.
Accordingly the independent SU(N) Fierz closure uses the **reversed** Wilson
word and yields the actual native color `+i CF/4`. Using the unreversed word
would instead produce the historical `-i CF/4`. No fitted family multiplier is
introduced to reconcile those words.

All native `C` components, in their declared multiplicative order, give

\[
C_{\mathrm{native}}=
\frac{C_F}{16}\left(\frac{\alpha_s}{\pi}\right)^3
e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)\Gamma(2-2\epsilon)
\pi^{3\epsilon}(I_0-2v_{14})^{3\epsilon}.
\]

The raw components are preserved and exactly one Euler component is added by
preparation. The accepted legacy a1 root has the corresponding `CA^2 CF`
projection. The full finite-N relation
`CF=CA^2 CF-2 CA CF^2` is checked independently of the scalar C product.
All 44 source C4 multiplicities are four. This control does not derive a new
generic rule replacing the accepted per-family C-to-I definitions.

Every endpoint bound may depend only on an endpoint already available in the
declared outer-to-inner order. The alpha measure remains the ordered simplex,
with `alpha3=1-alpha1-alpha2`. In the projective chart

\[
(\alpha_1,\alpha_2,\alpha_3)
=(a,(1-a)b,(1-a)(1-b)),\qquad |\det|=1-a,
\]

the checker integrates the probes `1+S` and `P^2`, at `z=2/3`, in two ways:
directly in the original dependent endpoint order and simplex, and after the
cube/projective maps. It compares exact rational values in each of the eight
families. These are deliberately **not** claimed to be full-J values. In
particular, substituting a positive integer epsilon does not make the full
source polynomial: the factor `(alpha1 alpha2 alpha3)^(-epsilon)` is still
singular.

## Full source relations

The source pairs are checked before comparing any Laurent vectors:

- T2: 05/04, 06/03, 07/02, 08/01, 20/19, 21/18, 31/17,
  16/32, 36/12, 37/11.
- T3: 39/09, 38/10, 30/26, 25/22, 24/23.
- T4: 35/13, 34/14, 33/15, 43/40, 29/27, 42/41.

Each check reconstructs the reflected edge chains, reverses their endpoint
order, checks every ordered gap, applies the stated alpha permutation and
calculates the endpoint and simplex determinants. It checks `S`, `P`, the pure
distance, tangent contraction and complete C. T2 has `J_target=J_source`.
For T3/T4 the z inversion gives the factors `S:1/z`, `P:1/z^2`, pure distance
and tangent contraction `1/z`, hence
`J_target(z)=z^(3 epsilon) J_source(1/z)`. The accepted crossing combination
then gives equal physical `I` for those pairs. Candidate coefficient equality
is a subsequent 147-position control, not the evidence for the source map.

## Independent Gamma and finite-coefficient controls

Define

\[
\begin{aligned}
G&=e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)\Gamma(1-2\epsilon),\\
A&=\Gamma(1+\epsilon)\Gamma(1-\epsilon),\qquad
B=B(1-\epsilon,1-\epsilon),\\
V&=\frac{\Gamma(1+\epsilon)\Gamma(1+2\epsilon)}{\Gamma(1+3\epsilon)},\qquad
R=\frac{\cosh(\epsilon\log z/2)}{\cosh(3\epsilon\log z/2)}.
\end{aligned}
\]

The original Gamma-route formulas are

\[
I_{01}=\frac{G(2A-B)}{576\epsilon^6},\qquad
I_{17}=\frac{GR(2A-B)(1-V)}{128\epsilon^6}.
\]

The already checked source reflections supply 08 and 31. The checker samples
`epsilon^6 I` on two complex epsilon circles, of radii 0.065 and 0.085, at
192 equally spaced nodes and 80 decimal digits. Fourier/Cauchy extraction gives
all seven coefficients. An independent log-polynomial HLOG series evaluator
uses 220 terms at `z=1/2`; the full-I comparison tolerance is `1e-52`.
This does not call the generation adapters' Gamma jets or HLOG algebra.

The same method independently compares raw `bareF0` through `bareF5` with

\[
F_{\rm bare}(\epsilon)
=\frac{(1-\epsilon)B-A}{\epsilon(1-2\epsilon)}.
\]

The apparent pole cancels. This primitive check is additional coverage, not a
new finite-value check of every graph using the primitive.

## T4 requirements and Omega6

The T4 numerator obeys an exact derivative identity with a nonzero recoil term:

\[
P=\frac{-zS_x-(1+z)S_y+S_w}{4}
+\frac{(1+z)\alpha_3}{2}
 [\alpha_1(1-x+zw)-z\alpha_2w].
\]

Dropping the recoil is rejected. Sixteen monomial moments check the
complementary 41/42 ordered domains. A separate symbolic factorization identity
checks their actual `1/8` weight.

For the empty T4 family, Q and its even combination must retain the absolute
powers epsilon^-2 through epsilon^2, and the prefactor through epsilon^4.
Three terms counted from a pole are insufficient. The checker uses an explicit
Laurent witness and rejects a relative-order budget. It checks the actual
candidate budget metadata for all eight empty T4 members. The current six
corrected roots 13/14/15/33/34/35 remain r002; this bookkeeping check is not
relabelled as a fresh independent integration of their finite term.

The checker parses the original raw Phi5, Tprimitive and Omega6 with its own
restricted linear-GPL reader. It proves exactly that every Tprimitive word
starts with `-1/2` and that removing that first letter yields Phi5, including
all period coordinates. Thus

\[
T_6'(z)=\Phi_5(z)/(z+1/2).
\]

An independent numerical evaluator uses `r=z/(1+z)` and rational differential
kernels. The letters `-1/2` and `-2` transform to `-1` and `2`; all real
positive-z probes converge in this chart. At 65 digits and 380 terms it checks
the explicit Omega derivative at `z=1/2,1,2`, crossing at `1/2` and `2`, and
the defining boundary `Omega6(1)=2 T6(1)`. The constants A4/A6 are evaluated from
their convergent real one-dimensional defining integrals, not fitted.

The candidate must retain the full original Omega, Phi and T coordinates and
the individual wide-letter sectors: coefficient `-1/64` for 27 and 29 and
`+1/32` for 28. These checks do not use the cancellation of Omega in the total
color channel. An actual mutant adds a wide-letter contribution to 27 and
subtracts it from 28; the sum is unchanged while the individual control rejects
the mutant. This establishes representation consistency of Omega and its
individual occurrence; it is not independent full-domain finite integration
of 27 or 28.

## Executed rejection cases and residual limitations

The report lists only rejection cases that actually execute. They include
duplicate Euler compensation, wrong orbit multiplicity, wrong assembly sign,
unavailable dependent bounds, omitted alpha swap, omitted endpoint reversal,
flipped numerator sign, dropped T4 recoil, omitted factor `1/8`, insufficient
relative Laurent budget and, when a frozen candidate is supplied, an individual
Omega mutation hidden by a cancelling pair.

A complete source bridge, raw-primitive replay, historical-table comparison,
representation check, independent finite control and formal proof are separate
pieces of evidence. Passing this checker must preserve those distinctions in
the publication status. The actual report supplies numerical residuals and
runtime; this document does not substitute predicted success for that report.
