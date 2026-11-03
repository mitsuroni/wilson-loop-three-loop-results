# Independent controls of the fifteen insertion entries

The batch distinguishes three questions: whether the original source has been interpreted correctly, whether its geometry coefficients have been calculated with the stated normalization, and whether the full insertion operator has been convolved to the required physical order. The control program answers them with different checks; one successful layer is not promoted to a proof of the others.

Run against a newly generated, frozen candidate from the package directory:

```text
python -B scripts/independent_controls.py --results PATH_TO_GENERATED_JSON --out PATH_TO_NEW_CONTROL_JSON
```

The program imports no batch evaluator and opens no historical final-coefficient table. It uses its own restricted scalar parser, exact rational HLOG/period arithmetic and numerical HLOG evaluator. Sources are the separately frozen native full-colour C/J packages under snapshots/benchmarks/current/prepared/3loop. The output must be a new file.

## Start from the complete source tensor

For every source, f is its distance polynomial, d its eikonal scalar product and u,w its two contractions with the separation. The complete one-loop insertion numerator is checked as

\[
N_{\Pi_1}=(1+2\epsilon)df+(1-2\epsilon)uw.
\]

The metric part is 2(1+epsilon)df and the longitudinal part is -df+(1-2epsilon)uw. Both are retained and their sum is compared with the original se_numerator, followed by reconstruction of the complete J integrand. For the common chain/Pi2 geometry the corresponding numerator is

\[
N_{(3)}=(1+4\epsilon)df+(1-3\epsilon)uw.
\]

Its geometry is shared, while the external operators remain distinct: a chain has two Pi1 insertions, whereas Pi2 is one two-loop 1PI operator. The checks also retain the original dependent order, orbit multiplicity, one Euler compensation and current assembly convention. There is no runtime family sign correction.

## Two ways to integrate the original polynomial J

At epsilon=1 and 2, the actual J integrands become polynomials on compact ordered parameter domains. For all fifteen entries the checker evaluates that polynomial by two separate methods at z=1/2:

1. Follow the original outer-to-inner order in reverse when integrating, retaining every dependent upper bound; integrate monomials exactly.
2. Construct a new triangular map s=x times its original parent bound, calculate the full determinant, and integrate the resulting cube polynomial using its monomial moments.

This gives thirty paired exact identities. At these regulator values the full Gamma prefactor can have poles; only J is evaluated. These tests validate source/domain manipulations at the stated epsilon values. They do not determine a divergent integral's Laurent finite coefficient.

The opposite-kernel algebra is checked independently for the two exponent rates2 and3, retaining the differentiated distance polynomial and all four endpoint values. This tests the full numerator derivative identity, rather than just its metric part. Family derivations still supply the actual endpoint/meromorphic-continuation argument.

## Restore the full external colour operator

Every geometry coefficient is compared with the actual coloured coefficient after an independent exact convolution. The one-insertion factor is

\[
B_1(\epsilon)=\frac{(10-6\epsilon)C_A-8(1-\epsilon)T_Fn_f}{3-2\epsilon}.
\]

The eleven Pi1 entries retain both AA and Af channels. The chain retains AA, Af and ff from B1 squared, including the cross term factor2. The Pi2 entries retain AA, Af and Ff from

\[
\phi_c(\epsilon)=4\epsilon^2\left(A_c(\epsilon)+B_c(\epsilon)/U(\epsilon)\right).
\]

The Pi2 A/B functions are read from the accepted original full-D operator, not from a stored phi coefficient table. The checker independently expands the defining log-Gamma expression for inverse U through order5 before the cubic source poles cancel. This provides phi through order4. It then derives each convolution budget from the actual geometry pole bound; the largest one-insertion rational colour budget reaches order5. Missing lower Laurent positions are not assumed zero.

The arithmetic uses rational coefficients indexed by HLOG word and powers of zeta2,zeta3,zeta5. Multiplication by an operator coefficient adds only period powers; it does not use the evaluator's algebra implementation. All seven physical positions are checked for each nonzero channel, and extra declared zero channels must actually vanish. This tests physical convolution given a geometry solution; a common error in the geometry and all its coloured products is not excluded by convolution alone.

## Check finite terms directly where a Gamma formula is available

The four Pi1 P1/P2 entries and the adjacent chain/Pi2 geometries have source-derived adjacent-chain Gamma formulas. Their complete C-derived rational prefactor, orbit factor and colour projection are reconstructed independently. The two distinct signs of the native chain/Pi2 geometry definitions are consumed from their own C objects, rather than inferred from sharing a J.

The checker extracts the Laurent coefficients of these six Gamma formulas by a discrete Cauchy/Fourier contour, using96 nodes at each of two radii, 1/20 and1/40. It compares them with the candidate coefficients evaluated by a separate log-polynomial Taylor solution of the 0,-1 HLOG differential kernels at z=1/2. The independent HLOG evaluator first checks a dilogarithm identity. All seven coefficient positions are retained; the contour controls reach the finite term, not only the leading pole.

These are observed high-precision numerical residuals. Discrete contour aliasing, numerical roundoff and HLOG series tails are not rigorous interval enclosures. The original family endpoint arguments remain necessary for the source-to-Gamma connection.

## Negative controls and limits

The script rejects omission of the longitudinal numerator, doubled Euler compensation, treatment of Pi2 as two chain insertions, an unavailable dependent integration order, an altered physical finite coefficient and deletion of the chain's flavour-squared channel.

Coverage is explicit: all fifteen entries receive original-source polynomial and operator/convolution checks; six adjacent Gamma entries receive an additional finite-coefficient Cauchy control. No new original multidimensional finite-integral numerical control is claimed for the other nine geometries. Their admitted analytic moment/HLOG evaluations and historical checks retain their original scope. The internal two-loop Pi2 literature operator is an adopted input, not independently reintegrated here. The batch and this control record both have formal=not_run.
